% ========================================================================
%  SECTION 6: TRACKING CONFIGURATION
% ========================================================================

cfg = struct();

%  6.1 System Parameters
% cfg.dt = 2.0; % hardcode ping interval [s]
cfg.dt = median(seconds(diff(pingTimes))); % assume constant ping interval





cfg.sigmaR = 10e-2;     % Range measurement std [m]
% cfg.sigmaTheta = sind(2*1.4);  
cfg.kappaTheta = 1.5e-3;  

%  6.2 Process Noise (anisotropic)
cfg.sigmaAccel = [2.5 2.5 3.7]*1e-3;  % Acceleration std [m/s^2] for [x, y, z]

%  6.4 Initial Covariance of position and velocity
cfg.P0_pos = [.2 .2 .2].^2;   % Initial position variance [m^2]
cfg.P0_vel = [.3 .3 .1].^2;   % Initial velocity variance [(m/s)^2]

%  6.3 Gating Parameters
%               Probability less than the critical value
% P_G        0.90      0.95      0.975     0.99      0.999
% 3          6.251     7.815     9.348    11.345    16.266
cfg.chi2Gate = 7.815;  % Chi-squared threshold
cfg.costOfNonAssignment_CMKF = 1.5; 

%  6.5 Track Lifecycle Parameters
cfg.maxMisses  = 5;    % Maximum consecutive misses before deletion
cfg.Mconf = 4;
cfg.Nconf = 6;

%  6.9 Measurement Covariance in Split-Beam Coordinates
%      z = [r, alpha, beta] where alpha = sin(theta_a), beta = sin(theta_b)
% R_sb = @(r) diag([cfg.sigmaR^2, ...
%                   cfg.sigmaTheta^2, ...
%                   cfg.sigmaTheta^2]);
R_sb = @(r) diag([ cfg.sigmaR^2, ...
                  (cfg.kappaTheta*r)^2, ...
                  (cfg.kappaTheta*r)^2]);






%  6.7 Bias Correction Flag
cfg.enableBiasCorrection = false;  % Enable z-coordinate bias correction

%  6.8 Build Kalman Filter Matrices
F = buildStateTransitionMatrix(cfg.dt);
Q = buildProcessNoiseMatrix(cfg.dt, cfg.sigmaAccel);

%  6.10 CMKF Observation Matrix (extracts position from state)
%       H = [I_3 | 0_3]  since z_cart = [x, y, z] = H * [x, y, z, vx, vy, vz]'
H = [eye(3) zeros(3)];


% ========================================================================
%  SECTION 7: MULTI-TARGET TRACKING (CMKF-D)
%  ========================================================================

% Initialize track storage
activeTracks = struct('id', {}, 'state', {}, 'P', {}, ...
                      'age', {}, 'hits', {}, 'misses', {}, ...
                      'innovHist', {}, 'posLog', {}, 'peakLog', {}, ...
                      'quality', {}, 'velFit', {}, 'speedFit', {}, 'headingFit', {}, ...
                      'velFitSE', {}, 'velFitR2', {}, 'posFit0', {}, ...
                      'nistAT', {}, 'logDetSAT', {}, 'costAT', {});

unassignedDetPerPing = zeros(numPings, 1);
archivedTracks = activeTracks;
nextTrackID = 1;

tic;
for pingIdx = 1:numPings
    
    %  7.1 Get Detections from Current Ping
    currentDetections = detections(detections.pingIdx == pingIdx, :);
    
    numDetections = height(currentDetections);
    numTracks = length(activeTracks);
    
    %  7.2 Predict State and Covariance for All Active Tracks
    predictedStates = zeros(numTracks, 6);
    predictedCovariances = cell(numTracks, 1);
    predictedSplitBeam = zeros(numTracks, 3);  % Predicted (r, alpha, beta) for Jacobian evaluation

    for trackIdx = 1:numTracks
        % State prediction: x_pred = F * x_prev
        predictedStates(trackIdx, :) = (F * activeTracks(trackIdx).state')';
        % Covariance prediction: P_pred = F * P * F' + Q
        predictedCovariances{trackIdx} = F * activeTracks(trackIdx).P * F' + Q;
        
        % Compute predicted split-beam coordinates for Jacobian evaluation
        % Note: Jacobian is evaluated at predicted position (DCMKF approach)
        predPos = predictedStates(trackIdx, 1:3);
        predictedSplitBeam(trackIdx, :) = cartesianToSplitBeam(predPos);
    end
    
    %  7.3 Handle No-Detection Case (rare)
    if numDetections == 0
        for trackIdx = 1:numTracks
            activeTracks(trackIdx).state = predictedStates(trackIdx, :);
            activeTracks(trackIdx).P = predictedCovariances{trackIdx};
            activeTracks(trackIdx).age = activeTracks(trackIdx).age + 1;
            activeTracks(trackIdx).misses = activeTracks(trackIdx).misses + 1;
        end
        [activeTracks, archivedTracks] = ...
            pruneAndArchiveTracks(activeTracks, archivedTracks, cfg, pingIdx, numPings);
        continue % to next ping
    end
    
    %  7.4 Extract Measurements and Convert to Cartesian (with bias correction)
    alongAngles   = currentDetections.AoA_csac(:, 1);  % theta_a [deg]
    athwartAngles = currentDetections.AoA_csac(:, 2);  % theta_b [deg]
        
    % Split-beam measurement: range and direction cosines
    ranges = currentDetections.Range;       % r [m]
    alpha  = sind(alongAngles);             % alpha = sin(theta_a)
    beta   = sind(athwartAngles);           % beta  = sin(theta_b)
    gamma  = sqrt(1 - alpha.^2 - beta.^2);  % gamma = cos(theta_o)
    
    % Convert to Cartesian with bias correction: z_cart = [x, y, z]
    x_meas = ranges .* alpha;
    y_meas = ranges .* beta;
    z_meas = ranges .* gamma;
    
    % Apply z-bias correction if enabled
    % if cfg.enableBiasCorrection
    %     bias_z = computeZBias(ranges, gamma, cfg.sigmaTheta);
    %     z_meas = z_meas + bias_z;
    % end
    
    measurementsCartesian = [x_meas, y_meas, z_meas];
    
    %  7.5 Build and Solve Cost Matrix for Assignment (Gating in Cartesian)
    costMatrix = inf(numTracks, numDetections);
    
    for trackIdx = 1:numTracks
        % Predicted position (Cartesian)
        x_pred_pos = predictedStates(trackIdx, 1:3)';
        P_pred = predictedCovariances{trackIdx};
        
        % Compute converted measurement covariance at predicted state
        % Note: Jacobian evaluated at predicted position (DCMKF approach)
        r_pred     = predictedSplitBeam(trackIdx, 1);
        alpha_pred = predictedSplitBeam(trackIdx, 2);
        beta_pred  = predictedSplitBeam(trackIdx, 3);
        R_cart = computeConvertedMeasurementCovariance(r_pred, alpha_pred, beta_pred, R_sb(r_pred));
        
        % Innovation covariance: S = H * P_pred * H' + R_cart
        % Since H = [I_3 | 0_3], this simplifies to: S = P_pred(1:3,1:3) + R_cart
        P_pred_pos = P_pred(1:3, 1:3);
        S = P_pred_pos + R_cart;
        
        for detIdx = 1:numDetections
            % Converted measurement (Cartesian)
            z_cart = measurementsCartesian(detIdx, :)';
            
            % Innovation: nu = z_cart - x_pred_pos
            innovation = z_cart - x_pred_pos;
            
            % Squared Mahalanobis distance: d^2 = nu' * S^{-1} * nu
            mahalDistSq = innovation' / S * innovation;
            logDetS = log(det(S));
            
            % Gating: accept if d^2 <= gamma
            if mahalDistSq <= cfg.chi2Gate
                costMatrix(trackIdx, detIdx) = mahalDistSq + logDetS;
            end
        end
    end
    
    %  7.6 Solve Assignment Problem (Hungarian Algorithm)
    [assignments, unassignedTracks, unassignedDetections] = matchpairs(costMatrix, cfg.costOfNonAssignment_CMKF);
    
    unassignedDetPerPing(pingIdx) = length(unassignedDetections);
    
    %  7.7 \\\\\\\\\UPDATE///////// Assigned Tracks
    for pairIdx = assignments.'
        trackIdx = pairIdx(1);
        detIdx = pairIdx(2);
        
        % Converted measurement (Cartesian, with bias correction)
        z_cart = measurementsCartesian(detIdx, :)';
        
        % Predicted state and covariance
        state_pred = predictedStates(trackIdx, :)';
        P_pred = predictedCovariances{trackIdx};
        
        % Compute converted measurement covariance at predicted state
        r_pred     = predictedSplitBeam(trackIdx, 1);
        alpha_pred = predictedSplitBeam(trackIdx, 2);
        beta_pred  = predictedSplitBeam(trackIdx, 3);
        % R_cart     = computeConvertedMeasurementCovariance(r_pred, alpha_pred, beta_pred, R_sb);
        R_cart     = computeConvertedMeasurementCovariance(r_pred, alpha_pred, beta_pred, R_sb(r_pred));
        
        % Innovation: nu = z_cart - H * x_pred = z_cart - x_pred(1:3)
        innovation = z_cart - state_pred(1:3);
        
        % Innovation covariance: S = H * P_pred * H' + R_cart = P_pred(1:3,1:3) + R_cart
        S = P_pred(1:3, 1:3) + R_cart;

        % Compute NIS (Normalized Innovation Squared)
        NIS = innovation' / S * innovation;  % Same as innovation' * inv(S) * innovation
        logDetS = log(det(S));
        cost = NIS + log(det(S));

        % Kalman gain: K = P_pred * H' * S^{-1}
        % Since H = [I_3 | 0_3], H' = [I_3; 0_3], so P_pred * H' = P_pred(:,1:3)
        K = P_pred(:, 1:3) / S;
        
        % State update: x_new = x_pred + K * nu
        state_new = state_pred + K * innovation;
        
        % Covariance update (Joseph form for numerical stability):
        % P_new = (I - K*H) * P_pred * (I - K*H)' + K * R_cart * K'
        IKH = eye(6) - K * H;
        P_new = IKH * P_pred * IKH' + K * R_cart * K';
        P_new = (P_new + P_new') / 2;  % Ensure symmetry
        
        % Store updated track
        activeTracks(trackIdx).state = state_new';
        activeTracks(trackIdx).P = P_new;
        activeTracks(trackIdx).age = activeTracks(trackIdx).age + 1;
        activeTracks(trackIdx).hits = activeTracks(trackIdx).hits + 1;
        activeTracks(trackIdx).misses = 0;
        

        
        activeTracks(trackIdx).innovHist(end+1) = norm(innovation);
        activeTracks(trackIdx).posLog(end+1, :) = {pingIdx, ...
                                                   measurementsCartesian(detIdx, :), ...
                                                   state_pred(1:3)', ...
                                                   state_new(1:3)'};
        activeTracks(trackIdx).peakLog(end+1, :) = currentDetections(detIdx, :);

        activeTracks(trackIdx).nistAT(end+1) = NIS;
        activeTracks(trackIdx).costAT(end+1) = cost;
        activeTracks(trackIdx).logDetSAT(end+1) = logDetS;
    end
    
    
    %  7.8 Update Unassigned Tracks (Missed Detection / Coasting)
    for trackIdx = unassignedTracks'
        activeTracks(trackIdx).state = predictedStates(trackIdx, :);
        activeTracks(trackIdx).P = predictedCovariances{trackIdx};
        activeTracks(trackIdx).age = activeTracks(trackIdx).age + 1;
        activeTracks(trackIdx).misses = activeTracks(trackIdx).misses + 1;
    end
    
    
    %  7.9 INITIALIZE New Tracks from Unassigned Detections
    for detIdx = unassignedDetections'
        newTrack = struct();
        newTrack.id = nextTrackID;
        nextTrackID = nextTrackID + 1;
        
        % Initial position from converted measurement (with bias correction)
        pos = measurementsCartesian(detIdx, :);
        % Initial velocity set to zero (no prior kinematic information)
        newTrack.state = [pos, 0, 0, 0];
        
        % Initial covariance from measurement uncertainty
        % Transform R_splitbeam to Cartesian using Jacobian at measured position
        % R_cart_init = computeConvertedMeasurementCovariance(ranges(detIdx), ...
        %                 alpha(detIdx), beta(detIdx), R_sb(ranges(detIdx)));
        % P_pos = R_cart_init/2; 
        P_pos = diag(cfg.P0_pos);
        P_vel = diag(cfg.P0_vel);
        newTrack.P = blkdiag(P_pos, P_vel);
        
        newTrack.age = 1;
        newTrack.hits = 1;
        newTrack.misses = 0;
        
        newTrack.innovHist = [];
        newTrack.posLog = {pingIdx, pos, pos, pos};
        newTrack.peakLog = currentDetections(detIdx, :);
        
        newTrack.quality = 0;
        newTrack.velFit = [NaN, NaN, NaN];
        newTrack.speedFit = NaN;
        newTrack.headingFit = [NaN, NaN];
        newTrack.velFitSE = [NaN, NaN, NaN];
        newTrack.velFitR2 = [NaN, NaN, NaN];
        newTrack.posFit0 = [NaN, NaN, NaN];

        newTrack.nistAT = [];
        newTrack.costAT = [];
        newTrack.logDetSAT = [];
        
        activeTracks = [activeTracks, newTrack]; %#ok<*AGROW>
    end
    
    
    %  7.10 Track Management (Confirmation and Termination)
    [activeTracks, archivedTracks] = ...
        pruneAndArchiveTracks(activeTracks, archivedTracks, cfg, pingIdx, numPings);
end
trackingTime = toc;

fprintf('Tracking completed in %.2f seconds\n', trackingTime);
fprintf('Active tracks: %d, Archived tracks: %d\n', length(activeTracks), length(archivedTracks));




visualize_ALL_Tracks
% visualize_ONE_Track


%% ========================================================================
%  HELPER FUNCTIONS: TRACK MANAGEMENT
%  ========================================================================

function [tracks, archivedTracks] = pruneAndArchiveTracks(tracks, archivedTracks, config, pingIdx, numPings)
    % Prune tracks based on hit ratio and miss count, archive qualified tracks
    
    toDelete = false(1, length(tracks));
    
    for trackIdx = 1:length(tracks)
        track = tracks(trackIdx);
        
        % % Delete young tracks with poor hit ratio
        % if track.age <= config.youngAge && track.age > 0
        %     hitRatio = track.hits / track.age;
        %     if hitRatio < config.minHitRatio
        %         toDelete(trackIdx) = true;
        %         continue
        %     end
        % end
        % 
        % % Delete tracks exceeding maximum misses
        % if track.misses >= config.maxMisses
        %     if track.hits >= config.minHits
        %         archivedTracks = [archivedTracks, track];
        %     end
        %     toDelete(trackIdx) = true;
        % end

        % --- Confirmation gate: delete tentative tracks that fail M/N ---
        if track.age == config.Nconf && track.hits < config.Mconf
            toDelete(trackIdx) = true;
            continue
        end

        % --- Deletion: remove tracks exceeding max consecutive misses ---
        if track.misses >= config.maxMisses
            if track.age >= config.Nconf   % confirmed → archive
                track.age = track.age - track.misses; % 3/8/2026
                archivedTracks = [archivedTracks, track];
            end
            toDelete(trackIdx) = true;
        end


    end
    
    tracks(toDelete) = [];
    
    % Archive remaining tracks at last ping
    if pingIdx == numPings && ~isempty(tracks)
        shouldArchive = [tracks.hits] >= config.Mconf;
        if any(shouldArchive)

            trackInds = find(shouldArchive);
            for i = 1:length(trackInds)
                trackIdx = trackInds(i);
                tracks(trackIdx).age = tracks(trackIdx).age - tracks(trackIdx).misses;
            end
            
            archivedTracks = [archivedTracks, tracks(shouldArchive)];
        end
    end


end



%% ========================================================================
%  CMKF MEASUREMENT CONVERSION FUNCTIONS
%  ========================================================================

function z_sb = cartesianToSplitBeam(position)
    % Convert Cartesian position to split-beam coordinates
    %
    % Input:
    %   position - [x, y, z] Cartesian position (1x3 or 3x1)
    %
    % Output:
    %   z_sb - [r, alpha, beta] split-beam coordinates
    %        r = range
    %        alpha = x/r (direction cosine, alongship)
    %        beta = y/r (direction cosine, athwartship)
    
    x = position(1);
    y = position(2);
    z = position(3);
    
    r = sqrt(x^2 + y^2 + z^2);
    
    if r < 1e-6
        z_sb = [0, 0, 0];
        return
    end
    
    alpha = x / r;
    beta = y / r;
    
    z_sb = [r, alpha, beta];
end


function bias_z = computeZBias(r, gamma, sigma_theta)
    % Compute the bias in the z-coordinate due to nonlinear transformation
    %
    % The z-coordinate is computed as z = r * gamma, where gamma = sqrt(1 - alpha^2 - beta^2).
    % Due to the concavity of the square root, measurement noise in alpha and beta
    % causes a negative bias in gamma (and thus z), making targets appear closer
    % to the transducer than they actually are.
    %
    % This function computes the additive bias correction:
    %   bias_z = r * sigma_theta^2 / (2 * gamma^3) * (1 + gamma^2)
    %
    % Input:
    %   r           - range [m] (scalar or vector)
    %   gamma       - sqrt(1 - alpha^2 - beta^2), cosine of off-axis angle (scalar or vector)
    %   sigma_theta - direction cosine measurement std [dimensionless] (scalar)
    %
    % Output:
    %   bias_z - bias correction to add to z [m] (same size as r)
    %
    % Reference:
    %   Lerro & Bar-Shalom (1993), adapted for split-beam coordinates
    
    % Ensure gamma is bounded away from zero to avoid singularity at beam edge
    % gamma = max(gamma, 0.1);
    
    % Bias formula: mu_z = r * sigma_theta^2 / (2 * gamma^3) * (1 + gamma^2)
    bias_z = r .* sigma_theta^2 ./ (2 * gamma.^3) .* (1 + gamma.^2);
end


function R_cart = computeConvertedMeasurementCovariance(r, alpha, beta, R_splitbeam)
    % Compute converted measurement covariance in Cartesian coordinates
    %
    % The measurement conversion is:
    %   x = r * alpha
    %   y = r * beta  
    %   z = r * gamma, where gamma = sqrt(1 - alpha^2 - beta^2)
    %
    % Input:
    %   r - range [m]
    %   alpha - direction cosine (alongship), alpha = sin(theta_a)
    %   beta - direction cosine (athwartship), beta = sin(theta_b)
    %   R_splitbeam - 3x3 measurement covariance in split-beam coordinates
    %                 diag([sigma_r^2, sigma_theta^2, sigma_theta^2])
    %
    % Output:
    %   R_cart - 3x3 converted measurement covariance in Cartesian coordinates
    %
    % Note: The Jacobian is evaluated at the position specified by (r, alpha, beta).
    %       For DCMKF (decorrelated CMKF), this should be the predicted position,
    %       not the measured position, to avoid correlation between the converted
    %       measurement covariance and the measurement noise.
    
    % gamma = sqrt(1 - alpha^2 - beta^2) = cos(off-axis angle)
    % gamma = sqrt(max(1e-6, 1 - alpha^2 - beta^2));
    gamma = sqrt(1 - alpha^2 - beta^2);
    
    % Jacobian of conversion: G = d(x,y,z) / d(r,alpha,beta)
    %
    % x = r * alpha      =>  dx/dr = alpha,   dx/dalpha = r,             dx/dbeta = 0
    % y = r * beta       =>  dy/dr = beta,    dy/dalpha = 0,             dy/dbeta = r
    % z = r * gamma      =>  dz/dr = gamma,   dz/dalpha = -r*alpha/gamma, dz/dbeta = -r*beta/gamma
    
    G = [alpha,  r,              0;
         beta,   0,              r;
         gamma, -r*alpha/gamma, -r*beta/gamma];
    
    R_cart = G * R_splitbeam * G';

    % Ensure symmetry
    R_cart = (R_cart + R_cart') / 2;
end


function F = buildStateTransitionMatrix(dt)
    % Build 6x6 state transition matrix for constant-velocity model
    %
    % State: x = [x, y, z, vx, vy, vz]'
    % Model: position(k) = position(k-1) + velocity(k-1) * dt
    %        velocity(k) = velocity(k-1)
    %
    % F = [I_3, dt*I_3]
    %     [0_3,    I_3]
    
    F = [  eye(3),  dt * eye(3);
         zeros(3),       eye(3)];
end


function Q = buildProcessNoiseMatrix(dt, sigmaAccel)
    % Build 6x6 process noise covariance for piecewise-constant white noise acceleration
    %
    % Input:
    %   dt - time step [s]
    %   sigmaAccel - sqrt[q_x, q_y, q_z] acceleration noise std [m/s^2] for each axis
    %
    % The process noise covariance is derived assuming piecewise constant
    % acceleration over each time interval:
    %
    %   Q = [dt^4/4 * Sigma_a,  dt^3/2 * Sigma_a]
    %       [dt^3/2 * Sigma_a,  dt^2   * Sigma_a]
    %
    % where Sigma_a = diag(q_x^2, q_y^2, q_z^2) is the acceleration noise
    % power spectral density matrix.
    %
    % Reference: Bar-Shalom et al. (2001), Estimation with Applications to 
    %            Tracking and Navigation
    
    Sigma_a = diag(sigmaAccel.^2);
    
    Q = [dt^4/4 * Sigma_a,  dt^3/2 * Sigma_a;
         dt^3/2 * Sigma_a,  dt^2   * Sigma_a];
end
