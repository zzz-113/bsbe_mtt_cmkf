function [tracks, stats] = trackTargetsCMKF(detections, env, cfg)
%TRACKTARGETSCMKF  Multi-target tracking with converted split-beam measurements.
%
%   [TRACKS, STATS] = trackTargetsCMKF(DETECTIONS, ENV, CFG) runs a converted
%   measurement Kalman filter over the detections produced by detectTargets and
%   returns the confirmed tracks.
%
%   Each detection is converted from split-beam coordinates [r, alpha, beta] to
%   Cartesian [x, y, z] before filtering, so the filter itself stays linear:
%   the observation matrix is H = [I3 0] and only the measurement covariance
%   carries the nonlinearity, through convertedMeasurementCovariance. The
%   conversion Jacobian is evaluated at the predicted position, which keeps the
%   converted covariance independent of the current measurement noise.
%
%   Per ping the filter
%     1. predicts every active track;
%     2. gates detections against each track on the Mahalanobis distance in
%        Cartesian space, at the cfg.tracking.chi2Gate threshold;
%     3. solves the resulting assignment problem with matchpairs;
%     4. updates assigned tracks (Joseph form), coasts unassigned tracks, and
%        starts a new tentative track from every unassigned detection;
%     5. confirms, deletes and archives tracks.
%
%   Inputs
%     detections - table from detectTargets
%     env        - geometry struct from detectTargets (uses pingInterval, numPings)
%     cfg        - configuration struct from defaultConfig
%
%   Outputs
%     tracks - struct array of confirmed tracks:
%        id         track identifier
%        state      final [x y z vx vy vz]
%        P          final 6x6 state covariance
%        age        number of pings the track spanned, excluding trailing coasts
%        hits       number of pings with an assigned detection
%        misses     consecutive missed updates at the point of archiving
%        posLog     {pingIndex, measured, predicted, updated} position per hit
%        peakLog    the detection table rows assigned to this track
%        innovHist  norm of the innovation per hit [m]
%        nisHist    normalised innovation squared per hit
%        logDetHist log(det(S)) per hit
%        costHist   assignment cost per hit, nisHist + logDetHist
%
%     stats - elapsedTime, numConfirmed, numStillActive and
%             unassignedPerPing, the count of detections that started a new
%             track in each ping
%
%   Only tracks that reach confirmation are returned. Tracks still alive at the
%   last ping are archived if they have more than cfg.tracking.Mconf hits.

arguments
    detections table
    env struct
    cfg struct
end

t = cfg.tracking;

dt = env.pingInterval;
[F, Q] = constantVelocityModel(dt, t.sigmaAccel);

% The converted measurement is the position itself, so the observation matrix
% is constant and the products below simplify: H*P*H' = P(1:3,1:3) and
% P*H' = P(:,1:3). H is kept explicitly for the Joseph-form update.
H = [eye(3) zeros(3)];

% Measurement covariance in split-beam coordinates. The angular terms scale
% with range because the direction-cosine error grows with distance.
Rsb = @(r) diag([t.sigmaR^2, (t.kappaTheta*r)^2, (t.kappaTheta*r)^2]);

TRACK_FIELDS = {'id', 'state', 'P', 'age', 'hits', 'misses', ...
                'posLog', 'peakLog', 'innovHist', 'nisHist', 'logDetHist', 'costHist'};
emptyTrack = [TRACK_FIELDS; repmat({{}}, 1, numel(TRACK_FIELDS))];
active   = struct(emptyTrack{:});
archived = active;

unassignedPerPing = zeros(env.numPings, 1);
nextTrackID = 1;

timer = tic;
for pingIdx = 1:env.numPings

    pingDetections = detections(detections.pingIndex == pingIdx, :);
    numDetections  = height(pingDetections);
    numTracks      = numel(active);

    %% Predict ----------------------------------------------------------
    predState = zeros(numTracks, 6);
    predCov   = cell(numTracks, 1);
    predSB    = zeros(numTracks, 3);   % predicted [r alpha beta], for the Jacobian

    for k = 1:numTracks
        predState(k, :) = (F * active(k).state')';
        predCov{k}      = F * active(k).P * F' + Q;
        predSB(k, :)    = cartesianToSplitBeam(predState(k, 1:3));
    end

    %% No detections in this ping: coast every track ---------------------
    if numDetections == 0
        for k = 1:numTracks
            active(k) = coast(active(k), predState(k, :), predCov{k});
        end
        [active, archived] = manageTracks(active, archived, t, pingIdx, env.numPings);
        continue
    end

    %% Convert measurements to Cartesian ---------------------------------
    ranges = pingDetections.range;
    alpha  = sind(pingDetections.aoa(:, 1));      % alongship direction cosine
    beta   = sind(pingDetections.aoa(:, 2));      % athwartship direction cosine
    gamma  = sqrt(1 - alpha.^2 - beta.^2);        % cosine of the off-axis angle

    measurements = [ranges .* alpha, ranges .* beta, ranges .* gamma];

    %% Gate and score every track/detection pair -------------------------
    % The converted covariance and the innovation covariance depend only on
    % the prediction, so both are formed once per track and reused by the
    % update below.
    costMatrix = inf(numTracks, numDetections);
    predRcart  = cell(numTracks, 1);
    predS      = cell(numTracks, 1);

    for k = 1:numTracks
        predRcart{k} = convertedMeasurementCovariance(predSB(k, 1), predSB(k, 2), ...
                                                      predSB(k, 3), Rsb(predSB(k, 1)));
        predS{k} = predCov{k}(1:3, 1:3) + predRcart{k};
        logDetS  = log(det(predS{k}));

        for j = 1:numDetections
            innovation = measurements(j, :)' - predState(k, 1:3)';
            nis = innovation' / predS{k} * innovation;

            % Adding log(det(S)) makes costs comparable between tracks with
            % different innovation covariances, so matchpairs is not biased
            % towards assigning to the track with the widest gate.
            if nis <= t.chi2Gate
                costMatrix(k, j) = nis + logDetS;
            end
        end
    end

    [assignments, unassignedTracks, unassignedDetections] = ...
        matchpairs(costMatrix, t.costOfNonAssignment);

    unassignedPerPing(pingIdx) = numel(unassignedDetections);

    %% Update assigned tracks --------------------------------------------
    for pair = assignments'
        k = pair(1);
        j = pair(2);

        statePred = predState(k, :)';
        P         = predCov{k};
        Rcart     = predRcart{k};
        S         = predS{k};

        innovation = measurements(j, :)' - statePred(1:3);

        nis     = innovation' / S * innovation;
        logDetS = log(det(S));

        K = P(:, 1:3) / S;
        stateNew = statePred + K * innovation;

        % Joseph form: stays positive semi-definite even when K is not the
        % exact optimal gain, which matters here because Rcart is itself an
        % approximation.
        IKH = eye(6) - K * H;
        Pnew = IKH * P * IKH' + K * Rcart * K';
        Pnew = (Pnew + Pnew') / 2;

        active(k).state  = stateNew';
        active(k).P      = Pnew;
        active(k).age    = active(k).age + 1;
        active(k).hits   = active(k).hits + 1;
        active(k).misses = 0;

        active(k).posLog(end+1, :) = {pingIdx, measurements(j, :), ...
                                      statePred(1:3)', stateNew(1:3)'};
        active(k).peakLog(end+1, :) = pingDetections(j, :);

        active(k).innovHist(end+1)  = norm(innovation);
        active(k).nisHist(end+1)    = nis;
        active(k).logDetHist(end+1) = logDetS;
        active(k).costHist(end+1)   = nis + logDetS;
    end

    %% Coast unassigned tracks -------------------------------------------
    for k = unassignedTracks'
        active(k) = coast(active(k), predState(k, :), predCov{k});
    end

    %% Start a tentative track from each unassigned detection -------------
    for j = unassignedDetections'
        newTrack = struct();
        newTrack.id = nextTrackID;
        nextTrackID = nextTrackID + 1;

        % Velocity is unobservable from a single detection, so it starts at
        % zero and P0vel is set wide enough to cover the expected speeds.
        pos = measurements(j, :);
        newTrack.state = [pos, 0, 0, 0];
        newTrack.P     = blkdiag(diag(t.P0pos), diag(t.P0vel));

        newTrack.age    = 1;
        newTrack.hits   = 1;
        newTrack.misses = 0;

        newTrack.posLog  = {pingIdx, pos, pos, pos};
        newTrack.peakLog = pingDetections(j, :);

        newTrack.innovHist  = [];
        newTrack.nisHist    = [];
        newTrack.logDetHist = [];
        newTrack.costHist   = [];

        active(end+1) = orderfields(newTrack, TRACK_FIELDS); %#ok<AGROW>
    end

    %% Confirm, delete, archive -------------------------------------------
    [active, archived] = manageTracks(active, archived, t, pingIdx, env.numPings);
end

tracks = archived;
stats = struct( ...
    'elapsedTime',       toc(timer), ...
    'numConfirmed',      numel(archived), ...
    'numStillActive',    numel(active), ...
    'unassignedPerPing', unassignedPerPing);
end


function track = coast(track, predictedState, predictedCov)
% Advance a track that received no detection: keep the prediction and count
% the miss. No innovation is logged, so the *Hist arrays index hits only.
track.state  = predictedState;
track.P      = predictedCov;
track.age    = track.age + 1;
track.misses = track.misses + 1;
end


function [tracks, archived] = manageTracks(tracks, archived, t, pingIdx, numPings)
% Apply the M/N confirmation rule and the deletion rule, and archive tracks
% that are worth keeping.

toDelete = false(1, numel(tracks));

for k = 1:numel(tracks)
    track = tracks(k);

    % Confirmation: a track that reaches age Nconf with fewer than Mconf hits
    % is discarded without being archived.
    if track.age == t.Nconf && track.hits < t.Mconf
        toDelete(k) = true;
        continue
    end

    % Deletion after maxMisses consecutive coasts. Archive only tracks that
    % lived long enough to be confirmed.
    if track.misses >= t.maxMisses
        if track.age >= t.Nconf
            % Trailing coasts carry no measurement, so roll age back to the
            % last hit before archiving.
            track.age = track.age - track.misses;
            archived(end+1) = track; %#ok<AGROW>
        end
        toDelete(k) = true;
    end
end

tracks(toDelete) = [];

% Tracks still alive when the recording ends were never given the chance to
% fail the deletion rule, so they are archived on hit count alone.
if pingIdx == numPings && ~isempty(tracks)
    keep = [tracks.hits] >= t.Mconf;
    for k = find(keep)
        tracks(k).age = tracks(k).age - tracks(k).misses;
    end
    archived = [archived, tracks(keep)];
end
end
