% ========================================================================
%  SECTION 8: POST-PROCESSING AND QUALITY ASSESSMENT
%  ========================================================================




%  8.1 Compute Track Quality Scores

for trackIdx = 1:length(archivedTracks)
    % archivedTracks(trackIdx).quality = computeTrackQuality(archivedTracks(trackIdx), cfg); %#ok<SAGROW>
    archivedTracks(trackIdx) = computeLinearFitVelocity(archivedTracks(trackIdx), cfg);
end

%  8.2 Filter Tracks by Quality
% highQualityTracks = archivedTracks([archivedTracks.quality] >= 0.65);
% highQualityTracks = archivedTracks([archivedTracks.hits] > 10);

highQualityTracks = archivedTracks([archivedTracks.hits] > 10);
% highQualityTracks(reviewResults.badIdx) = [];



% for i = 1:length(highQualityTracks)
%     nisMeanAlongTrack_ALL_Tracks(i) = std(highQualityTracks(i).nisHist);
%     costMeanAlongTrack_ALL_Tracks(i) = std(highQualityTracks(i).costHist);
% end
% 
% figure
% histogram(nisMeanAlongTrack_ALL_Tracks,50)
% figure
% histogram(costMeanAlongTrack_ALL_Tracks,50)
% % i = find(nisMeanAlongTrack_ALL_Tracks >= 3.9)
% % highQualityTracks(i) = [];


%  8.4 Print Summary
% disp(cfg);
fprintf('Total Detections:     %d\n', size(detections, 1));
fprintf('Confirmed Tracks:     %d\n', length(archivedTracks));
fprintf('High-Quality Tracks:  %d\n', length(highQualityTracks));
qualityScores = [archivedTracks.quality];
fprintf('Quality - Min: %.2f | Median: %.2f | Max: %.2f\n', min(qualityScores), median(qualityScores), max(qualityScores));

speeds = [highQualityTracks.speedFit];
r2Values = vertcat(highQualityTracks.velFitR2);
fprintf('Speed (fit) - Min: %.3f | Median: %.3f | Max: %.3f m/s\n', min(speeds), median(speeds), max(speeds));
fprintf('R² (fit) - Median [x,y,z]: [%.3f, %.3f, %.3f]\n', median(r2Values(:,1)), median(r2Values(:,2)), median(r2Values(:,3)));


%% ========================================================================
%  SECTION 9: VISUALIZATION - DETECTION AND TRACKING OVERVIEW
%  ========================================================================

figure('Units', 'inches', 'Position', [2, 2, 6.5, 3]*2)
tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact')


%  9.1 Echogram
ax1 = nexttile;
imagesc(pingTimes, rangeVector, db((rangeVector.^2) .* squeeze(mean(echoData, 2)) / 5.58^2 ));
axis xy
clim([-95 -45])
colormap(ax1, viridis)
xsecondarylabel('');
ylabel('Depth (m)')
title('(a) \itS_p')




%  9.2 Detection Scatter Plot
ax2 = nexttile;
scatter(detections.pingTime, detections.Range, detections.P2B/30, 'x','MarkerEdgeColor','#008080');
box on
title('(b) Detections')
xsecondarylabel('');
% title(sprintf('Detections (N=%d)', size(detections, 1)))
% ylabel('Range (m)')




%  9.3 3D Track Visualization
ax3 = nexttile;
colors = orderedcolors('gem');
nColors = size(colors,1);
hold on
for trackIdx = 1:length(highQualityTracks)

    pos = cell2mat(highQualityTracks(trackIdx).posLog(:, 4));
    plot3(pos(:,1), pos(:,2), pos(:,3), 'Color', colors(mod(trackIdx,nColors)+1,:));
    
end
daspect([1, 1, 0.7]);
axis([-10, 10, -10, 10, minDetectionRange, rangeVector(end)])
axis vis3d
% view([250, 10])
grid on; box on
xlabel('x (m)'); ylabel('y (m)'); 
% zlabel('Depth (m)')
zticks(0:10:60)
zticklabels(string(590 - (0:10:60)))
title('(c) 3D Tracks')




%  9.4 3D Track Velocity
ax4 = nexttile;
hold on
for trackIdx = 1:length(highQualityTracks)
    pos = cell2mat(highQualityTracks(trackIdx).posLog(:, 4));
    vel = highQualityTracks(trackIdx).velFit;

    pos = double(pos);
    vel = double(vel)*1e2;

    quiver3(0,0, median(pos(:,3)),...
        vel(1),vel(2),vel(3),'Color',colors(mod(trackIdx,nColors)+1,:),'AutoScale','off','ShowArrowHead','off');

    % z(trackIdx) = median(pos(:,3));
    % u_obs(trackIdx) = vel(1);
    % v_obs(trackIdx) = vel(2);
    
end
daspect([1, 1, 0.7]);
axis([-10, 10, -10, 10, minDetectionRange, rangeVector(end)])
set([ax3 ax4], 'Projection', 'Perspective', 'CameraViewAngle', 15,'CameraPosition',[-100 50 30])
axis vis3d
% view([250, 10])
grid on; box on
xlabel('u (cm/s)'); ylabel('v (cm/s)'); 
% zlabel('Depth (m)')
zticks(0:10:60)
zticklabels(string(590 - (0:10:60)))
title('(d) Velocity')


% addpath asl/
% [z, zi] = sort(z,'ascend');
% u_obs = u_obs(zi);
% v_obs = v_obs(zi);
% % [u_cv, v_cv, lam_opt] = smooth_uv_cv(z, u_obs, v_obs, 2);
% [u_cv, v_cv] = smooth_uv_profile(z, u_obs, v_obs, 1e3, 2);
% plot3(u_cv, v_cv, z, 'k', 'LineWidth', 1)




%  9.5 Track Validation
ax5 = nexttile;
imagesc(pingTimes, rangeVector, db((rangeVector.^2) .* squeeze(mean(echoData, 2))), 'AlphaData', 0.5);
axis xy; hold on
clim([-80, -10] + 36)
colormap(ax5, flipud(gray))

for trackIdx = 1:length(highQualityTracks)

    trackTime = highQualityTracks(trackIdx).peakLog.pingTime;
    trackRange = highQualityTracks(trackIdx).peakLog.Range;

    trackInnov = [0.1 highQualityTracks(trackIdx).innovHist];
    % trackNIS = [0.1 highQualityTracks(trackIdx).nistHist];
    % trackCost = [0.1 highQualityTracks(trackIdx).costHist];
    highQualityTracks(trackIdx).medianRange = median(trackRange);


    % meanNIS = mean(highQualityTracks(trackIdx).nistAT);
    % z = (meanNIS -1) / sqrt(6 / highQualityTracks(trackIdx).hits);
    % quality = 1 - normcdf(z);

    % q(trackIdx) = assessTrackQuality(highQualityTracks(trackIdx));

    % astd = sum(highQualityTracks(trackIdx).peakLog.AoA_csac_std(:,[1 2]),2);
    % scatter(trackTime, trackRange, astd*5,colors(mod(trackIdx,nColors)+1,:),'filled')

    

    plot(trackTime, trackRange, 'LineWidth', 0.88, ...
        'Color', colors(mod(trackIdx,nColors)+1,:))
end

% ylabel('Range (m)')
% title('Validation')
% Label with the parameters the filter ACTUALLY uses. R_sb is built from
% kappaTheta, not sigmaTheta, so report kappaTheta here -- otherwise the
% figure and the filename advertise a parameter that never enters R.
sigmaAccel = cfg.sigmaAccel;
title_name = sprintf(['N=%d, \\sigma_r=%.0f cm, \\kappa_\\theta=%.2g, ' ...
    'c_{na}=%.2f, \\sigma_a=[%.3g,%.3g,%.3g]'], ...
    length(highQualityTracks), ...
    cfg.sigmaR*1e2, ... % sigma_r [cm]
    cfg.kappaTheta, ...
    cfg.costOfNonAssignment_CMKF, ...
    sigmaAccel(1),sigmaAccel(2),sigmaAccel(3));

% title(title_name, 'Interpreter','tex')
title('(e) Validation')

xlim([ax1, ax2, ax5], pingTimes([1, end]))
ylim([ax1, ax2, ax5], [minDetectionRange, rangeVector(end)])

xsecondarylabel('');
yticks([ax1, ax2, ax5], 0:10:60)
yticklabels([ax1, ax2, ax5], string(590 - (0:10:60)))
% set([ax1, ax2, ax5], 'TickDir','out')
yticklabels([ax2 ax5],'')
zticklabels([ax3 ax4],'')



save_file_name = sprintf('N%d_r%.0f_kappa%.2g_cna%.2f_a%.3g-%.3g-%.3g', ...
    length(highQualityTracks), ...
    cfg.sigmaR*1e2, ... % sigma_r [cm]
    cfg.kappaTheta, ...
    cfg.costOfNonAssignment_CMKF, ...
    sigmaAccel(1),sigmaAccel(2),sigmaAccel(3))



set(findall(gcf, '-property', 'FontSize'), 'FontName','Arial', 'FontWeight','Normal');
set(findall(gcf, 'type', 'axes'), 'TitleHorizontalAlignment', 'left','FontSize', 18);



%%
function quality = computeTrackQuality(track, config)
    % Compute track quality score based on length and hit ratio
    LENGTH_SCALE = 10;
    lengthScore = min(1, max(0, (track.hits - config.Mconf) / LENGTH_SCALE));
    hitRatioScore = min(1, max(0, (track.hits / track.age - 0.5) / 0.5));
    quality = max(0, min(1, 0.5 * lengthScore + 0.5 * hitRatioScore));
end

% function quality = computeTrackQuality(track, config)
%     L = length(track.nistAT);
%     if L < 2
%         quality = 0;
%         return
%     end
% 
%     meanNIS = mean(track.nistAT);
% 
%     % Under correct model: E[NIS] = 3, Var[mean NIS] = 6/L
%     % z-score measures how many standard deviations from expected
%     z = (meanNIS - 3) / sqrt(6 / L);
% 
%     % Convert to quality: z=0 is perfect, large |z| is poor
%     % One-sided: we mainly care about z >> 0 (model mismatch)
%     % Map to [0, 1] via complementary Gaussian CDF
%     quality = 1 - normcdf(z);
% end






function track = computeLinearFitVelocity(track, config)
    % Compute velocity using linear least-squares fit to position history
    
    numPoints = size(track.posLog, 1);
    
    if numPoints < 2
        track.velFit = [NaN, NaN, NaN];
        track.speedFit = NaN;
        track.headingFit = [NaN, NaN];
        track.velFitSE = [NaN, NaN, NaN];
        track.velFitR2 = [NaN, NaN, NaN];
        track.posFit0 = [NaN, NaN, NaN];
        return
    end
    
    pingNumbers = cell2mat(track.posLog(:, 1));
    positions = cell2mat(track.posLog(:, 2));
    timeRelative = (pingNumbers - pingNumbers(1)) * config.dt;
    
    velocity = zeros(1, 3);
    position0 = zeros(1, 3);
    standardError = zeros(1, 3);
    rSquared = zeros(1, 3);
    
    for dim = 1:3
        y = positions(:, dim);
        X = [timeRelative(:), ones(numPoints, 1)];
        coefficients = X \ y;
        velocity(dim) = coefficients(1);
        position0(dim) = coefficients(2);
        
        yPredicted = X * coefficients;
        residuals = y - yPredicted;
        ssResidual = sum(residuals.^2);
        ssTotal = sum((y - mean(y)).^2);
        
        if ssTotal > 0
            rSquared(dim) = 1 - ssResidual / ssTotal;
        else
            rSquared(dim) = NaN;
        end
        
        if numPoints > 2
            mse = ssResidual / (numPoints - 2);
            XtXinv = inv(X' * X);
            standardError(dim) = sqrt(mse * XtXinv(1, 1));
        else
            standardError(dim) = NaN;
        end
    end
    
    track.velFit = velocity;
    track.speedFit = norm(velocity);
    track.velFitSE = standardError;
    track.velFitR2 = rSquared;
    track.posFit0 = position0;
    
    horizontalSpeed = norm(velocity(1:2));
    if track.speedFit > 1e-6
        azimuth = atan2d(velocity(2), velocity(1));
        elevation = atan2d(-velocity(3), horizontalSpeed);
        track.headingFit = [azimuth, elevation];
    else
        track.headingFit = [NaN, NaN];
    end
end


%% 


function q = assessTrackQuality(track)
    % --- Layer 1: Target strength consistency (filter-independent) ---
    ts_band_avg = cellfun(@(x) mean(db(abs(x))), track.peakLog.ts_f);
    q.ts_std = std(ts_band_avg);
    
    % Max ping-to-ping TS jump (detects ID swaps)
    q.ts_max_jump = max(abs(diff(ts_band_avg)));
    
    % --- Layer 2: Echo shape consistency (filter-independent) ---
    pulses = cellfun(@(x) abs(mean(x, 2)), ...
        track.peakLog.ts_data, 'UniformOutput', false);
    rho_echo = zeros(length(pulses)-1, 1);
    for k = 1:length(pulses)-1
        rho_echo(k) = corr(pulses{k}, pulses{k+1});
    end
    q.echo_corr_mean = mean(rho_echo);
    q.echo_corr_min  = min(rho_echo);   % single worst transition
    
    % --- Layer 3: Spatial smoothness (filter-independent) ---
    pos = cell2mat(track.posLog(:, 2));  % raw measurements
    pingIdxs = cell2mat(track.posLog(:, 1));
    t = (pingIdxs - pingIdxs(1));
    residuals = zeros(size(pos));
    for ax = 1:3
        p = polyfit(t, pos(:,ax), 1);
        residuals(:,ax) = pos(:,ax) - polyval(p, t);
    end
    % Per-axis lag-1 autocorrelation of residuals
    for ax = 1:3
        e = residuals(:,ax);
        q.rho1(ax) = sum(e(1:end-1).*e(2:end)) / sum(e.^2);
    end
    q.rho1_min = min(q.rho1);  % strongly negative = suspect
    
    % --- Layer 4: NIS consistency (filter-dependent, use after tuning) ---
    if ~isempty(track.nistAT)
        L = length(track.nistAT);
        q.nis_mean = mean(track.nistAT);
        q.nis_z = (q.nis_mean - 3) / sqrt(6/L);
    end
    
    % --- Summary flags ---
    q.suspect_id_swap   = q.ts_max_jump > 6 || q.echo_corr_min < 0.3;
    q.suspect_false     = q.ts_std > 5;
    q.suspect_merged    = q.rho1_min < -0.5;
end
