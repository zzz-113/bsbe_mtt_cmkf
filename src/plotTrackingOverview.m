function fig = plotTrackingOverview(detections, tracks, env, cfg)
%PLOTTRACKINGOVERVIEW  Five-panel detection and tracking summary figure.
%
%   FIG = plotTrackingOverview(DETECTIONS, TRACKS, ENV, CFG) draws, left to
%   right:
%
%     (a) S_p echogram, sector-averaged and range-compensated
%     (b) detections, marker area proportional to peak-to-background ratio
%     (c) filtered 3-D track positions
%     (d) fitted track velocity, drawn as a vector at the track's median depth
%     (e) tracks overlaid on a faded echogram, for visual validation
%
%   TRACKS must have been through fitTrackVelocity. Panels (a), (b) and (e)
%   share a time axis; panels (c) and (d) share a camera.
%
%   Vertical axes are labelled as depth below the surface, using
%   cfg.analysis.instrumentDepth for the upward-looking transducer.

arguments
    detections table
    tracks struct
    env struct
    cfg struct
end

a = cfg.analysis;

% Range in metres above the transducer -> depth below the surface.
rangeTicks = 0:10:ceil(env.rangeVector(end) / 10) * 10;
depthLabels = string(a.instrumentDepth - rangeTicks);
rangeLimits = [cfg.detection.minDetectionRange, env.rangeVector(end)];

% Range-compensated echo amplitude, shared by panels (a) and (e).
spEchogram = (env.rangeVector.^2) .* env.echoMeanSector;

trackColors = orderedcolors('gem');
numColors   = size(trackColors, 1);
colorOf = @(n) trackColors(mod(n, numColors) + 1, :);

fig = figure('Units', 'inches', 'Position', [2, 2, 6.5, 3] * 2);
tiledlayout(1, 5, 'TileSpacing', 'compact', 'Padding', 'compact')

%% (a) Echogram ---------------------------------------------------------
ax1 = nexttile;
imagesc(env.pingTime, env.rangeVector, db(spEchogram / a.spDisplayGain^2));
axis xy
clim([-95 -45])
colormap(ax1, viridis)
xsecondarylabel('');
ylabel('Depth (m)')
title('(a) \itS_p')

%% (b) Detections -------------------------------------------------------
ax2 = nexttile;
scatter(detections.pingTime, detections.range, detections.peakToBackground / 30, ...
        'x', 'MarkerEdgeColor', '#008080');
box on
xsecondarylabel('');
title('(b) Detections')

%% (c) Filtered 3-D tracks ----------------------------------------------
ax3 = nexttile;
hold on
for n = 1:numel(tracks)
    pos = cell2mat(tracks(n).posLog(:, 4));   % filtered position
    plot3(pos(:, 1), pos(:, 2), pos(:, 3), 'Color', colorOf(n));
end
daspect([1, 1, 0.7]);
axis([-a.horizontalSpan, a.horizontalSpan, -a.horizontalSpan, a.horizontalSpan, rangeLimits])
axis vis3d
grid on; box on
xlabel('x (m)'); ylabel('y (m)');
zticks(rangeTicks); zticklabels(depthLabels)
title('(c) 3D Tracks')

%% (d) Fitted velocity --------------------------------------------------
ax4 = nexttile;
hold on
for n = 1:numel(tracks)
    pos = double(cell2mat(tracks(n).posLog(:, 4)));
    vel = double(tracks(n).velFit) * a.velocityScale;

    % Each track contributes one vector, drawn from the axis origin at the
    % track's median depth, so the panel reads as a velocity profile.
    quiver3(0, 0, median(pos(:, 3)), vel(1), vel(2), vel(3), ...
            'Color', colorOf(n), 'AutoScale', 'off', 'ShowArrowHead', 'off');
end
daspect([1, 1, 0.7]);
axis([-a.horizontalSpan, a.horizontalSpan, -a.horizontalSpan, a.horizontalSpan, rangeLimits])
axis vis3d
grid on; box on
xlabel('u (cm/s)'); ylabel('v (cm/s)');
zticks(rangeTicks); zticklabels(depthLabels)
title('(d) Velocity')

set([ax3 ax4], 'Projection', 'Perspective', 'CameraViewAngle', 15, ...
    'CameraPosition', [-100 50 30])

%% (e) Tracks over the echogram -----------------------------------------
ax5 = nexttile;
imagesc(env.pingTime, env.rangeVector, db(spEchogram), 'AlphaData', 0.5);
axis xy; hold on
clim([-44 26])
colormap(ax5, flipud(gray))

for n = 1:numel(tracks)
    plot(tracks(n).peakLog.pingTime, tracks(n).peakLog.range, ...
         'LineWidth', 0.88, 'Color', colorOf(n))
end
title('(e) Validation')

%% Shared axis formatting -----------------------------------------------
timeAxes = [ax1, ax2, ax5];
xlim(timeAxes, env.pingTime([1, end]))
ylim(timeAxes, rangeLimits)
yticks(timeAxes, rangeTicks)
yticklabels(ax1, depthLabels)
yticklabels([ax2 ax5], '')
zticklabels([ax3 ax4], '')
xsecondarylabel('');

set(findall(fig, '-property', 'FontSize'), 'FontName', 'Arial', 'FontWeight', 'Normal');
set(findall(fig, 'type', 'axes'), 'TitleHorizontalAlignment', 'left', 'FontSize', 18);
end
