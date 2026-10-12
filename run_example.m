%RUN_EXAMPLE  Detect and track single targets in the example WBAT recording.
%
%   Run this script to reproduce the detection and tracking figure from the
%   paper on the five-minute recording shipped in data/. It is also the
%   template for running the pipeline on your own data: point pingFile at your
%   own export and adjust the values in defaultConfig.
%
%   Requires MATLAB R2023b or later and the Signal Processing Toolbox.

clear
% close all
% clc

repoDir = fileparts(mfilename('fullpath'));
addpath(fullfile(repoDir, 'src'), fullfile(repoDir, 'external'))

pingFile     = fullfile(repoDir, 'data', 'example_pings_20211105.mat');
transmitFile = fullfile(repoDir, 'data', 'transmit_signal.mat');

cfg = defaultConfig();

%% Load ------------------------------------------------------------------
pings = loadPings(pingFile, transmitFile);

% Sound speed from the moored CTD, averaged over the recording. With no CTD
% record, replace this with a constant, e.g. soundSpeed = 1490.
soundSpeed = mean(soundSpeedFromCTD(pings.time));

fprintf('%d pings, %s to %s, sound speed %.1f m/s\n', ...
    numel(pings.time), string(pings.time(1)), string(pings.time(end)), soundSpeed);

%% Detect ----------------------------------------------------------------
[detections, env] = detectTargets(pings, soundSpeed, cfg);
fprintf('Detections: %d\n', height(detections));

%% Track -----------------------------------------------------------------
[tracks, stats] = trackTargetsCMKF(detections, env, cfg);
tracks = fitTrackVelocity(tracks, env);

fprintf('Tracking completed in %.2f s\n', stats.elapsedTime);
fprintf('Confirmed tracks: %d (%d still active at the last ping)\n', ...
    stats.numConfirmed, stats.numStillActive);

%% Select and summarise ---------------------------------------------------
longTracks = tracks([tracks.hits] > cfg.analysis.minHits);

speeds = [longTracks.speedFit];
r2     = vertcat(longTracks.velFitR2);

fprintf('Tracks with more than %d hits: %d\n', cfg.analysis.minHits, numel(longTracks));
fprintf('Fitted speed [m/s] - min %.3f | median %.3f | max %.3f\n', ...
    min(speeds), median(speeds), max(speeds));
fprintf('Fitted R^2 median [x y z] - [%.3f %.3f %.3f]\n', ...
    median(r2(:, 1)), median(r2(:, 2)), median(r2(:, 3)));

%% Figure -----------------------------------------------------------------
plotTrackingOverview(detections, longTracks, env, cfg);

% exportgraphics(gcf, 'tracking_overview.pdf', 'ContentType', 'vector')
