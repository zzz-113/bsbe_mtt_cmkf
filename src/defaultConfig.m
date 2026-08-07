function cfg = defaultConfig()
%DEFAULTCONFIG  Parameter set used for the results in the paper.
%
%   cfg = defaultConfig() returns a struct with three sections:
%
%     cfg.detection - peak detection, angle estimation and target strength
%     cfg.tracking  - CMKF motion/measurement model, gating, track lifecycle
%     cfg.analysis  - post-processing and figure settings
%
%   Copy the struct and overwrite individual fields to run a different
%   configuration; every field is read at most once by the pipeline, so no
%   value is silently duplicated elsewhere.
%
%   Several detection parameters are expressed at a nominal sound speed and
%   are rescaled by the actual sound speed at run time (see detectTargets),
%   so that windows defined in metres keep the same extent in samples.

%% Detection ------------------------------------------------------------
d = struct();

d.nominalSoundSpeed = 1488;     % sound speed the length parameters below refer to [m/s]
d.maxAnalysisRange  = 45;       % process samples out to this range [m]
d.minDetectionRange = 4;        % discard peaks closer than this [m]
d.interpFactor      = 4;        % integer upsampling factor applied before detection

% Peak finding on the channel-averaged power envelope.
d.powerThreshold    = 10^(-93/10);  % minimum peak power [linear]
d.peakWidthLimits   = [0.09 0.15];  % [min max] width at half height [m]
d.minPeakSeparation = 0;            % minimum spacing between peaks [m]; 0 disables

% Background reference cells, as offsets either side of the peak. Used to
% form the peak-to-background ratio; the peak itself is excluded.
d.refCellLimits         = [-0.26 -0.13 0.13 0.26];  % [m]
d.peakToBackgroundLimits = [2 120];                 % [min max] accepted ratio

% Angle-of-arrival rejection. maxAngleStd is disabled by default: the paper
% uses only the off-axis limit.
d.offAxisAngleMax = 15;    % reject peaks beyond this off-axis angle [deg]
d.maxAngleStd     = Inf;   % reject if along+athwart angle spread exceeds this [deg]

% Window used for the target-strength spectrum, as offsets around the peak.
d.tsWindowLimits    = [-0.8 0.8];  % [m]
d.tsWindowTaperRatio = 0.3;        % Tukey window taper fraction

% Transducer geometry (Simrad ES38-18DK split-beam, three sectors).
d.transducerRadius = 130e-3 / 2;   % effective aperture radius [m]

cfg.detection = d;

%% Tracking -------------------------------------------------------------
t = struct();

% Measurement noise in split-beam coordinates z = [r, alpha, beta], where
% alpha = sin(theta_along) and beta = sin(theta_athwart). The direction-cosine
% variance grows with range, so it is parameterised as (kappaTheta * r)^2.
t.sigmaR     = 10e-2;    % range measurement std [m]
t.kappaTheta = 1.5e-3;   % direction-cosine std per metre of range [1/m]

% Process noise: piecewise-constant white acceleration, one std per axis.
t.sigmaAccel = [2.5 2.5 3.7] * 1e-3;   % [x y z] acceleration std [m/s^2]

% Covariance of a newly initialised track. Velocity is initialised to zero,
% so P0vel encodes the expected target speed rather than an estimate error.
t.P0pos = [0.2 0.2 0.2].^2;   % initial position variance [m^2]
t.P0vel = [0.3 0.3 0.1].^2;   % initial velocity variance [(m/s)^2]

% Gating and assignment. chi2Gate is the chi-squared quantile for 3 degrees
% of freedom:
%     P_G   0.90    0.95    0.975   0.99    0.999
%           6.251   7.815   9.348   11.345  16.266
t.chi2Gate            = 7.815;
t.costOfNonAssignment = 1.5;   % cost matchpairs charges for leaving a pair unmatched

% Track lifecycle. A track is confirmed once it survives to age Nconf with at
% least Mconf hits, and is deleted after maxMisses consecutive missed updates.
t.Mconf     = 4;
t.Nconf     = 6;
t.maxMisses = 5;

cfg.tracking = t;

%% Analysis and figures -------------------------------------------------
a = struct();

a.minHits        = 10;    % keep confirmed tracks with more hits than this
a.instrumentDepth = 590;  % depth of the upward-looking transducer [m], for depth axes
a.horizontalSpan = 10;    % half-width of the x and y axes in the 3-D panels [m]
a.velocityScale  = 1e2;   % m/s -> cm/s for the velocity quiver panel
a.spDisplayGain  = 5.58;  % amplitude scale applied to the Sp echogram only

cfg.analysis = a;
end
