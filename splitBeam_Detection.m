%% ========================================================================
%  ACOUSTIC TARGET DETECTION AND TRACKING ALGORITHM 
%  Converted Measurement Kalman Filter
%  ========================================================================
%  Description:
%      This script performs acoustic target detection, tracking, and target 
%      strength measurement using broadband split-beam acoustic data. The 
%      algorithm includes:
%        - Peak-based target detection
%        - Kalman filter tracking with converted measurements (CMKF)
%        - Jacobian evaluated at predicted position
%
%  References:
%      - Lerro & Bar-Shalom (1993): Tracking with consistent 
%        converted measurements versus EKF
%
%  Input Data:
%      - Pulse-compressed split-beam data (.mat file with 'data' structure)
%      - Transmit signal ('WBAT_tx.mat')
%
%  Output:
%      - Detections with peak properties and measurements
%      - Tracked targets with position, velocity, and quality metrics
%      - Visualization figures
%
%  Author:   Zhaozhong Zhuang
%  Date:     Jan. 2026
%  ========================================================================

clear
close all
clc

% ========================================================================
%  SECTION 1: DATA LOADING AND CONFIGURATION
% ========================================================================

% % load('AR59_2-Phase3-D20211004-T153012-0.mat')
% % load('AR59_2-Phase3-D20211004-T173012-0.mat')
% % load('AR59_2-Phase3-D20211009-T050012-0.mat')


% load('AR59_2-Phase3-D20211018-T110012-0.mat') % ***
% load('AR59_2-Phase4-D20211107-T220012-0.mat')
% load('AR59_2-Phase3-D20211006-T223012-0.mat')

load('AR59_2-Phase4-D20211105-T103012-0.mat') % in paper
% load('AR59_2-Phase3-D20211005-T103012-0.mat')

% load('July2022-Phase3-D20220915-T103012-0.mat')
% load('AR59_2-Phase5-D20211207-T203012-0.mat')

% load('July2022-Phase0-D20220715-T160012-0.mat')

% test close-range false target
% load('July2022-Phase4-D20220927-T070012-0.mat')
% load('July2022-Phase6-D20221031-T183012-0.mat')

% calculate sound speed for the 5 min dataset
soundSpeed = mean(getWBAT_SoundSpeed(datetime([data.echodata.timestamp],'ConvertFrom','datenum')));

%  1.2 Load Transmit Signal
% -- Outputs: 'fs_dec' (sampling frequency after decimation) and 'ytx' (transmit signal)
load('WBAT_tx');

%  1.3 Extract System Parameters
NOMINAL_SOUND_SPEED = 1488;  % Nominal sound speed for parameter scaling [m/s]

% -- Pulse and frequency parameters
pulseDuration   = data.param(1).PulseDuration;
frequencyStart  = data.param(1).FrequencyStart;
frequencyEnd    = data.param(1).FrequencyEnd;
frequencyCenter = (frequencyStart + frequencyEnd) / 2;

% -- Data dimensions
numPings          = size(data.echodata, 2);
baseSamplingFreq  = 1 / data.param(1).SampleInterval;

% -- Range-time conversion factors
samplesPerMeter = 2 / soundSpeed * baseSamplingFreq;
rangePerSample  = 1 / samplesPerMeter;

% -- Restrict analysis to depths < 45 m
maxAnalysisRange = 45;
pingLength = round(maxAnalysisRange * samplesPerMeter);

%  1.4 Construct Range and Time Vectors
rangeVector = (0:pingLength-1).' * rangePerSample - pulseDuration * soundSpeed / 2;
timeVector  = rangeVector * 2 / soundSpeed;
pingTimes   = datetime([data.echodata.timestamp], 'ConvertFrom', 'datenum');

%  1.5 Extract and Reshape Echo Data
% -- Pulse-compressed echo data: reshape to [range x channel x ping]
echoData = [data.echodata.compressed];
echoData = echoData(1:pingLength, :);
echoData = reshape(echoData, pingLength, 3, numPings);


% ========================================================================
%  SECTION 2: PREPROCESSING
%  ========================================================================

%  2.1 Upsampling Configuration
INTERPOLATION_FACTOR = 4;
samplingFreq = baseSamplingFreq * INTERPOLATION_FACTOR;

%  2.2 Resample Echo Data and Transmit Signal
[echoData, timeVector] = resample(echoData, timeVector, samplingFreq);

% -- Update derived parameters after interpolation
rangeVector     = timeVector * soundSpeed / 2;
samplesPerMeter = samplesPerMeter * INTERPOLATION_FACTOR;
rangePerSample  = rangePerSample / INTERPOLATION_FACTOR;
pingLength      = length(timeVector);

% -- Interpolate transmit signal
transmitSignal = resample(ytx, INTERPOLATION_FACTOR, 1);


% ========================================================================
%  SECTION 3: TARGET DETECTION CONFIGURATION
%  ========================================================================

%  3.1 Detection Thresholds
% Note: Parameters scale with sound speed ratio for consistency

speedRatio = soundSpeed / NOMINAL_SOUND_SPEED;

minPeakSeparation  = 0 * speedRatio;         % Minimum distance between peaks [m]
minDetectionRange  = 4;         % Ignore detections closer than this [m]
peakWidthLimits    = [0.09, 0.15] * speedRatio;  % [min, max] peak width at half-height [m], [9 15]
peakToBackgroundLimits = [2, 120];           % [min, max] peak-to-background ratio, [2, 120]
powerThreshold     = 10^(-93/10);            % Minimum peak power threshold [linear]

% Angle quality thresholds (currently disabled)
maxAngleStd      = inf;
offAxisAngleMax  = 15;


%  3.2 Background Power Reference Window
% Define reference cells for computing background power (excludes peak region)
refCellLimits  = [-0.26, -0.13, 0.13, 0.26] * speedRatio;  % [m]
refCellIndices = [round(refCellLimits(1) * samplesPerMeter) : round(refCellLimits(2) * samplesPerMeter), ...
                  round(refCellLimits(3) * samplesPerMeter) : round(refCellLimits(4) * samplesPerMeter)];


%  3.3 Target Strength (TS) Measurement Window
tsWindowLimits = [-0.8, 0.8] * speedRatio;  % Window extent [m]
tsWindowIndices = round(tsWindowLimits(1) * samplesPerMeter) : ...
                  round(tsWindowLimits(2) * samplesPerMeter);
tsWindowLength = length(tsWindowIndices);
tsWindowTaperRatio = 0.3;
tsWindow = tukeywin(tsWindowLength, tsWindowTaperRatio);

%  3.4 Frequency Axis Construction (handles bandpass sampling aliasing)
nfft = 2^(1 + nextpow2(max([length(transmitSignal), tsWindowLength])));

% -- Compute aliasing index
aliasingIndex = floor(frequencyEnd / baseSamplingFreq) - floor(frequencyStart / baseSamplingFreq);

% -- Construct frequency axis accounting for aliasing
freqAxis = ((0:nfft-1) / nfft - floor(2 * (0:nfft-1) / nfft) * aliasingIndex) * samplingFreq ...
           + floor(frequencyEnd / baseSamplingFreq) * baseSamplingFreq;
[sortedFreq, freqSortMap] = sort(freqAxis.', 'ascend'); %#ok<UDIM>

% -- Find indices corresponding to signal bandwidth
[~, freqLimitIndices] = min(abs(sortedFreq - [frequencyStart, frequencyEnd]));
frequencyBandIndices = freqLimitIndices(1) : freqLimitIndices(2);

%  3.5 Pre-compute Transmit Signal FFT
transmitFFT = fft(transmitSignal, nfft);
transmitFFT = transmitFFT(freqSortMap);

%  3.6 Transducer Geometry Parameters
TRANSDUCER_EFFECTIVE_RADIUS = 130e-3 / 2;  % Effective radius [m]
sectorCentroidRadius = 2 * TRANSDUCER_EFFECTIVE_RADIUS * sin(pi/3) / pi;
sectorSeparation = sectorCentroidRadius * sqrt(3);  % Distance between sector centroids [m]


% ========================================================================
%  SECTION 4: TARGET DETECTION PROCESSING
%  ========================================================================

detections = struct;
detectionCount = 0;

for pingIdx = 1:numPings
    
    %  4.1 Extract Current Ping Data
    channelData = echoData(:, :, pingIdx);
    averagePower = abs(mean(channelData, 2)).^2;
    
    %  4.2 Peak Finding
    [peakPowers, peakLocations, peakWidths] = findpeaks(averagePower, ...
        'MinPeakDistance', minPeakSeparation * samplesPerMeter, ...
        'MinPeakHeight',   powerThreshold, ...
        'MinPeakWidth',    peakWidthLimits(1) * samplesPerMeter, ...
        'MaxPeakWidth',    peakWidthLimits(2) * samplesPerMeter, ...
        'WidthReference',  'halfheight');
    
    %  4.3 Process Each Detected Peak
    for peakIdx = 1:length(peakLocations)
        
        peakPower    = peakPowers(peakIdx);
        peakLocation = peakLocations(peakIdx);
        peakWidth    = peakWidths(peakIdx);
        
        %  4.3.1 Boundary Check - Skip peaks too close to edges
        minWindowOffset = min([refCellIndices(1), tsWindowIndices(1)]);
        maxWindowOffset = max([refCellIndices(end), tsWindowIndices(end)]);
        
        if (minWindowOffset + peakLocation < 1) || ...
           (maxWindowOffset + peakLocation > pingLength) || ...
           (rangeVector(peakLocation) < minDetectionRange)
            continue
        end
        
        %  4.3.2 Compute Peak-to-Background Ratio
        % -- Integrate background power over reference cells
        backgroundPower = mean(averagePower(refCellIndices + peakLocation));
        peakToBackground = peakPower / backgroundPower;
        
        % -- Apply P2B filter
        if (peakToBackground < peakToBackgroundLimits(1)) || ...
           (peakToBackground > peakToBackgroundLimits(2))
            continue
        end
        
        %%%%  4.3.3.i Time-Domain AoA Estimation (Split-Aperture Correlator)
        % -- Used for tracking - single sample at peak
        peakSamples = channelData(peakLocation, :);
        
        % -- Sensitivity factor at center frequency
        sensitivityCenterFreq = soundSpeed / (2 * pi * frequencyCenter) / sectorSeparation;
        
        % -- Sector pairs phase differences
        dphi_LF = angle(peakSamples(1)/peakSamples(3));  % Left-Front
        dphi_RF = angle(peakSamples(2)/peakSamples(3));  % Right-Front
        
        % -- Convert phase to angles
        alongAngleCSAC   = asind((dphi_LF + dphi_RF) * sensitivityCenterFreq / sqrt(3));
        athwartAngleCSAC = asind((dphi_LF - dphi_RF) * sensitivityCenterFreq);
        offAxisAngleCSAC = asind(sqrt(sind(alongAngleCSAC)^2 + sind(athwartAngleCSAC)^2));

        % -- Optional: Reject large-angle detections
        if offAxisAngleCSAC > offAxisAngleMax
            continue
        end

        %%%%  4.3.3.ii Time-Domain AoA STD (Split-Aperture Correlator)
        % -- Used for tracking - single sample at peak
        peakSamples = channelData(peakLocation+(-5:5), :);
        
        % -- Sensitivity factor at center frequency
        sensitivityCenterFreq = soundSpeed / (2 * pi * frequencyCenter) / sectorSeparation;
        
        % -- Sector pairs phase differences
        dphi_LF = angle(peakSamples(:,1)./peakSamples(:,3));  % Left-Front
        dphi_RF = angle(peakSamples(:,2)./peakSamples(:,3));  % Right-Front
        
        % -- Convert phase to angles
        alongAngleCSAC_std   = std(asind((dphi_LF + dphi_RF) * sensitivityCenterFreq / sqrt(3)));
        athwartAngleCSAC_std = std(asind((dphi_LF - dphi_RF) * sensitivityCenterFreq));
        offAxisAngleCSAC_std = std(asind(sqrt(sind(alongAngleCSAC).^2 + sind(athwartAngleCSAC).^2)));

        % -- Optional: Reject large-angle STD
        if alongAngleCSAC_std + athwartAngleCSAC_std > maxAngleStd
            continue
        end


        
        %  4.3.4 Frequency-Domain AoA Estimation
        tsData = channelData(tsWindowIndices + peakLocation, :);
        tsFFT = fft(tsData .* tsWindow, nfft);
        tsFFT = tsFFT(freqSortMap, :);
        
        %  4.3.5 Target Strength Estimation
        targetStrength = mean(tsFFT, 2) ./ abs(transmitFFT).^2 * INTERPOLATION_FACTOR;
        targetStrength = targetStrength .* rangeVector(peakLocation)^2;
        
        %  4.3.6 Beam Pattern Compensation
        beamPattern = 2 * besselj(1, 2*pi*sortedFreq/soundSpeed * TRANSDUCER_EFFECTIVE_RADIUS * sind(offAxisAngleCSAC)) ./ ...
                                    (2*pi*sortedFreq/soundSpeed * TRANSDUCER_EFFECTIVE_RADIUS * sind(offAxisAngleCSAC) + eps);
        targetStrength = targetStrength ./ beamPattern.^2;
        
        %  4.3.7 Frequency-Domain Angle Estimation
        sensitivityFreq = soundSpeed ./ (2 * pi * sortedFreq) / sectorSeparation;
        
        crossCorrLF_freq = tsFFT(:, 1) .* conj(tsFFT(:, 3));
        crossCorrRF_freq = tsFFT(:, 2) .* conj(tsFFT(:, 3));
        
        alongAngleFreq   = asind((angle(crossCorrLF_freq) + angle(crossCorrRF_freq)) .* sensitivityFreq / sqrt(3));
        athwartAngleFreq = asind((angle(crossCorrLF_freq) - angle(crossCorrRF_freq)) .* sensitivityFreq);
        
        %  4.3.8 Compute AoA Statistics
        alongAngleStd   = std(alongAngleFreq(frequencyBandIndices));
        athwartAngleStd = std(athwartAngleFreq(frequencyBandIndices));

        alongAngleAvg   = mean(alongAngleFreq(frequencyBandIndices));
        athwartAngleAvg = mean(athwartAngleFreq(frequencyBandIndices));
        offAxisAngleAvg = asind(sqrt(sind(alongAngleAvg)^2 + sind(athwartAngleAvg)^2));
        

        %  4.3.9 Store Detection
        detectionCount = detectionCount + 1;
        
        % -- Metadata
        detections(detectionCount).pingIdx           = pingIdx;
        detections(detectionCount).pingTime          = pingTimes(pingIdx);
        detections(detectionCount).samplingFreq0     = baseSamplingFreq;
        detections(detectionCount).interpFactor      = INTERPOLATION_FACTOR;
        detections(detectionCount).SoundSpeed        = soundSpeed;
        detections(detectionCount).SoundSpeedNominal = NOMINAL_SOUND_SPEED;
        
        % -- Peak properties
        detections(detectionCount).Pwr   = peakPower;
        detections(detectionCount).Loc   = peakLocation;
        detections(detectionCount).Range = rangeVector(peakLocation);
        detections(detectionCount).Width = peakWidth * rangePerSample;
        
        % -- Background and P2B ratio
        detections(detectionCount).Bkg  = backgroundPower;
        detections(detectionCount).P2B  = peakToBackground;
        
        % -- Angle estimates
        detections(detectionCount).AoA_csac  = [alongAngleCSAC, athwartAngleCSAC, offAxisAngleCSAC];
        detections(detectionCount).AoA_csac_std  = [alongAngleCSAC_std, athwartAngleCSAC_std, offAxisAngleCSAC_std];
        detections(detectionCount).AoA_f     = [alongAngleAvg, athwartAngleAvg, offAxisAngleAvg];
        detections(detectionCount).AoA_f_std = [alongAngleStd, athwartAngleStd];
        
        % -- Target strength data
        detections(detectionCount).ts_data = {tsData};
        detections(detectionCount).ts_ind  = {tsWindowIndices + peakLocation};
        detections(detectionCount).ts_f    = {targetStrength(frequencyBandIndices)};
        detections(detectionCount).f       = {sortedFreq(frequencyBandIndices)};
    end
end

% Convert to table for easier manipulation
detections = struct2table(detections);
size(detections,1)



splitBeam_Tracking_CMKF
% exportgraphics(gcf,[save_file_name '_CMKF.pdf'],'ContentType','vector')


% cfg.gatingSpace = 'cartesian';  % 'observation' or 'cartesian'
% splitBeam_Tracking_gh
% exportgraphics(gcf,[save_file_name '_conv_CART.pdf'],'ContentType','vector')


% cfg.gatingSpace = 'observation';  % 'observation' or 'cartesian'
% splitBeam_Tracking_gh
% exportgraphics(gcf,[save_file_name '_conv_SB.pdf'],'ContentType','vector')