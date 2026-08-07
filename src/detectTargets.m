function [detections, env] = detectTargets(pings, soundSpeed, cfg)
%DETECTTARGETS  Detect single targets in broadband split-beam echoes.
%
%   [DETECTIONS, ENV] = detectTargets(PINGS, SOUNDSPEED, CFG) runs, for every
%   ping in PINGS:
%
%     1. upsampling of the pulse-compressed echoes by cfg.detection.interpFactor;
%     2. peak picking on the sector-averaged power envelope;
%     3. rejection of peaks by range, width, peak-to-background ratio and
%        off-axis angle;
%     4. angle of arrival from the split-aperture correlator, both at the peak
%        sample and across the signal band;
%     5. an uncalibrated target-strength spectrum, compensated for spherical
%        spreading and for the circular-piston beam pattern.
%
%   Inputs
%     pings      - struct from loadPings
%     soundSpeed - sound speed for the recording [m/s], scalar
%     cfg        - configuration struct from defaultConfig
%
%   Outputs
%     detections - table, one row per accepted peak:
%        pingIndex        index of the ping the peak was found in
%        pingTime         datetime of that ping
%        peakPower        peak of the sector-averaged power envelope [linear]
%        peakSample       sample index of the peak
%        range            range to the peak [m]
%        peakWidth        peak width at half height [m]
%        backgroundPower  mean power over the reference cells [linear]
%        peakToBackground peakPower / backgroundPower
%        aoa              [along athwart offAxis] at the peak sample [deg]
%        aoaStd           std of the above over +/-5 samples [deg]
%        aoaSpectral      [along athwart offAxis] averaged over the band [deg]
%        aoaSpectralStd   [along athwart] std across the band [deg]
%        ts               {complex TS spectrum over the band}
%        frequency        {frequency axis of ts [Hz]}
%        echoWindow       {[nWindow x 3] complex echo used for ts}
%
%     env - acquisition geometry shared with the tracker and the figures:
%        soundSpeed, samplingFreq, baseSamplingFreq, interpFactor,
%        samplesPerMetre, rangePerSample, rangeVector, pingTime, numPings,
%        pingInterval, echoMeanSector
%
%   Angles are signed, in degrees, in the transducer frame: 'along' is the
%   alongship axis and 'athwart' the athwartship axis.
%
%   Requires the Signal Processing Toolbox (resample, findpeaks, tukeywin).

arguments
    pings struct
    soundSpeed (1,1) double
    cfg struct
end

d = cfg.detection;

% Parameters below are quoted at a nominal sound speed. Rescaling them keeps
% window lengths fixed in metres as the in situ sound speed changes.
speedRatio = soundSpeed / d.nominalSoundSpeed;

%% Range axis and echo data --------------------------------------------
baseSamplingFreq = 1 / pings.sampleInterval;
samplesPerMetre  = 2 / soundSpeed * baseSamplingFreq;
rangePerSample   = 1 / samplesPerMetre;

numPings   = numel(pings.time);
pingLength = round(d.maxAnalysisRange * samplesPerMetre);
if pingLength > size(pings.compressed, 1)
    error('detectTargets:rangeTooLarge', ...
        ['cfg.detection.maxAnalysisRange (%g m) exceeds the %g m stored in ' ...
         'the ping data.'], d.maxAnalysisRange, size(pings.compressed, 1) * rangePerSample);
end

% The pulse-compressed peak sits at the centre of the transmit pulse, so the
% range axis is shifted back by half a pulse length.
rangeVector = (0:pingLength-1).' * rangePerSample - pings.pulseDuration * soundSpeed / 2;
timeVector  = rangeVector * 2 / soundSpeed;

echoData = pings.compressed(1:pingLength, :, :);

%% Upsampling -----------------------------------------------------------
samplingFreq = baseSamplingFreq * d.interpFactor;

[echoData, timeVector] = resample(echoData, timeVector, samplingFreq);
transmitSignal = resample(pings.transmitSignal, d.interpFactor, 1);

rangeVector     = timeVector * soundSpeed / 2;
samplesPerMetre = samplesPerMetre * d.interpFactor;
rangePerSample  = rangePerSample / d.interpFactor;
pingLength      = numel(timeVector);

%% Analysis windows -----------------------------------------------------
% Background reference cells, as sample offsets either side of a peak.
refCellLimits  = d.refCellLimits * speedRatio;
refCellIndices = [round(refCellLimits(1) * samplesPerMetre) : round(refCellLimits(2) * samplesPerMetre), ...
                  round(refCellLimits(3) * samplesPerMetre) : round(refCellLimits(4) * samplesPerMetre)];

% Tapered window over which the target-strength spectrum is taken.
tsWindowLimits  = d.tsWindowLimits * speedRatio;
tsWindowIndices = round(tsWindowLimits(1) * samplesPerMetre) : round(tsWindowLimits(2) * samplesPerMetre);
tsWindow        = tukeywin(numel(tsWindowIndices), d.tsWindowTaperRatio);

% Samples either side of the peak used for the angle spread.
STD_HALF_WINDOW = 5;
stdIndices = -STD_HALF_WINDOW:STD_HALF_WINDOW;

% A peak needs all three windows to fit inside the ping.
minWindowOffset = min([refCellIndices, tsWindowIndices, stdIndices]);
maxWindowOffset = max([refCellIndices, tsWindowIndices, stdIndices]);

%% Frequency axis -------------------------------------------------------
% The receiver samples the 36-44 kHz band below Nyquist, so the spectrum is
% aliased. aliasingIndex counts how many Nyquist zones the band spans, and the
% axis is unwrapped back to true frequency before sorting.
nfft = 2^(1 + nextpow2(max(numel(transmitSignal), numel(tsWindowIndices))));

aliasingIndex = floor(pings.frequencyEnd / baseSamplingFreq) - ...
                floor(pings.frequencyStart / baseSamplingFreq);

freqAxis = ((0:nfft-1) / nfft - floor(2 * (0:nfft-1) / nfft) * aliasingIndex) * samplingFreq ...
           + floor(pings.frequencyEnd / baseSamplingFreq) * baseSamplingFreq;
[sortedFreq, freqSortMap] = sort(freqAxis.', 'ascend'); %#ok<UDIM>

[~, freqLimitIndices] = min(abs(sortedFreq - [pings.frequencyStart, pings.frequencyEnd]));
bandIndices = freqLimitIndices(1):freqLimitIndices(2);
bandFreq    = sortedFreq(bandIndices);

transmitFFT = fft(transmitSignal, nfft);
transmitFFT = transmitFFT(freqSortMap);

%% Split-aperture geometry ---------------------------------------------
% The three sectors are 120-degree pie slices. Their centroids lie at radius
% 2*R*sin(pi/3)/pi, and adjacent centroids are sqrt(3) times that apart.
sectorCentroidRadius = 2 * d.transducerRadius * sin(pi/3) / pi;
sectorSeparation     = sectorCentroidRadius * sqrt(3);

% Phase difference to angle, at the centre of the band.
frequencyCentre = (pings.frequencyStart + pings.frequencyEnd) / 2;
angleSensitivity = soundSpeed / (2 * pi * frequencyCentre) / sectorSeparation;

% Same conversion, resolved per frequency bin.
angleSensitivitySpectral = soundSpeed ./ (2 * pi * sortedFreq) / sectorSeparation;

%% Detection loop -------------------------------------------------------
detectionCells = cell(1, numPings);

for pingIdx = 1:numPings

    sectorData   = echoData(:, :, pingIdx);
    averagePower = abs(mean(sectorData, 2)).^2;

    [peakPowers, peakSamples, peakWidths] = findpeaks(averagePower, ...
        'MinPeakDistance', d.minPeakSeparation * speedRatio * samplesPerMetre, ...
        'MinPeakHeight',   d.powerThreshold, ...
        'MinPeakWidth',    d.peakWidthLimits(1) * speedRatio * samplesPerMetre, ...
        'MaxPeakWidth',    d.peakWidthLimits(2) * speedRatio * samplesPerMetre, ...
        'WidthReference',  'halfheight');

    pingDetections = cell(1, numel(peakSamples));
    nKept = 0;

    for peakIdx = 1:numel(peakSamples)

        peakSample = peakSamples(peakIdx);

        % -- Reject peaks whose analysis windows fall off the ping, and peaks
        %    inside the near-field blanking range.
        if peakSample + minWindowOffset < 1 || ...
           peakSample + maxWindowOffset > pingLength || ...
           rangeVector(peakSample) < d.minDetectionRange
            continue
        end

        % -- Peak-to-background ratio over the reference cells.
        backgroundPower  = mean(averagePower(refCellIndices + peakSample));
        peakToBackground = peakPowers(peakIdx) / backgroundPower;

        if peakToBackground < d.peakToBackgroundLimits(1) || ...
           peakToBackground > d.peakToBackgroundLimits(2)
            continue
        end

        % -- Angle of arrival at the peak sample. The split-aperture
        %    correlator takes the phase of each sector pair; the sum resolves
        %    the alongship angle and the difference the athwartship angle.
        peakSectors = sectorData(peakSample, :);
        dPhiLF = angle(peakSectors(1) / peakSectors(3));
        dPhiRF = angle(peakSectors(2) / peakSectors(3));

        alongAngle   = asind((dPhiLF + dPhiRF) * angleSensitivity / sqrt(3));
        athwartAngle = asind((dPhiLF - dPhiRF) * angleSensitivity);
        offAxisAngle = asind(sqrt(sind(alongAngle)^2 + sind(athwartAngle)^2));

        if offAxisAngle > d.offAxisAngleMax
            continue
        end

        % -- Spread of the same estimator over the samples spanning the peak.
        %    A single target gives a stable angle across the echo; a merged or
        %    noise-dominated peak does not.
        nearSectors = sectorData(peakSample + stdIndices, :);
        dPhiLFnear = angle(nearSectors(:, 1) ./ nearSectors(:, 3));
        dPhiRFnear = angle(nearSectors(:, 2) ./ nearSectors(:, 3));

        alongNear   = asind((dPhiLFnear + dPhiRFnear) * angleSensitivity / sqrt(3));
        athwartNear = asind((dPhiLFnear - dPhiRFnear) * angleSensitivity);
        offAxisNear = asind(sqrt(sind(alongNear).^2 + sind(athwartNear).^2));

        angleStd = [std(alongNear), std(athwartNear), std(offAxisNear)];

        if angleStd(1) + angleStd(2) > d.maxAngleStd
            continue
        end

        % -- Spectrum of the tapered echo window, one column per sector.
        echoWindow = sectorData(tsWindowIndices + peakSample, :);
        echoFFT    = fft(echoWindow .* tsWindow, nfft);
        echoFFT    = echoFFT(freqSortMap, :);

        % -- Target strength: sector-averaged echo spectrum divided by the
        %    transmit power spectrum, corrected for two-way spherical
        %    spreading. The interpFactor term undoes the gain the upsampling
        %    puts into the FFT. No system calibration is applied, so this is a
        %    relative measure.
        ts = mean(echoFFT, 2) ./ abs(transmitFFT).^2 * d.interpFactor;
        ts = ts .* rangeVector(peakSample)^2;

        % -- Circular-piston beam pattern at the estimated off-axis angle.
        %    eps keeps the on-axis limit 2*J1(u)/u -> 1 finite.
        u = 2 * pi * sortedFreq / soundSpeed * d.transducerRadius * sind(offAxisAngle);
        beamPattern = 2 * besselj(1, u) ./ (u + eps);
        ts = ts ./ beamPattern.^2;

        % -- Angle of arrival resolved per frequency bin, from the same
        %    cross-spectra. Its spread across the band is a second single-target
        %    indicator, independent of the time-domain one above.
        crossLF = echoFFT(:, 1) .* conj(echoFFT(:, 3));
        crossRF = echoFFT(:, 2) .* conj(echoFFT(:, 3));

        alongSpectral   = asind((angle(crossLF) + angle(crossRF)) .* angleSensitivitySpectral / sqrt(3));
        athwartSpectral = asind((angle(crossLF) - angle(crossRF)) .* angleSensitivitySpectral);

        alongSpectralMean   = mean(alongSpectral(bandIndices));
        athwartSpectralMean = mean(athwartSpectral(bandIndices));
        offAxisSpectralMean = asind(sqrt(sind(alongSpectralMean)^2 + sind(athwartSpectralMean)^2));

        % -- Accept.
        nKept = nKept + 1;
        pingDetections{nKept} = struct( ...
            'pingIndex',        pingIdx, ...
            'pingTime',         pings.time(pingIdx), ...
            'peakPower',        peakPowers(peakIdx), ...
            'peakSample',       peakSample, ...
            'range',            rangeVector(peakSample), ...
            'peakWidth',        peakWidths(peakIdx) * rangePerSample, ...
            'backgroundPower',  backgroundPower, ...
            'peakToBackground', peakToBackground, ...
            'aoa',              [alongAngle, athwartAngle, offAxisAngle], ...
            'aoaStd',           angleStd, ...
            'aoaSpectral',      [alongSpectralMean, athwartSpectralMean, offAxisSpectralMean], ...
            'aoaSpectralStd',   [std(alongSpectral(bandIndices)), std(athwartSpectral(bandIndices))], ...
            'ts',               {{ts(bandIndices)}}, ...
            'frequency',        {{bandFreq}}, ...
            'echoWindow',       {{echoWindow}});
    end

    detectionCells{pingIdx} = [pingDetections{1:nKept}];
end

detections = struct2table([detectionCells{:}]);

%% Geometry shared with downstream stages -------------------------------
env = struct( ...
    'soundSpeed',       soundSpeed, ...
    'baseSamplingFreq', baseSamplingFreq, ...
    'samplingFreq',     samplingFreq, ...
    'interpFactor',     d.interpFactor, ...
    'samplesPerMetre',  samplesPerMetre, ...
    'rangePerSample',   rangePerSample, ...
    'rangeVector',      rangeVector, ...
    'pingTime',         pings.time, ...
    'numPings',         numPings, ...
    'pingInterval',     median(seconds(diff(pings.time))), ...
    'echoMeanSector',   squeeze(mean(echoData, 2)));
end
