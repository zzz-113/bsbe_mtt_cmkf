function makeExampleDataset(rawDir, outDir)
%MAKEEXAMPLEDATASET  Build the slimmed datasets shipped in data/.
%
%   makeExampleDataset(RAWDIR, OUTDIR) reads the original WBAT exports from
%   RAWDIR and writes reduced copies to OUTDIR. It is kept in the repository
%   so the provenance of the shipped example data is reproducible; running it
%   is not required to reproduce the paper figure.
%
%   Two reductions are applied, both purely to keep the repository small:
%
%     1. Ping data. The raw export stores both the raw complex samples
%        ('complexsamples') and the pulse-compressed samples ('compressed').
%        Only the latter is used, so 'complexsamples' is dropped. Samples are
%        then truncated to MAXRANGE metres; the analysis in run_example.m
%        reaches 45 m.
%
%     2. CTD data. Only the four fields consumed by soundSpeedFromCTD are
%        kept (mtime, temp, salt, press); the redundant broken-out date fields
%        are dropped.
%
%   Inputs
%     rawDir - folder holding the original exports. Default: <repo>/raw, which
%              is git-ignored: the unreduced exports are not distributed.
%     outDir - destination folder. Default: <repo>/data.
%
%   Original -> shipped file names:
%     AR59_2-Phase4-D20211105-T103012-0.mat -> example_pings_20211105.mat
%     WBAT_tx.mat                           -> transmit_signal.mat
%     WBAT_CTD_2021.mat                     -> ctd_2021.mat
%     WBAT_CTD_2022.mat                     -> ctd_2022.mat

arguments
    rawDir (1,:) char = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'raw')
    outDir (1,:) char = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data')
end

MAX_RANGE = 64;          % keep samples out to this range [m]
NOMINAL_SOUND_SPEED = 1500;   % upper bound, so the cut is conservative [m/s]

if ~isfolder(outDir)
    mkdir(outDir)
end

%% Ping data ------------------------------------------------------------
rawPingFile = 'AR59_2-Phase4-D20211105-T103012-0.mat';
raw = load(fullfile(rawDir, rawPingFile));
data = raw.data;

sampleInterval  = data.param(1).SampleInterval;
samplesPerMetre = 2 / NOMINAL_SOUND_SPEED / sampleInterval;
nKeep = min(round(MAX_RANGE * samplesPerMetre), data.echodata(1).count);

for k = 1:numel(data.echodata)
    data.echodata(k).compressed = data.echodata(k).compressed(1:nKeep, :);
    data.echodata(k).range      = data.echodata(k).range(1:nKeep);
    data.echodata(k).count      = nKeep;
    data.echodata(k).maxRange   = double(data.echodata(k).range(end));
end
data.echodata = rmfield(data.echodata, 'complexsamples');

data.provenance = struct( ...
    'originalFile',   rawPingFile, ...
    'trimmedToRange', data.echodata(1).maxRange, ...
    'droppedFields',  {{'echodata.complexsamples'}}, ...
    'note', ['Reduced copy for repository distribution: raw (uncompressed) ' ...
             'samples removed and pulse-compressed samples truncated in range.']);

save(fullfile(outDir, 'example_pings_20211105.mat'), 'data', '-v7.3')

%% Transmit signal ------------------------------------------------------
tx = load(fullfile(rawDir, 'WBAT_tx.mat'));   % fs_dec, ytx
save(fullfile(outDir, 'transmit_signal.mat'), '-struct', 'tx', 'fs_dec', 'ytx')

%% CTD casts ------------------------------------------------------------
KEEP = {'mtime', 'temp', 'salt', 'press'};
for year = [2021 2022]
    ctd = load(fullfile(rawDir, sprintf('WBAT_CTD_%d.mat', year)));
    Dat = rmfield(ctd.Dat, setdiff(fieldnames(ctd.Dat), KEEP));
    save(fullfile(outDir, sprintf('ctd_%d.mat', year)), 'Dat', '-v7.3')
end

fprintf('Wrote reduced datasets to %s\n', outDir);
end
