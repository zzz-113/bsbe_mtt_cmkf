function pings = loadPings(matFile, transmitFile)
%LOADPINGS  Read a pulse-compressed WBAT/EK80 export into a compact struct.
%
%   pings = loadPings(MATFILE, TRANSMITFILE) reads MATFILE, which must contain
%   a variable 'data' in the layout produced by the WBAT export used in this
%   study, together with the transmit signal in TRANSMITFILE, and returns only
%   the fields the pipeline needs.
%
%   Expected layout of the stored 'data' variable:
%     data.echodata(k).compressed  [nSample x 3] complex pulse-compressed
%                                  samples for the three transducer sectors
%     data.echodata(k).timestamp   scalar MATLAB datenum
%     data.param(1).SampleInterval sample period [s]
%     data.param(1).PulseDuration  transmit pulse length [s]
%     data.param(1).FrequencyStart start of the transmit sweep [Hz]
%     data.param(1).FrequencyEnd   end of the transmit sweep [Hz]
%
%   TRANSMITFILE must contain 'ytx', the complex transmit signal at the same
%   sample rate as the echoes, as written by the WBAT export.
%
%   Output fields:
%     .compressed     [nSample x 3 x nPing] complex, sectors in columns
%     .time           [1 x nPing] datetime of each ping
%     .sampleInterval sample period [s]
%     .pulseDuration  transmit pulse length [s]
%     .frequencyStart start of the transmit sweep [Hz]
%     .frequencyEnd   end of the transmit sweep [Hz]
%     .transmitSignal complex transmit signal, column vector
%
%   Ping parameters are read from the first ping and assumed constant over
%   the file, which holds for a single WBAT recording phase.
%
%   To run the pipeline on data in a different format, write a replacement
%   for this function that returns the same fields; nothing downstream reads
%   the raw file.

arguments
    matFile (1,:) char
    transmitFile (1,:) char
end

raw = load(matFile, 'data');
if ~isfield(raw, 'data')
    error('loadPings:missingVariable', ...
        '%s does not contain a variable named ''data''.', matFile);
end
data = raw.data;

nPing   = numel(data.echodata);
nSample = size(data.echodata(1).compressed, 1);

% [data.echodata.compressed] concatenates ping-major as
% [ping1sector1 ping1sector2 ping1sector3 ping2sector1 ...], which reshapes
% directly to [sample x sector x ping].
compressed = [data.echodata.compressed];
pings.compressed = reshape(compressed, nSample, 3, nPing);

pings.time = datetime([data.echodata.timestamp], 'ConvertFrom', 'datenum');

p = data.param(1);
pings.sampleInterval = p.SampleInterval;
pings.pulseDuration  = p.PulseDuration;
pings.frequencyStart = p.FrequencyStart;
pings.frequencyEnd   = p.FrequencyEnd;

tx = load(transmitFile, 'ytx');
pings.transmitSignal = tx.ytx(:);
end
