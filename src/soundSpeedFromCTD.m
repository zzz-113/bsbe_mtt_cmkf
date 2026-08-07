function c = soundSpeedFromCTD(t, ctdFiles)
%SOUNDSPEEDFROMCTD  Sound speed at given times from moored CTD records.
%
%   c = soundSpeedFromCTD(T) returns the sound speed [m/s] at the datetimes T,
%   linearly interpolated from the CTD casts shipped in data/. Times outside
%   the CTD record return NaN.
%
%   c = soundSpeedFromCTD(T, CTDFILES) uses the cell array of file paths
%   CTDFILES instead. Each file must contain a struct 'Dat' with fields
%   mtime (MATLAB datenum), temp [degC], salt and press [dbar]; records from
%   all files are concatenated in the order given.
%
%   Sound speed is evaluated with gsw_sound_speed from the GSW Oceanographic
%   Toolbox (see external/). The CTD salinity and temperature are passed
%   through as recorded, which treats them as Absolute Salinity and
%   Conservative Temperature.

arguments
    t datetime
    ctdFiles cell = defaultCtdFiles()
end

mtime = [];
temp  = [];
salt  = [];
press = [];

for k = 1:numel(ctdFiles)
    d = load(ctdFiles{k}, 'Dat');
    mtime = [mtime, d.Dat.mtime]; %#ok<AGROW>
    temp  = [temp,  d.Dat.temp];  %#ok<AGROW>
    salt  = [salt,  d.Dat.salt];  %#ok<AGROW>
    press = [press, d.Dat.press]; %#ok<AGROW>
end

ctdTime = datetime(mtime, 'ConvertFrom', 'datenum');
c = interp1(ctdTime, gsw_sound_speed(salt, temp, press), t);
end


function files = defaultCtdFiles()
% Locate the CTD records shipped alongside the code, in data/.
dataDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'data');
files = {fullfile(dataDir, 'ctd_2021.mat'), fullfile(dataDir, 'ctd_2022.mat')};
end
