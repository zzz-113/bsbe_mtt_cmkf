function splitBeam = cartesianToSplitBeam(position)
%CARTESIANTOSPLITBEAM  Cartesian position to split-beam coordinates.
%
%   SPLITBEAM = cartesianToSplitBeam(POSITION) maps POSITION = [x y z] to
%   SPLITBEAM = [r alpha beta], where
%
%       r     = sqrt(x^2 + y^2 + z^2)   range [m]
%       alpha = x / r                   alongship direction cosine
%       beta  = y / r                   athwartship direction cosine
%
%   This is the inverse of the measurement conversion the tracker applies to
%   each detection, and is used to evaluate the conversion Jacobian at the
%   predicted position rather than at the measurement.
%
%   Positions within 1e-6 m of the transducer return [0 0 0], since the
%   direction cosines are undefined there.

x = position(1);
y = position(2);
z = position(3);

r = sqrt(x^2 + y^2 + z^2);

if r < 1e-6
    splitBeam = [0 0 0];
    return
end

splitBeam = [r, x / r, y / r];
end
