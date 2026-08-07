function Rcart = convertedMeasurementCovariance(r, alpha, beta, Rsb)
%CONVERTEDMEASUREMENTCOVARIANCE  Split-beam measurement covariance in Cartesian.
%
%   RCART = convertedMeasurementCovariance(R, ALPHA, BETA, RSB) propagates the
%   3x3 measurement covariance RSB, given in split-beam coordinates
%   [r, alpha, beta], through the conversion
%
%       x = r * alpha
%       y = r * beta
%       z = r * gamma,      gamma = sqrt(1 - alpha^2 - beta^2)
%
%   using the first-order (Jacobian) approximation RCART = G * RSB * G', with
%
%       G = d(x,y,z) / d(r,alpha,beta)
%         = [ alpha    r                0
%             beta     0                r
%             gamma   -r*alpha/gamma   -r*beta/gamma ].
%
%   The point (R, ALPHA, BETA) at which G is evaluated should be the predicted
%   position, not the measured one. Evaluating at the measurement would make
%   RCART a function of the measurement noise and correlate the two, which is
%   the error the debiased/decorrelated CMKF is formulated to avoid.
%
%   gamma appears in the denominator, so RCART grows without bound as the
%   target approaches the 90-degree off-axis edge of the beam. In practice
%   detectTargets rejects detections well inside that limit.
%
%   Reference
%     Lerro, D. and Bar-Shalom, Y. (1993). Tracking with debiased consistent
%     converted measurements versus EKF. IEEE Trans. Aerospace and Electronic
%     Systems, 29(3), 1015-1022.

gamma = sqrt(1 - alpha^2 - beta^2);

G = [alpha,  r,               0;
     beta,   0,               r;
     gamma, -r*alpha/gamma,  -r*beta/gamma];

Rcart = G * Rsb * G';
Rcart = (Rcart + Rcart') / 2;   % discard asymmetry from round-off
end
