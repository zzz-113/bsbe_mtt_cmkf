function [F, Q] = constantVelocityModel(dt, sigmaAccel)
%CONSTANTVELOCITYMODEL  State transition and process noise for a 3-D CV model.
%
%   [F, Q] = constantVelocityModel(DT, SIGMAACCEL) returns the 6x6 matrices of
%   a constant-velocity model with state x = [x y z vx vy vz]' over a time
%   step DT:
%
%       F = [ I3  dt*I3 ]      position(k) = position(k-1) + dt*velocity(k-1)
%           [ 0   I3    ]      velocity(k) = velocity(k-1)
%
%   SIGMAACCEL is a 1x3 vector of per-axis acceleration standard deviations
%   [m/s^2]. Treating the acceleration as piecewise constant and white over
%   each step gives
%
%       Q = [ dt^4/4 * Sa   dt^3/2 * Sa ],     Sa = diag(sigmaAccel.^2)
%           [ dt^3/2 * Sa   dt^2   * Sa ]
%
%   The axes are kept separate so that vertical and horizontal manoeuvring can
%   be given different weights.
%
%   Reference
%     Bar-Shalom, Y., Li, X.-R. and Kirubarajan, T. (2001). Estimation with
%     Applications to Tracking and Navigation. Wiley, section 6.2.

arguments
    dt (1,1) double
    sigmaAccel (1,3) double
end

F = [  eye(3), dt * eye(3);
     zeros(3),      eye(3)];

Sa = diag(sigmaAccel.^2);

Q = [dt^4/4 * Sa, dt^3/2 * Sa;
     dt^3/2 * Sa, dt^2   * Sa];
end
