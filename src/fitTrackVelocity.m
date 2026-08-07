function tracks = fitTrackVelocity(tracks, env)
%FITTRACKVELOCITY  Least-squares velocity of each track from its measurements.
%
%   TRACKS = fitTrackVelocity(TRACKS, ENV) fits an independent straight line
%   in time to the x, y and z measurements logged for every track and adds:
%
%     velFit     [vx vy vz] fitted velocity [m/s]
%     speedFit   norm of velFit [m/s]
%     headingFit [azimuth elevation] of velFit [deg]; elevation is positive
%                upwards, i.e. towards decreasing z
%     posFit0    fitted position at the track's first ping [m]
%     velFitSE   standard error of each velocity component [m/s]
%     velFitR2   coefficient of determination per axis
%
%   The fit uses the converted measurements (column 2 of posLog), not the
%   filter output, so the result is an independent estimate rather than a
%   restatement of the filter's velocity state.
%
%   Tracks with fewer than two logged positions get NaN, as do the standard
%   errors of tracks with exactly two, which have no residual degrees of
%   freedom.

arguments
    tracks struct
    env struct
end

for n = 1:numel(tracks)
    posLog    = tracks(n).posLog;
    numPoints = size(posLog, 1);

    velocity      = nan(1, 3);
    position0     = nan(1, 3);
    standardError = nan(1, 3);
    rSquared      = nan(1, 3);
    heading       = nan(1, 2);

    if numPoints >= 2
        pingIndices  = cell2mat(posLog(:, 1));
        positions    = cell2mat(posLog(:, 2));
        timeRelative = (pingIndices - pingIndices(1)) * env.pingInterval;

        X = [timeRelative(:), ones(numPoints, 1)];

        for axis = 1:3
            y = positions(:, axis);
            coefficients = X \ y;
            velocity(axis)  = coefficients(1);
            position0(axis) = coefficients(2);

            residuals  = y - X * coefficients;
            ssResidual = sum(residuals.^2);
            ssTotal    = sum((y - mean(y)).^2);

            if ssTotal > 0
                rSquared(axis) = 1 - ssResidual / ssTotal;
            end

            if numPoints > 2
                % The (1,1) element of inv(X'X) scales the slope variance, so
                % the inverse itself is needed rather than a solve.
                XtXinv = inv(X' * X);
                standardError(axis) = sqrt(ssResidual / (numPoints - 2) * XtXinv(1, 1));
            end
        end

        if norm(velocity) > 1e-6
            heading = [atan2d(velocity(2), velocity(1)), ...
                       atan2d(-velocity(3), norm(velocity(1:2)))];
        end
    end

    tracks(n).velFit     = velocity;
    tracks(n).speedFit   = norm(velocity);
    tracks(n).headingFit = heading;
    tracks(n).posFit0    = position0;
    tracks(n).velFitSE   = standardError;
    tracks(n).velFitR2   = rSquared;
end
end
