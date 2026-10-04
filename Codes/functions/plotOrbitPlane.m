function [x, y] = plotOrbitPlane(a, e, om, th0, thf)
%PLOTORBITPLANE Sample a Keplerian trajectory in its orbital plane.
%
% [x, y] = plotOrbitPlane(a, e, om, th0, thf)
%
% Inputs:
%   a      - Semi-major axis [km].
%   e      - Eccentricity [-].
%   om     - Argument of periapsis measured from the x-axis [rad].
%   th0    - Initial true anomaly [rad].
%   thf    - Final true anomaly [rad].
%
% Outputs:
%   x, y   - Sampled Cartesian coordinates in the orbital plane [km].
%
% The x-axis is aligned with the ascending-node direction, and the positive
% out-of-plane direction follows the orbital angular momentum vector.

dth = 0.01; % angular sampling step [rad]

p = (1-e^2)*a;

% Preallocate the coordinate arrays
x = zeros( size(th0:dth:thf+dth,2), 1);
y = zeros( size(th0:dth:thf+dth,2), 1);

k = 1;
for th = th0 : dth: thf+dth
    % Polar equation of the conic and in-plane Cartesian coordinates
    r = p/(1 + e*cos(th));
    x(k) = r * cos(om + th);
    y(k) = r * sin(om + th);
    k = k+1;
end

end
