function [x, y, z] = plotOrbitSpace(a, e, i, OM, om, th0, thf)
%PLOTORBITSPACE Sample a Keplerian trajectory in three-dimensional space.
%
% [x, y, z] = plotOrbitSpace(a, e, i, OM, om, th0, thf)
%
% Inputs:
%   a      - Semi-major axis [km].
%   e      - Eccentricity [-].
%   i      - Inclination [rad].
%   OM     - Right ascension of the ascending node (RAAN) [rad].
%   om     - Argument of periapsis [rad].
%   th0    - Initial true anomaly [rad].
%   thf    - Final true anomaly [rad].
%
% Outputs:
%   x, y, z - Sampled inertial Cartesian coordinates [km].

dth = 0.01; % angular sampling step [rad]
 
p = (1-e^2)*a;

% Rotation from the perifocal frame to the inertial frame
R_OM    = rotz(OM *180/pi);
R_i     = rotx(i  *180/pi);
R_om    = rotz(om *180/pi);
R = R_OM*R_i*R_om;

% If necessary, continue the arc through 2*pi
thf = thf + (thf<th0)*2*pi;

% Preallocate the coordinate arrays
x = zeros( size(th0:dth:thf+dth,2), 1);
y = zeros( size(th0:dth:thf+dth,2), 1);
z = zeros( size(th0:dth:thf+dth,2), 1);

k = 1;
for th = th0 : dth: thf+dth
    % Position in the perifocal and inertial frames
    r = p/(1 + e*cos(th));
    r_per = r*[cos(th); sin(th); 0];
    r_cel = R * r_per;
    
    x(k) = r_cel(1);
    y(k) = r_cel(2);
    z(k) = r_cel(3);
    k = k+1;
end

end
