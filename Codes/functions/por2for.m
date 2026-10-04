function [a, e, th] = por2for(r, rA, i, OM, om)
%POR2FOR Determine orbit shape from a position and prescribed orientation.
%
% [a, e, th] = por2for(r, rA, i, OM, om)
%
% Determines the elliptic orbit that passes through r, has apoapsis radius
% rA, and has the prescribed spatial orientation.
%
% Inputs:
%   r      - Position vector on the orbit [km].
%   rA     - Apoapsis radius [km].
%   i      - Inclination [rad].
%   OM     - Right ascension of the ascending node (RAAN) [rad].
%   om     - Argument of periapsis [rad].
%
% Outputs:
%   a      - Semi-major axis [km].
%   e      - Eccentricity [-].
%   th     - True anomaly corresponding to r [rad].

r_mod = norm(r);

% Rotation from the perifocal frame to the inertial frame
R_OM    = rotz(OM *180/pi);
R_i     = rotx(i  *180/pi);
R_om    = rotz(om *180/pi);
R = R_OM*R_i*R_om;

% Transform the position to the perifocal frame and recover true anomaly
r_per = R' * r;
th = atan2(r_per(2),r_per(1));
th = th + 2*pi*(th<0);

% Determine eccentricity from the conic equation and apoapsis constraint
e = (rA - r_mod) / (r_mod*cos(th) + rA);

% Semi-major axis from the apoapsis radius
a = rA / (1+e);

end
