function [out1,out2] = par2car(a, e, i, OM, om, th, mu, q)
%PAR2CAR Convert Keplerian orbital elements to Cartesian state vectors.
%
% [r, v] = par2car(a, e, i, OM, om, th, mu)
% [v, r] = par2car(a, e, i, OM, om, th, mu, 1)
%
% Inputs:
%   a      - Semi-major axis [km].
%   e      - Eccentricity [-].
%   i      - Inclination [rad].
%   OM     - Right ascension of the ascending node (RAAN) [rad].
%   om     - Argument of periapsis [rad].
%   th     - True anomaly [rad].
%   mu     - Gravitational parameter [km^3/s^2].
%   q      - Optional output-order flag. If omitted or zero, position is
%            returned first; if one, velocity is returned first.
%
% Outputs:
%   r      - Position vector [km].
%   v      - Velocity vector [km/s].

% Semi-latus rectum and orbital radius
p = (1-e^2)*a;
r_mod = p/(1 + e*cos(th));

% Position and velocity in the perifocal frame
r_per = r_mod*[cos(th); sin(th); 0];
v_per = sqrt(mu/p) * [-sin(th); e + cos(th); 0];

% Rotation from the perifocal frame to the inertial frame. The Aerospace
% Toolbox rotation functions use degrees.
R_OM    = rotz(OM *180/pi);
R_i     = rotx(i  *180/pi);
R_om    = rotz(om *180/pi);

R = R_OM*R_i*R_om;

% Position and velocity in the inertial frame
r = R * r_per;
v = R * v_per;

% Select the requested output order
if nargin == 7
    q = 0;
end
out1 = ~q*r +  q*v;
out2 =  q*r + ~q*v;


end
