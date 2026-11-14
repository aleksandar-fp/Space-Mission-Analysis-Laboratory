function [a, e, i, OM, om, th] = car2par(R, V, mu)
%CAR2PAR Convert Cartesian state vectors to Keplerian orbital elements.
%
% [a, e, i, OM, om, th] = car2par(R, V, mu)
%
% Inputs:
%   R      - Position vector [km].
%   V      - Velocity vector [km/s].
%   mu     - Gravitational parameter [km^3/s^2].
%
% Outputs:
%   a      - Semi-major axis [km].
%   e      - Eccentricity [-].
%   i      - Inclination [rad].
%   OM     - Right ascension of the ascending node (RAAN) [rad].
%   om     - Argument of periapsis [rad].
%   th     - True anomaly [rad].
%
% The classical-element formulation is singular for circular or equatorial
% orbits and is therefore intended for nondegenerate cases.

% Position and velocity magnitudes
r = norm(R);
v = norm(V);

% Semi-major axis from the specific orbital energy
a = (2/r - v^2/mu)^(-1);

% Specific angular momentum
H = cross(R, V);
h = norm(H);

% Eccentricity vector and magnitude
E = cross(V, H) / mu - R/r;
e = norm(E);

% Inclination and unit vector along the line of nodes
i = acos(H(3)/h);
N = cross([0;0;1], H) / norm(cross([0;0;1], H));

% RAAN, argument of periapsis, and true anomaly with quadrant corrections
OM = (N(2)<0)*2*pi + (-1)^(N(2)<0)*acos(N(1));
om = (E(3)<0)*2*pi + (-1)^(E(3)<0)*acos(dot(N,E)/e);
th = (dot(V,R)<0)*2*pi + (-1)^(dot(V,R)<0)*acos(dot(R,E)/(r*e));
th = th + 2*pi*(th<0);

end
