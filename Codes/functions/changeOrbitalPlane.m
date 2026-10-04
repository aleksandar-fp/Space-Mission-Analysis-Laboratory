function [DeltaV, om_f, theta] = changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)
%CHANGEORBITALPLANE Compute an impulsive orbital-plane change.
%
% [DeltaV, om_f, theta] = changeOrbitalPlane(a, e, i_i, OM_i, om_i, i_f, OM_f, mu)
%
% The orbit shape is preserved while inclination and RAAN are changed. The
% maneuver is performed at the intersection of the two orbital planes.
%
% Inputs:
%   a, e       - Semi-major axis [km] and eccentricity [-].
%   i_i, OM_i  - Initial inclination and RAAN [rad].
%   om_i       - Initial argument of periapsis [rad].
%   i_f, OM_f  - Final inclination and RAAN [rad].
%   mu         - Gravitational parameter [km^3/s^2].
%
% Outputs:
%   DeltaV     - Plane-change impulse magnitude [km/s].
%   om_f       - Final argument of periapsis [rad].
%   theta      - Initial-orbit true anomaly at the maneuver point [rad].

% Changes in inclination and RAAN
delta_i = i_f - i_i;
delta_OM = OM_f - OM_i;

% Angle between the initial and final orbital planes and the associated
% arguments of latitude
alfa = acos( cos(i_i)*cos(i_f) + sin(i_i)*sin(i_f)*cos(delta_OM) );
sin_ui = sin(delta_OM)/sin(alfa)*sin(i_f);
sin_uf = sin(delta_OM)/sin(alfa)*sin(i_i);
cos_ui = (-cos(i_f)+cos(alfa)*cos(i_i)) / (sin(alfa)*sin(i_i));
cos_uf = (cos(i_i)-cos(alfa)*cos(i_f)) / (sin(alfa)*sin(i_f));

% Quadrant-selection flags
n = (delta_OM >= 0 && delta_i < 0) || (delta_OM < 0 && delta_i < 0);
m = (delta_OM >= 0 && delta_i < 0) || (delta_OM < 0 && delta_i >= 0);

% Maneuver geometry
u_i = atan2(sin_ui, (-1)^n * cos_ui);
u_f = atan2(sin_uf, (-1)^n * cos_uf);

theta = m * 2*pi +(-1)^m * u_i -om_i;
theta = theta + (theta<0)*2*pi;
om_f = m * 2*pi +(-1)^m * u_f -theta;
om_f = om_f + (om_f<0)*2*pi;

% Plane-change impulse at the maneuver point
DeltaV = 2*sqrt(mu/(a*(1-e^2))) * (1+e*cos(theta)) * sin(alfa/2);

end
