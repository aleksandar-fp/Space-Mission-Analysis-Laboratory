function deltat = timeOfFlight(a, e, th1, th2, mu)
%TIMEOFFLIGHT Compute the forward time of flight along an elliptic orbit.
%
% deltat = timeOfFlight(a, e, th1, th2, mu)
%
% Returns the time required to move from th1 to th2 in the direction of
% increasing true anomaly. If th2 < th1, one passage through 2*pi is added.
%
% Inputs:
%   a, e   - Semi-major axis [km] and eccentricity [-].
%   th1    - Initial true anomaly [rad].
%   th2    - Final true anomaly [rad].
%   mu     - Gravitational parameter [km^3/s^2].
%
% Output:
%   deltat - Forward time of flight [s].
%
% This function is intended for elliptic orbits (0 <= e < 1).

% Convert true anomalies to eccentric anomalies
E1 = 2*atan2(sqrt((1-e)/(1+e)) * sin(th1/2), cos(th1/2));
E2 = 2*atan2(sqrt((1-e)/(1+e)) * sin(th2/2), cos(th2/2));

% Compute elapsed times from periapsis using Kepler's equation
t1 = sqrt(a^3/mu) * (E1 - e*sin(E1));
t2 = sqrt(a^3/mu) * (E2 - e*sin(E2) + (th2<th1)*2*pi);

% Time of flight between the two positions
deltat = t2 - t1;

end
