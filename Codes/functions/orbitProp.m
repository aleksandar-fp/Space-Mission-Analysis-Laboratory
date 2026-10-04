function th_out = orbitProp(a, e, Dt, th_in, mu, type)
%ORBITPROP Propagate the true anomaly along an elliptic orbit.
%
% th1 = orbitProp(a, e, Dt, th2, mu, 1)
% th2 = orbitProp(a, e, Dt, th1, mu, 2)
%
% For type = 1, the function propagates backward by Dt and returns the
% anomaly from which th_in is reached after Dt. For type = 2, it propagates
% forward by Dt from th_in.
%
% Inputs:
%   a, e   - Semi-major axis [km] and eccentricity [-].
%   Dt     - Propagation interval [s].
%   th_in  - Known true anomaly [rad].
%   mu     - Gravitational parameter [km^3/s^2].
%   type   - Propagation direction: 1 backward, 2 forward.
%
% Output:
%   th_out - Propagated true anomaly in [0, 2*pi) [rad].
%
% This function is intended for elliptic orbits (0 <= e < 1).

if type == 1
    % Convert the known final anomaly to time since periapsis
    th2 = th_in;
    E2 = 2*atan2(sqrt((1-e)/(1+e)) * sin(th2/2), cos(th2/2));
    t2 = sqrt(a^3/mu) * (E2 - e*sin(E2));

    % Propagate backward and solve Kepler's equation
    t1 = t2 - Dt;
    f = @(E1) sqrt(a^3/mu) * (E1 - e*sin(E1)) - t1;
    E1 = fsolve(f,1);

    th1 = 2*atan2(sqrt((1+e)/(1-e)) * sin(E1/2), cos(E1/2));
    th1 = th1 +2*pi*(th1<0) -(th1>2*pi)*2*pi;
    
    th_out = th1;


elseif type == 2
    % Convert the known initial anomaly to time since periapsis
    th1 = th_in;
    E1 = 2*atan2(sqrt((1-e)/(1+e)) * sin(th1/2), cos(th1/2));
    t1 = sqrt(a^3/mu) * (E1 - e*sin(E1));
    
    % Propagate forward and solve Kepler's equation
    t2 = t1 + Dt;
    f = @(E2) sqrt(a^3/mu) * (E2 - e*sin(E2)) - t2;
    E2 = fsolve(f,1);

    th2 = 2*atan2(sqrt((1+e)/(1-e)) * sin(E2/2), cos(E2/2));
    th2 = th2 +2*pi*(th2<0) -(th2>2*pi)*2*pi;

    th_out = th2;
end

th_out = th_out + 2*pi*(th_out<0);

end
