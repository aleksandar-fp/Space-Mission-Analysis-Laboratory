function [DeltaV, thi, thf] = changePericenterArg(a, e, om_i, om_f, mu)
%CHANGEPERICENTERARG Compute an argument-of-periapsis change maneuver.
%
% [DeltaV, thi, thf] = changePericenterArg(a, e, om_i, om_f, mu)
%
% The initial and final orbits have the same shape and orbital plane but
% different arguments of periapsis. The two possible intersection points
% produce the same impulse magnitude.
%
% Inputs:
%   a, e       - Semi-major axis [km] and eccentricity [-].
%   om_i, om_f - Initial and final arguments of periapsis [rad].
%   mu         - Gravitational parameter [km^3/s^2].
%
% Outputs:
%   DeltaV     - Maneuver impulse magnitude [km/s].
%   thi        - True anomalies on the initial orbit [rad].
%   thf        - Corresponding true anomalies on the final orbit [rad].
%


delta_om = om_f - om_i;

% True anomalies of the two orbit-intersection points
thi = [        delta_om/2;     pi + delta_om/2];
thf = [ 2*pi - delta_om/2;     pi - delta_om/2];

% Express all anomalies in the interval [0, 2*pi]
for k=1:2
    thi(k) = thi(k) +2*pi*(thi(k)<0) -2*pi*(thi(k)>2*pi);
    thf(k) = thf(k) +2*pi*(thf(k)<0) -2*pi*(thf(k)>2*pi);
end

% Impulse magnitude at either intersection point
p = a*(1-e^2);
DeltaV = 2*sqrt(mu/p) * e * sin(abs(delta_om/2));

end
