function tt = transferTime(th1_e, th2_a, om_guess)
%TRANSFERTIME Compute the flight time of the minimum-cost transfer orbit.
%
% tt = transferTime(th1_e, th2_a, om_guess)
%
% For prescribed Earth-departure and asteroid-arrival true anomalies, this
% function minimizes the total endpoint delta-v with respect to the transfer
% orbit argument of periapsis and returns the corresponding flight time.
%
% Inputs:
%   th1_e     - Earth true anomaly at departure [rad].
%   th2_a     - Asteroid true anomaly at rendezvous [rad].
%   om_guess  - Initial guess for the transfer-orbit argument of periapsis
%               [rad].
%
% Output:
%   tt        - Transfer time for the minimum-cost orbit [s].
%
% This scenario-specific function uses the fixed heliocentric orbital data
% defined in Section 2.

% Solar gravitational parameter
mu = 132712440018;  % [km^3/s^2]

% Earth orbital elements
a_e = 1.4946e+08;
e_e = 0.016;
i_e = 9.1920e-05;
OM_e = 2.7847;
om_e = 5.2643;

% Asteroid orbital elements
a_a = 1.893919*149597871;
e_a = 0.444290;
i_a = 4.19*pi/180;
OM_a = 298.91*pi/180;
om_a = 168.65*pi/180;

% Endpoint positions and velocities
[r1, v1_e]  = par2car(a_e, e_e, i_e, OM_e, om_e, th1_e, mu);
[r2, v2_a]  = par2car(a_a, e_a, i_a, OM_a, om_a, th2_a, mu);

% Transfer-orbit parameters as functions of the argument of periapsis
a_t     = @(om)  direct(r1,r2,om, 'a');
e_t     = @(om)  direct(r1,r2,om, 'e');
i_t     = @(om)  direct(r1,r2,om, 'i');
OM_t    = @(om)  direct(r1,r2,om, 'OM');
th1_t   = @(om)  direct(r1,r2,om, 'th1');
th2_t   = @(om)  direct(r1,r2,om, 'th2');

% Spacecraft velocities at the transfer endpoints
v1_t = @(om) par2car(a_t(om), e_t(om), i_t(om), OM_t(om), om, th1_t(om), mu, 1);
v2_t = @(om) par2car(a_t(om), e_t(om), i_t(om), OM_t(om), om, th2_t(om), mu, 1);

% Departure, arrival, and total delta-v
DV1 = @(om) norm(v1_t(om) - v1_e);
DV2 = @(om) norm(v2_t(om) - v2_a);
DV =  @(om) DV1(om) + DV2(om);

% Minimum-cost argument of periapsis
omS_t = fminsearch(DV, om_guess);

% Flight time along the resulting transfer orbit
tt = timeOfFlight(a_t(omS_t), e_t(omS_t), th1_t(omS_t), th2_t(omS_t), mu);

end
