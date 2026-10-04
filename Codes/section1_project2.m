
% ======================================================================= %
%                                                                         %
%                    SPACE MISSION ANALYSIS LABORATORY                    %
%                                                                         %
%           Final Report - Aerospace Engineering Degree Program           %
%                                                                         %
%                        Politecnico di Milano                            %
%                                                                         %
%                    Aleksandar Franchini Petrovic                        %
%                                                                         %
%              ------------------------------------------                 %
%									                                      %
%              Section 1: Transfer to the Parking Orbit                   %
%                                                                         %
%                Project 2: Near-Single-Impulse Maneuver                  %
%                                                                         %
% ======================================================================= %

% PURPOSE
% Design an alternative transfer from the initial geocentric orbit to the
% parking orbit. The maneuver sequence consists of a small correction at
% the initial-orbit apoapsis followed by the main impulse at the point where
% the transfer orbit intersects the final orbit.
%
% UNITS
% Distances: [km]
% Velocities and delta-v costs: [km/s]
% Times: [s]
% Angles: [rad]
%
% NOTATION
% a, e, i, OM, om, and th denote the semi-major axis, eccentricity,
% inclination, RAAN, argument of periapsis, and true anomaly, respectively.
% Subscripts i, f, and t denote the initial, final, and transfer orbits.

mu = 398600.4;          % Earth's gravitational parameter [km^3/s^2]

%% INITIAL ORBIT DATA
% Orbital parameters
a_i      = 24400.00;         % [km]  semi-major axis
e_i      = 0.728300;         % [-]   eccentricity
i_i      = 0.104700;         % [rad] inclination
OM_i     = 1.400000;         % [rad] right ascension of the ascending node (RAAN)
om_i     = 3.107000;         % [rad] argument of periapsis
th_start = 1.487000;         % [rad] true anomaly at the initial time

% Position and velocity at the initial time
[r_start, v_start]= par2car(a_i, e_i, i_i, OM_i, om_i, th_start, mu);

%% FINAL ORBIT DATA
% Position and velocity at the final time
r_end = [-5985.618600;  -7930.756300;  6091.925600];
v_end = [    4.532000;     -5.596000;    -2.832000];

% Final-orbit parameters and true anomaly at the final time
[a_f, e_f, i_f, OM_f, om_f, th_end] = car2par(r_end, v_end, mu);


%% INITIAL AND FINAL ORBITS REPRESENTATION
[xi, yi, zi] = plotOrbitSpace(a_i, e_i, i_i, OM_i, om_i, 0, 2*pi);
[xf, yf, zf] = plotOrbitSpace(a_f, e_f, i_f, OM_f, om_f, 0, 2*pi);

figure(1);
plot3(xi,yi,zi, linewidth=2,color="blue");
hold on
plot3(xf,yf,zf, linewidth=2,color="red");
plot3(0,0,0, marker='*',MarkerSize=10,LineWidth=2);

plot3(r_start(1),r_start(2),r_start(3), marker=".",MarkerSize=20,color="blue");
plot3(r_end(1),r_end(2),r_end(3), marker=".",MarkerSize=20,color="red");
hold off

legend("Initial Orbit", "Final Orbit", "Earth", ...
        "Start Position", "End Position");
title("Initial and Final Orbits");
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), axis equal, view([46,30]);
axis ([-3e4 7e4 -3e4 8e4 -6e4 2e4]), grid on;


%% MAIN-IMPULSE POINT DETERMINATION
% ----------------------------------------------------------------------- %

% Find the point on the final orbit that lies in the initial orbital plane.
% Two points satisfy this condition; select the one closest to the point of
% minimum distance between the two orbits. At the desired point, the final
% orbit position vector must be perpendicular to the initial-orbit angular
% momentum vector; therefore, their dot product is set equal to zero.

% Final-orbit position vector as a function of true anomaly
rF = @(th_f) par2car(a_f, e_f, i_f, OM_f, om_f, th_f, mu);

% Initial-orbit unit specific angular momentum vector
h = cross(r_start, v_start) / norm(cross(r_start, v_start));

% Dot product between the final-orbit position vector and the initial-orbit
% unit specific angular momentum vector
F = @(th_f) dot(rF(th_f), h);

% Plot F to identify the initial guess for solving F = 0
n = 50;
Z = zeros(n,1);
th_f = linspace(0,2*pi,n);

for k=1:n
    Z(k) = F(th_f(k));
end

figure(2);
plot(th_f, Z, LineWidth=3);
xlabel("\theta_f [rad]"), ylabel("F(\theta_f) [km]"), grid on;
title("Main-Impulse Point Determination");

% The solution near the point of minimum distance between the two orbits
% corresponds to the second root.

thI_f0 = 4.2;  % initial guess [rad]

% Solution: final-orbit true anomaly at the main-impulse point
thI_f = fsolve(F, thI_f0);

% Corresponding main-impulse position
rI = par2car(a_f, e_f, i_f, OM_f, om_f, thI_f, mu);



%% TRANSFER ORBIT CHARACTERIZATION
% --------------------------------------------------------------------- %

% The transfer orbit lies in the initial orbital plane and is apsidally
% aligned with the initial orbit; in particular, they share the same
% apoapsis, where the correction maneuver is performed. The transfer
% orbit is determined by requiring it to pass through the point
% identified in the previous section. 


% Initial-orbit apoapsis position, where the correction is performed
rA_i = par2car(a_i,e_i,i_i,OM_i,om_i,pi,mu);

rA_t = norm(rA_i);  % transfer-orbit apoapsis radius
i_t = i_i;          % inclination
OM_t = OM_i;        % right ascension of the ascending node
om_t = om_i;        % argument of periapsis

% Semi-major axis, eccentricity, and transfer-orbit true anomaly at the
% main-impulse point
[a_t, e_t, thI_t] = por2for(rI, rA_t, i_t, OM_t, om_t);

%% ORBIT REPRESENTATION
% --------------------------------------------------------------------- %

[xt, yt, zt] = plotOrbitSpace(a_t, e_t, i_t, OM_t, om_t, pi, thI_t);

figure(3);
plot3(xi,yi,zi, linewidth=2,color="blue");                              
hold on
plot3(xf,yf,zf, linewidth=2,color="red");
plot3(0,0,0,marker='*', MarkerSize=10, LineWidth=2);
plot3(xt,yt,zt, linewidth=2,color="cyan", LineStyle=":");

plot3(r_start(1),r_start(2),r_start(3), marker=".",MarkerSize=25,color="blue");
plot3(r_end(1),r_end(2),r_end(3), marker=".",MarkerSize=25,color="red");

plot3(rA_i(1),rA_i(2),rA_i(3),marker="square",MarkerSize=15, color="cyan",LineWidth=3);
plot3(rI(1),rI(2),rI(3), marker="o",MarkerSize=15,color="cyan",LineWidth=3);
hold off

legend("Initial Orbit", "Final Orbit", "Earth", "Transfer Orbit",...
        "Start Position", "End Position", "Correction Maneuver", ...
        "Main Maneuver");

title("Near-Single-Impulse Transfer");
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), axis equal;
axis ([-3e4 7e4 -3e4 8e4 -6e4 2e4]), view([46,30]), grid on;

%% TOTAL TRANSFER TIME AND COST
% ----------------------------------------------------------------------- %

% Time spent on the transfer orbit [s]
DT = timeOfFlight(a_t, e_t, pi, thI_t, mu);

% Correction-maneuver cost [km/s]
[~, vA_i] = par2car(a_i, e_i, i_i, OM_i, om_i, pi, mu);
[~, vA_t] = par2car(a_t, e_t, i_t, OM_t, om_t, pi, mu);
DV1 = norm(vA_t - vA_i);

% Main-maneuver cost [km/s]
[~, vI_t] = par2car(a_t, e_t, i_t, OM_t, om_t, thI_t, mu);
[~, vI_f] = par2car(a_f, e_f, i_f, OM_f, om_f, thI_f, mu);
DV2 = norm(vI_f - vI_t);

DV = DV1 + DV2;     % total delta-v cost [km/s]

% ============================= END =====================================
