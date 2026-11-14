
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
%                 Project 1: Standard Maneuver Sequence                   %
%                                                                         %
% ======================================================================= %

% PURPOSE
% Design a transfer from the initial geocentric orbit to the parking orbit
% using a standard sequence of four impulsive maneuvers:
% 1) orbital plane change;
% 2) argument-of-periapsis change;
% 3-4) apoapsis-to-periapsis bitangent transfer.
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
% Subscripts i, f, and t1-t3 denote the initial, final, and transfer orbits.

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


%% ORBITAL PLANE-CHANGE MANEUVER (MANEUVER 1)
% ----------------------------------------------------------------------- %

% Characterization of transfer orbit 1: same shape as the initial orbit
% and same orbital plane as the final orbit
a_t1  = a_i;         % semi-major axis    
e_t1  = e_i;         % eccentricity
i_t1  = i_f;         % inclination
OM_t1 = OM_f;        % right ascension of the ascending node (RAAN)

% Maneuver cost, argument of periapsis of transfer orbit 1,
% and true anomaly on the transfer orbit at maneuver 1
[DV1, om_t1, th1] = changeOrbitalPlane(a_i,  e_i,  i_i,   OM_i, ...
                                       om_i, i_t1, OM_t1, mu);

% Position of maneuver 1
r1 = par2car(a_i, e_i, i_i, OM_i, om_i, th1, mu);

% Initial orbit and transfer orbit 1 representation
[xt1, yt1, zt1] = plotOrbitSpace(a_t1, e_t1, i_t1, OM_t1, om_t1, 0, 2*pi);

figure(2);
plot3(xi,yi,zi, linewidth=2,color="blue");
hold on
plot3(xt1,yt1,zt1, linewidth=2,color="cyan");
plot3(0,0,0, marker='*',MarkerSize=10,LineWidth=2);

plot3(r_start(1),r_start(2),r_start(3), marker=".",MarkerSize=20,color="blue");
plot3(r1(1),r1(2),r1(3), marker="o",markersize=10,linewidth=2,color="cyan");
hold off

legend("Initial Orbit", "Transfer Orbit 1", "Earth", ...
        "Start Position", "Maneuver 1 Position");
title("Orbital Plane-Change Maneuver");
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), axis equal, view([46,30]);
axis ([-3e4 7e4 -3e4 8e4 -6e4 2e4]), grid on;


%% ARGUMENT-OF-PERIAPSIS CHANGE MANEUVER (MANEUVER 2)
% ----------------------------------------------------------------------- %

% Characterization of transfer orbit 2: same shape as transfer orbit 1
% and same argument of periapsis as the final orbit
a_t2  = a_t1;       % semi-major axis
e_t2  = e_t1;       % eccentricity
i_t2  = i_t1;       % inclination
OM_t2 = OM_t1;      % right ascension of the ascending node (RAAN)
om_t2 = om_f;       % argument of periapsis

% Maneuver cost, true anomaly on transfer orbit 1 at maneuver 2,
% and true anomaly on transfer orbit 2 at maneuver 2
[DV2, th2_t1, th2_t2] = changePericenterArg(a_t1, e_t1, om_t1, om_t2, mu);

% Two pairs of true anomalies satisfy the maneuver geometry. Select the pair
% corresponding to the maneuver point closest to periapsis in order to
% minimize the transfer time.
th2_t1 = th2_t1(1);
th2_t2 = th2_t2(1);

% Maneuver 2 position
r2 = par2car(a_t1, e_t1, i_t1, OM_t1, om_t1, th2_t1, mu);

% Transfer orbit 1 and transfer orbit 2 representation
[xt2, yt2, zt2] = plotOrbitSpace(a_t2, e_t2, i_t2, OM_t2, om_t2, 0, 2*pi);

figure(3);
plot3(xt1,yt1,zt1, linewidth=2,color="cyan");
hold on
plot3(xt2,yt2,zt2, linewidth=2,color="#FF8C00");
plot3(0,0,0, marker='*',MarkerSize=10,LineWidth=2);

plot3(r1(1),r1(2),r1(3), marker="o",markersize=10,linewidth=2,color="cyan");
plot3(r2(1),r2(2),r2(3), marker="o",MarkerSize=10,linewidth=2,color="#FF8C00");
hold off

legend("Transfer Orbit 1", "Transfer Orbit 2", "Earth", ...
       "Maneuver 1 Position", "Maneuver 2 Position");
title("Argument-of-Periapsis Change Maneuver");
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), axis equal, view([46,30]);
axis ([-3e4 7e4 -3e4 8e4 -6e4 2e4]), grid on;


%% ORBITAL SHAPE-CHANGE MANEUVER (MANEUVERS 3 AND 4)
% ----------------------------------------------------------------------- %

% Apoapsis-to-periapsis bitangent transfer.
% Maneuver 3: apoapsis, maneuver 4: periapsis.
% Characterization of transfer orbit 3: same orbital plane and argument
% of periapsis as the final orbit, same apoapsis as transfer orbit 2,
% and same periapsis as the final orbit

rA_t3 = a_t2 * (1+e_t2);             % apoapsis radius
rP_t3 = a_f  * (1-e_f );             % periapsis radius
a_t3  = (rA_t3+rP_t3)/2;             % semi-major axis
e_t3  = (rA_t3-rP_t3)/(rA_t3+rP_t3); % eccentricity
i_t3  = i_t2;                        % inclination
OM_t3 = OM_t2;                       % right ascension of the ascending node (RAAN)
om_t3 = om_t2;                       % argument of periapsis
 
% Costs of maneuvers 3 and 4
[DV3, DV4] = bitangentTransfer(a_t2, e_t2, a_f, e_f, 'ap', mu);

% Positions of maneuvers 3 and 4
r3 = par2car(a_t3, e_t3, i_t3, OM_t3, om_t3, pi, mu);
r4 = par2car(a_t3, e_t3, i_t3, OM_t3, om_t3, 0, mu);

% Transfer orbit 2 and transfer orbit 3 representation
[xt3, yt3, zt3] = plotOrbitSpace(a_t3, e_t3, i_t3, OM_t3, om_t3, 0, 2*pi);

figure(4);
plot3(xt2,yt2,zt2, linewidth=2,color="#FF8C00");
hold on
plot3(xt3,yt3,zt3, linewidth=2,color="#006400");
plot3(0,0,0, marker='*',MarkerSize=10,LineWidth=2);

plot3(r2(1),r2(2),r2(3), marker="o",MarkerSize=10,linewidth=2,color="#FF8C00");
plot3(r3(1),r3(2),r3(3), marker="o",MarkerSize=10,linewidth=2,color="#006400");
hold off

legend("Transfer Orbit 2", "Transfer Orbit 3", "Earth", ...
        "Maneuver 2 Position", "Maneuver 3 Position");
title("Bitangent Transfer: Apoapsis Maneuver");
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), axis equal, view([46,30]);
axis ([-3e4 7e4 -3e4 8e4 -6e4 2e4]), grid on;

% Transfer orbit 3 and final orbit representation
figure(5)
hold on
plot3(xt3,yt3,zt3, linewidth=2,color="#006400");
plot3(xf,yf,zf, linewidth=2,color="red");
plot3(0,0,0, marker='*',MarkerSize=10,LineWidth=2);

plot3(r3(1),r3(2),r3(3), marker="o",MarkerSize=10,linewidth=2,color="#006400");
plot3(r4(1),r4(2),r4(3), marker="o",MarkerSize=10,linewidth=2,color="black");
plot3(r_end(1),r_end(2),r_end(3), marker=".",MarkerSize=20,color="red");
hold off

legend("Transfer Orbit 3", "Final Orbit", "Earth", ...
       "Maneuver 3 Position", "Maneuver 4 Position", "End Position");
title("Bitangent Transfer: Periapsis Maneuver");
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), axis equal, view([46,30]);
axis ([-3e4 7e4 -3e4 8e4 -6e4 2e4]), grid on;

%% TOTAL TRANSFER TIME AND COST
% ----------------------------------------------------------------------- %

DV = DV1 + DV2 + DV3 + DV4;     % total delta-v cost [km/s]

DT1 = timeOfFlight(a_t1, e_t1, th1, th2_t1, mu); % time on transfer orbit 1 [s]
DT2 = timeOfFlight(a_t2, e_t2, th2_t2, pi, mu);  % time on transfer orbit 2 [s]
DT3 = timeOfFlight(a_t3, e_t3, pi, 0, mu);       % time on transfer orbit 3 [s]

DT = DT1 + DT2 + DT3;           % total transfer time [s]

% ============================= END =====================================
