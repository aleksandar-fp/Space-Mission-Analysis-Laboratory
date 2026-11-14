
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
%                       Section 2: Direct Transfer                       %
%                                                                         %
% ======================================================================= %

% PURPOSE
% Design a direct heliocentric transfer from Earth to the asteroid. For a
% prescribed Earth departure anomaly, determine the transfer orbit that
% minimizes the total delta-v while ensuring rendezvous with the asteroid.
% The departure anomaly is then varied over discrete values to identify the
% lowest-cost solution.
%
% UNITS
% Distances: [km]
% Velocities and delta-v costs: [km/s]
% Times: [s]
% Angles: [rad], unless otherwise specified
%
% NOTATION
% a, e, i, OM, om, and th denote the semi-major axis, eccentricity,
% inclination, RAAN, argument of periapsis, and true anomaly, respectively.
% Subscripts e, a, and t denote Earth, asteroid, and transfer orbit.
% Subscripts 1 and 2 denote departure and arrival, respectively.

mu = 132712440018;      % Sun's gravitational parameter [km^3/s^2]

%% EARTH ORBIT DATA
% Orbital parameters
a_e  = 1.4946e+08;       % [km] semi-major axis
e_e  = 0.016;            % [-]   eccentricity
i_e  = 9.1920e-05;       % [rad] inclination
OM_e = 2.7847;           % [rad] right ascension of the ascending node (RAAN)
om_e = 5.2643;           % [rad] argument of periapsis

%% ASTEROID ORBIT DATA
% Orbital parameters
a_a  = 1.893919*149597871;   % [km] semi-major axis
e_a  = 0.444290;             % [-]   eccentricity
i_a  = 4.19*pi/180;          % [rad] inclination
OM_a = 298.91*pi/180;        % [rad] right ascension of the ascending node (RAAN)
om_a = 168.65*pi/180;        % [rad] argument of periapsis


%% POINT OF MINIMUM DISTANCE BETWEEN THE ORBITS
% Position along Earth's orbit as a function of true anomaly
rE = @(th_e) par2car(a_e, e_e, i_e, OM_e, om_e, th_e, mu);

% Position along the asteroid's orbit as a function of true anomaly
rA = @(th_a) par2car(a_a, e_a, i_a, OM_a, om_a, th_a, mu);

% Distance between the orbits as a function of the two true anomalies
Dr = @(TH) norm(rE(TH(1)) -  rA(TH(2)));

% Plot the distance function to identify an initial guess for minimization
n = 50;
m = 50;

th_e = linspace(-pi,pi,n);
th_a = linspace(-pi,pi,m);

Z = zeros(n,m);

for k=1:n
    for l=1:m
        Z(k,l) = Dr([th_e(k); th_a(l)]);
    end
end

figure(1)
surf(th_e,th_a, Z');
xlabel("\theta_e [rad]"), ylabel("\theta_a [rad]"), zlabel("||\Delta r|| [km]");
title("Distance Between the Orbits");


TH0 = [0;0];                    % initial guess
THS = fminsearch(Dr, TH0);      % solution

thMD_e = THS(1);         % Earth's true anomaly at minimum distance
thMD_a = THS(2);         % asteroid's true anomaly at minimum distance

% Express both true anomalies in the interval [0, 2*pi)
thMD_e = thMD_e + 2*pi*(thMD_e<0);
thMD_a = thMD_a + 2*pi*(thMD_a<0);

% Earth-orbit position at the point of minimum distance
rMD_e = par2car(a_e, e_e, i_e, OM_e, om_e, thMD_e, mu);
% asteroid-orbit position at the point of minimum distance
rMD_a = par2car(a_a, e_a, i_a, OM_a, om_a, thMD_a, mu);

%% INITIAL ORBIT REPRESENTATION
[xe, ye, ze] = plotOrbitSpace(a_e, e_e, i_e, OM_e, om_e, 0, 2*pi);
[xa, ya, za] = plotOrbitSpace(a_a, e_a, i_a, OM_a, om_a, 0, 2*pi);

figure(2)
plot3(xe,ye,ze, linewidth=2, color="blue");
hold on
plot3(xa,ya,za, linewidth=2,color="red");
plot3(0,0,0,marker='*', MarkerSize=10, LineWidth=2);

plot3(rMD_e(1),rMD_e(2),rMD_e(3),marker='.', MarkerSize=20, LineWidth=2, color="blue");
plot3(rMD_a(1), rMD_a(2), rMD_a(3),marker='.', MarkerSize=20, LineWidth=2, color="red");
hold off

legend("Earth Orbit", "Asteroid Orbit", "Sun", ...
        "Earth on 9 Jan 2092", "Asteroid on 9 Jan 2092")
title("Reference Configuration at Minimum Orbital Distance");
axis equal, axis ([-4e8 4e8 -5e8 3e8 -3e8 3e8])
xlabel("x [km]"), ylabel("y [km]"), zlabel("z [km]"), grid on;


%% ORBITAL PROPAGATION
% --------------------------------------------------------------------- %

% Treat Earth's true anomaly at departure as a prescribed parameter and
% vary it over discrete values to identify the minimum transfer cost.
th1_e = 340*pi/180;  % prescribed Earth departure true anomaly [rad]

                                                                
% Earth's position and velocity at departure
[r1_e, v1_e]  = par2car(a_e, e_e, i_e, OM_e, om_e, th1_e, mu);

% Departure lead time relative to the reference date (9 January 2092)
tAdv = timeOfFlight(a_e,e_e,th1_e,thMD_e,mu);

% Asteroid true anomaly at spacecraft departure from Earth
th1_a = orbitProp(a_a,e_a,tAdv,thMD_a,mu,1);

% Asteroid position at spacecraft departure from Earth
r1_a  = par2car(a_a, e_a, i_a, OM_a, om_a, th1_a, mu);


%% TRANSFER ORBIT CHARACTERIZATION
% --------------------------------------------------------------------- %

% Once Earth's departure position is fixed, the free variables are the
% asteroid's true anomaly at arrival (th2_a) and the transfer-orbit
% argument of periapsis (om).


% Asteroid position and velocity vectors at arrival
r2_a = @(th2_a) par2car(a_a, e_a, i_a, OM_a, om_a, th2_a, mu);
v2_a = @(th2_a) par2car(a_a, e_a, i_a, OM_a, om_a, th2_a, mu, 1);

% Transfer-orbit parameters
a_t     = @(th2_a,om)  direct(r1_e,r2_a(th2_a),om, 'a');   % semi-major axis
e_t     = @(th2_a,om)  direct(r1_e,r2_a(th2_a),om, 'e');   % eccentricity
i_t     = @(th2_a,om)  direct(r1_e,r2_a(th2_a),om, 'i');   % inclination
OM_t    = @(th2_a,om)  direct(r1_e,r2_a(th2_a),om, 'OM');  % RAAN

% spacecraft true anomaly on the transfer orbit at departure
th1_t   = @(th2_a,om)  direct(r1_e,r2_a(th2_a),om, 'th1');

% spacecraft true anomaly on the transfer orbit at arrival
th2_t   = @(th2_a,om)  direct(r1_e,r2_a(th2_a),om, 'th2');
 
%% TOTAL TRANSFER COST
% Spacecraft velocity at departure
v1_t = @(th2_a,om) par2car(a_t(th2_a,om), e_t(th2_a,om), i_t(th2_a,om), ...
                             OM_t(th2_a,om), om, th1_t(th2_a,om), mu, 1);

% Spacecraft velocity at arrival
v2_t = @(th2_a,om) par2car(a_t(th2_a,om), e_t(th2_a,om), i_t(th2_a,om), ...
                             OM_t(th2_a,om), om, th2_t(th2_a,om), mu, 1);

DV1 = @(th2_a,om) norm(v1_t(th2_a,om) - v1_e       );  % departure delta-v
DV2 = @(th2_a,om) norm(v2_t(th2_a,om) - v2_a(th2_a));  % arrival delta-v
DV  =  @(th2_a,om) DV1(th2_a,om) + DV2(th2_a,om);      % total delta-v

%% ARRIVAL-TIME DELAY
% Time required for the asteroid to travel from its position at
% departure to the designated rendezvous point
tA = @(th2_a)     timeOfFlight(a_a, e_a, th1_a, th2_a, mu);

% Time required for the spacecraft to travel from its departure
% position to the designated rendezvous point
tT = @(th2_a,om)  timeOfFlight(a_t(th2_a,om), e_t(th2_a,om), ...
                      th1_t(th2_a,om), th2_t(th2_a,om), mu);

% Spacecraft arrival delay relative to the asteroid
DT = @(th2_a,om)  tT(th2_a,om) - tA(th2_a);


% The following sections implement the procedure described in detail in
% the report.


%% STEPS 1-2: ADMISSIBLE RANGE OF om

th2_a_ = pi/4;              % prescribed trial value of th2_a [rad]

% For a prescribed th2_a, use the eccentricity plot as a function of om
% to determine the interval in which om defines an elliptic transfer orbit.
% The eccentricity must satisfy 0 < e_t < 1.

e_t0 = @(om) e_t(th2_a_, om); 
n = 500;
om = linspace(0, 2*pi, n);
Z = zeros(n,1);
for k=1:n
    Z(k) = e_t0(om(k));
end

figure(3);
plot(om, Z, om, zeros(n,1), om, ones(n,1), linewidth=5);
grid on, xlabel("\omega_t [rad]"),ylabel("e_t [-]"), axis([0,2*pi,-1,2])
title("Admissible Argument-of-Periapsis Interval");


%% STEP 3: INITIAL GUESS FOR om

% The admissible range of om is obtained from the plot defined above
om_min      = 0.7;
om_max      = 3.5;

% Find the value of om that minimizes the cost, using the plot to
% determine the initial guess

DV0 = @(om) DV(th2_a_, om);
n = 200;
om = linspace(om_min, om_max, n);
Z = zeros(n,1);
for k=1:n
    Z(k) = DV0(om(k));
end

figure(4);
plot(om, Z, LineWidth=5);
xlim([om_min;om_max]); grid on,
xlabel("\omega_t [rad]"),ylabel("\Delta v [km/s]")
title("Transfer Cost at the Prescribed Arrival Anomaly");
        

%% STEPS 4-7: SOLUTION OF THE OPTIMIZATION PROBLEM

% Determine the initial guess for om from the preceding plot
om_guess    = 3.2;    

% Minimize the transfer cost with respect to om for each th2_a. The transfer
% time and the resulting delay therefore become functions of th2_a alone.

tT_ = @(th2_a) transferTime(th1_e, th2_a, om_guess);
DT_ = @(th2_a) tT_(th2_a) - tA(th2_a);

Dth = pi/12;
TH = [th2_a_-Dth; th2_a_+Dth];

n = 100;
th2_a = linspace(TH(1),TH(2), n);
Z = zeros(n,1);

for k=1:n
    Z(k) = DT_(th2_a(k));
end

figure(5);
plot(th2_a, Z/3600, th2_a, zeros(size(th2_a)), LineWidth=5);
grid on, xlabel("\theta_{a,2} [rad]"), ylabel("\Delta T [h]"),
xlim([TH(1),TH(2)])
title("Arrival-Time Delay");


%% COMPLETE THE TRANSFER ORBIT CHARACTERIZATION
% ---------------------------------------------------------------------- %

% System solution
th2_a0 = 0.75; % initial guess for th2_a [rad]

% Numerically determine the zero of the delay by minimizing its magnitude
DT_ = @(th2_a) abs(DT_(th2_a));
th2_a_s = fminsearch(DT_, th2_a0);           % th2_a solution

% Minimize the total delta-v at the rendezvous anomaly
DV_s = @(om) DV(th2_a_s, om);
om_s = fminsearch(DV_s,om_guess);            % om solution

% Transfer-orbit parameters
a_t   =  a_t   (th2_a_s,om_s);      % semi-major axis
e_t   =  e_t   (th2_a_s,om_s);      % eccentricity
i_t   =  i_t   (th2_a_s,om_s);      % inclination
OM_t  =  OM_t  (th2_a_s,om_s);      % RAAN

% spacecraft true anomaly on the transfer orbit at departure
th1_t =  th1_t (th2_a_s,om_s);
% spacecraft true anomaly on the transfer orbit at arrival
th2_t =  th2_t (th2_a_s,om_s);
% Spacecraft position on the transfer orbit at departure
r1_t = par2car(a_t,e_t,i_t,OM_t,om_s,th1_t,mu);

% Spacecraft position on the transfer orbit at arrival
r2_t = par2car(a_t,e_t,i_t,OM_t,om_s,th2_t,mu);

DV     =  DV    (th2_a_s,om_s);         % total delta-v [km/s]
tT     =  tT    (th2_a_s,om_s);         % total transfer time [s]

r2_a  =  r2_a  (th2_a_s);               % rendezvous position

om_t = om_s;
th2_a = th2_a_s;

% Earth's true anomaly at spacecraft--asteroid rendezvous
th2_e = orbitProp(a_e,e_e,tT,th1_e,mu,2);
% Earth's position at spacecraft--asteroid rendezvous
r2_e = par2car(a_e,e_e,i_e,OM_e,om_e,th2_e,mu);

%% MANEUVER REPRESENTATION
% --------------------------------------------------------------------- %

[xe, ye, ze] = plotOrbitSpace(a_e, e_e, i_e, OM_e, om_e, 0, 2*pi);
[xa, ya, za] = plotOrbitSpace(a_a, e_a, i_a, OM_a, om_a, 0, 2*pi);
[xt, yt, zt] = plotOrbitSpace(a_t, e_t, i_t, OM_t, om_t, 0, 2*pi);

figure(6)
plot3(xe,ye,ze, linewidth=2, color="blue");
hold on
plot3(xa,ya,za, linewidth=2,color="red");
plot3(0,0,0,marker='*', MarkerSize=10, LineWidth=2);

plot3(r1_e(1),r1_e(2),r1_e(3),marker='.', MarkerSize=20, LineWidth=2, color="blue");
plot3(r2_e(1),r2_e(2),r2_e(3), marker='.', MarkerSize=20, LineWidth=2, color="blue");

plot3(r1_a(1),r1_a(2),r1_a(3),marker='.', MarkerSize=20, LineWidth=2, color="red");
plot3(r2_a(1),r2_a(2),r2_a(3), marker='.', MarkerSize=20, LineWidth=2, color="red");

plot3(xt,yt,zt, linewidth=2, LineStyle=":", color="cyan");
plot3(r1_t(1),r1_t(2),r1_t(3),marker='o', MarkerSize=10, LineWidth=2, color="cyan");
plot3(r2_t(1),r2_t(2),r2_t(3),marker='o', MarkerSize=10, LineWidth=2, color="cyan");
hold off

legend("Earth Orbit", "Asteroid Orbit", "Sun", ...
       "Earth", "", ...
       "", "Asteroid", ...
       "Transfer Orbit", "", ...
       "Spacecraft")
title("Direct Interplanetary Transfer");
axis equal, axis ([-4e8 4e8 -5e8 3e8 -3e8 3e8])
xlabel("x [km]");ylabel("y [km]");zlabel("z [km]");
xlim([-1.8974e+08, 6.5163e+07]),ylim([-3.5877e+07, 2.1903e+08])
zlim([-7.9188e+07, 1.1199e+08]),view([-54 28]), grid on


%% DATA FOR DISCRETE DEPARTURE TRUE ANOMALIES
% ---------------------------------------------------------------------- %
%
%       th1_e = 0     deg       ->       DV = 7.62 km/s
%       th1_e = 355   deg       ->       DV = 7.46 km/s
%       th1_e = 350   deg       ->       DV = 7.33 km/s
%       th1_e = 345   deg       ->       DV = 7.23 km/s
%       th1_e = 342.5 deg       ->       DV = 7.20 km/s
%       th1_e = 340   deg       ->       DV = 7.19 km/s       <- minimum
%       th1_e = 337.5 deg       ->       DV = 7.22 km/s 
%       th1_e = 335   deg       ->       DV = 7.28 km/s 
%       th1_e = 330   deg       ->       DV = 7.49 km/s

A = [330; 335; 337.5; 340; 342.5; 345; 350; 355; 360];
B = [7.49;7.28;7.22;7.19;7.20;7.23;7.33;7.46;7.62];
figure(7);
plot(A,B, 'o', linewidth=10);
grid on, xlabel("\theta_{e,1} [deg]"), ylabel("\Delta v [km/s]")
title("Transfer Cost Versus Earth Departure Anomaly");
axis([327.5,362.5,7.15,7.65])

% =============================== END ===================================
