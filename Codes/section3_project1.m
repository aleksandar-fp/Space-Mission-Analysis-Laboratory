
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
%                   Section 3: Hyperbolic Trajectories                    %
%                                                                         %
%                  Project 1: Earth Departure Hyperbola                   %
%                                                                         %
% ======================================================================= %

% PURPOSE
% Characterize the geocentric departure hyperbola connecting the parking
% orbit defined in Section 1 to the heliocentric transfer trajectory defined
% in Section 2. The maneuver point on the parking orbit is varied over
% discrete values to minimize the required injection delta-v.
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
% Subscripts p and h denote the parking orbit and departure hyperbola.
% In variable names, I denotes the asymptotic condition.

mu = 398600.4;           % Earth's gravitational parameter [km^3/s^2]

%% PARKING ORBIT DATA
% Position and velocity at the reference time
rS = [-5985.618600; -7930.756300; 6091.925600];
vS = [4.532000; -5.596000; -2.832000];

% Parking-orbit parameters and true anomaly at the reference time
[a_p, e_p, i_p, OM_p, om_p, thS_p] = car2par(rS, vS, mu);

%% ASYMPTOTIC VELOCITY
% Hyperbolic excess velocity in the heliocentric ecliptic reference frame
vI_elio = [ - 3.635205978285434  ;
              0.738135418064654  ;
            - 1.405891388922091 ];

% Rotate through the obliquity of the ecliptic to express the velocity in
% the geocentric equatorial reference frame
vI = rotx(23.45)*vI_elio;   % hyperbolic excess velocity vector [km/s]
rI = vI/norm(vI);           % unit vector of the asymptotic direction [-]


%% DEPARTURE HYPERBOLA
% ----------------------------------------------------------------------- %

a_h = -mu/norm(vI)^2;       % hyperbola semi-major axis [km]

% Prescribe the true anomaly of the injection point on the parking orbit.
% This parameter is varied over discrete values to minimize the maneuver cost.
thH_p = pi/180*177.5;  % parking-orbit true anomaly at injection [rad]

% Corresponding position and velocity on the parking orbit
[rH, vH_p] = par2car(a_p, e_p, i_p, OM_p, om_p, thH_p, mu);


%% DEFINITION OF THE NONLINEAR SYSTEM

% Define F, whose root provides the hyperbolic eccentricity e_h, the true
% anomaly thH_h at injection, and the asymptotic true anomaly thI.

rH_mod = norm(rH);                         % injection radius [km]
alfa = acos(dot(rH/norm(rH),rI));          % angle between rH and rI [rad]

% Unknown vector: X = [e_h, thH_h, thI]
f1 = @(X) rH_mod/a_h * ( 1 + X(1)*cos(X(2)) ) - (1 - X(1)^2);
f2 = @(X) cos(X(3)) + 1/X(1);
f3 = @(X) X(3) - alfa - X(2);
F = @(X) [f1(X); f2(X); f3(X)];


% Set e_h = e_h0 to plot g(thH_h, thI) = F(e_h0, thH_h, thI) and determine
% a suitable initial guess for solving the nonlinear system.

% Plotting variables: Y = [thH_h, thI]
e_h_ = 4;
g = @(Y) norm(F([e_h_; Y(1); Y(2)]));

n = 50;
m = 50;

thH_h = linspace(0,2*pi,n);
thI   = linspace(pi/2,pi,m);

Z = zeros(n,m);

for k = 1:n
    for l = 1:m
        Z(k,l) = g([thH_h(k); thI(l)]);
    end
end

figure(1)
surf(thH_h,thI,Z');
xlabel('\theta_{h,H} [rad]'), ylabel('\theta_{\infty} [rad]'), zlabel('g [-]');
title('Initial-Guess Search for the Nonlinear System');

%% SOLUTION OF THE NONLINEAR SYSTEM
% Initial guess
thH_h0 = 1.0;
thI0 = 2.5;             

X0 = [e_h_, thH_h0, thI0];
XS = fsolve(F, X0);

% Solution
e_h    = XS(1);
thH_h  = XS(2);
thI    = XS(3);

% Express the angular solutions in the interval [0, 2*pi)
thH_h = thH_h + 2*pi*(thH_h<0);
thI = thI + 2*pi*(thI<0);

%% REMAINING HYPERBOLA PARAMETERS
% Unit specific angular momentum vector and inclination
h_h = cross(rH/norm(rH), rI) / norm(cross(rH/norm(rH), rI));
i_h = acos(h_h(3));

% Unit vector along the line of nodes and RAAN
N_h = cross([0;0;1],h_h)/norm(cross([0;0;1],h_h));
OM_h = (N_h(2)<0)*2*pi + (-1)^(N_h(2)<0) * acos(N_h(1));

% Argument of periapsis obtained from the argument of latitude at injection
sinB = dot(cross(N_h, rH/norm(rH)), h_h);
cosB = dot(N_h, rH/norm(rH));
om_h = atan2(sinB,cosB) - thH_h;
om_h = om_h + 2*pi*(om_h<0);

% Turning angle, impact parameter, and periapsis radius
del = 2*asin(1/e_h);                     % turning angle [rad]
DEL = -a_h*e_h*cos(del/2);               % impact parameter [km]
rP = mu/norm(vI)^2 * (e_h-1);            % periapsis radius [km]

% Injection delta-v [km/s]
[~, vH_h] = par2car(a_h, e_h, i_h, OM_h, om_h, thH_h, mu);
DV = norm(vH_h - vH_p);


%% ORBIT REPRESENTATION
% --------------------------------------------------------------------- %
[xp, yp, zp] = plotOrbitSpace(a_p, e_p, i_p, OM_p, om_p, 0, 2*pi);
[xh, yh, zh] = plotOrbitSpace(a_h, e_h, i_h, OM_h, om_h, -thI*0.4, thI*0.7);

figure(2);
plot3(xp,yp,zp, linewidth=2,color="blue");
hold on
plot3(xh,yh,zh, linewidth=2,color="red");
plot3(0,0,0, marker='*',MarkerSize=10,LineWidth=2);
% Scale the unit asymptotic direction for visualization only
quiver3(0,0,0,rI(1)*1e5,rI(2)*1e5,rI(3)*1e5, LineWidth=1.5,color="cyan");
plot3(rH(1),rH(2),rH(3), marker="o",markersize=10,LineWidth=2,color="black");
hold off

legend("Parking Orbit", "Departure Hyperbola", ...
       "Earth", "Asymptotic Departure Direction", "Maneuver Point");
title("Earth Departure Hyperbola");
axis equal, grid on,view([113 21]);
xlabel("x [km]"),ylabel("y [km]"),zlabel("z [km]");


%% DATA FOR DISCRETE DEPARTURE TRUE ANOMALIES
% --------------------------------------------------------------------- %
%
% thH_p = 170    deg       ->            DV = 4.456 km/s
% thH_p = 172.5  deg       ->            DV = 4.415 km/s
% thH_p = 175    deg       ->            DV = 4.388 km/s
% thH_p = 176.25 deg       ->            DV = 4.381 km/s
% thH_p = 177.5  deg       ->            DV = 4.377 km/s  <- minimum
% thH_p = 178.75 deg       ->            DV = 4.378 km/s
% thH_p = 180    deg       ->            DV = 4.382 km/s
% thH_p = 182.5  deg       ->            DV = 4.404 km/s
% thH_p = 185    deg       ->            DV = 4.442 km/s
% thH_p = 187.5  deg       ->            DV = 4.496 km/s
% thH_p = 190    deg       ->            DV = 4.567 km/s

A = [  170, 172.5,   175, 176.25, 177.5, 178.75, 180,   182.5, 185,   187.5, 190  ]';
B = [4.456, 4.415, 4.388,  4.381, 4.377, 4.378,  4.382, 4.404, 4.442, 4.496, 4.567]';

figure(3);
plot(A,B,'o', linewidth=10);
grid on, axis([167.5,192.5,4.35,4.6])
xlabel("\theta_{p,H} [deg]"), ylabel("\Delta v [km/s]"),
title("Injection Cost Versus Parking-Orbit Departure Anomaly");

% ============================== END ====================================
