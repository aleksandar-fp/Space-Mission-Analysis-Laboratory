
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
%                  Project 2: Asteroid Arrival Hyperbola                 %
%                                                                         %
% ======================================================================= %

% PURPOSE
% Characterize the asteroid-centered arrival hyperbola and the circular
% capture orbit. Because the asteroid gravitational parameter is extremely
% small, the relevant quantities are interpreted in the limit mu/rP -> 0.
%
% UNITS
% Distances: [km]
% Velocities and delta-v costs: [km/s]
% Masses: [kg]
% Angles: [rad], unless otherwise specified
%
% NOTATION
% Subscripts h and c denote the arrival hyperbola and circular capture orbit.
% In variable names, I and P denote asymptotic and periapsis conditions.

M = 6.7021e10;          % asteroid mass [kg]
G = 6.6743e-20;         % gravitational constant [km^3/(kg*s^2)]
mu = G*M;               % asteroid gravitational parameter [km^3/s^2]

%% ARRIVAL HYPERBOLA
% --------------------------------------------------------------------- %

vI = 3.226237965872079;     % hyperbolic excess speed relative to the asteroid [km/s]

% Hyperbola semi-major axis: a_h -> 0 as mu -> 0
a_h = -mu/vI^2;

% Prescribed periapsis radius, equal to the capture-orbit radius
rP = 2.5;                   % [km]

% Hyperbola eccentricity: e_h -> infinity as mu/rP -> 0
e_h = 1-rP/a_h;

% Asymptotic true anomaly, turning angle, and impact parameter
thI_h = acos(1/e_h);        % thI_h -> 90 deg as mu/rP -> 0
del = 2*asin(1/e_h);        % turning angle -> 0
DEL = -a_h/tan(del/2);      % impact parameter -> rP [km]
  
%% CAPTURE MANEUVER
% Periapsis speed on the arrival hyperbola and speed on the circular
% capture orbit
vP_h = sqrt(vI^2 + 2*mu/rP);    % hyperbolic periapsis speed [km/s]
vP_c = sqrt(mu/rP);             % circular-orbit speed [km/s]

% Capture delta-v: DV -> vI as mu/rP -> 0
DV = abs(vP_h-vP_c);


%% ORBIT REPRESENTATION
% --------------------------------------------------------------------- %

[xh, yh] = plotOrbitPlane(a_h, e_h,  pi/2, -thI_h*0.9,  thI_h*0.9);
[xc, yc] = plotOrbitPlane(rP, 0, pi/2, 0, 2*pi);

figure(1);
plot(xh,yh, linewidth=2,color="blue");
hold on
plot(xc,yc, linewidth=2,color="red");
plot(0,0, marker="*",linewidth=2,markersize=10);
plot(0, rP, marker="o",linewidth=2, markersize=20, color="cyan");
hold off

grid on,axis equal,axis ([-8,8,-3,3]),xlabel('x [km]'),ylabel('y [km]'),
legend('Arrival Hyperbola', 'Capture Orbit', 'Asteroid', 'Maneuver Point');
title('Asteroid Arrival and Circular Capture');

% ============================== END =====================================
