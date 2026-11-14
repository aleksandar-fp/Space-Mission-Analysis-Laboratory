function OUT = direct(r1, r2, om, type)
%DIRECT Determine an orbit through two prescribed position vectors.
%
% OUT = direct(r1, r2, om, type)
%
% The two positions and the prescribed argument of periapsis define the
% orbital plane, shape, and true anomalies at the two points.
%
% Inputs:
%   r1, r2 - Position vectors [km].
%   om     - Prescribed argument of periapsis [rad].
%   type   - Requested output: 'a', 'e', 'i', 'OM', 'th1', or 'th2'.
%
% Output:
%   OUT    - Requested orbital parameter. The semi-major axis is returned
%            in km, eccentricity is dimensionless, and angles are in rad.
%
% The formulation assumes non-collinear position vectors and a
% non-equatorial orbital plane.

I = [1;0;0];
J = [0;1;0];
K = [0;0;1];

% Orbital-plane normal and inclination
h    = cross(r1,r2)/norm(cross(r1,r2));
i         = acos(h(3));

% Line of nodes and RAAN
N         = cross(K,h)/norm(cross(K,h));
OM        = (N(2)<0)* 2*pi + (-1)^(N(2)<0)  * acos(N(1));

% Transform the positions to the perifocal frame
R_OM    =    rotz(OM *180/pi); 
R_i     =    rotx(i  *180/pi);
R_om    =    rotz(om *180/pi);
R       =    R_OM * R_i * R_om;

r1_pf   =    R'*r1;
r2_pf   =    R'*r2;

% True anomalies at the prescribed positions
sin1    =    dot(r1_pf,J)/norm(r1_pf);
cos1    =    dot(r1_pf,I)/norm(r1_pf);
sin2    =    dot(r2_pf,J)/norm(r2_pf);
cos2    =    dot(r2_pf,I)/norm(r2_pf);
th1     =    atan2(sin1,cos1);
th2     =    atan2(sin2,cos2);

% Eccentricity and semi-major axis
R1  =  norm(r1);
R2  =  norm(r2);

e =   (R1 - R2) / (R2*cos2 - R1*cos1);
a =    R1 * (1 + e*cos1) / (1 - e^2);

% Return the requested orbital parameter
switch type
    case "a"
        OUT = a;
    case "e"
        OUT = e;
    case "i"
        OUT = i;
    case "OM"
        OUT = OM;
    case "th1"
        OUT = th1;
    case "th2"
        OUT = th2;
end

end
