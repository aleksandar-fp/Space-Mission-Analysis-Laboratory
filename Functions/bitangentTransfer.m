function [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)
%BITANGENTTRANSFER Compute a bitangent transfer between two elliptic orbits.
%
% [DeltaV1, DeltaV2, Deltat] = bitangentTransfer(a_i, e_i, a_f, e_f, type, mu)
%
% Computes the two tangential velocity increments and the transfer time
% for coplanar, apsidally aligned elliptic orbits.
%
% Inputs:
%   a_i, e_i - Initial semi-major axis [km] and eccentricity [-].
%   a_f, e_f - Final semi-major axis [km] and eccentricity [-].
%   type     - Transfer type: 'pa', 'ap', 'pp', or 'aa', where 'p' and
%              'a' denote periapsis and apoapsis.
%   mu       - Gravitational parameter [km^3/s^2].
%
% Outputs:
%   DeltaV1  - First tangential velocity increment [km/s].
%   DeltaV2  - Second tangential velocity increment [km/s].
%   Deltat   - Time along the transfer arc [s].

% Initial- and final-orbit apsidal radii
Rp_i = a_i * (1 - e_i);      % initial periapsis radius
Ra_i = a_i * (1 + e_i);      % initial apoapsis radius
Rp_f = a_f * (1 - e_f);      % final periapsis radius
Ra_f = a_f * (1 + e_f);      % final apoapsis radius


switch type

    % Periapsis to apoapsis
    case "pa"
        % Transfer-orbit geometry
        Rp_T = Rp_i;
        Ra_T = Ra_f;
        a_T = (Rp_T + Ra_T)/2;
        
        % Tangential velocity increments
        DeltaV1 =   sqrt(mu) * sqrt(2/Rp_T - 1/a_T) + ...
                  - sqrt(mu) * sqrt(2/Rp_i - 1/a_i);
        DeltaV2 =   sqrt(mu) * sqrt(2/Ra_f - 1/a_f) + ...
                  - sqrt(mu) * sqrt(2/Ra_T - 1/a_T);


    % Apoapsis to periapsis
    case "ap"
        % Transfer-orbit geometry
        Rp_T = Rp_f;
        Ra_T = Ra_i;
        a_T = (Rp_T + Ra_T)/2;

        % Tangential velocity increments
        DeltaV1 =   sqrt(mu) * sqrt(2/Ra_T - 1/a_T) + ...
                  - sqrt(mu) * sqrt(2/Ra_i - 1/a_i);
        DeltaV2 =   sqrt(mu) * sqrt(2/Rp_f - 1/a_f) + ...
                  - sqrt(mu) * sqrt(2/Rp_T - 1/a_T);

    % Periapsis to periapsis
     case "pp"
        % Transfer-orbit geometry
        Rp_T = Rp_i;
        Ra_T = Rp_f;
        a_T = (Rp_T + Ra_T)/2;
        
        % Tangential velocity increments
        DeltaV1 =   sqrt(mu) * sqrt(2/Rp_T - 1/a_T) + ...
                  - sqrt(mu) * sqrt(2/Rp_i - 1/a_i);
        DeltaV2 =   sqrt(mu) * sqrt(2/Rp_f - 1/a_f) + ...
                  - sqrt(mu) * sqrt(2/Ra_T - 1/a_T);

    % Apoapsis to apoapsis
     case "aa"
        % Transfer-orbit geometry
        Rp_T = Ra_i;
        Ra_T = Ra_f;
        a_T = (Rp_T + Ra_T)/2;
        
        % Tangential velocity increments
        DeltaV1 =   sqrt(mu) * sqrt(2/Rp_T - 1/a_T) + ...
                  - sqrt(mu) * sqrt(2/Ra_i - 1/a_i);
        DeltaV2 =   sqrt(mu) * sqrt(2/Ra_f - 1/a_f) + ...
                  - sqrt(mu) * sqrt(2/Ra_T - 1/a_T);
end

% Each transfer arc spans half of the transfer-orbit period
Deltat = pi * sqrt(a_T^3 / mu);

end
