function S_TO = calculate_takeoff_distance_full(TOGW_lbm, T_static_lbf, CL_max, S_ref_ft2, CD0, mu_r)
% calculate_takeoff_distance_full: Compute takeoff distance with drag and friction
%
% INPUTS:
%   TOGW_lbm     - Takeoff gross weight [lbm]
%   T_static_lbf - Static thrust at sea level with max afterburner [lbf]
%   CL_max       - Maximum lift coefficient during takeoff (with flaps)
%   S_ref_ft2    - Reference wing area [ft^2]
%   CD0          - Zero-lift drag coefficient (optional, default: 0.025)
%   mu_r         - Rolling friction coefficient (optional, default: 0.025 for concrete)
%
% OUTPUTS:
%   S_TO         - Takeoff distance [ft]
%
% EQUATION (Embry-Riddle, Eq. 23):
%   S_TO = 1.44 * W^2 / (g * rho * S * CL_max * (T - D_avg - mu_r*(W - L_avg)))
%
% Where drag and lift are evaluated at 0.7 * V_LO:
%   V_LO = 1.2 * V_stall = 1.2 * sqrt(2*W/(rho*S*CL_max))
%
% REFERENCE:
%   https://www.aerodynamics4students.com/propulsion/takeoff-and-landing-performance.php
%   Embry-Riddle Aeronautical University
%
% Author: Dominic Larin (Chief Engineer)
% Date: 2025-10-08

    % Constants - Sea Level Standard Atmosphere
    g = 32.174;                      % ft/s^2
    rho_SL_slug = 0.002377 / 32.174; % slugs/ft^3 (convert from lbm/ft^3)

    % Default parameters
    if nargin < 5 || isempty(CD0)
        CD0 = 0.025;      % Typical for clean fighter configuration
    end
    if nargin < 6 || isempty(mu_r)
        mu_r = 0.025;     % Concrete runway
    end

    % Convert TOGW to weight in lbf (numerically equal to lbm at Earth surface)
    W_lbf = TOGW_lbm;  % Weight in lbf
    W_slug = TOGW_lbm / 32.174;  % Mass in slugs

    % Calculate liftoff speed (Eq. 25)
    % V_LO = 1.2 * sqrt(2*W/(rho*S*CL_max))
    V_stall = sqrt(2 * W_lbf / (rho_SL_slug * g * S_ref_ft2 * CL_max));  % ft/s
    V_LO = 1.2 * V_stall;  % ft/s

    % Evaluate at 0.7 * V_LO (Eq. 24)
    V_eval = 0.7 * V_LO;  % ft/s
    q_eval = 0.5 * rho_SL_slug * V_eval^2;  % Dynamic pressure [lbf/ft^2] (slug-based)

    % Calculate lift at evaluation point
    % At 0.7*V_LO, velocity is 70% of liftoff, so dynamic pressure is 0.49x
    % To maintain same lift coefficient relationship:
    CL_eval = CL_max * 0.49;  % Approximate CL during ground roll
    L_avg = q_eval * S_ref_ft2 * CL_eval;  % lbf

    % Calculate drag at evaluation point
    % Ground roll: induced drag is small, mostly zero-lift drag
    K = 1 / (pi * 7.5 * 0.85);  % Assume AR~7.5, e~0.85 for fighter
    CD_eval = CD0 + K * CL_eval^2;
    D_avg = q_eval * S_ref_ft2 * CD_eval;  % lbf

    % Calculate rolling friction
    FR_avg = mu_r * (W_lbf - L_avg);  % lbf

    % Net accelerating force
    F_net = T_static_lbf - D_avg - FR_avg;  % lbf

    % Check for negative net force (can't take off!)
    if F_net <= 0
        S_TO = Inf;  % Cannot take off
        return;
    end

    % Calculate takeoff distance (Eq. 23)
    % S_TO = 1.44 * W^2 / (g * rho * S * CL_max * F_net)
    S_TO = 1.44 * W_lbf^2 / (g * rho_SL_slug * g * S_ref_ft2 * CL_max * F_net);

end
