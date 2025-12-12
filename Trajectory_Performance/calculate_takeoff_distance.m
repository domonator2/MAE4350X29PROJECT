function S_TO = calculate_takeoff_distance(TOGW_lbm, T_static_lbf, CL_max, S_ref_ft2)
% calculate_takeoff_distance: Compute takeoff distance for screening
%
% INPUTS:
%   TOGW_lbm     - Takeoff gross weight [lbm]
%   T_static_lbf - Static thrust at sea level with max afterburner [lbf]
%   CL_max       - Maximum lift coefficient during takeoff (with flaps)
%   S_ref_ft2    - Reference wing area [ft^2] (default: 185 ft^2 for X-29)
%
% OUTPUTS:
%   S_TO         - Takeoff distance [ft]
%
% EQUATION:
%   S_TO = 1.44 * W^2 / (g * rho * S_ref * CL_max * T)
%
% CONSTANTS (Sea Level, Standard Atmosphere):
%   g = 32.174 ft/s^2
%   rho = 0.002377 lbm/ft^3
%   S_ref = 185 ft^2 (X-29 wing area)
%
% SCREENING CRITERIA:
%   Edwards AFB runway length = 2 miles = 10,560 ft
%   Vehicles with S_TO > 10,560 ft are unfeasible
%
% Author: Dominic Larin (Chief Engineer)
% Date: 2025-10-08

    % Constants - Sea Level Standard Atmosphere
    g = 32.174;           % ft/s^2
    rho_SL = 0.002377;    % lbm/ft^3

    % Default S_ref to X-29 value if not provided
    if nargin < 4 || isempty(S_ref_ft2)
        S_ref_ft2 = 185;  % X-29 wing area
    end

    % Calculate takeoff distance
    S_TO = 1.44 * TOGW_lbm^2 / (g * rho_SL * S_ref_ft2 * CL_max * T_static_lbf);

end
