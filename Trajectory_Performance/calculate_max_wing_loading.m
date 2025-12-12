function WL_max_kg_m2 = calculate_max_wing_loading(S_ref_m2, S_TO_max_ft, T_static_lbf, CL_max)
% calculate_max_wing_loading: Back-calculate maximum wing loading from takeoff constraint
%
% INPUTS:
%   S_ref_m2     - Reference wing area [m^2]
%   S_TO_max_ft  - Maximum allowable takeoff distance [ft]
%   T_static_lbf - Static thrust at sea level [lbf]
%   CL_max       - Maximum lift coefficient during takeoff
%
% OUTPUTS:
%   WL_max_kg_m2 - Maximum allowable wing loading [kg/m^2]
%
% DERIVATION:
%   From: S_TO = 1.44 * W^2 / (g * rho * CL_max * T)
%   With: W = WL * S_ref
%   Then: S_TO = 1.44 * (WL * S_ref)^2 / (g * rho * CL_max * T)
%   Solve for WL:
%   WL^2 = S_TO * g * rho * CL_max * T / (1.44 * S_ref^2)
%   WL = sqrt(S_TO * g * rho * CL_max * T / (1.44 * S_ref^2))
%
% Author: Dominic Larin (Chief Engineer)
% Date: 2025-10-08

    % Constants
    g = 32.174;           % ft/s^2
    rho_SL = 0.002377;    % lbm/ft^3

    % Unit conversions
    m2_to_ft2 = 10.7639;
    lbm_ft2_to_kg_m2 = 4.88243;

    % Convert S_ref to ft^2
    S_ref_ft2 = S_ref_m2 * m2_to_ft2;

    % Calculate maximum wing loading in lbm/ft^2
    WL_max_lbm_ft2 = sqrt(S_TO_max_ft * g * rho_SL * CL_max * T_static_lbf / (1.44 * S_ref_ft2^2));

    % Convert to kg/m^2
    WL_max_kg_m2 = WL_max_lbm_ft2 * lbm_ft2_to_kg_m2;

end
