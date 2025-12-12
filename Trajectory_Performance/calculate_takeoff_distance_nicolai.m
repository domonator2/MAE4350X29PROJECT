function S_TO = calculate_takeoff_distance_nicolai(TOGW_lbf, S_ref_ft2, T_static_lbf, CL_max)
% calculate_takeoff_distance_nicolai: Compute takeoff distance using Nicolai's equation
%
% INPUTS:
%   TOGW_lbf     - Takeoff gross weight [lbf]
%   S_ref_ft2    - Reference wing area [ft^2]
%   T_static_lbf - Static thrust at sea level with max afterburner [lbf]
%   CL_max       - Maximum lift coefficient during takeoff (with flaps)
%
% OUTPUTS:
%   S_TO         - Takeoff distance [ft]
%
% EQUATION (Nicolai):
%   S_TO = 20.9 * (W/S) / (CL_max * T/W) + 69.6 * sqrt((W/S) / CL_max)
%
% This is an empirical correlation for jet aircraft takeoff distance.
%
% REFERENCE:
%   Nicolai, L. M., & Carichner, G. E. (2010).
%   Fundamentals of Aircraft and Airship Design, Volume 1: Aircraft Design.
%   AIAA Education Series.
%
% Author: Dominic Larin (Chief Engineer)
% Date: 2025-10-08

    % Calculate wing loading [lbf/ft^2]
    WS = TOGW_lbf / S_ref_ft2;

    % Calculate thrust-to-weight ratio
    TW = T_static_lbf / TOGW_lbf;

    % Nicolai's takeoff distance equation
    % Term 1: Acceleration distance (depends on T/W)
    term1 = 20.9 * WS / (CL_max * TW);

    % Term 2: Rotation and liftoff distance (depends on stall speed)
    term2 = 69.6 * sqrt(WS / CL_max);

    % Total takeoff distance
    S_TO = term1 + term2;

    % Sanity check
    if S_TO < 0 || ~isfinite(S_TO)
        S_TO = Inf;  % Invalid configuration
    end

end
