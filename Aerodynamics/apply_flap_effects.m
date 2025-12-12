function aero_flapped = apply_flap_effects(aero_clean, delta_flap_deg, geom, delta_strake_deg)
% APPLY_FLAP_EFFECTS Add flaperon and strake effects to clean aerodynamics
%
% This function adds incremental lift and drag effects due to flaperon and
% strake deflection to the clean configuration aerodynamics.
%
% Inputs:
%   aero_clean      - Clean configuration aero structure from raymer_aero_model_elliptical
%   delta_flap_deg  - Flaperon deflection [deg], positive = trailing edge down
%   geom            - Geometry structure (must contain flaperon data)
%   delta_strake_deg - Strake deflection [deg] (optional, default = 0)
%
% Outputs:
%   aero_flapped    - Modified aero structure with flap/strake effects
%
% X-29 Flaperon Data:
%   Sf = 2.4903 m² (flaperon area)
%   Sf/S = 0.1427 (14.27% of wing area)
%   δf_max = 24.75° (max deflection, trailing edge down)
%
% Typical Deflections:
%   Takeoff:  Flaperon 15°, Strake 6°
%   Landing:  Flaperon 20°, Strake 10°
%
% Reference: Raymer Ch. 12, Roskam Part VI

% =========================================================================
% DEFAULT STRAKE DEFLECTION
% =========================================================================
if nargin < 4 || isempty(delta_strake_deg)
    delta_strake_deg = 0;
end

% =========================================================================
% X-29 CONTROL SURFACE LIFT EFFECTIVENESS (FROM ANALYSIS)
% =========================================================================
% These values are derived from Raymer's method for control surface
% lift effectiveness on the X-29 configuration
%
% CL_delta_flaperon = 0.0132 /deg (trailing edge flaperon)
% CL_delta_strake   = 0.0066 /deg (forward strake/canard)

CL_delta_flaperon = 0.0132;  % [1/deg]
CL_delta_strake = 0.0066;    % [1/deg]

% =========================================================================
% LIFT COEFFICIENT INCREMENT
% =========================================================================
% ΔCL = CL_delta_flaperon * δ_flaperon + CL_delta_strake * δ_strake

delta_CL_flaperon = CL_delta_flaperon * delta_flap_deg;
delta_CL_strake = CL_delta_strake * delta_strake_deg;
delta_CL = delta_CL_flaperon + delta_CL_strake;

% =========================================================================
% DRAG COEFFICIENT INCREMENT
% =========================================================================
% Flap deflection increases parasite drag (pressure drag, interference)
% Typical values:
%   Takeoff flaps (15°): ΔCD0 ≈ 0.010-0.015
%   Landing flaps (20°): ΔCD0 ≈ 0.015-0.025
%
% Use quadratic model to capture increasing drag at higher deflections:
%   ΔCD0 = K1 * δf + K2 * δf² (for flaperon)
%   ΔCD0_strake = K1_s * δs + K2_s * δs² (for strake)

K1_CD_flap = 0.0005;    % Linear term [1/deg] for flaperon
K2_CD_flap = 0.00002;   % Quadratic term [1/deg²] for flaperon
K1_CD_strake = 0.0003;  % Linear term [1/deg] for strake (smaller surface)
K2_CD_strake = 0.00001; % Quadratic term [1/deg²] for strake

delta_CD0_flap = K1_CD_flap * delta_flap_deg + K2_CD_flap * delta_flap_deg^2;
delta_CD0_strake = K1_CD_strake * delta_strake_deg + K2_CD_strake * delta_strake_deg^2;
delta_CD0 = delta_CD0_flap + delta_CD0_strake;

% =========================================================================
% APPLY INCREMENTS TO CLEAN AERO
% =========================================================================
aero_flapped = aero_clean;

% Update CL
aero_flapped.CL = aero_clean.CL + delta_CL;

% Update CD0 and total CD
aero_flapped.CD0 = aero_clean.CD0 + delta_CD0;
aero_flapped.CDi = aero_clean.K * aero_flapped.CL^2;  % Induced drag with new CL
aero_flapped.CD = aero_flapped.CD0 + aero_flapped.CDi;

% Note: Forces (L, D) are not recalculated here since they depend on
% flight conditions. The calling function should recalculate forces from
% the updated coefficients (CL, CD) if needed.

% =========================================================================
% VALIDATION
% =========================================================================
% Check for reasonable values
if delta_flap_deg > 30
    warning('Flap deflection %.1f° exceeds typical range. X-29 max = %.1f°', ...
            delta_flap_deg, 24.75);
end

if aero_flapped.CL > 2.5
    warning('CL = %.2f may be unrealistic (post-stall)', aero_flapped.CL);
end

end
