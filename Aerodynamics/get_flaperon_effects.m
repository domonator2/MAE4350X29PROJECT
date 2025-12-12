function [delta_CL, delta_CD, delta_Cm] = get_flaperon_effects(delta_flap_deg, M, geom)
%GET_FLAPERON_EFFECTS Calculate aerodynamic increments from flaperon deflection
%
% Models the X-29 flaperon (combined flap/aileron) aerodynamic effects using
% Raymer's plain flap methodology with corrections for swept wings.
%
% The X-29 has full-span flaperons on the forward-swept wing trailing edge.
% In flap mode, both surfaces deflect symmetrically (trailing edge down = positive).
%
% Inputs:
%   delta_flap_deg  - Flaperon deflection angle [deg] (positive = TE down)
%                     Typical range: 0° (clean) to 25° (full flaps)
%   M               - Mach number [-]
%   geom            - Geometry structure (optional, uses defaults if empty)
%
% Outputs:
%   delta_CL        - Lift coefficient increment [-]
%   delta_CD        - Drag coefficient increment [-]
%   delta_Cm        - Pitching moment increment [-] (positive = nose up)
%
% References:
%   [1] Raymer, D.P., "Aircraft Design: A Conceptual Approach", Ch. 12
%   [2] NASA TM-4140, "Flight Test Results for the F-29 Forward-Swept Wing"
%   [3] ESDU 74011, "Rate of change of lift coefficient with control deflection"
%
% X-29 Flaperon Specifications (measured):
%   - Area (one side): Sf = 2.49 m²
%   - Span (one side): bf = 2.594 m
%   - Chord ratio: cf/c = 0.25 (25% of local wing chord)
%   - Max deflection: 24.75° TED
%   - Derived: Sf_total/Sref ≈ 0.284, bf_total/b ≈ 0.69
%
% Author: Dominic Larin
% Date: December 2025

%% X-29 Flaperon Geometry (measured data)
% Source: Project team measurements
Sf_one_side = 2.49;     % Flaperon area, one side [m²]
bf_one_side = 2.594;    % Flaperon span, one side [m]
cf_c = 0.25;            % Chord ratio (flaperon chord / local wing chord)
delta_max = 24.75;      % Maximum deflection [deg] TED

% Total flaperon (both sides)
Sf_total = 2 * Sf_one_side;  % = 4.98 m²
bf_total = 2 * bf_one_side;  % = 5.188 m

% Reference geometry
if nargin >= 3 && ~isempty(geom) && isfield(geom, 'wing')
    Sref = geom.wing.S;
    b_wing = sqrt(geom.wing.AR * Sref);  % Wing span from AR
    Lambda_c4 = geom.wing.Lambda_c25_deg;
else
    Sref = 17.54;       % Wing reference area [m²]
    b_wing = 7.5;       % Wing span [m] (approx from AR=4)
    Lambda_c4 = -30;    % Quarter-chord sweep [deg]
end

% Derived ratios
Sf_Sref = Sf_total / Sref;   % ≈ 0.284 (28.4% of wing area)
bf_b = bf_total / b_wing;    % ≈ 0.69 (69% span coverage)

%% Clamp deflection to physical limits
delta_flap_deg = max(0, min(delta_flap_deg, delta_max));

%% =========================================================================
%  LIFT INCREMENT (Raymer Eq. 12.22 with sweep correction)
%  =========================================================================

% Base flap effectiveness (from ESDU 74011 for plain flaps)
% CL_delta_flap = 2*pi * tau_f where tau_f = f(cf/c)
% For cf/c = 0.25: tau_f ≈ 0.58
% CL_delta_2D ≈ 2*pi * 0.58 = 3.65 per radian ≈ 0.064 per degree

% Empirical fit for tau_f vs cf/c (from Raymer Fig. 12.23)
tau_f = 0.9 * sqrt(cf_c);  % Simplified approximation

% 2D lift effectiveness [per deg]
CL_delta_2D = (2 * pi * tau_f) * (pi/180);  % Convert to per degree

% Sweep correction (Raymer Eq. 12.25)
% CLalpha_swept / CLalpha_2D = cos(Lambda)
cos_sweep = cosd(abs(Lambda_c4));

% Compressibility correction (Prandtl-Glauert)
if M < 0.85
    beta_PG = sqrt(1 - M^2);
else
    beta_PG = sqrt(abs(M^2 - 1));  % Approximate for transonic
end
compress_factor = 1 / beta_PG;
compress_factor = min(compress_factor, 1.5);  % Limit for transonic

% 3D lift effectiveness with corrections
CL_delta_3D = CL_delta_2D * cos_sweep * compress_factor * Sf_Sref / cf_c;

% Account for partial span - use simple spanwise integration factor
spanwise_factor = bf_b;  % Proportional to span coverage

% Final lift increment
delta_CL = CL_delta_3D * spanwise_factor * delta_flap_deg;

% Empirical calibration to match X-29 CL_max = 1.6
% Clean CL_max ≈ 0.8-0.9 (from aero model at alpha = 12°)
% Need delta_CL ≈ 0.7-0.8 at delta = 20-25°
% Calibration factor
cal_CL = 1.15;
delta_CL = delta_CL * cal_CL;

%% =========================================================================
%  DRAG INCREMENT (Raymer Eq. 12.37)
%  =========================================================================

% Flap drag has three components:
% 1. Profile drag increase (boundary layer/pressure drag)
% 2. Induced drag from additional lift
% 3. Interference drag

% Profile drag increment (empirical, Raymer Fig. 12.36)
% CD_flap_profile ≈ 0.0023 * (delta/10)^1.5 for plain flaps
CD_profile_base = 0.0023;
delta_CD_profile = CD_profile_base * (delta_flap_deg / 10)^1.5;

% Induced drag increment
% delta_CDi = delta_CL^2 / (pi * AR * e)
AR = 4.0;  % X-29 aspect ratio
e_flap = 0.70;  % Reduced Oswald efficiency with flaps
delta_CDi = delta_CL^2 / (pi * AR * e_flap);

% Interference drag (flap-wing junction)
CD_interference = 0.0005 * (delta_flap_deg / delta_max);

% Total drag increment
delta_CD = delta_CD_profile + delta_CDi + CD_interference;

% Calibration to match landing config delta_CD ≈ 0.010 at delta = 20°
% Scale to expected value
if delta_flap_deg > 0
    expected_delta_CD_20deg = 0.010;
    actual_delta_CD_20deg = CD_profile_base * (20/10)^1.5 + ...
                           (CL_delta_3D * spanwise_factor * 20 * cal_CL)^2 / (pi * AR * e_flap) + ...
                           0.0005 * (20/delta_max);
    cal_CD = expected_delta_CD_20deg / actual_delta_CD_20deg;
    cal_CD = max(0.5, min(cal_CD, 2.0));  % Limit calibration range
else
    cal_CD = 1.0;
end
delta_CD = delta_CD * cal_CD;

%% =========================================================================
%  PITCHING MOMENT INCREMENT
%  =========================================================================

% Trailing edge down flap creates nose-down pitching moment
% Cm_delta_flap ≈ -0.25 * CL_delta (typical for c_flap/c = 0.25)

% Moment arm: flap is aft of CG, so TE-down creates nose-down moment
moment_factor = -0.30;  % Empirical (depends on flap location relative to CG)

delta_Cm = moment_factor * delta_CL;

%% =========================================================================
%  SPECIAL CASE: TAKEOFF vs LANDING
%  =========================================================================

% Note: This function returns the raw aerodynamic increments.
% The calling function should apply these based on flight phase:
%
% TAKEOFF:
%   - Partial flap (10-15°) for increased CL without excessive drag
%   - Flaperon provides CL_max ≈ 1.3-1.4
%
% LANDING:
%   - Full flap (20-25°) for maximum CL and drag
%   - Flaperon provides CL_max ≈ 1.6
%   - Higher drag is beneficial for approach angle control

end
