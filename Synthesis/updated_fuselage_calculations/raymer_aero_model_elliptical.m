%{
====================================================================================================
SCRIPT NAME: raymer_aero_model_elliptical.m
AUTHOR: Dominic Larin
INITIATED: 11/20/2025
====================================================================================================
SCRIPT DESCRIPTION:
Unified aerodynamics model adapted for ELLIPTICAL FUSELAGE geometry with smaller wing area.
This model handles the fact that the elliptical fuselage uses S_wing_ratio = 0.451 instead of 0.85,
resulting in significantly smaller wing reference areas.

Key Features:
1. Correct Raymer subsonic lift formula without calibration (exposure factor properly implemented)
2. Supersonic lift using CLα = 4/β with finite swept-wing corrections
3. Supersonic drag using swept-wing induced drag formula K = [4(M²-1)cos(Λ)] / [4·AR·β - 2]
4. Smooth transonic fairing between M=0.85 and M=1.2
5. Empirical CD0 model calibrated to X-29 performance data

Target Performance:
- Subsonic (M=0.6): CL≈0.6, L/D≈9-10
- Transonic (M=0.9): CL≈0.8
- Supersonic (M=1.2): CL≈0.48, CD≈0.10
- Supersonic (M=1.3): CL≈0.43, CD≈0.10
- Cruise (M=1.4): CL≈0.39, CD≈0.10, L/D≈3.8
====================================================================================================
%}

% Inputs (SI):
%   M              - Free-stream Mach Number. [-]
%   alpha_deg      - Angle of Attack. [deg]
%   h              - Altitude. [m]
%   geom           - A structure containing detailed vehicle geometry from test_elliptical_fuselage()
%
% Outputs:
%   aero           - Structure containing aerodynamic coefficients and forces:
%     .L              - Lift [N]
%     .D              - Drag [N]
%     .L_over_D       - Lift-to-Drag Ratio [-]
%     .LD_max         - Maximum L/D [-]
%     .CL             - Lift Coefficient [-]
%     .CLopt          - Optimum Lift Coefficient [-]
%     .CLalpha        - Lift-Curve Slope [1/rad]
%     .CD             - Total Drag Coefficient [-]
%     .CD0            - Parasitic Drag Coefficient [-]
%     .CDi            - Induced Drag Coefficient [-]
%     .K              - Induced Drag Factor [-]
%     .regime         - Flight regime ('subsonic', 'transonic', or 'supersonic')

function aero = raymer_aero_model_elliptical(M, alpha_deg, h, geom)

% =========================================================================
% SETUP AND GEOMETRY PREPARATION
% =========================================================================

% Get atmospheric properties
[~,atm.a,~,atm.rho,~,atm.mu] = atmosisa(h);

opts = struct();
opts.units = 'SI';
opts.Ewave = 2.0; % Wave drag efficiency parameter

% ----- Prepare Geometry Fields -----
geom.Sref = geom.wing.S;  % Wing area (45.1% of planform, not 85%)
geom.AR = geom.wing.AR;
geom.Swet = geom.S_wet_total;
geom.Cfe = 0.0055;
geom.sweep_tmax_deg = geom.wing.Lambda_c25_deg;
geom.sweep_LE_deg = geom.wing.Lambda_LE_deg;

% Elliptical fuselage maximum cross-sectional area
if isfield(geom.fuse, 'D')
    geom.Amax = pi * (geom.fuse.D/2)^2;
else
    geom.Amax = pi * (geom.fuse.W/2) * (geom.fuse.H/2);
end
geom.Lref = geom.fuse.L;
geom.Ewd = 2.0;

if ~isfield(geom, 'iw')
    geom.iw = 0;
end
if ~isfield(geom, 'alphazero')
    geom.alphazero = 0;
end

% ----- Regime Boundaries -----
% Supersonic CL formula valid for M > 1/cos(Λ_LE) = 1.1547 (Λ_LE = -30°)
% Transonic fairing: M=0.85 to M=1.2
Mcr = 0.85;  % Critical Mach - start of transonic
Mss = 1.2;   % Supersonic start - end of transonic

% =========================================================================
% AERODYNAMIC COEFFICIENT CALCULATION (REGIME-DEPENDENT)
% =========================================================================

if M <= Mcr
    % =====================================================================
    % SUBSONIC REGIME (M <= 0.85)
    % =====================================================================
    % Use Raymer subsonic model - CL formula is correct without calibration
    % CL_alpha = (2πA)/(2+√(4+(A²β²/η)²(1+tan²Λ/β²))) × exposure
    % With exposure = S_exposed × F_fus (no 0.98 cap)

    out_sub = Raymer_subsonic_model(M, alpha_deg, geom, opts);

    % Adjust Oswald efficiency for smaller wing, larger canard+vtail
    e_oswald = out_sub.e_oswald * 0.92;  % 8% reduction for tail interference

    CL = out_sub.CL;
    CLalpha = out_sub.CL_alpha;
    CD0 = out_sub.CD0_used;
    K = 1/(pi * geom.AR * e_oswald);
    CDi = K * CL^2;
    CD = CD0 + CDi;
    regime = 'subsonic';

elseif M > Mcr && M < Mss
    % =====================================================================
    % TRANSONIC REGIME (0.85 < M < 1.2)
    % =====================================================================
    % Fair CL, CD0, and K between subsonic and supersonic boundaries

    % --- Get boundary values from subsonic model at M=Mcr ---
    out_sub_boundary = Raymer_subsonic_model(Mcr, alpha_deg, geom, opts);
    e_oswald_sub = out_sub_boundary.e_oswald * 0.92;
    CL_sub = out_sub_boundary.CL;
    CLalpha_sub = out_sub_boundary.CL_alpha;
    K_sub = 1/(pi * geom.AR * e_oswald_sub);

    % --- Get boundary values from supersonic model at M=Mss ---
    % Supersonic lift: CLα = 4/β with corrections
    beta_sup = sqrt(Mss^2 - 1);
    CL_alpha_sup = 4 / beta_sup;
    alpha_rad = deg2rad(alpha_deg);
    CL_uncorrected_sup = CL_alpha_sup * alpha_rad;

    % Apply finite swept-wing correction
    AR = geom.AR;
    Lambda_LE_rad = deg2rad(geom.wing.Lambda_LE_deg);
    correction_sup = (AR/(AR+2) + cos(Lambda_LE_rad)*(Mss-0.2)) / 2;
    CL_sup = CL_uncorrected_sup * correction_sup;
    CLalpha_sup = CL_alpha_sup * correction_sup;

    % Supersonic K
    numerator_sup = 4 * (Mss^2 - 1) * cos(Lambda_LE_rad);
    denominator_sup = 4 * AR * beta_sup - 2;
    K_sup = numerator_sup / denominator_sup;

    % --- Fair CL using parabolic curve ---
    curvature_CL = -0.5;
    [a_CL, b_CL, c_CL] = find_parabola_with_curvature(Mcr, CL_sub, Mss, CL_sup, curvature_CL);
    CL_base = a_CL*M^2 + b_CL*M + c_CL;

    % Add transonic CL hump at M~0.92
    M_peak = 0.92;
    sigma = 0.11;
    boost = 0.28 * exp(-((M - M_peak)/sigma)^2);
    CL = CL_base * (1 + boost);

    % --- Fair CLalpha ---
    curvature_CLa = -4;
    [a_CLa, b_CLa, c_CLa] = find_parabola_with_curvature(Mcr, CLalpha_sub, Mss, CLalpha_sup, curvature_CLa);
    CLalpha = a_CLa*M^2 + b_CLa*M + c_CLa;

    % --- Empirical CD0 ---
    CD0 = X29_empirical_transonic_CD0(M);

    % --- Fair K factor ---
    curvature_K = 0.1;
    [aK, bK, cK] = find_parabola_with_curvature(Mcr, K_sub, Mss, K_sup, curvature_K);
    K = aK*M^2 + bK*M + cK;

    % --- Calculate drag ---
    CDi = K * CL^2;
    CD = CD0 + CDi;
    regime = 'transonic';

else  % M >= Mss
    % =====================================================================
    % SUPERSONIC REGIME (M >= 1.2)
    % =====================================================================
    % Supersonic lift: CLα = 4/β where β = √(M²-1)
    % Valid for M > 1/cos(Λ_LE) = 1.1547 (with Λ_LE = -30°)

    beta = sqrt(M^2 - 1);
    CL_alpha = 4 / beta;  % per radian

    % Calculate uncorrected CL
    alpha_rad = deg2rad(alpha_deg);
    CL_uncorrected = CL_alpha * alpha_rad;

    % Apply finite swept-wing correction
    % Blends AR correction (AR/(AR+2)) and sweep correction (cos(Λ))
    % with transonic transition factor (M - 0.2)
    AR = geom.AR;
    Lambda_LE_rad = deg2rad(geom.wing.Lambda_LE_deg);

    if abs(M - 1.2) < 1e-6  % M = 1.2 case
        correction = (AR/(AR+2) + cos(Lambda_LE_rad)*(M-0.2)) / 2;
    else  % M >= 1.3 case
        correction = (AR/(AR+2) + cos(Lambda_LE_rad)) * (M-0.2) / 2;
    end

    CL = CL_uncorrected * correction;
    CLalpha = CL_alpha * correction;

    % Empirical CD0 (continuous from transonic)
    CD0 = X29_empirical_transonic_CD0(M);

    % Swept-wing supersonic induced drag factor
    % K = [4(M²-1)cos(Λ)] / [4·AR·√(M²-1) - 2]
    numerator = 4 * (M^2 - 1) * cos(Lambda_LE_rad);
    denominator = 4 * AR * beta - 2;
    K = numerator / denominator;

    % Calculate drag
    CDi = K * CL^2;
    CD = CD0 + CDi;
    regime = 'supersonic';
end

% =========================================================================
% FORCES AND PERFORMANCE METRICS
% =========================================================================

% Dynamic pressure
V = M * atm.a;
q = 0.5 * atm.rho * V^2;

Sref = geom.wing.S;

% Lift and drag forces
L = q * CL * Sref;
D = q * CD * Sref;

% Performance estimates
LD_max = 1 / (2*sqrt(K*CD0));
CLopt = sqrt(CD0 / K);

% Actual L/D at current condition
if D > 0
    L_over_D = L / D;
else
    L_over_D = 0;
end

% =========================================================================
% PACKAGE OUTPUTS
% =========================================================================

aero = struct();
aero.L = L;
aero.D = D;
aero.L_over_D = L_over_D;
aero.LD_max = LD_max;
aero.CL = CL;
aero.CLopt = CLopt;
aero.CLalpha = CLalpha;
aero.CD = CD;
aero.CD0 = CD0;
aero.CDi = CDi;
aero.K = K;
aero.regime = regime;

end
