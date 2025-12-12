function prop = get_F404_thrust_physics(M, h_m, PLA)
%GET_F404_THRUST_PHYSICS Physics-based F404-GE-400 engine model
%
% Pure thermodynamic cycle model with NO empirical corrections.
% Based on station-by-station analysis: Inlet -> Fan -> HPC -> Combustor ->
% HPT -> LPT -> Mixer -> Afterburner -> Nozzle
%
% INPUTS:
%   M     - Mach number
%   h_m   - Altitude in METERS (uses atmosisa for consistency)
%   PLA   - Power Lever Angle (30=idle, 87=MIL, 130=MAX AB)
%
% OUTPUTS:
%   prop.thrust_N        - Net thrust [N]
%   prop.thrust_lbf      - Net thrust [lbf]
%   prop.fuel_flow_kg_s  - Fuel flow [kg/s]
%   prop.fuel_flow_lb_hr - Fuel flow [lb/hr]
%   prop.TSFC_SI         - TSFC [kg/N/s]
%   prop.TSFC_imperial   - TSFC [lb/hr/lbf]
%   prop.mode            - 'DRY' or 'AB'
%
% Author: Dominic Larin (converted from Python physics model)
% Date: December 2025

%% Design Point Constants (SLS)
ref_airflow = 144.0;   % lb/s
ref_BPR = 0.34;
ref_OPR = 25.0;
ref_FPR = 4.0;
ref_TIT_R = 3010.0;    % ~1670K
ref_eta_fan = 0.88;
ref_eta_hpc = 0.87;
ref_eta_turb = 0.95;

% Thermodynamic properties
LHV = 18400.0 * 778.169;      % BTU/lb to ft-lbf/lbm
Cp_c = 0.24 * 778.169;        % Cold section Cp
Cp_h = 0.27 * 778.169;        % Hot section Cp
gamma_c = 1.4;
gamma_h = 1.3;
g0 = 32.174;

%% 1. Atmosphere Model - Use MATLAB's atmosisa (same as trajectory sim)
[T_K, ~, P_Pa, ~] = atmosisa(h_m);

% Convert to Imperial units for cycle calculations
T_amb = T_K * 1.8;           % Kelvin to Rankine
P_amb = P_Pa * 0.020885;     % Pascals to lbf/ft²

%% 2. Inlet Recovery (X-29 Fixed Geometry, eta_r = 0.96 baseline)
if M <= 1.0
    eta_inlet = 0.96;
else
    eta_inlet = 0.96 * (1.0 - 0.10 * (M - 1.0)^1.5);
end

%% 3. Inlet Ram Compression
tr_fac = 1.0 + 0.5 * (gamma_c - 1) * M^2;
Tt2 = T_amb * tr_fac;
Pt2 = P_amb * (tr_fac)^(gamma_c / (gamma_c - 1)) * eta_inlet;

%% 4. Control Schedule (PLA -> RPM/TIT mapping)
PLA = max(30, min(130, PLA));

if PLA <= 87
    throt = (PLA - 30) / (87 - 30);
    ab_ratio = 0.0;
else
    throt = 1.0;
    ab_ratio = (PLA - 87) / (130 - 87);
end

% Map throttle to corrected speed
N2_corr = 0.6 + 0.4 * throt;  % 60% to 100%

% Off-design scaling
flow_scalar = N2_corr;
pr_scalar = N2_corr^2.5;

%% 5. Fan Compression
Wa = ref_airflow * flow_scalar * (Pt2/2116.22) / sqrt(Tt2/518.67);
FPR = 1.0 + (ref_FPR - 1.0) * pr_scalar;
Pt13 = Pt2 * FPR;
Tt13 = Tt2 * (1 + (FPR^((gamma_c-1)/gamma_c) - 1)/ref_eta_fan);

% Splitter
BPR = ref_BPR;
Wa_core = Wa / (1 + BPR);
Wa_byp = Wa - Wa_core;

%% 5b. Bleed and Cooling Extractions
% Bleed: Customer bleed, ECS, anti-ice (extracted from HPC, does not return)
% Cooling: Turbine cooling air (extracted from HPC, rejoins after HPT)
% Values based on F404 data and Carlos v4.10 calibration
if PLA <= 87
    % Dry power - higher bleed fraction at lower power
    bleed_frac = 0.05;   % 5% bleed at MIL
    cool_frac = 0.022;   % 2.2% cooling at MIL
else
    % Afterburner - reduced bleed/cooling
    bleed_frac = 0.04;   % 4% bleed at AB
    cool_frac = 0.020;   % 2.0% cooling at AB
end

Wa_bleed = Wa_core * bleed_frac;   % Bleed air (lost from cycle)
Wa_cool = Wa_core * cool_frac;     % Cooling air (rejoins after HPT)
Wa_core_eff = Wa_core - Wa_bleed - Wa_cool;  % Effective core flow through combustor

%% 6. High Pressure Compressor
CPR_des = ref_OPR / ref_FPR;
CPR = 1.0 + (CPR_des - 1.0) * pr_scalar;
Pt3 = Pt13 * CPR;
Tt3 = Tt13 * (1 + (CPR^((gamma_c-1)/gamma_c) - 1)/ref_eta_hpc);

%% 7. Combustor (Main)
TIT_target = Tt2 + (ref_TIT_R - 518.67) * throt;
if TIT_target < Tt3
    TIT_target = Tt3 + 100;
end

% Fuel flow from energy balance (uses effective core flow, not total)
Wf_main = (Wa_core_eff * Cp_h * (TIT_target - Tt3)) / (LHV * 0.98);
Pt4 = Pt3 * 0.95;  % Burner pressure loss

%% 8. HPT & LPT (Work Balance)
% HPT drives HPC (full core flow compresses, but only effective flow through combustor)
dHt_hpc = Wa_core * Cp_c * (Tt3 - Tt13);
dHt_hpt = dHt_hpc / 0.99;  % 1% mechanical loss

% HPT expansion (only effective core + fuel passes through HPT)
Wa_hpt = Wa_core_eff + Wf_main/32.2;
Tt45 = TIT_target - dHt_hpt / (Wa_hpt * Cp_h);
Pt45 = Pt4 * (Tt45/TIT_target)^(gamma_h/(gamma_h-1));

% Cooling air rejoins after HPT (mixes at Tt3 temperature)
% Energy balance for mixing: (Wa_hpt * Tt45 + Wa_cool * Tt3) / (Wa_hpt + Wa_cool)
Wa_post_cool = Wa_hpt + Wa_cool;
Tt45_mixed = (Wa_hpt * Cp_h * Tt45 + Wa_cool * Cp_c * Tt3) / (Wa_post_cool * Cp_h);
% Pressure unchanged (simplified - cooling injection has ~same Pt)

% LPT drives Fan
dHt_fan = Wa * Cp_c * (Tt13 - Tt2);
dHt_lpt = dHt_fan / 0.99;
Tt5 = Tt45_mixed - dHt_lpt / (Wa_post_cool * Cp_h);
Pt5 = Pt45 * (Tt5/Tt45_mixed)^(gamma_h/(gamma_h-1));

%% 9. Mixer
% Core flow at mixer = effective + cooling + fuel (bleed is lost)
Wa_core_mixer = Wa_core_eff + Wa_cool + Wf_main/32.2;  % = Wa_post_cool

Pt16 = Pt13 * 0.98;  % Bypass duct loss
Hm_core = Wa_core_mixer * Cp_h * Tt5;
Hm_byp = Wa_byp * Cp_c * Tt13;
Wm_tot = Wa_core_mixer + Wa_byp;  % Total mixed flow (excludes bleed)
Tt_mix = (Hm_core + Hm_byp) / (Wm_tot * Cp_h);
Pt_mix = (Pt5 * (Wa_core - Wa_bleed) + Pt16 * Wa_byp) / (Wa - Wa_bleed);

%% 10. Afterburner
Wf_ab = 0.0;
Tt7 = Tt_mix;
Pt7 = Pt_mix * 0.97;  % Mixer/AB duct loss

if ab_ratio > 0
    % Max AB Temp: 4000 R (2220 K) - realistic for F404
    Max_AB_T = 4000.0;
    Tt7_target = Tt_mix + (Max_AB_T - Tt_mix) * ab_ratio;
    Q_ab = Wm_tot * Cp_h * (Tt7_target - Tt_mix);
    Wf_ab = Q_ab / (LHV * 0.90);  % 90% AB efficiency
    Tt7 = Tt7_target;
    Pt7 = Pt7 * (1.0 - 0.05 * ab_ratio);  % Rayleigh loss
end

%% 11. Nozzle Expansion
NPR = Pt7 / P_amb;
NPR = max(NPR, 1.0);

% Nozzle Coefficients (realistic C-D nozzle losses)
% Cd = discharge coefficient (mass flow capacity)
% Cv = velocity coefficient (kinetic energy efficiency)
% Cfg = gross thrust coefficient = Cd * Cv
Cd_noz = 0.98;   % Typical for well-designed C-D nozzle
Cv_noz = 0.98;   % Accounts for boundary layer and non-ideal expansion
Cfg = Cd_noz * Cv_noz;  % ~0.96 combined

if NPR > 1.89  % Choked
    V9 = sqrt(2 * g0 * Cp_h * Tt7 * (1 - (1/NPR)^((gamma_h-1)/gamma_h)));
else
    M9 = sqrt(2/(gamma_h-1) * (NPR^((gamma_h-1)/gamma_h) - 1));
    V9 = M9 * sqrt(gamma_h * 32.174 * 53.35 * Tt7);
end

%% 12. Thrust Calculation
% Gross thrust (with nozzle efficiency)
Fg = Cfg * (Wm_tot + Wf_ab/32.2) / 32.2 * V9;

% Ram drag
V_flight = M * sqrt(1.4 * 32.174 * 53.35 * T_amb);
Ram_Drag = (Wa / 32.2) * V_flight;

% Net thrust
Fn_lbf = Fg - Ram_Drag;
Fn_lbf = max(Fn_lbf, 0);  % Can't be negative

%% 13. Fuel Flow and TSFC
Total_Wf_lb_s = Wf_main + Wf_ab;

% Minimum idle fuel flow for combustor stability
% F404 ground idle is typically 500-700 lb/hr (~600 lb/hr nominal)
% This accounts for real combustor fuel scheduling at low power
min_idle_fuel_lb_hr = 600;  % lb/hr (typical F404 ground idle)
min_idle_fuel_lb_s = min_idle_fuel_lb_hr / 3600;

if PLA <= 40  % At or near idle
    Total_Wf_lb_s = max(Total_Wf_lb_s, min_idle_fuel_lb_s);
end

Total_Wf_lb_hr = Total_Wf_lb_s * 3600.0;

if Fn_lbf > 0
    TSFC_imperial = Total_Wf_lb_hr / Fn_lbf;
else
    TSFC_imperial = 0;
end

%% Package Output
prop = struct();
prop.thrust_lbf = Fn_lbf;
prop.thrust_N = Fn_lbf * 4.44822;
prop.fuel_flow_lb_hr = Total_Wf_lb_hr;
prop.fuel_flow_kg_s = Total_Wf_lb_s * 0.453592;
prop.TSFC_imperial = TSFC_imperial;
prop.TSFC_SI = TSFC_imperial * 2.832545e-5;  % Convert to kg/N/s

if ab_ratio > 0
    prop.mode = 'AB';
else
    prop.mode = 'DRY';
end

% Additional outputs for debugging
prop.Wa = Wa;
prop.Wa_core = Wa_core;
prop.Wa_core_eff = Wa_core_eff;
prop.Wa_bleed = Wa_bleed;
prop.Wa_cool = Wa_cool;
prop.bleed_frac = bleed_frac;
prop.cool_frac = cool_frac;
prop.OPR = Pt3/Pt2;
prop.TIT_R = TIT_target;
prop.NPR = NPR;
prop.eta_inlet = eta_inlet;
prop.Cfg = Cfg;           % Nozzle gross thrust coefficient
prop.eta_AB = 0.90;       % AB efficiency (when active)

end
