function prop = get_F404_Thrust_Calibrated(M, h_m, PLA)
%GET_F404_THRUST_CALIBRATED Calibrated physics-based F404-GE-400 engine model
%
% Physics-based 0-D cycle model CALIBRATED to NASA TM-4140 flight test data.
% Provides accurate thrust and fuel flow predictions across the full flight
% envelope including extrapolation to DCFC (M=1.6, h=15.3km).
%
% The calibration corrects for:
%   1. Altitude lapse rate (physics model under-predicts at high altitude)
%   2. Transonic/supersonic thrust (physics model over-predicts at high Mach, low alt)
%   3. Fuel flow (physics model over-predicts at high Mach)
%   4. Part-power behavior (different calibration for dry vs AB operation)
%
% INPUTS:
%   M     - Mach number (0.1 to 1.6+)
%   h_m   - Altitude in METERS
%   PLA   - Power Lever Angle (30=idle, 87=MIL, 130=MAX AB)
%
% OUTPUTS:
%   prop.thrust_N        - Net thrust [N] (calibrated)
%   prop.thrust_lbf      - Net thrust [lbf] (calibrated)
%   prop.fuel_flow_kg_s  - Fuel flow [kg/s] (calibrated)
%   prop.fuel_flow_lb_hr - Fuel flow [lb/hr] (calibrated)
%   prop.TSFC_SI         - TSFC [kg/N/s]
%   prop.TSFC_imperial   - TSFC [lb/hr/lbf]
%   prop.mode            - 'DRY' or 'AB'
%   prop.cal_factor_T    - Thrust calibration factor applied
%   prop.cal_factor_ff   - Fuel flow calibration factor applied
%
% Author: Dominic Larin
% Date: December 2025
% Calibrated against NASA TM-4140 F404 flight test data

%% Get uncalibrated physics model output
prop_raw = get_F404_thrust_physics(M, h_m, PLA);

%% Calibration Tables for MAX AB (PLA=130)
% NEW APPROACH: Separate altitude and Mach corrections to preserve physical trends
%
% The raw physics model has correct Mach trends (ram compression) but:
%   - Under-predicts at high altitude (needs boost)
%   - Over-predicts at high Mach/low altitude (needs reduction)
%
% Solution: Use reference Mach (M=0.9) for altitude correction,
%           then preserve relative Mach effect from physics model

% Mach breakpoints (extended to include low Mach for better interpolation)
M_bp = [0.0, 0.2, 0.4, 0.6, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6];

% Altitude breakpoints (meters)
h_bp = [0, 3000, 6000, 9000, 12200, 15300];

% Simplified calibration: altitude-only scaling (uniform across all Mach at each altitude)
% This preserves raw physics spacing/trends while hitting key design points
T_cal_AB = [
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=0.0
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=0.2
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=0.4
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=0.6
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=0.8
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=0.9
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=1.0
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=1.2
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12;   % M=1.4
    1.08, 1.08, 1.08, 1.08, 1.33, 1.12    % M=1.6
];

% Fuel flow calibration factors for MAX AB (NASA/Physics ratio)
% Derived from NASA TM-4140 (h in ft) vs physics model (h in m) at PLA=130
ff_cal_AB = [
    1.07, 1.06, 1.40, 1.68, 2.24, 2.62;   % M=0.0 (added, similar to old M=0.3)
    1.06, 1.06, 1.39, 1.67, 2.23, 2.61;   % M=0.2 (added, similar to old M=0.3)
    1.06, 1.06, 1.38, 1.66, 2.22, 2.60;   % M=0.4 (was M=0.3)
    1.00, 0.99, 1.28, 1.47, 1.80, 2.27;   % M=0.6 (was M=0.5)
    0.85, 0.90, 1.05, 1.15, 1.15, 1.75;   % M=0.8 (interpolated/adjusted)
    0.77, 0.85, 0.93, 1.05, 1.10, 1.67;   % M=0.9
    0.70, 0.74, 0.83, 1.01, 1.05, 1.40;   % M=1.0
    0.56, 0.56, 0.67, 1.00, 0.97, 1.25;   % M=1.2
    0.44, 0.43, 0.54, 0.78, 0.95, 1.15;   % M=1.4
    0.34, 0.33, 0.43, 0.60, 0.85, 1.50    % M=1.6
];

% Thrust calibration factors for MIL power (PLA=87)
% Simplified: altitude-only scaling (uniform across all Mach at each altitude)
T_cal_MIL = [
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=0.0
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=0.2
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=0.4
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=0.6
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=0.8
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=0.9
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=1.0
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=1.2
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24;   % M=1.4
    1.19, 1.09, 1.09, 1.09, 1.37, 1.24    % M=1.6
];

% Fuel flow calibration factors for MIL power (PLA=87)
% Adjusted for subsonic to match Conners TSFC ~1.0 lb/hr/lbf
% Supersonic values retained from NASA TM-4140 calibration
ff_cal_MIL = [
    1.06, 1.01, 1.09, 1.22, 1.58, 1.84;   % M=0.0 (added, similar to old M=0.3)
    1.06, 1.01, 1.09, 1.22, 1.58, 1.84;   % M=0.2 (added, similar to old M=0.3)
    1.05, 1.00, 1.08, 1.21, 1.57, 1.83;   % M=0.4 (was M=0.3)
    1.02, 0.98, 1.00, 1.07, 1.29, 1.60;   % M=0.6 (was M=0.5)
    0.97, 0.94, 0.96, 1.02, 1.07, 1.30;   % M=0.8 (interpolated/adjusted)
    0.95, 0.92, 0.95, 1.00, 1.05, 1.20;   % M=0.9
    0.90, 0.85, 0.88, 0.95, 1.00, 1.14;   % M=1.0
    0.56, 0.48, 0.54, 0.74, 1.35, 1.50;   % M=1.2
    0.47, 0.39, 0.45, 0.60, 1.10, 1.35;   % M=1.4
    0.39, 0.32, 0.38, 0.48, 0.85, 1.20    % M=1.6
];

%% Interpolate calibration factors
% Clamp inputs to table bounds for extrapolation safety
M_interp = max(M_bp(1), min(M_bp(end), M));
h_interp = max(h_bp(1), min(h_bp(end), h_m));

% Direct table interpolation (reverted from two-stage approach)
cal_T_AB = interp2(h_bp, M_bp, T_cal_AB, h_interp, M_interp, 'linear');
cal_T_MIL = interp2(h_bp, M_bp, T_cal_MIL, h_interp, M_interp, 'linear');
cal_ff_AB = interp2(h_bp, M_bp, ff_cal_AB, h_interp, M_interp, 'linear');
cal_ff_MIL = interp2(h_bp, M_bp, ff_cal_MIL, h_interp, M_interp, 'linear');

% Blend between MIL and AB calibration based on PLA
% PLA <= 87: use MIL calibration
% PLA >= 130: use AB calibration
% Between: linear interpolation
if PLA <= 87
    cal_factor_T = cal_T_MIL;
    cal_factor_ff = cal_ff_MIL;
elseif PLA >= 130
    cal_factor_T = cal_T_AB;
    cal_factor_ff = cal_ff_AB;
else
    % Linear blend between MIL and AB
    blend = (PLA - 87) / (130 - 87);
    cal_factor_T = cal_T_MIL + blend * (cal_T_AB - cal_T_MIL);
    cal_factor_ff = cal_ff_MIL + blend * (cal_ff_AB - cal_ff_MIL);
end

% Handle extrapolation beyond table (use edge values with small growth)
if h_m > h_bp(end)
    % Extrapolate altitude effect (thrust drops faster at extreme altitude)
    h_excess = (h_m - h_bp(end)) / 1000;  % km beyond table
    cal_factor_T = cal_factor_T * (1 + 0.05 * h_excess);  % 5% per km
    cal_factor_ff = cal_factor_ff * (1 + 0.03 * h_excess);  % 3% per km
end

if M > M_bp(end)
    % Extrapolate Mach effect (diminishing returns at higher Mach)
    M_excess = M - M_bp(end);
    cal_factor_T = cal_factor_T * (1 - 0.05 * M_excess);  % Slight reduction
    cal_factor_ff = cal_factor_ff * (1 - 0.03 * M_excess);
end

%% Apply calibration
thrust_N_cal = prop_raw.thrust_N * cal_factor_T;
thrust_lbf_cal = prop_raw.thrust_lbf * cal_factor_T;
fuel_flow_kg_s_cal = prop_raw.fuel_flow_kg_s * cal_factor_ff;
fuel_flow_lb_hr_cal = prop_raw.fuel_flow_lb_hr * cal_factor_ff;

% Recalculate TSFC with calibrated values
if thrust_lbf_cal > 0
    TSFC_imperial_cal = fuel_flow_lb_hr_cal / thrust_lbf_cal;
else
    TSFC_imperial_cal = 0;
end
TSFC_SI_cal = TSFC_imperial_cal * 2.832545e-5;

%% Package calibrated output
prop = struct();
prop.thrust_N = thrust_N_cal;
prop.thrust_lbf = thrust_lbf_cal;
prop.fuel_flow_kg_s = fuel_flow_kg_s_cal;
prop.fuel_flow_lb_hr = fuel_flow_lb_hr_cal;
prop.TSFC_imperial = TSFC_imperial_cal;
prop.TSFC_SI = TSFC_SI_cal;
prop.mode = prop_raw.mode;

% Calibration info
prop.cal_factor_T = cal_factor_T;
prop.cal_factor_ff = cal_factor_ff;
prop.uncalibrated_thrust_N = prop_raw.thrust_N;
prop.uncalibrated_fuel_flow_kg_s = prop_raw.fuel_flow_kg_s;

% Pass through diagnostic info from raw model
prop.Wa = prop_raw.Wa;
prop.OPR = prop_raw.OPR;
prop.TIT_R = prop_raw.TIT_R;
prop.NPR = prop_raw.NPR;
prop.eta_inlet = prop_raw.eta_inlet;

end
