function [cruise] = sim_cruise_canards(params, geom)
%SIM_CRUISE_CANARDS Simulate X-29 cruise with canard/stability tracking
%
% Simulates steady level supersonic cruise at constant altitude and Mach.
% Tracks CG movement due to fuel burn and corresponding trim changes.
%
% Uses CALIBRATED physics-based propulsion model via get_F404_thrust_calibrated.
%
% Inputs:
%   params.W_initial      - Initial weight [kg]
%   params.W_fuel_initial - Initial fuel weight [kg]
%   params.h              - Cruise altitude [m]
%   params.M              - Cruise Mach number
%   params.duration_min   - Cruise duration [minutes]
%   params.verbose        - Print detailed output
%   geom                  - Geometry structure for Raymer aero model
%
% Outputs:
%   cruise - Structure with time histories and final state
%
% Author: Dominic Larin
% Date: December 2025

%% Default parameters
if nargin < 1, params = struct(); end
if nargin < 2, geom = []; end

W_initial = get_param(params, 'W_initial', 7500);
W_fuel_initial = get_param(params, 'W_fuel_initial', 1500);
h = get_param(params, 'h', 12200);
M_cruise = get_param(params, 'M', 1.4);
duration_min = get_param(params, 'duration_min', 5);
verbose = get_param(params, 'verbose', true);

%% Constants
g = 9.81;
ft2m = 0.3048;

%% Aircraft parameters
S_ref = 17.54;              % Wing area [m²]
S_canard = 3.4923;          % Canard area [m²]
c_bar = 8.0 * ft2m;         % MAC [m]

CL_alpha_wing = 0.08;       % [1/deg]
CL_0_wing = 0.02;           % Zero-alpha lift (less at supersonic)
CL_alpha_canard = 0.0263;   % [1/deg]

%% Atmosphere at cruise
[T_atm, a, P, rho] = atmosisa(h);
V_cruise = M_cruise * a;
q = 0.5 * rho * V_cruise^2;
h_ft = h / ft2m;

%% Aerodynamic model - using Raymer via get_X29_aero()
% CD will be computed at each timestep using the full aero model

%% Simulation setup
dt = 1.0;  % Time step [s] - same as original
duration_s = duration_min * 60;
n_steps = round(duration_s / dt) + 1;

% Initialize arrays
t = zeros(n_steps, 1);
h_arr = h * ones(n_steps, 1);
V = V_cruise * ones(n_steps, 1);
M_arr = M_cruise * ones(n_steps, 1);
W = zeros(n_steps, 1);
W_fuel = zeros(n_steps, 1);
x_cg = zeros(n_steps, 1);
alpha = zeros(n_steps, 1);
delta_canard = zeros(n_steps, 1);
CL = zeros(n_steps, 1);
CD = zeros(n_steps, 1);
LD = zeros(n_steps, 1);
T_arr = zeros(n_steps, 1);
D_arr = zeros(n_steps, 1);
PLA_arr = zeros(n_steps, 1);
fuel_flow_arr = zeros(n_steps, 1);
Cm_alpha_arr = zeros(n_steps, 1);

% Initial conditions
W(1) = W_initial;
W_fuel(1) = W_fuel_initial;
[x_cg(1), ~] = calculate_X29_CG(W_initial, W_fuel_initial, 1877);

% Get stability derivatives at cruise Mach
[Cm_alpha, Cm_delta_c, Cm_delta_s, Cm_delta_f] = get_Cm_derivatives(M_cruise);

if verbose
    fprintf('  Cruise conditions: h=%.0f m (%.0f ft), M=%.2f, V=%.0f m/s\n', h, h_ft, M_cruise, V_cruise);
    fprintf('  Dynamic pressure: %.2f kPa\n', q/1000);
    fprintf('  Using Raymer aero model via get_X29_aero()\n');
    fprintf('  Cm_alpha = %.3f /rad (unstable > 0)\n', Cm_alpha);
end

%% Cruise simulation with proper fuel flow
for i = 1:n_steps
    t(i) = (i-1) * dt;

    if i > 1
        % Update weight from previous fuel burn
        W(i) = W(i-1) - fuel_flow_arr(i-1) * dt;
        W_fuel(i) = W_fuel(i-1) - fuel_flow_arr(i-1) * dt;

        % Ensure fuel doesn't go negative
        if W_fuel(i) < 0
            W_fuel(i) = 0;
            W(i) = W_initial - W_fuel_initial;  % OEW
        end
    end

    % CG position
    [x_cg(i), cg_data] = calculate_X29_CG(W(i), W_fuel(i), 1877);

    % Required CL for level flight (L = W)
    CL_required = W(i) * g / (q * S_ref);
    CL(i) = CL_required;

    % Estimate alpha from CL for aero model - use Raymer CL_alpha
    if ~isempty(geom)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M_cruise, 5, h, geom);
        alpha_for_aero = CL_required / CL_alpha;
    else
        alpha_for_aero = CL_required / 0.08;  % Fallback
    end
    alpha_for_aero = max(0, min(alpha_for_aero, 15));

    % Get drag from Raymer aero model
    [~, CD_val, ~, ~, ~, ~] = get_X29_aero(M_cruise, alpha_for_aero, h, geom);
    CD(i) = CD_val;
    D_val = q * S_ref * CD(i);
    D_arr(i) = D_val;

    % L/D
    LD(i) = CL(i) / CD(i);

    % Solve for PLA such that T = D (using physics model)
    [PLA_sol, T_sol, fuel_flow_sol, ~] = solve_for_PLA(D_val, M_cruise, h);
    PLA_arr(i) = PLA_sol;
    T_arr(i) = T_sol;
    fuel_flow_arr(i) = fuel_flow_sol;

    % Trim calculation for alpha and canard
    % Estimate alpha from CL
    beta = sqrt(M_cruise^2 - 1);
    CLalpha = 4 / beta * (4.0 / (4.0 + 2)) * 0.8;  % Supersonic, finite wing
    alpha_est = CL_required / CLalpha * (180/pi);  % degrees
    alpha_est = max(0, min(alpha_est, 10));
    alpha(i) = alpha_est;

    % Canard deflection for moment balance
    alpha_rad = alpha_est * pi/180;
    if abs(Cm_delta_c) > 0.01
        delta_rad = -Cm_alpha * alpha_rad / Cm_delta_c;
        delta_deg = delta_rad * 180/pi;
    else
        delta_deg = 5;
    end
    delta_canard(i) = max(-32, min(58, delta_deg));

    Cm_alpha_arr(i) = Cm_alpha;
end

%% Package output
cruise = struct();
cruise.t = t;
cruise.h = h_arr;
cruise.V = V;
cruise.M = M_arr;
cruise.W = W;
cruise.W_fuel = W_fuel;
cruise.x_cg = x_cg;
cruise.alpha = alpha;
cruise.delta_canard = delta_canard;
cruise.CL = CL;
cruise.CD = CD;
cruise.LD = LD;
cruise.T = T_arr;
cruise.D = D_arr;
cruise.PLA = PLA_arr;
cruise.fuel_flow = fuel_flow_arr;
cruise.Cm_alpha = Cm_alpha_arr;

cruise.W_final = W(end);
cruise.W_fuel_final = W_fuel(end);
cruise.fuel_burned = W_fuel(1) - W_fuel(end);
cruise.avg_LD = mean(LD);
cruise.avg_fuel_flow = mean(fuel_flow_arr);

% CG travel during cruise
cruise.cg_travel = x_cg(end) - x_cg(1);

if verbose
    fprintf('  Cruise duration: %.1f min\n', duration_min);
    fprintf('  Fuel burned: %.1f kg (%.1f kg/min)\n', cruise.fuel_burned, cruise.fuel_burned/duration_min);
    fprintf('  Average fuel flow: %.3f kg/s\n', cruise.avg_fuel_flow);
    fprintf('  CG travel: %.3f m (%.2f in)\n', cruise.cg_travel, cruise.cg_travel/0.0254);
    fprintf('  Alpha range: %.2f° to %.2f°\n', min(alpha), max(alpha));
    fprintf('  Canard range: %.1f° to %.1f°\n', min(delta_canard), max(delta_canard));
    fprintf('  PLA range: %.1f° to %.1f°\n', min(PLA_arr), max(PLA_arr));
end

end

%% Helper functions
function val = get_param(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

%% Solve for PLA such that T = D (using calibrated physics model)
function [PLA_sol, T_sol, fuel_flow_sol, TSFC_sol] = solve_for_PLA(D_target, M, h_m)
% Binary search for PLA using calibrated physics propulsion model
% Note: h_m is altitude in METERS (calibrated model expects meters)
PLA_min = 70;
PLA_max = 130;

% Check bounds
prop_min = get_F404_thrust_calibrated(M, h_m, PLA_min);
prop_max = get_F404_thrust_calibrated(M, h_m, PLA_max);

if D_target < prop_min.thrust_N
    % Use minimum power
    PLA_sol = PLA_min;
    prop_sol = prop_min;
elseif D_target > prop_max.thrust_N
    % Use maximum power (can't maintain condition)
    PLA_sol = PLA_max;
    prop_sol = prop_max;
else
    % Binary search
    for iter = 1:50
        PLA_mid = (PLA_min + PLA_max) / 2;
        prop_mid = get_F404_thrust_calibrated(M, h_m, PLA_mid);

        if prop_mid.thrust_N < D_target
            PLA_min = PLA_mid;
        else
            PLA_max = PLA_mid;
        end

        if abs(prop_mid.thrust_N - D_target) < 10  % Within 10 N
            break;
        end
    end
    PLA_sol = PLA_mid;
    prop_sol = prop_mid;
end

T_sol = prop_sol.thrust_N;
fuel_flow_sol = prop_sol.fuel_flow_kg_s;
TSFC_sol = prop_sol.TSFC_imperial;
end
