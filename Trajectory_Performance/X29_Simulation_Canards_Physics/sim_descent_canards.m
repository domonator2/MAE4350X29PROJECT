function [descent] = sim_descent_canards(params, geom)
%SIM_DESCENT_CANARDS Simulate X-29 descent with canard/stability tracking
%
% Descends from cruise altitude to pattern altitude while decelerating.
% Uses idle/partial throttle with speed brakes as needed.
%
% Uses Raymer-based aerodynamics model via get_X29_aero().
%
% Inputs:
%   params.W_initial      - Initial weight [kg]
%   params.W_fuel_initial - Initial fuel weight [kg]
%   params.h_initial      - Starting altitude [m]
%   params.h_final        - Target altitude [m] (pattern altitude)
%   params.verbose        - Print detailed output
%   geom                  - Geometry structure for Raymer aero model
%
% Outputs:
%   descent - Structure with time histories and final state
%
% Author: Dominic Larin
% Date: December 3 2025

%% Default parameters
if nargin < 1, params = struct(); end
if nargin < 2, geom = []; end

W_initial = get_param(params, 'W_initial', 7200);
W_fuel_initial = get_param(params, 'W_fuel_initial', 1200);
h_initial = get_param(params, 'h_initial', 12200);
h_final = get_param(params, 'h_final', 1000);
verbose = get_param(params, 'verbose', true);

%% Constants
g = 9.81;
ft2m = 0.3048;

%% Aircraft parameters
S_ref = 17.45;              % Wing area [m²]
S_canard = 3.4923;          % Canard area [m²]
c_bar = 8.0 * ft2m;         % MAC [m]

CL_alpha_wing = 0.08;       % [1/deg]
CL_0_wing = 0.05;
CL_alpha_canard = 0.0263;   % [1/deg]
CL_max_landing = 1.6;       % Max CL with flaps for landing

%% Calculate target approach speed (must match landing phase)
% Landing approach speed = 1.2 * V_stall (Raymer methodology)
rho_sl = 1.225;  % Sea level density [kg/m³]
V_stall = sqrt(2 * W_initial * g / (rho_sl * S_ref * CL_max_landing));
V_approach = 1.2 * V_stall;  % Target speed at end of descent

%% Descent profile
% Target: descend at 4200 ft/min while decelerating (faster descent)
descent_rate_target = -21.4;  % m/s (~4200 ft/min)

%% Simulation setup
dt = 1.0;  % Time step [s]
t_max = 585;  % Max descent time [s]
n_max = round(t_max / dt) + 1;

% Initialize arrays
t = zeros(n_max, 1);
h = zeros(n_max, 1);
V = zeros(n_max, 1);
M = zeros(n_max, 1);
W = zeros(n_max, 1);
W_fuel = zeros(n_max, 1);
x_cg = zeros(n_max, 1);
alpha = zeros(n_max, 1);
delta_canard = zeros(n_max, 1);
CL = zeros(n_max, 1);
CD = zeros(n_max, 1);
LD = zeros(n_max, 1);
T_arr = zeros(n_max, 1);
D_arr = zeros(n_max, 1);

% Initial conditions
h(1) = h_initial;
[~, a_init, ~, ~] = atmosisa(h_initial);
M_init = 1.4;  % Start from cruise Mach
V(1) = M_init * a_init;
M(1) = M_init;
W(1) = W_initial;
W_fuel(1) = W_fuel_initial;
[x_cg(1), ~] = calculate_X29_CG(W_initial, W_fuel_initial, 1877);

%% Descent simulation
i = 1;
reached_target = false;

while i < n_max && ~reached_target
    i = i + 1;
    t(i) = t(i-1) + dt;

    % Current state
    h_curr = h(i-1);
    V_curr = V(i-1);
    W_curr = W(i-1);
    W_fuel_curr = W_fuel(i-1);

    % Atmosphere
    [T_atm, a, P, rho] = atmosisa(h_curr);
    M_curr = V_curr / a;
    q = 0.5 * rho * V_curr^2;

    % Target speed (decelerate as we descend)
    % Linear interpolation from cruise speed to approach speed
    h_frac = (h_curr - h_final) / (h_initial - h_final);
    h_frac = max(0, min(1, h_frac));

    % Target velocity: cruise speed at top, approach speed at bottom
    % This ensures smooth handoff to landing phase
    V_cruise = M_init * a_init;  % Initial cruise velocity
    V_target = V_approach + h_frac * (V_cruise - V_approach);

    % Required CL for current descent
    % At descent, L ≈ W * cos(gamma)
    gamma_deg = asind(descent_rate_target / V_curr);
    CL_required = W_curr * g * cosd(gamma_deg) / (q * S_ref);
    CL_required = min(CL_required, 1.0);
    CL(i) = CL_required;

    % Trim calculation - use CL_alpha from Raymer model
    if ~isempty(geom)
        [~, ~, ~, ~, ~, CL_alpha_raymer] = get_X29_aero(M_curr, 5, h_curr, geom);
        alpha_est = CL_required / CL_alpha_raymer;
    else
        alpha_est = (CL_required - CL_0_wing) / CL_alpha_wing;  % Fallback
    end
    alpha_est = max(0, min(alpha_est, 10));

    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_curr);

    if abs(Cm_delta_c) > 0.01
        alpha_rad = alpha_est * pi/180;
        delta_rad = -Cm_alpha * alpha_rad / Cm_delta_c;
        delta_deg = delta_rad * 180/pi;
    else
        delta_deg = 5;
    end
    delta_deg = max(-32, min(58, delta_deg));

    alpha(i) = alpha_est;
    delta_canard(i) = delta_deg;

    % Drag from Raymer aero model
    [~, CD_clean, ~, ~, ~, ~] = get_X29_aero(M_curr, alpha_est, h_curr, geom);
    CD(i) = CD_clean;

    % Add speed brake drag if decelerating too slowly
    if V_curr > V_target * 1.1
        CD(i) = CD(i) + 0.02;  % Speed brake increment
    end

    D = q * S_ref * CD(i);
    D_arr(i) = D;
    LD(i) = CL(i) / CD(i);

    % Thrust at idle - physics-based model with PLA=35 (idle)
    prop_idle = get_F404_thrust_calibrated(max(0.3, M_curr), h_curr, 35);
    T_idle = prop_idle.thrust_N;
    fuel_flow_idle = prop_idle.fuel_flow_kg_s;
    T_arr(i) = T_idle;

    % Descent dynamics - controlled descent at target rate
    % Use target descent rate directly
    dh_dt = descent_rate_target;

    % Deceleration
    if V_curr > V_target
        dV_dt = -2;  % Decelerate at 2 m/s²
    else
        dV_dt = 0;
    end

    % Update state
    % Don't go below approach speed (which could be ~70 m/s)
    V(i) = max(V_curr + dV_dt * dt, V_approach * 0.95);
    h(i) = h_curr + dh_dt * dt;

    % Mach
    [~, a_new, ~, ~] = atmosisa(h(i));
    M(i) = V(i) / a_new;

    % Fuel burn at idle
    W_fuel(i) = W_fuel_curr - fuel_flow_idle * dt;
    W(i) = W_curr - fuel_flow_idle * dt;

    % CG
    [x_cg(i), ~] = calculate_X29_CG(W(i), W_fuel(i), 1877);

    % Check if reached pattern altitude
    if h(i) <= h_final
        reached_target = true;
        h(i) = h_final;
    end
end

%% Trim arrays
n_actual = i;
t = t(1:n_actual);
h = h(1:n_actual);
V = V(1:n_actual);
M = M(1:n_actual);
W = W(1:n_actual);
W_fuel = W_fuel(1:n_actual);
x_cg = x_cg(1:n_actual);
alpha = alpha(1:n_actual);
delta_canard = delta_canard(1:n_actual);
CL = CL(1:n_actual);
CD = CD(1:n_actual);
LD = LD(1:n_actual);
T_arr = T_arr(1:n_actual);
D_arr = D_arr(1:n_actual);

%% Package output
descent = struct();
descent.t = t;
descent.h = h;
descent.V = V;
descent.M = M;
descent.W = W;
descent.W_fuel = W_fuel;
descent.x_cg = x_cg;
descent.alpha = alpha;
descent.delta_canard = delta_canard;
descent.CL = CL;
descent.CD = CD;
descent.LD = LD;
descent.T = T_arr;
descent.D = D_arr;

descent.W_final = W(end);
descent.W_fuel_final = W_fuel(end);
descent.descent_time = t(end);
descent.V_final = V(end);
descent.V_approach_target = V_approach;

if verbose
    fprintf('  Descent time: %.1f s (%.1f min)\n', t(end), t(end)/60);
    fprintf('  Altitude: %.0f m → %.0f m\n', h(1), h(end));
    fprintf('  Speed: %.0f m/s → %.0f m/s (target approach: %.0f m/s)\n', V(1), V(end), V_approach);
    fprintf('  Mach: %.2f → %.2f\n', M(1), M(end));
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
