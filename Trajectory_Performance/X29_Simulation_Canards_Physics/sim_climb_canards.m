function [climb] = sim_climb_canards(params, geom)
%SIM_CLIMB_CANARDS Simulate X-29 climb with canard/trim tracking
%
% Uses Hull's equations of motion with ODE45 numerical integration.
% Same 3-phase climb schedule as sim_climb.m, but tracks canard deflection.
%
% Hull's Equations (vertical plane):
%   x' = V cos(gamma)
%   h' = V sin(gamma)
%   V' = (g/W)[T cos(alpha) - D - W sin(gamma)]
%   gamma' = (g/WV)[T sin(alpha) + L - W cos(gamma)]
%   W' = -fuel_flow * g
%
% Climb Schedule:
%   Phase 1: Accelerate at low altitude (h~300m) to q = 21 kPa
%   Phase 2: Maintain constant q = 21 kPa, climb until M = 0.9
%   Phase 3: Maintain constant M = 0.9, climb until h = h_cruise
%
% Inputs:
%   params.W_initial      - Initial weight [kg]
%   params.W_fuel_initial - Initial fuel weight [kg]
%   params.h_initial      - Starting altitude [m]
%   params.V_initial      - Starting velocity [m/s]
%   params.h_target       - Target cruise altitude [m]
%   params.M_target       - Target Mach at end of climb (default 0.9)
%   params.verbose        - Print detailed output
%   geom                  - Geometry structure for Raymer aero model
%
% Outputs:
%   climb - Structure with time histories and final state
%
% Author: Dominic Larin
% Date: December 2025

%% Default parameters
if nargin < 1, params = struct(); end
if nargin < 2, geom = []; end

W0 = get_param(params, 'W_initial', 7900);
W_fuel_initial = get_param(params, 'W_fuel_initial', 1800);
h0 = get_param(params, 'h_initial', 15);
V0 = get_param(params, 'V_initial', 100);
h_cruise = get_param(params, 'h_target', 12200);
M_target = get_param(params, 'M_target', 0.9);
verbose = get_param(params, 'verbose', true);

%% Constants and targets
g = 9.81;
q_target = 21000;       % Pa - target dynamic pressure
M_transition = 0.9;     % Mach for phase 2->3 transition

% Altitude for level acceleration (Phase 1)
h_accel = max(300, h0 + 100);  % m

%% Aerodynamic parameters
S_ref = 17.54;  % Wing reference area [m²]
aero_params.geom = geom;
aero_params.use_raymer = ~isempty(geom);
aero_params.S_ref = S_ref;
aero_params.CD0 = 0.025;
aero_params.AR = 4.0;
aero_params.e = 0.80;
aero_params.K = 1 / (pi * aero_params.AR * aero_params.e);
aero_params.CL_alpha = 0.08;  % per degree (fallback - use Raymer when available)

%% Header
if verbose
    fprintf('==================================================================\n');
    fprintf('    X-29 CLIMB SIMULATION (ODE45 + HULL EQUATIONS + CANARDS)     \n');
    fprintf('==================================================================\n\n');

    fprintf('Hull''s Equations of Motion:\n');
    fprintf('  x'' = V*cos(gamma)\n');
    fprintf('  h'' = V*sin(gamma)\n');
    fprintf('  V'' = (g/W)[T*cos(alpha) - D - W*sin(gamma)]\n');
    fprintf('  W'' = -fuel_flow*g\n\n');

    fprintf('CLIMB SCHEDULE:\n');
    fprintf('  Phase 1: Climb to h=%.0fm, then accel to q = %.0f kPa\n', h_accel, q_target/1000);
    fprintf('  Phase 2: Constant q = %.0f kPa climb until M = %.2f\n', q_target/1000, M_transition);
    fprintf('  Phase 3: Constant M = %.2f climb until h = %.1f km\n\n', M_transition, h_cruise/1000);

    fprintf('INITIAL CONDITIONS:\n');
    fprintf('  Weight:    %.0f kg\n', W0);
    fprintf('  Altitude:  %.0f m\n', h0);
    fprintf('  Velocity:  %.1f m/s (%.0f kts)\n\n', V0, V0*1.944);
end

%% ODE options
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1.0);

% Storage for all phases
all_t = [];
all_Y = [];
all_phase = [];

%% PHASE 1a: Initial climb to acceleration altitude
if verbose
    fprintf('--- PHASE 1a: Climb to h = %.0f m ---\n', h_accel);
end

% State: [x, h, V, W] - W in Newtons
y0_1a = [0; h0; V0; W0*g];
gamma_1a = deg2rad(10);  % Fixed 10 deg climb

% Event: stop when h >= h_accel
opts_1a = odeset(ode_opts, 'Events', @(t,y) event_altitude(t, y, h_accel));

[t_1a, Y_1a] = ode45(@(t,y) climb_EOM_fixed_gamma(t, y, gamma_1a, aero_params, g), ...
    [0, 120], y0_1a, opts_1a);

if verbose
    fprintf('  Completed: t=%.1fs, h=%.0fm, V=%.0fm/s\n', t_1a(end), Y_1a(end,2), Y_1a(end,3));
end

% Store
all_t = [all_t; t_1a];
all_Y = [all_Y; Y_1a];
all_phase = [all_phase; ones(length(t_1a), 1)];

%% PHASE 1b: Level acceleration to q = 21 kPa
if verbose
    fprintf('--- PHASE 1b: Level accel to q = %.0f kPa ---\n', q_target/1000);
end

y0_1b = Y_1a(end, :)';
t_offset_1b = t_1a(end);

% Event: stop when q >= q_target
opts_1b = odeset(ode_opts, 'Events', @(t,y) event_dynamic_pressure(t, y, q_target));

[t_1b, Y_1b] = ode45(@(t,y) climb_EOM_level(t, y, aero_params, g), ...
    [0, 120], y0_1b, opts_1b);

[~, a_1b, ~, rho_1b] = atmosisa(Y_1b(end,2));
M_1b = Y_1b(end,3) / a_1b;
q_1b = 0.5 * rho_1b * Y_1b(end,3)^2;

if verbose
    fprintf('  Completed: t=%.1fs, M=%.2f, q=%.1fkPa\n', t_offset_1b + t_1b(end), M_1b, q_1b/1000);
end

% Store (offset time)
all_t = [all_t; t_1b(2:end) + t_offset_1b];
all_Y = [all_Y; Y_1b(2:end, :)];
all_phase = [all_phase; ones(length(t_1b)-1, 1)];

t_phase1 = t_offset_1b + t_1b(end);
fuel_phase1 = (W0*g - Y_1b(end,4)) / g;

%% PHASE 2: Constant q climb until M = 0.9
if verbose
    fprintf('--- PHASE 2: Constant q = %.0f kPa climb until M = %.2f ---\n', q_target/1000, M_transition);
end

y0_2 = Y_1b(end, :)';
t_offset_2 = t_phase1;

% Event: stop when M >= M_transition
opts_2 = odeset(ode_opts, 'Events', @(t,y) event_mach(t, y, M_transition));

[t_2, Y_2] = ode45(@(t,y) climb_EOM_const_q(t, y, q_target, aero_params, g), ...
    [0, 600], y0_2, opts_2);

[~, a_2, ~, ~] = atmosisa(Y_2(end,2));
M_2 = Y_2(end,3) / a_2;

if verbose
    fprintf('  Completed: t=%.1fs, h=%.0fm, M=%.2f\n', t_offset_2 + t_2(end), Y_2(end,2), M_2);
end

% Store
all_t = [all_t; t_2(2:end) + t_offset_2];
all_Y = [all_Y; Y_2(2:end, :)];
all_phase = [all_phase; 2*ones(length(t_2)-1, 1)];

t_phase2 = t_offset_2 + t_2(end);
fuel_phase2 = (W0*g - Y_2(end,4)) / g;

%% PHASE 3: Constant M = 0.9 climb to cruise altitude
if verbose
    fprintf('--- PHASE 3: Constant M = %.2f climb to h = %.1f km ---\n', M_transition, h_cruise/1000);
end

y0_3 = Y_2(end, :)';
t_offset_3 = t_phase2;

% Event: stop when h >= h_cruise
opts_3 = odeset(ode_opts, 'Events', @(t,y) event_altitude(t, y, h_cruise));

[t_3, Y_3] = ode45(@(t,y) climb_EOM_const_M(t, y, M_transition, aero_params, g), ...
    [0, 600], y0_3, opts_3);

if verbose
    fprintf('  Completed: t=%.1fs, h=%.0fm\n', t_offset_3 + t_3(end), Y_3(end,2));
end

% Store
all_t = [all_t; t_3(2:end) + t_offset_3];
all_Y = [all_Y; Y_3(2:end, :)];
all_phase = [all_phase; 3*ones(length(t_3)-1, 1)];

%% Final state
t_final = t_offset_3 + t_3(end);
h_final = Y_3(end, 2);
V_final = Y_3(end, 3);
W_final = Y_3(end, 4) / g;  % Convert back to kg
fuel_burned = W0 - W_final;

[~, a_final, ~, ~] = atmosisa(h_final);
M_final = V_final / a_final;

%% Compute derived quantities (alpha, canard, CG)
n_points = length(all_t);
M_hist = zeros(n_points, 1);
alpha_hist = zeros(n_points, 1);
delta_canard_hist = zeros(n_points, 1);
CL_hist = zeros(n_points, 1);
CD_hist = zeros(n_points, 1);
x_cg_hist = zeros(n_points, 1);

for i = 1:n_points
    h_i = all_Y(i, 2);
    V_i = all_Y(i, 3);
    W_i = all_Y(i, 4);  % Weight in Newtons
    W_kg = W_i / g;

    [~, a_i, ~, rho_i] = atmosisa(h_i);
    M_hist(i) = V_i / a_i;
    q_i = 0.5 * rho_i * V_i^2;

    % CL required for level flight approximation
    CL_i = W_i / (q_i * S_ref);
    CL_i = max(0.1, min(1.2, CL_i));
    CL_hist(i) = CL_i;

    % Alpha estimate - use CL_alpha from Raymer model
    if aero_params.use_raymer && ~isempty(geom)
        [~, ~, ~, ~, ~, CL_alpha_i] = get_X29_aero(M_hist(i), 5, h_i, geom);
        alpha_hist(i) = CL_i / CL_alpha_i;  % deg
    else
        alpha_hist(i) = CL_i / 0.08;  % deg (fallback)
    end
    alpha_hist(i) = max(0, min(15, alpha_hist(i)));

    % CD from Raymer model
    if aero_params.use_raymer && ~isempty(geom)
        [~, CD_hist(i), ~, ~, ~, ~] = get_X29_aero(M_hist(i), alpha_hist(i), h_i, geom);
    else
        CD_hist(i) = aero_params.CD0 + aero_params.K * CL_i^2;
    end

    % Canard deflection for trim
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_hist(i));
    if abs(Cm_delta_c) > 0.01
        delta_rad = -Cm_alpha * (alpha_hist(i) * pi/180) / Cm_delta_c;
        delta_canard_hist(i) = delta_rad * 180/pi;
    else
        delta_canard_hist(i) = 5;
    end
    delta_canard_hist(i) = max(-32, min(58, delta_canard_hist(i)));

    % CG position
    W_fuel_i = W_fuel_initial - (W0 - W_kg);
    W_fuel_i = max(0, W_fuel_i);
    [x_cg_hist(i), ~] = calculate_X29_CG(W_kg, W_fuel_i, W_fuel_initial);
end

%% Summary
if verbose
    fprintf('\n==================================================================\n');
    fprintf('CLIMB SUMMARY (ODE45 + Hull''s Equations)\n');
    fprintf('==================================================================\n');
    fprintf('Phase 1 (Accel to q=21kPa):\n');
    fprintf('  Time:        %.1f s\n', t_phase1);
    fprintf('  Fuel:        %.1f kg\n', fuel_phase1);
    fprintf('\nPhase 2 (Const q climb):\n');
    fprintf('  Time:        %.1f s (%.1f s cumulative)\n', t_phase2-t_phase1, t_phase2);
    fprintf('  Fuel:        %.1f kg (%.1f kg cumulative)\n', fuel_phase2-fuel_phase1, fuel_phase2);
    fprintf('\nPhase 3 (Const M climb):\n');
    fprintf('  Time:        %.1f s (%.1f s cumulative)\n', t_final-t_phase2, t_final);
    fprintf('  Fuel:        %.1f kg (%.1f kg cumulative)\n', fuel_burned-fuel_phase2, fuel_burned);
    fprintf('\n==================================================================\n');
    fprintf('TOTALS:\n');
    fprintf('  Total time:        %.1f s (%.2f min)\n', t_final, t_final/60);
    fprintf('  Total fuel:        %.1f kg\n', fuel_burned);
    fprintf('  Final altitude:    %.0f m (%.0f ft)\n', h_final, h_final*3.281);
    fprintf('  Final Mach:        %.2f\n', M_final);
    fprintf('  Final weight:      %.0f kg\n', W_final);
    fprintf('  Canard range:      %.1f deg to %.1f deg\n', min(delta_canard_hist), max(delta_canard_hist));
    fprintf('==================================================================\n');
end

%% Package results
climb = struct();
climb.t = all_t;
climb.x = all_Y(:, 1);
climb.h = all_Y(:, 2);
climb.V = all_Y(:, 3);
climb.W = all_Y(:, 4) / g;  % Convert to kg
climb.M = M_hist;
climb.alpha = alpha_hist;
climb.delta_canard = delta_canard_hist;
climb.CL = CL_hist;
climb.CD = CD_hist;
climb.x_cg = x_cg_hist;
climb.phase = all_phase;

climb.W_final = W_final;
climb.W_fuel_final = W_fuel_initial - fuel_burned;
climb.total_time = t_final;
climb.total_fuel = fuel_burned;
climb.climb_time = t_final;

end

%% Helper function
function val = get_param(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end

%% Event functions
function [value, isterminal, direction] = event_altitude(t, y, h_target)
    value = y(2) - h_target;
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = event_dynamic_pressure(t, y, q_target)
    h = y(2);
    V = y(3);
    [~, ~, ~, rho] = atmosisa(h);
    q = 0.5 * rho * V^2;
    value = q - q_target;
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = event_mach(t, y, M_target)
    h = y(2);
    V = y(3);
    [~, a, ~, ~] = atmosisa(h);
    M = V / a;
    value = M - M_target;
    isterminal = 1;
    direction = 1;
end

%% Equations of Motion - Fixed gamma climb (Phase 1a)
function dydt = climb_EOM_fixed_gamma(t, y, gamma, aero, g)
    h = max(0, y(2));
    V = max(10, y(3));
    W = y(4);  % Weight in Newtons

    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);
    M = V / a;
    q = 0.5 * rho * V^2;

    % Propulsion (max AB for acceleration) - physics-based model
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Aerodynamics
    CL = W * cos(gamma) / (q * aero.S_ref);
    CL = max(0.1, min(1.2, CL));

    if aero.use_raymer && ~isempty(aero.geom)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL / CL_alpha;
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, aero.geom);
    else
        CD = aero.CD0 + aero.K * CL^2;
    end
    D = q * aero.S_ref * CD;

    % Hull's equations (fixed gamma)
    dxdt = V * cos(gamma);
    dhdt = V * sin(gamma);
    dVdt = (g / W) * (T - D - W * sin(gamma));
    dWdt = -fuel_flow * g;

    dydt = [dxdt; dhdt; dVdt; dWdt];
end

%% Equations of Motion - Level acceleration (Phase 1b)
function dydt = climb_EOM_level(t, y, aero, g)
    h = max(0, y(2));
    V = max(10, y(3));
    W = y(4);

    gamma = 0;  % Level flight

    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);
    M = V / a;
    q = 0.5 * rho * V^2;

    % Propulsion (max AB for acceleration) - physics-based model
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Aerodynamics (level flight: L = W)
    CL = W / (q * aero.S_ref);
    CL = max(0.1, min(1.2, CL));

    if aero.use_raymer && ~isempty(aero.geom)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL / CL_alpha;
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, aero.geom);
    else
        CD = aero.CD0 + aero.K * CL^2;
    end
    D = q * aero.S_ref * CD;

    % Hull's equations (level flight)
    dxdt = V;
    dhdt = 0;
    dVdt = (g / W) * (T - D);
    dWdt = -fuel_flow * g;

    dydt = [dxdt; dhdt; dVdt; dWdt];
end

%% Equations of Motion - Constant q climb (Phase 2)
function dydt = climb_EOM_const_q(t, y, q_target, aero, g)
    h = max(0, y(2));
    V = y(3);
    W = y(4);

    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);

    % For constant q: V = sqrt(2*q_target/rho)
    V_target = sqrt(2 * q_target / rho);
    M = V_target / a;
    q = q_target;

    % Constant q climb rate calculation
    H_scale = 8500;
    dV_dh = V_target / (2 * H_scale);

    % Aerodynamics
    CL_est = W / (q * aero.S_ref);
    CL_est = max(0.1, min(1.0, CL_est));

    if aero.use_raymer && ~isempty(aero.geom)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL_est / CL_alpha;
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, aero.geom);
    else
        CD = aero.CD0 + aero.K * CL_est^2;
    end
    D = q * aero.S_ref * CD;

    % Propulsion - MAX AB for max rate climb
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Solve for actual gamma based on available thrust
    T_minus_D = T - D;
    denominator = (V_target^2 / (2 * H_scale)) + g;
    sin_gamma = (g / W) * T_minus_D / denominator;
    sin_gamma = max(0.01, min(0.6, sin_gamma));
    gamma = asin(sin_gamma);

    % Kinematics
    dxdt = V_target * cos(gamma);
    dhdt = V_target * sin(gamma);
    dVdt = dV_dh * dhdt;
    dWdt = -fuel_flow * g;

    dydt = [dxdt; dhdt; dVdt; dWdt];
end

%% Equations of Motion - Constant M climb (Phase 3)
function dydt = climb_EOM_const_M(t, y, M_target, aero, g)
    h = max(0, y(2));
    V = y(3);
    W = y(4);

    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);

    % Maintain constant Mach
    V_target = M_target * a;
    q = 0.5 * rho * V_target^2;

    % For constant M: dV/dh = M * da/dh
    da_dh = -0.0033;
    dV_dh = M_target * da_dh;

    % Aerodynamics
    CL_est = W / (q * aero.S_ref);
    CL_est = max(0.1, min(1.2, CL_est));

    if aero.use_raymer && ~isempty(aero.geom)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M_target, 5, h, aero.geom);
        alpha_deg = CL_est / CL_alpha;
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M_target, alpha_deg, h, aero.geom);
    else
        CD = aero.CD0 + aero.K * CL_est^2;
    end
    D = q * aero.S_ref * CD;

    % Propulsion - MAX AB for max rate climb
    prop = get_F404_thrust_calibrated(M_target, h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Solve for actual gamma based on available thrust
    sin_gamma = (T - D) / W;
    sin_gamma = max(0, min(0.5, sin_gamma));
    gamma = asin(sin_gamma);

    % Kinematics
    dxdt = V_target * cos(gamma);
    dhdt = V_target * sin(gamma);
    dVdt = dV_dh * dhdt;
    dWdt = -fuel_flow * g;

    dydt = [dxdt; dhdt; dVdt; dWdt];
end

%% Helper function - Solve for PLA to achieve required thrust
function [PLA_sol, fuel_flow] = solve_PLA_for_thrust(T_required, M, h)
    % Binary search for PLA that gives required thrust
    PLA_min = 30;   % Idle
    PLA_max = 130;  % Max AB

    % Check bounds
    prop_min = get_F404_thrust_calibrated(max(0.3, M), h, PLA_min);
    prop_max = get_F404_thrust_calibrated(max(0.3, M), h, PLA_max);

    if T_required <= prop_min.thrust_N
        PLA_sol = PLA_min;
        fuel_flow = prop_min.fuel_flow_kg_s;
        return;
    elseif T_required >= prop_max.thrust_N
        PLA_sol = PLA_max;
        fuel_flow = prop_max.fuel_flow_kg_s;
        return;
    end

    % Binary search
    for iter = 1:20
        PLA_mid = (PLA_min + PLA_max) / 2;
        prop = get_F404_thrust_calibrated(max(0.3, M), h, PLA_mid);

        if prop.thrust_N < T_required
            PLA_min = PLA_mid;
        else
            PLA_max = PLA_mid;
        end

        if abs(prop.thrust_N - T_required) < 100  % 100 N tolerance
            break;
        end
    end

    PLA_sol = PLA_mid;
    fuel_flow = prop.fuel_flow_kg_s;
end
