function climb = sim_climb_elevated(W_initial, S_ref, h_initial, V_initial, h_target, M_target, geom)
%SIM_CLIMB_ELEVATED Simulate X-29 climb from elevated airfield using Hull's equations
%
% Uses Hull's equations of motion with ODE45 numerical integration.
% 3-phase climb schedule optimized for fuel efficiency.
%
% Hull's Equations (vertical plane):
%   x' = V cos(gamma)
%   h' = V sin(gamma)
%   V' = (g/W)[T cos(alpha) - D - W sin(gamma)]
%   W' = -fuel_flow * g
%
% Climb Schedule:
%   Phase 1: Climb to safe altitude, then level accel to q = 21 kPa
%   Phase 2: Maintain constant q = 21 kPa, climb until M = 0.9
%   Phase 3: Maintain constant M = 0.9, climb until h = h_target
%
% Inputs:
%   W_initial - Initial weight [kg]
%   S_ref     - Wing reference area [m²]
%   h_initial - Starting altitude [m] (field elevation + 50ft AGL)
%   V_initial - Starting velocity [m/s] (liftoff speed)
%   h_target  - Target cruise altitude [m]
%   M_target  - Target cruise Mach number
%   geom      - Geometry structure for Raymer aero model
%
% Outputs:
%   climb     - Structure with climb results
%
% Author: Dominic Larin
% Date: December 2025

g = 9.81;

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║   X-29 CLIMB (Hull''s Equations + ODE45 + Elevated Start)      ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

%% Constants and targets
q_target = 21000;       % Pa - target dynamic pressure for efficient climb
M_transition = 0.9;     % Mach for phase 2->3 transition

% Altitude for level acceleration (Phase 1) - at least 500m above start
h_accel = max(h_initial + 500, 1000);

%% Aerodynamic parameters
aero_params.geom = geom;
aero_params.use_raymer = ~isempty(geom);
aero_params.S_ref = S_ref;
aero_params.CD0 = 0.025;
aero_params.AR = 4.0;
aero_params.e = 0.80;
aero_params.K = 1 / (pi * aero_params.AR * aero_params.e);

%% Initial conditions
[~, a_init, ~, ~] = atmosisa(h_initial);
M_init = V_initial / a_init;

fprintf('Hull''s Equations of Motion:\n');
fprintf('  x'' = V*cos(gamma)\n');
fprintf('  h'' = V*sin(gamma)\n');
fprintf('  V'' = (g/W)[T - D - W*sin(gamma)]\n');
fprintf('  W'' = -fuel_flow*g\n\n');

fprintf('CLIMB SCHEDULE:\n');
fprintf('  Phase 1: Climb to h=%.0fm, then accel to q = %.0f kPa\n', h_accel, q_target/1000);
fprintf('  Phase 2: Constant q = %.0f kPa climb until M = %.2f\n', q_target/1000, M_transition);
fprintf('  Phase 3: Constant M = %.2f climb until h = %.1f km\n\n', M_transition, h_target/1000);

fprintf('INITIAL CONDITIONS (Elevated Airfield):\n');
fprintf('  Weight:    %.0f kg\n', W_initial);
fprintf('  Altitude:  %.0f m (%.0f ft)\n', h_initial, h_initial/0.3048);
fprintf('  Velocity:  %.1f m/s (%.0f kts), M=%.2f\n\n', V_initial, V_initial*1.944, M_init);

%% ODE options
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1.0);

% Storage for all phases
all_t = [];
all_Y = [];

%% PHASE 1a: Initial climb to acceleration altitude
fprintf('--- PHASE 1a: Climb to h = %.0f m ---\n', h_accel);

% State: [x, h, V, W] - W in Newtons
y0_1a = [0; h_initial; V_initial; W_initial*g];
gamma_1a = deg2rad(10);  % Fixed 10 deg climb

% Event: stop when h >= h_accel
opts_1a = odeset(ode_opts, 'Events', @(t,y) event_altitude(t, y, h_accel));

[t_1a, Y_1a] = ode45(@(t,y) climb_EOM_fixed_gamma(t, y, gamma_1a, aero_params, g), ...
    [0, 300], y0_1a, opts_1a);

fprintf('  Completed: t=%.1fs, h=%.0fm, V=%.0fm/s\n', t_1a(end), Y_1a(end,2), Y_1a(end,3));

% Store
all_t = [all_t; t_1a];
all_Y = [all_Y; Y_1a];

%% PHASE 1b: Level acceleration to q = 21 kPa
fprintf('--- PHASE 1b: Level accel to q = %.0f kPa ---\n', q_target/1000);

y0_1b = Y_1a(end, :)';
t_offset_1b = t_1a(end);

% Event: stop when q >= q_target
opts_1b = odeset(ode_opts, 'Events', @(t,y) event_dynamic_pressure(t, y, q_target));

[t_1b, Y_1b] = ode45(@(t,y) climb_EOM_level(t, y, aero_params, g), ...
    [0, 180], y0_1b, opts_1b);

[~, a_1b, ~, rho_1b] = atmosisa(Y_1b(end,2));
M_1b = Y_1b(end,3) / a_1b;
q_1b = 0.5 * rho_1b * Y_1b(end,3)^2;

fprintf('  Completed: t=%.1fs, M=%.2f, q=%.1fkPa\n', t_offset_1b + t_1b(end), M_1b, q_1b/1000);

% Store (offset time)
all_t = [all_t; t_1b(2:end) + t_offset_1b];
all_Y = [all_Y; Y_1b(2:end, :)];

t_phase1 = t_offset_1b + t_1b(end);
fuel_phase1 = (W_initial*g - Y_1b(end,4)) / g;

%% PHASE 2: Constant q climb until M = 0.9
fprintf('--- PHASE 2: Constant q = %.0f kPa climb until M = %.2f ---\n', q_target/1000, M_transition);

y0_2 = Y_1b(end, :)';
t_offset_2 = t_phase1;

% Event: stop when M >= M_transition
opts_2 = odeset(ode_opts, 'Events', @(t,y) event_mach(t, y, M_transition));

[t_2, Y_2] = ode45(@(t,y) climb_EOM_const_q(t, y, q_target, aero_params, g), ...
    [0, 600], y0_2, opts_2);

[~, a_2, ~, ~] = atmosisa(Y_2(end,2));
M_2 = Y_2(end,3) / a_2;

fprintf('  Completed: t=%.1fs, h=%.0fm, M=%.2f\n', t_offset_2 + t_2(end), Y_2(end,2), M_2);

% Store
all_t = [all_t; t_2(2:end) + t_offset_2];
all_Y = [all_Y; Y_2(2:end, :)];

t_phase2 = t_offset_2 + t_2(end);
fuel_phase2 = (W_initial*g - Y_2(end,4)) / g;

%% PHASE 3: Constant M = 0.9 climb to cruise altitude
fprintf('--- PHASE 3: Constant M = %.2f climb to h = %.1f km ---\n', M_transition, h_target/1000);

y0_3 = Y_2(end, :)';
t_offset_3 = t_phase2;

% Event: stop when h >= h_target
opts_3 = odeset(ode_opts, 'Events', @(t,y) event_altitude(t, y, h_target));

[t_3, Y_3] = ode45(@(t,y) climb_EOM_const_M(t, y, M_transition, aero_params, g), ...
    [0, 600], y0_3, opts_3);

fprintf('  Completed: t=%.1fs, h=%.0fm\n', t_offset_3 + t_3(end), Y_3(end,2));

% Store
all_t = [all_t; t_3(2:end) + t_offset_3];
all_Y = [all_Y; Y_3(2:end, :)];

%% Final state
t_final = t_offset_3 + t_3(end);
h_final = Y_3(end, 2);
V_final = Y_3(end, 3);
W_final = Y_3(end, 4) / g;  % Convert back to kg
fuel_burned = W_initial - W_final;

[~, a_final, ~, ~] = atmosisa(h_final);
M_final = V_final / a_final;

%% Compute Mach history
n_points = length(all_t);
M_hist = zeros(n_points, 1);
gamma_hist = zeros(n_points, 1);
for i = 1:n_points
    h_i = all_Y(i, 2);
    V_i = all_Y(i, 3);
    [~, a_i, ~, ~] = atmosisa(h_i);
    M_hist(i) = V_i / a_i;

    % Estimate gamma from altitude rate
    if i > 1
        dh = all_Y(i,2) - all_Y(i-1,2);
        dt = all_t(i) - all_t(i-1);
        if dt > 0 && V_i > 10
            gamma_hist(i) = asind(min(1, max(-1, (dh/dt)/V_i)));
        end
    end
end

%% Summary
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('CLIMB SUMMARY (Hull''s Equations + ODE45)\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('Phase 1 (Climb + Accel to q=21kPa):\n');
fprintf('  Time:        %.1f s\n', t_phase1);
fprintf('  Fuel:        %.1f kg\n', fuel_phase1);
fprintf('\nPhase 2 (Constant q climb to M=0.9):\n');
fprintf('  Time:        %.1f s (%.1f s cumulative)\n', t_phase2-t_phase1, t_phase2);
fprintf('  Fuel:        %.1f kg (%.1f kg cumulative)\n', fuel_phase2-fuel_phase1, fuel_phase2);
fprintf('\nPhase 3 (Constant M climb to altitude):\n');
fprintf('  Time:        %.1f s (%.1f s cumulative)\n', t_final-t_phase2, t_final);
fprintf('  Fuel:        %.1f kg (%.1f kg cumulative)\n', fuel_burned-fuel_phase2, fuel_burned);
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('TOTALS:\n');
fprintf('  Duration:      %.1f s (%.1f min)\n', t_final, t_final/60);
fprintf('  Altitude:      %.0f m → %.0f m\n', h_initial, h_final);
fprintf('  Mach:          %.2f → %.2f\n', M_init, M_final);
fprintf('  Distance:      %.1f km (%.1f NM)\n', all_Y(end,1)/1000, all_Y(end,1)/1852);
fprintf('  Fuel burned:   %.1f kg\n', fuel_burned);
fprintf('  Final weight:  %.0f kg\n', W_final);
fprintf('════════════════════════════════════════════════════════════════\n');

%% Package output
climb = struct();
climb.t = all_t(:);
climb.h = all_Y(:, 2);
climb.V = all_Y(:, 3);
climb.M = M_hist(:);
climb.W = all_Y(:, 4) / g;  % Convert to kg
climb.x = all_Y(:, 1);
climb.gamma = gamma_hist(:);
climb.total_time = t_final;
climb.total_fuel = fuel_burned;
climb.W_initial = W_initial;
climb.W_final = W_final;
climb.distance = all_Y(end, 1);
climb.h_final = h_final;
climb.V_final = V_final;
climb.M_final = M_final;

end

%% Event functions
function [value, isterminal, direction] = event_altitude(~, y, h_target)
    value = y(2) - h_target;
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = event_dynamic_pressure(~, y, q_target)
    h = y(2);
    V = y(3);
    [~, ~, ~, rho] = atmosisa(h);
    q = 0.5 * rho * V^2;
    value = q - q_target;
    isterminal = 1;
    direction = 1;
end

function [value, isterminal, direction] = event_mach(~, y, M_target)
    h = y(2);
    V = y(3);
    [~, a, ~, ~] = atmosisa(h);
    M = V / a;
    value = M - M_target;
    isterminal = 1;
    direction = 1;
end

%% Equations of Motion - Fixed gamma climb (Phase 1a)
function dydt = climb_EOM_fixed_gamma(~, y, gamma, aero, g)
    h = max(0, y(2));
    V = max(10, y(3));
    W = y(4);  % Weight in Newtons

    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);
    M = V / a;
    q = 0.5 * rho * V^2;

    % Aerodynamics
    CL = W * cos(gamma) / (q * aero.S_ref);
    CL = max(0.1, min(1.2, CL));

    if aero.use_raymer && ~isempty(aero.geom)
        % Get CL_alpha from Raymer model (use nominal alpha=5 deg)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL / CL_alpha;  % Use proper CL_alpha from model
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, aero.geom);
    else
        CD = aero.CD0 + aero.K * CL^2;
    end
    D = q * aero.S_ref * CD;

    % Solve for required thrust: T = D + W*sin(gamma) (for steady climb)
    % Add small acceleration margin to ensure climb rate is achieved
    T_required = D + W * sin(gamma) + 0.1 * W;  % 0.1g acceleration margin
    [~, fuel_flow] = solve_PLA_for_thrust(T_required, max(0.3, M), h);

    % Hull's equations (fixed gamma) - use required thrust
    dxdt = V * cos(gamma);
    dhdt = V * sin(gamma);
    dVdt = (g / W) * (T_required - D - W * sin(gamma));
    dWdt = -fuel_flow * g;

    dydt = [dxdt; dhdt; dVdt; dWdt];
end

%% Equations of Motion - Level acceleration (Phase 1b)
function dydt = climb_EOM_level(~, y, aero, g)
    h = max(0, y(2));
    V = max(10, y(3));
    W = y(4);

    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);
    M = V / a;
    q = 0.5 * rho * V^2;

    % Propulsion (max AB for acceleration) - NASA TM-4140 calibrated model
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Aerodynamics (level flight: L = W)
    CL = W / (q * aero.S_ref);
    CL = max(0.1, min(1.2, CL));

    if aero.use_raymer && ~isempty(aero.geom)
        % Get CL_alpha from Raymer model (use nominal alpha=5 deg)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL / CL_alpha;  % Use proper CL_alpha from model
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
function dydt = climb_EOM_const_q(~, y, q_target, aero, g)
    h = max(0, y(2));
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
        % Get CL_alpha from Raymer model (use nominal alpha=5 deg)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL_est / CL_alpha;  % Use proper CL_alpha from model
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, aero.geom);
    else
        CD = aero.CD0 + aero.K * CL_est^2;
    end
    D = q * aero.S_ref * CD;

    % Propulsion - MAX AB for max rate climb (NASA TM-4140 calibrated)
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Solve for actual gamma based on available thrust
    % From T = D + W*sin(gamma): sin(gamma) = (T - D) / W
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
function dydt = climb_EOM_const_M(~, y, M_target, aero, g)
    h = max(0, y(2));
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
        % Get CL_alpha from Raymer model (use nominal alpha=5 deg)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M_target, 5, h, aero.geom);
        alpha_deg = CL_est / CL_alpha;  % Use proper CL_alpha from model
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

%% Throttle Solver - Binary search for PLA that gives required thrust
function [PLA_sol, fuel_flow] = solve_PLA_for_thrust(T_required, M, h)
%SOLVE_PLA_FOR_THRUST Binary search to find PLA that gives required thrust
% Extended range 30-130 to include afterburner for climb

    PLA_min = 30;   % Idle
    PLA_max = 130;  % Max AB

    % Handle edge cases
    prop_min = get_F404_thrust_calibrated(M, h, PLA_min);
    prop_max = get_F404_thrust_calibrated(M, h, PLA_max);

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
        prop = get_F404_thrust_calibrated(M, h, PLA_mid);
        T_mid = prop.thrust_N;

        if T_mid < T_required
            PLA_min = PLA_mid;
        else
            PLA_max = PLA_mid;
        end

        if abs(T_mid - T_required) < 100  % Within 100 N
            break;
        end
    end

    PLA_sol = PLA_mid;
    fuel_flow = prop.fuel_flow_kg_s;
end
