function [takeoff] = sim_takeoff_canards(params, geom)
%SIM_TAKEOFF_CANARDS Simulate X-29 takeoff with canard/rotation dynamics
%
% Simulates:
%   1. Ground roll (acceleration to rotation speed)
%   2. Rotation (pitch-up using canard)
%   3. Liftoff and initial climb to screen height (50 ft / 15 m)
%
% Tracks canard deflection, CG, and trim throughout.
%
% Inputs:
%   params.W_total  - Total aircraft weight [kg]
%   params.W_fuel   - Fuel weight [kg]
%   params.verbose  - Print detailed output
%   geom            - Geometry structure for Raymer aero model
%
% Outputs:
%   takeoff - Structure with time histories and final state
%
% Author: Dominic Larin
% Date: December 2025

%% Default parameters
if nargin < 1, params = struct(); end
if nargin < 2, geom = []; end

W_total = get_param(params, 'W_total', 8074);
W_fuel = get_param(params, 'W_fuel', 1877);
verbose = get_param(params, 'verbose', true);

%% Constants
g = 9.81;
rho_sl = 1.225;     % Sea level density [kg/m³]
ft2m = 0.3048;

%% Aircraft parameters
S_ref = 17.54;              % Wing area [m²]
S_canard = 3.4923;          % Canard area [m²]
c_bar = 8.0 * ft2m;         % MAC [m]

% Aerodynamic parameters - get CL_alpha from Raymer model at takeoff conditions
% Use M=0.25 (typical takeoff speed ~85 m/s) and alpha=5 deg as reference
if ~isempty(geom)
    [~, ~, ~, ~, ~, CL_alpha_raymer] = get_X29_aero(0.25, 5, 0, geom);
    CL_alpha_wing = CL_alpha_raymer;  % [1/deg] from Raymer model
    fprintf('  Using Raymer CL_alpha = %.4f [1/deg] for takeoff\n', CL_alpha_wing);
else
    CL_alpha_wing = 0.08;   % Fallback [1/deg]
    fprintf('  Using fallback CL_alpha = %.4f [1/deg] (no geometry provided)\n', CL_alpha_wing);
end
CL_0_wing = 0.0;             % Zero-alpha lift (Raymer model already accounts for this)
CL_alpha_canard = 0.0263;   % [1/deg] - canard contribution (for pitch control)
delta_CD_gear = 0.015;      % Gear-down drag increment (added to Raymer model)

% Flaperon and strake configuration for takeoff
% Uses apply_flap_effects() for consistent lift/drag increments
delta_flaperon_takeoff = 15;  % Takeoff flaperon setting [deg]
delta_strake_takeoff = 6;     % Takeoff strake setting [deg]

% Ground friction
mu_roll = 0.02;             % Rolling friction coefficient

%% Takeoff speeds
V_stall = 75;               % Approximate stall speed [m/s]
V_rotate = 1.1 * V_stall;   % Rotation speed [m/s] (~82 m/s)
V_liftoff = 1.15 * V_stall; % Liftoff speed [m/s] (~86 m/s)

%% Phase 1: Ground Roll
dt = 0.1;  % Time step [s]
t_max_ground = 60;  % Max ground roll time [s]

% Initialize arrays
n_max = round(t_max_ground / dt) + 1;
t_ground = zeros(n_max, 1);
V_ground = zeros(n_max, 1);
x_ground = zeros(n_max, 1);
W_ground = zeros(n_max, 1);
W_fuel_ground = zeros(n_max, 1);
delta_canard_ground = zeros(n_max, 1);

% Initial state
V_ground(1) = 0;
x_ground(1) = 0;
W_ground(1) = W_total;
W_fuel_ground(1) = W_fuel;
delta_canard_ground(1) = 0;  % Canard neutral during initial roll

% Get initial CG
[x_cg_init, ~] = calculate_X29_CG(W_total, W_fuel, 1877);

% Thrust (afterburner for takeoff) - calibrated physics model
prop_max = get_F404_thrust_calibrated(0, 0, 130);  % Sea level, M=0, max AB
T_max = prop_max.thrust_N;

i = 1;
reached_rotate = false;

while i < n_max && ~reached_rotate
    i = i + 1;
    t_ground(i) = t_ground(i-1) + dt;

    V = V_ground(i-1);
    W = W_ground(i-1);

    % Dynamic pressure and Mach
    q = 0.5 * rho_sl * V^2;
    M = V / 340.3;  % Approximate speed of sound at sea level

    % Aerodynamics (ground roll, alpha ~ 2 deg)
    alpha_ground = 2;

    % Get CL and CD from Raymer aero model, then apply flap effects
    if ~isempty(geom) && M > 0.05
        aero_clean = raymer_aero_model_elliptical(max(0.1, M), alpha_ground, 0, geom);
        % Apply flaperon + strake effects using apply_flap_effects
        aero_flapped = apply_flap_effects(aero_clean, delta_flaperon_takeoff, geom, delta_strake_takeoff);
        CL = aero_flapped.CL;
        CD = aero_flapped.CD + delta_CD_gear;  % Add gear drag
    else
        CL_clean = CL_0_wing + CL_alpha_wing * alpha_ground;
        CD_clean = 0.025 + 0.05 * CL_clean^2;  % Fallback drag estimate
        CL = CL_clean + 0.238;  % Fallback flap increment
        CD = CD_clean + delta_CD_gear + 0.015;  % Fallback flap drag
    end

    L = q * S_ref * CL;
    D = q * S_ref * CD;

    % Ground reaction
    N = max(0, W * g - L);

    % Friction force
    F_friction = mu_roll * N;

    % Thrust - physics-based model
    prop = get_F404_thrust_calibrated(max(0.3, M), 0, 130);  % Max AB for takeoff
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Net acceleration
    a = (T - D - F_friction) / (W);  % m/s²

    % Update velocity and position
    V_ground(i) = V + a * dt;
    x_ground(i) = x_ground(i-1) + V * dt + 0.5 * a * dt^2;

    % Fuel burn (fuel_flow is in kg/s from engine model)
    W_fuel_ground(i) = W_fuel_ground(i-1) - fuel_flow * dt;
    W_ground(i) = W_ground(i-1) - fuel_flow * dt;

    % Canard starts deflecting as we approach rotation
    if V_ground(i) > 0.9 * V_rotate
        % Ramp up canard deflection
        ramp = (V_ground(i) - 0.9*V_rotate) / (0.1*V_rotate);
        delta_canard_ground(i) = ramp * 25;  % Target ~25 deg at rotation
    else
        delta_canard_ground(i) = 0;
    end

    % Check for rotation speed
    if V_ground(i) >= V_rotate
        reached_rotate = true;
    end
end

n_ground = i;
t_ground = t_ground(1:n_ground);
V_ground = V_ground(1:n_ground);
x_ground = x_ground(1:n_ground);
W_ground = W_ground(1:n_ground);
W_fuel_ground = W_fuel_ground(1:n_ground);
delta_canard_ground = delta_canard_ground(1:n_ground);

%% Phase 2: Rotation (pitch-up to liftoff)
t_rotate_max = 5;  % Max rotation time [s]
n_rotate = round(t_rotate_max / dt);

t_rotate = zeros(n_rotate, 1);
V_rotate_arr = zeros(n_rotate, 1);
h_rotate = zeros(n_rotate, 1);
alpha_rotate = zeros(n_rotate, 1);
delta_canard_rotate = zeros(n_rotate, 1);
W_rotate = zeros(n_rotate, 1);
W_fuel_rotate = zeros(n_rotate, 1);

% Initial state from ground roll
V_rotate_arr(1) = V_ground(end);
h_rotate(1) = 0;
alpha_rotate(1) = 2;
delta_canard_rotate(1) = delta_canard_ground(end);
W_rotate(1) = W_ground(end);
W_fuel_rotate(1) = W_fuel_ground(end);

% Pitch rate during rotation
pitch_rate = 5;  % deg/s (typical rotation rate)

i = 1;
lifted_off = false;

while i < n_rotate && ~lifted_off
    i = i + 1;
    t_rotate(i) = t_rotate(i-1) + dt;

    V = V_rotate_arr(i-1);
    W = W_rotate(i-1);
    alpha = alpha_rotate(i-1);

    q = 0.5 * rho_sl * V^2;
    M = V / 340.3;

    % Increase alpha (rotation)
    alpha_new = min(alpha + pitch_rate * dt, 12);  % Limit to 12 deg

    % Get trim canard deflection for this alpha
    % Simplified: we're rotating, so canard provides nose-up moment
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M);

    % Canard deflection for rotation
    delta_canard_rotate(i) = 20 + 2 * (alpha_new - 2);  % Ramps up with alpha

    % Total lift - use Raymer aero model + apply_flap_effects
    if ~isempty(geom)
        aero_clean = raymer_aero_model_elliptical(max(0.1, M), alpha_new, 0, geom);
        aero_flapped = apply_flap_effects(aero_clean, delta_flaperon_takeoff, geom, delta_strake_takeoff);
        CL_wing = aero_flapped.CL;
        CL_wing_clean = aero_clean.CL;
        CD_flapped = aero_flapped.CD;
    else
        CL_wing_clean = CL_0_wing + CL_alpha_wing * alpha_new;
        CL_wing = CL_wing_clean + 0.238;  % Fallback flap increment
        CD_flapped = 0.025 + 0.05 * CL_wing^2 + 0.015;
    end

    CL_canard = CL_alpha_canard * delta_canard_rotate(i);
    L_wing = q * S_ref * CL_wing;
    L_canard = q * S_canard * CL_canard;
    L_total = L_wing + L_canard;

    % Drag - use flapped aero + gear
    CD = CD_flapped + delta_CD_gear;
    D = q * S_ref * CD;

    % Thrust - physics-based model
    prop = get_F404_thrust_calibrated(max(0.3, M), 0, 130);  % Max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Check for liftoff (require minimum alpha of 10 deg AND L >= W)
    alpha_min_liftoff = 10.0;  % Minimum rotation angle before liftoff [deg]
    if L_total >= W * g && alpha_new >= alpha_min_liftoff
        lifted_off = true;
        % Save liftoff conditions
        alpha_liftoff = alpha_new;
        V_liftoff_actual = V;
        M_liftoff = M;
        CL_liftoff = CL_total;
        CL_wing_liftoff = CL_wing;
        CL_clean_liftoff = CL_wing_clean;
        CD_liftoff = CD;
    end

    % Continue accelerating
    a = (T - D) / W - g * sind(0);  % Assume level for now
    V_rotate_arr(i) = V + a * dt;

    % Update state
    alpha_rotate(i) = alpha_new;
    W_fuel_rotate(i) = W_fuel_rotate(i-1) - fuel_flow * dt;
    W_rotate(i) = W_rotate(i-1) - fuel_flow * dt;
end

n_rotate_actual = i;
t_rotate = t_rotate(1:n_rotate_actual);
V_rotate_arr = V_rotate_arr(1:n_rotate_actual);
h_rotate = h_rotate(1:n_rotate_actual);
alpha_rotate = alpha_rotate(1:n_rotate_actual);
delta_canard_rotate = delta_canard_rotate(1:n_rotate_actual);
W_rotate = W_rotate(1:n_rotate_actual);
W_fuel_rotate = W_fuel_rotate(1:n_rotate_actual);

%% Phase 3: Initial Climb to Screen Height (50 ft = 15 m)
h_screen = 15;  % Screen height [m]
t_climb_max = 30;  % Max initial climb time [s]
n_climb = round(t_climb_max / dt);

t_climb = zeros(n_climb, 1);
V_climb = zeros(n_climb, 1);
h_climb = zeros(n_climb, 1);
gamma_climb = zeros(n_climb, 1);  % Flight path angle
alpha_climb = zeros(n_climb, 1);
delta_canard_climb = zeros(n_climb, 1);
W_climb = zeros(n_climb, 1);
W_fuel_climb = zeros(n_climb, 1);
CL_climb = zeros(n_climb, 1);
CD_climb = zeros(n_climb, 1);
T_climb = zeros(n_climb, 1);
D_climb = zeros(n_climb, 1);

% Initial state from rotation
V_climb(1) = V_rotate_arr(end);
h_climb(1) = 0;
gamma_climb(1) = 5;  % Initial climb angle [deg]
alpha_climb(1) = alpha_rotate(end);
delta_canard_climb(1) = delta_canard_rotate(end);
W_climb(1) = W_rotate(end);
W_fuel_climb(1) = W_fuel_rotate(end);

i = 1;
reached_screen = false;

while i < n_climb && ~reached_screen
    i = i + 1;
    t_climb(i) = t_climb(i-1) + dt;

    V = V_climb(i-1);
    h = h_climb(i-1);
    W = W_climb(i-1);
    gamma = gamma_climb(i-1);

    % Atmosphere (still essentially sea level)
    [~, a, ~, rho] = atmosisa(h);
    q = 0.5 * rho * V^2;
    M = V / a;

    % Get trim condition
    trim = calculate_trim_simple(V, h, W, W_fuel_climb(i-1));

    alpha_climb(i) = trim.alpha;
    delta_canard_climb(i) = trim.delta_canard;
    CL_climb(i) = trim.CL;
    CD_climb(i) = trim.CD;

    L = q * S_ref * CL_climb(i);
    D = q * S_ref * CD_climb(i);
    D_climb(i) = D;

    % Thrust - physics-based model
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % Max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;
    T_climb(i) = T;

    % Flight path dynamics
    % L = W * cos(gamma) approximately for small gamma
    % T - D = W * sin(gamma) + W/g * dV/dt

    % Climb angle
    sin_gamma = (T - D) / (W * g);
    sin_gamma = min(max(sin_gamma, 0), 0.5);  % Limit climb angle
    gamma_climb(i) = asind(sin_gamma);

    % Update position
    V_climb(i) = V + ((T - D) / W - g * sin_gamma) * dt;
    V_climb(i) = max(V_climb(i), V_liftoff);  % Don't slow below liftoff speed

    h_climb(i) = h + V * sind(gamma) * dt;

    % Fuel burn (fuel_flow is in kg/s from engine model)
    W_fuel_climb(i) = W_fuel_climb(i-1) - fuel_flow * dt;
    W_climb(i) = W_climb(i-1) - fuel_flow * dt;

    % Check for screen height
    if h_climb(i) >= h_screen
        reached_screen = true;
    end
end

n_climb_actual = i;
t_climb = t_climb(1:n_climb_actual);
V_climb = V_climb(1:n_climb_actual);
h_climb = h_climb(1:n_climb_actual);
gamma_climb = gamma_climb(1:n_climb_actual);
alpha_climb = alpha_climb(1:n_climb_actual);
delta_canard_climb = delta_canard_climb(1:n_climb_actual);
W_climb = W_climb(1:n_climb_actual);
W_fuel_climb = W_fuel_climb(1:n_climb_actual);
CL_climb = CL_climb(1:n_climb_actual);
CD_climb = CD_climb(1:n_climb_actual);
T_climb = T_climb(1:n_climb_actual);
D_climb = D_climb(1:n_climb_actual);

%% Combine all phases
% Time
t_total = [t_ground; t_ground(end) + t_rotate; t_ground(end) + t_rotate(end) + t_climb];

% Build output arrays with proper CG calculation
n_total = length(t_total);
takeoff = struct();
takeoff.t = t_total;
takeoff.h = [zeros(n_ground, 1); h_rotate; h_climb];
takeoff.V = [V_ground; V_rotate_arr; V_climb];

% Mach number
takeoff.M = zeros(n_total, 1);
for j = 1:n_total
    [~, a_local, ~, ~] = atmosisa(takeoff.h(j));
    takeoff.M(j) = takeoff.V(j) / a_local;
end

takeoff.W = [W_ground; W_rotate; W_climb];
takeoff.W_fuel = [W_fuel_ground; W_fuel_rotate; W_fuel_climb];

% CG for each point
takeoff.x_cg = zeros(n_total, 1);
for j = 1:n_total
    [takeoff.x_cg(j), ~] = calculate_X29_CG(takeoff.W(j), takeoff.W_fuel(j), 1877);
end

% Alpha (ground roll is fixed at ~2 deg)
takeoff.alpha = [2*ones(n_ground, 1); alpha_rotate; alpha_climb];

% Canard deflection
takeoff.delta_canard = [delta_canard_ground; delta_canard_rotate; delta_canard_climb];

% Aerodynamics (calculate for each point)
takeoff.CL = zeros(n_total, 1);
takeoff.CD = zeros(n_total, 1);
takeoff.LD = zeros(n_total, 1);
takeoff.T = zeros(n_total, 1);
takeoff.D = zeros(n_total, 1);

for j = 1:n_total
    [~, a_local, ~, rho] = atmosisa(takeoff.h(j));
    q = 0.5 * rho * takeoff.V(j)^2;
    M = takeoff.V(j) / a_local;

    % Lift and drag - use Raymer aero model + apply_flap_effects
    if ~isempty(geom) && M > 0.05
        aero_clean = raymer_aero_model_elliptical(max(0.1, M), takeoff.alpha(j), takeoff.h(j), geom);
        aero_flapped = apply_flap_effects(aero_clean, delta_flaperon_takeoff, geom, delta_strake_takeoff);
        CL_wing = aero_flapped.CL;
        CD_flapped = aero_flapped.CD;
    else
        CL_wing_clean = CL_0_wing + CL_alpha_wing * takeoff.alpha(j);
        CL_wing = CL_wing_clean + 0.238;  % Fallback flap increment
        CD_flapped = 0.025 + 0.05 * CL_wing^2 + 0.015;
    end

    CL_canard = CL_alpha_canard * takeoff.delta_canard(j);
    L_wing = q * S_ref * CL_wing;
    L_canard = q * S_canard * CL_canard;
    L_total = L_wing + L_canard;
    takeoff.CL(j) = L_total / (q * S_ref);

    % Drag - use flapped aero + gear
    takeoff.CD(j) = CD_flapped + delta_CD_gear;
    takeoff.D(j) = q * S_ref * takeoff.CD(j);

    % L/D
    if takeoff.D(j) > 0
        takeoff.LD(j) = L_total / takeoff.D(j);
    else
        takeoff.LD(j) = 0;
    end

    % Thrust - physics-based model
    prop = get_F404_thrust_calibrated(max(0.3, M), takeoff.h(j), 130);
    takeoff.T(j) = prop.thrust_N;
end

% Final state
takeoff.W_final = takeoff.W(end);
takeoff.W_fuel_final = takeoff.W_fuel(end);
takeoff.h_final = takeoff.h(end);
takeoff.V_final = takeoff.V(end);

takeoff.ground_roll_distance = x_ground(end);
takeoff.total_time = t_total(end);

% Flaperon and strake deflection (constant during takeoff)
takeoff.delta_flaperon = delta_flaperon_takeoff * ones(n_total, 1);
takeoff.delta_strake = delta_strake_takeoff * ones(n_total, 1);

% Liftoff aerodynamics (calculate from apply_flap_effects)
takeoff.alpha_liftoff = alpha_liftoff;
takeoff.V_liftoff = V_liftoff_actual;
takeoff.M_liftoff = M_liftoff;
% Calculate liftoff aero using apply_flap_effects for consistency
if ~isempty(geom)
    aero_liftoff_clean = raymer_aero_model_elliptical(max(0.1, M_liftoff), alpha_liftoff, 0, geom);
    aero_liftoff_flapped = apply_flap_effects(aero_liftoff_clean, delta_flaperon_takeoff, geom, delta_strake_takeoff);
    takeoff.CL_wing_liftoff = aero_liftoff_flapped.CL;
    takeoff.CL_clean_liftoff = aero_liftoff_clean.CL;
    takeoff.CD_liftoff = aero_liftoff_flapped.CD + delta_CD_gear;
    takeoff.delta_CL_flaperon = aero_liftoff_flapped.CL - aero_liftoff_clean.CL;  % Total flap increment
else
    takeoff.CL_wing_liftoff = CL_wing_liftoff;
    takeoff.CL_clean_liftoff = CL_clean_liftoff;
    takeoff.CD_liftoff = CD_liftoff;
    takeoff.delta_CL_flaperon = 0.238;
end
takeoff.CL_total_liftoff = CL_liftoff;          % Total CL including canard
takeoff.delta_CL_strake = 0;  % Included in delta_CL_flaperon via apply_flap_effects

if verbose
    fprintf('  Ground roll distance: %.0f m\n', takeoff.ground_roll_distance);
    fprintf('  Rotation speed: %.1f m/s (%.0f kts)\n', V_rotate, V_rotate*1.944);
    fprintf('  Liftoff speed: %.1f m/s (%.0f kts)\n', V_liftoff_actual, V_liftoff_actual*1.944);
    fprintf('  Max canard during rotation: %.1f deg\n', max(takeoff.delta_canard));
    fprintf('  Flaperon deflection: %.1f deg\n', delta_flaperon_takeoff);
    fprintf('  Strake deflection: %.1f deg\n', delta_strake_takeoff);
    fprintf('\n');
    fprintf('  LIFTOFF AERODYNAMICS (using apply_flap_effects):\n');
    fprintf('    Angle of attack:   %.1f deg\n', takeoff.alpha_liftoff);
    fprintf('    CL_clean (α):      %.3f\n', takeoff.CL_clean_liftoff);
    fprintf('    ΔCL (flap+strake): %.3f\n', takeoff.delta_CL_flaperon);
    fprintf('    ─────────────────────────\n');
    fprintf('    CL_wing:           %.3f\n', takeoff.CL_wing_liftoff);
    fprintf('    CD_total:          %.4f\n', takeoff.CD_liftoff);
    fprintf('    L/D at liftoff:    %.2f\n', takeoff.CL_wing_liftoff/takeoff.CD_liftoff);
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

function trim = calculate_trim_simple(V, h, W, W_fuel)
    % Simplified trim calculation for takeoff/landing
    % Returns alpha and delta_canard for L=W

    [~, a, ~, rho] = atmosisa(h);
    q = 0.5 * rho * V^2;
    M = V / a;
    g = 9.81;

    S_ref = 17.54;
    S_canard = 3.4923;
    CL_alpha_wing = 0.08;
    CL_0_wing = 0.1;
    CL_alpha_canard = 0.0263;

    % Required CL for level flight
    CL_required = W * g / (q * S_ref);

    % Simplified: assume wing provides most lift
    % alpha = (CL_required - CL_0) / CL_alpha
    alpha = (CL_required - CL_0_wing) / CL_alpha_wing;
    alpha = min(max(alpha, 0), 15);  % Limit

    % Canard deflection for trim (simplified)
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M);

    % Very rough: canard deflection proportional to alpha
    delta_canard = 3 * alpha;  % Empirical relationship
    delta_canard = min(max(delta_canard, 0), 58);

    % Output
    trim.alpha = alpha;
    trim.delta_canard = delta_canard;
    trim.CL = CL_required;
    trim.CD = 0.035 + 0.10 * CL_required^2;  % Rough
end
