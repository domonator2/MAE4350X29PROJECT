function [landing] = sim_landing_canards(params, geom)
%SIM_LANDING_CANARDS Simulate X-29 landing using Raymer methodology
%
% Landing trajectory analysis following Raymer's methodology [44]:
%   Phase 1: Approach - From pattern altitude along 3° glideslope
%   Phase 2: Flare - Circular arc transition to touchdown
%   Phase 3: Ground Roll - Free roll + braking to full stop
%
% Uses Raymer-based aerodynamics model via get_X29_aero() with gear/flap increments.
%
% Key equations from Raymer:
%   V_stall = sqrt(2W / (rho * S_ref * CL_max))           [Eq. 52]
%   V_approach = 1.2 * V_stall (military)                 [Eq. 53]
%   V_touchdown = 1.1 * V_stall
%   R = V_F^2 / (g * (n-1)), n = 1.2                      [Eq. 55]
%   h_F = R * (1 - cos(gamma_a))                          [Eq. 56]
%   S_F = R * sin(gamma_a)                                [Eq. 57]
%   S_a = (h_obs - h_F) / tan(gamma_a)                    [Eq. 54]
%   S_G = (1/2gK_A) * ln((K_T + K_A*Vi^2)/(K_T + K_A*Vf^2)) [Eq. 58]
%
% Inputs:
%   params.W_total       - Total aircraft weight [kg]
%   params.W_fuel        - Fuel weight [kg]
%   params.h_pattern     - Pattern altitude [m] (default: 610 m = 2000 ft)
%   params.CL_max        - Maximum lift coefficient with flaps (default: 1.6)
%   params.verbose       - Print detailed output
%   geom                 - Geometry structure for Raymer aero model
%
% Outputs:
%   landing - Structure with time histories, distances, and final state
%
% Reference: Raymer, D.P., "Aircraft Design: A Conceptual Approach" [44]
%
% Author: Dominic Larin
% Date: December 2025

%% Default parameters
if nargin < 1, params = struct(); end
if nargin < 2, geom = []; end

W_total = get_param(params, 'W_total', 6800);       % Landing weight [kg]
W_fuel = get_param(params, 'W_fuel', 500);          % Remaining fuel [kg]
h_pattern = get_param(params, 'h_pattern', 610);    % Pattern altitude [m] (2000 ft)
CL_max = get_param(params, 'CL_max', 1.6);          % With flaperons deflected
verbose = get_param(params, 'verbose', true);

%% Physical Constants
g = 9.81;                   % Gravitational acceleration [m/s²]
rho_sl = 1.225;             % Sea level density [kg/m³]
ft2m = 0.3048;
m2ft = 1/ft2m;
ms2kts = 1.94384;

%% Landing Configuration Parameters (Table 30)
h_obstacle = 15.24;         % Obstacle height [m] (50 ft)
gamma_a = 3;                % Approach angle [deg] (standard ILS)
mu_brake = 0.3;             % Braking coefficient (military, dry concrete)
t_free_roll = 2;            % Free roll time before braking [s]
n_flare = 1.2;              % Load factor during flare

%% Aircraft Geometry
S_ref = 17.54;              % Wing reference area [m²]
S_canard = 3.4923;          % Canard area [m²]
c_bar = 8.0 * ft2m;         % Mean aerodynamic chord [m]

%% Aerodynamic Parameters (landing configuration)
CL_alpha_wing = 0.08;       % Lift curve slope [1/deg]
CL_0_landing = 0.3;         % Zero-alpha CL with flaps

% Drag increments for landing configuration
delta_CD_gear = 0.015;      % Landing gear drag increment

% Flaperon and strake configuration for landing
delta_flaperon_landing = 20;  % Landing flaperon setting [deg]
delta_strake_landing = 10;    % Landing strake setting [deg]

%% =========================================================================
%  PHASE 0: CALCULATE REFERENCE SPEEDS (Raymer Eqs. 52-53)
%  =========================================================================

% Calculate CL_max using apply_flap_effects for consistency
% Get clean aero at stall conditions
alpha_stall = 12;   % Stall angle [deg]
M_approach_est = 0.25;  % Estimated approach Mach

if ~isempty(geom)
    aero_stall_clean = raymer_aero_model_elliptical(M_approach_est, alpha_stall, 0, geom);
    aero_stall_flapped = apply_flap_effects(aero_stall_clean, delta_flaperon_landing, geom, delta_strake_landing);
    CL_max_clean = aero_stall_clean.CL;
    CL_max_used = aero_stall_flapped.CL;
    delta_CL_total = CL_max_used - CL_max_clean;
else
    CL_alpha = 0.0835;  % [1/deg] fallback
    CL_max_clean = CL_alpha * alpha_stall;
    delta_CL_total = 0.330;  % Fallback flap increment
    CL_max_used = CL_max_clean + delta_CL_total;
end

fprintf('  CL_max calculation (using apply_flap_effects):\n');
fprintf('    CL_clean (α=12°):  %.3f\n', CL_max_clean);
fprintf('    ΔCL (flap+strake): %.3f (flap=%.0f°, strake=%.0f°)\n', delta_CL_total, delta_flaperon_landing, delta_strake_landing);
fprintf('    CL_max:            %.3f\n', CL_max_used);

% Stall speed (Eq. 52)
V_stall = sqrt(2 * W_total * g / (rho_sl * S_ref * CL_max_used));

% Approach speed - 1.2 * Vstall for military (Eq. 53)
V_approach = 1.2 * V_stall;

% Touchdown speed - 1.1 * Vstall for military [44]
V_touchdown = 1.1 * V_stall;

% Average flare velocity - 1.15 * Vstall [44]
V_flare_avg = 1.15 * V_stall;

%% =========================================================================
%  PHASE 0B: CALCULATE FLARE GEOMETRY (Raymer Eqs. 55-57)
%  =========================================================================

% Flare radius from circular arc approximation (Eq. 55)
R_flare = V_flare_avg^2 / (g * (n_flare - 1));

% Flare initiation height (Eq. 56)
h_flare_init = R_flare * (1 - cosd(gamma_a));

% Horizontal flare distance (Eq. 57)
S_flare = R_flare * sind(gamma_a);

% Approach distance from obstacle to flare (Eq. 54)
S_approach = (h_obstacle - h_flare_init) / tand(gamma_a);

%% =========================================================================
%  PHASE 1: APPROACH (Pattern altitude to flare initiation)
%  =========================================================================

dt = 0.5;  % Time step [s]

% Approach time and distance from pattern altitude
h_approach_start = h_pattern;
h_approach_end = h_flare_init;

% Distance along glideslope from pattern to flare
dist_pattern_to_flare = (h_approach_start - h_approach_end) / sind(gamma_a);

% Time for approach phase
t_approach_total = dist_pattern_to_flare / V_approach;
n_approach = ceil(t_approach_total / dt) + 5;

% Initialize approach arrays
t_app = zeros(n_approach, 1);
h_app = zeros(n_approach, 1);
x_app = zeros(n_approach, 1);  % Horizontal position
V_app = zeros(n_approach, 1);
W_app = zeros(n_approach, 1);
W_fuel_app = zeros(n_approach, 1);
delta_canard_app = zeros(n_approach, 1);
alpha_app = zeros(n_approach, 1);
fuel_flow_app = zeros(n_approach, 1);

% Initial conditions
h_app(1) = h_approach_start;
x_app(1) = 0;
V_app(1) = V_approach;
W_app(1) = W_total;
W_fuel_app(1) = W_fuel;

i = 1;
reached_flare = false;

while i < n_approach && ~reached_flare
    i = i + 1;
    t_app(i) = t_app(i-1) + dt;

    h_curr = h_app(i-1);
    V_curr = V_app(i-1);
    W_curr = W_app(i-1);
    W_fuel_curr = W_fuel_app(i-1);

    % Atmosphere at current altitude
    [~, a_local, ~, rho] = atmosisa(h_curr);
    q = 0.5 * rho * V_curr^2;
    M_curr = V_curr / a_local;

    % Required CL for steady glide
    CL_required = W_curr * g * cosd(gamma_a) / (q * S_ref);

    % Estimate angle of attack - use CL_alpha from Raymer model if available
    if ~isempty(geom)
        [~, ~, ~, ~, ~, CL_alpha_raymer] = get_X29_aero(M_curr, 5, h_curr, geom);
        alpha_est = CL_required / CL_alpha_raymer;
    else
        alpha_est = (CL_required - CL_0_landing) / CL_alpha_wing;  % Fallback
    end
    alpha_est = max(0, min(alpha_est, 12));
    alpha_app(i) = alpha_est;

    % Calculate canard deflection for trim
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_curr);
    if abs(Cm_delta_c) > 0.01
        delta_deg = -Cm_alpha * (alpha_est * pi/180) / Cm_delta_c * 180/pi;
    else
        delta_deg = 15;
    end
    delta_canard_app(i) = max(-32, min(58, delta_deg));

    % Get idle thrust fuel flow from physics model
    % PLA = 30 for idle
    h_ft = h_curr * m2ft;
    if h_ft > 0
        prop = get_F404_thrust_calibrated(M_curr, h_curr, 30);
        fuel_flow_app(i) = prop.fuel_flow_kg_s;
    else
        fuel_flow_app(i) = 0.15;  % Approximate idle fuel flow [kg/s]
    end

    % Update position along glideslope
    dh = -V_curr * sind(gamma_a) * dt;
    dx = V_curr * cosd(gamma_a) * dt;

    h_app(i) = h_curr + dh;
    x_app(i) = x_app(i-1) + dx;
    V_app(i) = V_approach;  % Constant approach speed

    % Fuel consumption
    dm_fuel = fuel_flow_app(i) * dt;
    W_fuel_app(i) = max(0, W_fuel_curr - dm_fuel);
    W_app(i) = W_curr - dm_fuel;

    if h_app(i) <= h_flare_init
        reached_flare = true;
        h_app(i) = h_flare_init;
    end
end

% Trim arrays
n_app_actual = i;
t_app = t_app(1:n_app_actual);
h_app = h_app(1:n_app_actual);
x_app = x_app(1:n_app_actual);
V_app = V_app(1:n_app_actual);
W_app = W_app(1:n_app_actual);
W_fuel_app = W_fuel_app(1:n_app_actual);
alpha_app = alpha_app(1:n_app_actual);
delta_canard_app = delta_canard_app(1:n_app_actual);
fuel_flow_app = fuel_flow_app(1:n_app_actual);

%% =========================================================================
%  PHASE 2: FLARE (Circular arc transition to touchdown)
%  =========================================================================

% Time for flare (estimate from distance and average velocity)
t_flare_est = S_flare / V_flare_avg;
n_flare_pts = ceil(t_flare_est / dt) + 10;

% Initialize flare arrays
t_flare = zeros(n_flare_pts, 1);
h_flare_arr = zeros(n_flare_pts, 1);
x_flare = zeros(n_flare_pts, 1);
V_flare = zeros(n_flare_pts, 1);
alpha_flare = zeros(n_flare_pts, 1);
delta_canard_flare = zeros(n_flare_pts, 1);
W_flare = zeros(n_flare_pts, 1);
W_fuel_flare = zeros(n_flare_pts, 1);

% Initial conditions from approach end
h_flare_arr(1) = h_app(end);
x_flare(1) = x_app(end);
V_flare(1) = V_app(end);
W_flare(1) = W_app(end);
W_fuel_flare(1) = W_fuel_app(end);
alpha_flare(1) = alpha_app(end);
delta_canard_flare(1) = delta_canard_app(end);

% Flare along circular arc
% theta goes from gamma_a to 0
theta_flare = gamma_a;  % Current flight path angle

i = 1;
touched_down = false;

while i < n_flare_pts && ~touched_down
    i = i + 1;
    t_flare(i) = t_flare(i-1) + dt;

    h_curr = h_flare_arr(i-1);
    V_curr = V_flare(i-1);
    W_curr = W_flare(i-1);

    % Decelerate from approach to touchdown speed during flare
    decel_rate = (V_approach - V_touchdown) / t_flare_est;
    V_flare(i) = max(V_touchdown, V_curr - decel_rate * dt);

    % Flight path angle decreases along circular arc
    % Rate of change: d(theta)/dt = V/R
    d_theta = (V_curr / R_flare) * (180/pi) * dt;
    theta_flare = max(0, theta_flare - d_theta);

    % Update position
    dh = -V_curr * sind(theta_flare) * dt;
    dx = V_curr * cosd(theta_flare) * dt;

    h_flare_arr(i) = max(0, h_curr + dh);
    x_flare(i) = x_flare(i-1) + dx;

    % Alpha increases during flare (pitch-up maneuver)
    alpha_flare(i) = min(alpha_flare(i-1) + 0.5 * dt, 10);

    % Canard deflection for trim
    M_curr = V_curr / 340;
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_curr);
    if abs(Cm_delta_c) > 0.01
        delta_deg = -Cm_alpha * (alpha_flare(i) * pi/180) / Cm_delta_c * 180/pi;
    else
        delta_deg = 20;
    end
    delta_canard_flare(i) = max(-32, min(58, delta_deg));

    % Minimal fuel burn during flare (idle thrust)
    W_fuel_flare(i) = W_fuel_flare(i-1) - 0.15 * dt;
    W_flare(i) = W_flare(i-1) - 0.15 * dt;

    if h_flare_arr(i) <= 0
        touched_down = true;
        h_flare_arr(i) = 0;
    end
end

% Trim arrays
n_flare_actual = i;
t_flare = t_flare(1:n_flare_actual);
h_flare_arr = h_flare_arr(1:n_flare_actual);
x_flare = x_flare(1:n_flare_actual);
V_flare = V_flare(1:n_flare_actual);
alpha_flare = alpha_flare(1:n_flare_actual);
delta_canard_flare = delta_canard_flare(1:n_flare_actual);
W_flare = W_flare(1:n_flare_actual);
W_fuel_flare = W_fuel_flare(1:n_flare_actual);

%% =========================================================================
%  PHASE 3: GROUND ROLL (Free roll + Braking)
%  =========================================================================

% Landing weight at touchdown
W_TD = W_flare(end);
V_TD = V_flare(end);

% Ground roll parameters - use apply_flap_effects for consistency
CL_ground = 0.3;            % Reduced CL on ground
% Get base CD from Raymer model at low speed + gear + flaperon using apply_flap_effects
if ~isempty(geom)
    aero_ground_clean = raymer_aero_model_elliptical(0.15, 2, 0, geom);
    aero_ground_flapped = apply_flap_effects(aero_ground_clean, delta_flaperon_landing, geom, delta_strake_landing);
    CD_ground = aero_ground_flapped.CD + delta_CD_gear;  % CD with gear + flaperon
else
    CD_ground = 0.035 + delta_CD_gear + 0.015;  % Fallback
end

%% Phase 3A: Free Roll (before brakes applied)
S_free_roll = V_TD * t_free_roll;

%% Phase 3B: Braking Distance using Raymer's method (Eqs. 58-60)

% Calculate K_T and K_A coefficients
% K_T = T/W - mu (Eq. 59)
% At idle, assume T/W ≈ 0.02 (very low idle thrust)
T_W_idle = 0.02;
K_T = T_W_idle - mu_brake;  % Will be negative (braking dominates)

% K_A = (rho/2) * (1/(W/S)) * (mu*CL - CD) (Eq. 60)
% Note: CD_ground already includes induced drag from Raymer model
W_S = W_TD * g / S_ref;     % Wing loading [N/m²]
K_A = (rho_sl / 2) * (1 / W_S) * (mu_brake * CL_ground - CD_ground);

% Braking distance from V_TD to 0 using analytical formula (Eq. 58)
% S_G = (1 / 2gK_A) * ln((K_T + K_A*Vi^2) / (K_T + K_A*Vf^2))
if abs(K_A) > 1e-10
    V_i = V_TD;     % Initial velocity (after free roll)
    V_f = 0;        % Final velocity (full stop)

    numerator = K_T + K_A * V_i^2;
    denominator = K_T + K_A * V_f^2;

    if numerator > 0 && denominator > 0
        S_brake = (1 / (2 * g * K_A)) * log(numerator / denominator);
    else
        % Fallback to simplified method if K_A method fails
        a_avg = mu_brake * g;
        S_brake = V_TD^2 / (2 * a_avg);
    end
else
    % Simplified braking distance (Eq. 70 in report)
    a_avg = mu_brake * g;
    S_brake = V_TD^2 / (2 * a_avg);
end

% Total ground roll
S_ground_roll = S_free_roll + S_brake;

%% Simulate ground roll with time history

% Estimate total ground roll time
t_free = t_free_roll;
t_brake_est = V_TD / (mu_brake * g);  % Rough estimate
t_roll_total = t_free + t_brake_est + 5;
n_roll = ceil(t_roll_total / dt);

% Initialize ground roll arrays
t_roll = zeros(n_roll, 1);
h_roll = zeros(n_roll, 1);
x_roll = zeros(n_roll, 1);
V_roll = zeros(n_roll, 1);
alpha_roll = zeros(n_roll, 1);
delta_canard_roll = zeros(n_roll, 1);
W_roll = zeros(n_roll, 1);
W_fuel_roll = zeros(n_roll, 1);

% Initial conditions
V_roll(1) = V_TD;
x_roll(1) = x_flare(end);
W_roll(1) = W_TD;
W_fuel_roll(1) = W_fuel_flare(end);
alpha_roll(1) = 2;          % Ground attitude
delta_canard_roll(1) = 0;   % Neutral on ground

i = 1;
stopped = false;
brake_applied = false;
t_since_TD = 0;

while i < n_roll && ~stopped
    i = i + 1;
    t_roll(i) = t_roll(i-1) + dt;
    t_since_TD = t_roll(i);

    V_curr = V_roll(i-1);
    W_curr = W_roll(i-1);
    x_curr = x_roll(i-1);

    % Dynamic pressure
    q = 0.5 * rho_sl * V_curr^2;

    % Aerodynamic forces (CD_ground already includes full drag from Raymer model)
    L = q * S_ref * CL_ground;
    D = q * S_ref * CD_ground;

    % Normal force (weight minus lift)
    N = max(0, W_curr * g - L);

    % Check if brakes applied (after free roll time)
    if t_since_TD >= t_free_roll
        brake_applied = true;
    end

    % Calculate deceleration
    if brake_applied
        F_brake = mu_brake * N;
    else
        F_brake = 0;  % Free rolling
    end

    % Total retarding force
    F_retard = F_brake + D;

    % Deceleration
    a_decel = -F_retard / W_curr;

    % Update velocity and position
    V_roll(i) = max(0, V_curr + a_decel * dt);
    x_roll(i) = x_curr + V_curr * dt + 0.5 * a_decel * dt^2;

    % Ground attitude
    alpha_roll(i) = 2;
    delta_canard_roll(i) = 0;

    % No fuel burn on ground (engines at idle, negligible)
    W_fuel_roll(i) = W_fuel_roll(i-1);
    W_roll(i) = W_roll(i-1);

    if V_roll(i) <= 0.1
        stopped = true;
        V_roll(i) = 0;
    end
end

% Trim arrays
n_roll_actual = i;
t_roll = t_roll(1:n_roll_actual);
h_roll = h_roll(1:n_roll_actual);
x_roll = x_roll(1:n_roll_actual);
V_roll = V_roll(1:n_roll_actual);
alpha_roll = alpha_roll(1:n_roll_actual);
delta_canard_roll = delta_canard_roll(1:n_roll_actual);
W_roll = W_roll(1:n_roll_actual);
W_fuel_roll = W_fuel_roll(1:n_roll_actual);

%% =========================================================================
%  COMBINE ALL PHASES AND PACKAGE OUTPUT
%  =========================================================================

% Time history (continuous time)
t_total = [t_app;
           t_app(end) + t_flare;
           t_app(end) + t_flare(end) + t_roll];

h_total = [h_app; h_flare_arr; h_roll];
x_total = [x_app; x_flare; x_roll];
V_total = [V_app; V_flare; V_roll];
alpha_total = [alpha_app; alpha_flare; alpha_roll];
delta_total = [delta_canard_app; delta_canard_flare; delta_canard_roll];
W_total_arr = [W_app; W_flare; W_roll];
W_fuel_total = [W_fuel_app; W_fuel_flare; W_fuel_roll];

n_total = length(t_total);

%% Package output structure
landing = struct();

% Time histories
landing.t = t_total;
landing.h = h_total;
landing.x = x_total;  % Horizontal position [m]
landing.V = V_total;

% Calculate Mach for each point
landing.M = zeros(n_total, 1);
for j = 1:n_total
    [~, a_local, ~, ~] = atmosisa(max(0, landing.h(j)));
    landing.M(j) = landing.V(j) / a_local;
end

landing.W = W_total_arr;
landing.W_fuel = W_fuel_total;

% CG tracking
landing.x_cg = zeros(n_total, 1);
for j = 1:n_total
    [landing.x_cg(j), ~] = calculate_X29_CG(landing.W(j), landing.W_fuel(j), 1877);
end

landing.alpha = alpha_total;
landing.delta_canard = delta_total;

% Aerodynamics
landing.CL = zeros(n_total, 1);
landing.CD = zeros(n_total, 1);
landing.LD = zeros(n_total, 1);
landing.T = zeros(n_total, 1);
landing.D = zeros(n_total, 1);

for j = 1:n_total
    [~, a_local, ~, rho] = atmosisa(max(0, landing.h(j)));
    q = 0.5 * rho * landing.V(j)^2;
    M_curr = landing.V(j) / a_local;

    if landing.h(j) > 0 && ~isempty(geom)
        % Use apply_flap_effects for consistency
        aero_clean = raymer_aero_model_elliptical(max(0.1, M_curr), landing.alpha(j), landing.h(j), geom);
        aero_flapped = apply_flap_effects(aero_clean, delta_flaperon_landing, geom, delta_strake_landing);
        CL = aero_flapped.CL;
        CD = aero_flapped.CD + delta_CD_gear;
    else
        CL = CL_ground;
        CD = CD_ground;  % Already includes gear + flaperon drag
    end

    landing.CL(j) = CL;
    landing.CD(j) = CD;
    landing.D(j) = q * S_ref * CD;
    if landing.D(j) > 0
        landing.LD(j) = q * S_ref * CL / landing.D(j);
    end
    landing.T(j) = 0;  % Idle thrust
end

% Flaperon deflection (constant during landing)
landing.delta_flaperon = delta_flaperon_landing * ones(n_total, 1);
landing.CL_max_used = CL_max_used;

%% Landing Performance Summary (Raymer methodology)

% Reference speeds
landing.V_stall = V_stall;
landing.V_approach = V_approach;
landing.V_touchdown = V_touchdown;
landing.V_flare_avg = V_flare_avg;

% Flare geometry
landing.R_flare = R_flare;
landing.h_flare_init = h_flare_init;

% Distances (analytical from Raymer)
landing.S_approach_from_obstacle = S_approach;  % From 50 ft obstacle to flare
landing.S_flare = S_flare;
landing.S_free_roll = S_free_roll;
landing.S_brake = S_brake;
landing.S_ground_roll = S_ground_roll;
landing.S_total_from_obstacle = S_approach + S_flare + S_ground_roll;

% Distance from pattern altitude (simulated)
landing.S_approach_from_pattern = x_app(end);
landing.S_total_from_pattern = x_roll(end);

% Times
landing.t_approach = t_app(end);
landing.t_flare = t_flare(end);
landing.t_ground_roll = t_roll(end);
landing.total_time = t_total(end);

% Fuel consumed
landing.fuel_consumed_approach = W_fuel - W_fuel_app(end);
landing.fuel_consumed_total = W_fuel - W_fuel_roll(end);

% Final weights
landing.W_final = landing.W(end);
landing.W_fuel_final = landing.W_fuel(end);

%% Display Results (matching Table 31 format)
if verbose
    fprintf('\n');
    fprintf('==========================================================\n');
    fprintf('         X-29 LANDING PERFORMANCE (Raymer Method)         \n');
    fprintf('==========================================================\n');
    fprintf('\n--- Reference Speeds ---\n');
    fprintf('  Stall Speed:        %6.1f m/s  (%5.0f kts)\n', V_stall, V_stall*ms2kts);
    fprintf('  Approach Speed:     %6.1f m/s  (%5.0f kts)  [1.2 Vstall]\n', V_approach, V_approach*ms2kts);
    fprintf('  Touchdown Speed:    %6.1f m/s  (%5.0f kts)  [1.1 Vstall]\n', V_touchdown, V_touchdown*ms2kts);
    fprintf('  Flare Avg Speed:    %6.1f m/s  (%5.0f kts)  [1.15 Vstall]\n', V_flare_avg, V_flare_avg*ms2kts);

    fprintf('\n--- Flare Geometry ---\n');
    fprintf('  Flare Radius:       %6.0f m   (%6.0f ft)\n', R_flare, R_flare*m2ft);
    fprintf('  Flare Height:       %6.1f m   (%6.1f ft)\n', h_flare_init, h_flare_init*m2ft);

    fprintf('\n--- Landing Distances (from 50 ft obstacle) ---\n');
    fprintf('  Approach Distance:  %6.0f m   (%6.0f ft)\n', S_approach, S_approach*m2ft);
    fprintf('  Flare Distance:     %6.0f m   (%6.0f ft)\n', S_flare, S_flare*m2ft);
    fprintf('  Free Roll Distance: %6.0f m   (%6.0f ft)\n', S_free_roll, S_free_roll*m2ft);
    fprintf('  Braking Distance:   %6.0f m   (%6.0f ft)\n', S_brake, S_brake*m2ft);
    fprintf('  Total Ground Roll:  %6.0f m   (%6.0f ft)\n', S_ground_roll, S_ground_roll*m2ft);
    fprintf('  ------------------------------------------\n');
    fprintf('  TOTAL LANDING DIST: %6.0f m   (%6.0f ft)\n', ...
        landing.S_total_from_obstacle, landing.S_total_from_obstacle*m2ft);

    fprintf('\n--- Times ---\n');
    fprintf('  Approach Time:      %6.1f s   (%5.2f min)\n', landing.t_approach, landing.t_approach/60);
    fprintf('  Flare Time:         %6.1f s\n', landing.t_flare);
    fprintf('  Ground Roll Time:   %6.1f s\n', landing.t_ground_roll);
    fprintf('  Total Landing Time: %6.1f s   (%5.2f min)\n', landing.total_time, landing.total_time/60);

    fprintf('\n--- Fuel & Weight ---\n');
    fprintf('  Landing Weight:     %6.0f kg  (%6.0f lb)\n', W_total, W_total*2.205);
    fprintf('  Fuel at Touchdown:  %6.1f kg  (%6.0f lb)\n', W_fuel_flare(end), W_fuel_flare(end)*2.205);
    fprintf('  Fuel Consumed:      %6.1f kg  (%6.1f lb)\n', landing.fuel_consumed_total, landing.fuel_consumed_total*2.205);

    fprintf('\n--- Canard Deflection ---\n');
    fprintf('  Max During Approach: %5.1f deg\n', max(delta_canard_app));
    fprintf('  Max During Flare:    %5.1f deg\n', max(delta_canard_flare));

    fprintf('\n--- Flaperon ---\n');
    fprintf('  Flaperon deflection: %5.1f deg\n', delta_flaperon_landing);
    fprintf('  CL_max (with flap):  %5.2f\n', CL_max_used);

    fprintf('==========================================================\n');
end

end

%% =========================================================================
%  HELPER FUNCTIONS
%  =========================================================================

function val = get_param(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
