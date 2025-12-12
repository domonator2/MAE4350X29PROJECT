function takeoff = sim_takeoff_elevated(W_initial, S_ref, h_field, geom)
%SIM_TAKEOFF_ELEVATED Simulate takeoff from elevated airfield
%
% Inputs:
%   W_initial - Initial weight [kg]
%   S_ref     - Wing reference area [m²]
%   h_field   - Field elevation [m]
%   geom      - Geometry structure
%
% Outputs:
%   takeoff   - Structure with takeoff results

if nargin < 4
    geom = [];
end

g = 9.81;

% Get atmospheric conditions at field elevation
[T_atm, a, P, rho] = atmosisa(h_field);

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║         X-29 ELEVATED TAKEOFF SIMULATION (ODE45)               ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('FIELD CONDITIONS:\n');
fprintf('  Elevation:     %.0f m (%.0f ft)\n', h_field, h_field/0.3048);
fprintf('  Temperature:   %.1f K (%.1f °C)\n', T_atm, T_atm - 273.15);
fprintf('  Pressure:      %.0f Pa\n', P);
fprintf('  Density:       %.4f kg/m³ (%.1f%% of sea level)\n', rho, 100*rho/1.225);
fprintf('  Speed of sound: %.1f m/s\n', a);
fprintf('\n');

% Aerodynamic parameters
CL_ground = 0.3;      % At ground roll alpha
mu = 0.04;            % Rolling friction

% Flaperon configuration for takeoff
delta_flaperon_takeoff = 10;  % Takeoff flaperon setting [deg]

% Get CL_max from Raymer model at alpha=12 (clean) + flaperon increment
[CL_max_clean, ~, ~, ~, ~, ~] = get_X29_aero(0.25, 12, h_field, geom);
[delta_CL_flap, ~, ~] = get_flaperon_effects(delta_flaperon_takeoff, 0.25, geom);
CL_max = CL_max_clean + delta_CL_flap;


% Calculate critical speeds (adjusted for density altitude)
V_stall = sqrt(2 * W_initial * g / (rho * S_ref * CL_max));
V_rotate = 1.05 * V_stall;  % Rotation speed
V_LOF = 1.1 * V_stall;      % Liftoff speed

fprintf('CRITICAL SPEEDS (density adjusted):\n');
fprintf('  V_stall:   %.1f m/s (%.0f kts)\n', V_stall, V_stall * 1.944);
fprintf('  V_rotate:  %.1f m/s (%.0f kts)\n', V_rotate, V_rotate * 1.944);
fprintf('  V_LOF:     %.1f m/s (%.0f kts)\n', V_LOF, V_LOF * 1.944);
fprintf('\n');

% Get thrust at field conditions (MAX AB)
PLA = 130;  % Maximum afterburner
prop = get_F404_thrust_calibrated(0.1, h_field, PLA);
T_SL = prop.thrust_N;

fprintf('THRUST (at field elevation):\n');
fprintf('  T_available: %.1f kN (%.0f lbf)\n', T_SL/1000, T_SL/4.448);
fprintf('  T/W ratio:   %.3f\n', T_SL / (W_initial * g));
fprintf('\n');

% Ground roll simulation
dt = 0.1;  % Time step
t_max = 60;  % Max time

t = 0;
x = 0;
V = 0;
W = W_initial;
alpha = 2;  % Ground roll alpha

t_arr = [0];
x_arr = [0];
V_arr = [0];
h_arr = [h_field];

% Ground roll phase
fprintf('PHASE 1: GROUND ROLL\n');
while V < V_LOF && t < t_max
    % Atmosphere at current altitude
    [~, ~, ~, rho_curr] = atmosisa(h_field);

    % Dynamic pressure
    q = 0.5 * rho_curr * V^2;

    % Get Mach
    M = V / a;

    % Thrust (adjusted for velocity)
    prop = get_F404_thrust_calibrated(M, h_field, PLA);
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Lift and drag during ground roll
    if V < V_rotate
        CL = CL_ground;
    else
        % Rotation: increase alpha gradually
        alpha = 2 + 8 * (V - V_rotate) / (V_LOF - V_rotate);
        alpha = min(alpha, 10);
        % Get CL from Raymer aero model
        [CL_clean, ~, ~, ~, ~, ~] = get_X29_aero(max(0.1, M), alpha, h_field, geom);
        % Add flaperon increment
        [dCL_flap, dCD_flap, ~] = get_flaperon_effects(delta_flaperon_takeoff, M, geom);
        CL = CL_clean + dCL_flap;
    end

    L = q * S_ref * CL;

    % Get drag from Raymer aero model + flaperon increment
    [~, CD_clean, ~, ~, ~] = get_X29_aero(max(0.1, M), alpha, h_field, geom);
    [~, dCD_flap, ~] = get_flaperon_effects(delta_flaperon_takeoff, M, geom);
    CD = CD_clean + dCD_flap;
    D = q * S_ref * CD;

    % Ground reaction
    N = max(0, W * g - L);
    F_friction = mu * N;

    % Acceleration
    a_x = (T - D - F_friction) / (W);

    % Update state
    V = V + a_x * dt;
    x = x + V * dt;
    t = t + dt;
    W = W - fuel_flow * dt;

    % Store
    t_arr(end+1) = t;
    x_arr(end+1) = x;
    V_arr(end+1) = V;
    h_arr(end+1) = h_field;
end

V_liftoff = V;
x_ground_roll = x;
t_ground_roll = t;
fuel_ground = W_initial - W;

fprintf('  Ground roll:    %.0f m (%.0f ft)\n', x_ground_roll, x_ground_roll/0.3048);
fprintf('  Time:           %.1f s\n', t_ground_roll);
fprintf('  Liftoff speed:  %.1f m/s (%.0f kts)\n', V_liftoff, V_liftoff * 1.944);
fprintf('  Fuel burned:    %.1f kg\n', fuel_ground);

% Climb to 50 ft (15.24 m) above field
h_screen = h_field + 15.24;
fprintf('\nPHASE 2: CLIMB TO 50 FT AGL\n');

alpha = 8;  % Climb alpha
gamma = 8;  % Climb angle

while h_arr(end) < h_screen && t < t_max
    h_curr = h_arr(end);
    [~, a_curr, ~, rho_curr] = atmosisa(h_curr);
    M = V / a_curr;
    q = 0.5 * rho_curr * V^2;

    % Thrust
    prop = get_F404_thrust_calibrated(M, h_curr, PLA);
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Aero - use Raymer model + flaperon
    [CL_clean, CD_clean, ~, ~, ~, ~] = get_X29_aero(max(0.1, M), alpha, h_curr, geom);
    [dCL_flap, dCD_flap, ~] = get_flaperon_effects(delta_flaperon_takeoff, M, geom);
    CL = CL_clean + dCL_flap;
    CD = CD_clean + dCD_flap;
    L = q * S_ref * CL;
    D = q * S_ref * CD;

    % Climb equations
    dh = V * sind(gamma) * dt;
    dx = V * cosd(gamma) * dt;

    % Update
    t = t + dt;
    x = x + dx;
    W = W - fuel_flow * dt;

    t_arr(end+1) = t;
    x_arr(end+1) = x;
    V_arr(end+1) = V;
    h_arr(end+1) = h_curr + dh;
end

fuel_total = W_initial - W;

fprintf('  Climb distance: %.0f m\n', x - x_ground_roll);
fprintf('  Climb time:     %.1f s\n', t - t_ground_roll);
fprintf('  Final altitude: %.0f m (%.0f ft AGL)\n', h_arr(end), (h_arr(end) - h_field)/0.3048);

% Package output
takeoff = struct();
takeoff.t = t_arr(:);
takeoff.x = x_arr(:);
takeoff.V = V_arr(:);
takeoff.h = h_arr(:);
takeoff.W_initial = W_initial;
takeoff.W_final = W;
takeoff.fuel_burned = fuel_total;
takeoff.total_time = t;
takeoff.total_distance = x;
takeoff.V_LOF = V_liftoff;
takeoff.ground_roll = x_ground_roll;
takeoff.delta_flaperon = delta_flaperon_takeoff;
takeoff.CL_max = CL_max;

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('TAKEOFF SUMMARY:\n');
fprintf('  Total distance:  %.0f m (%.0f ft)\n', x, x/0.3048);
fprintf('  Total time:      %.1f s\n', t);
fprintf('  Fuel burned:     %.1f kg\n', fuel_total);
fprintf('  Final weight:    %.0f kg\n', W);
fprintf('  Flaperon:        %.1f deg\n', delta_flaperon_takeoff);
fprintf('  CL_max:          %.2f\n', CL_max);
fprintf('════════════════════════════════════════════════════════════════\n');

end
