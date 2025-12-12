function landing = sim_landing_elevated(W_initial, S_ref, h_pattern, V_pattern, h_field, geom)
%SIM_LANDING_ELEVATED Simulate landing at elevated airfield
%
% Uses Raymer methodology for approach, flare, and ground roll
%
% Inputs:
%   W_initial - Initial weight [kg]
%   S_ref     - Wing reference area [m²]
%   h_pattern - Pattern altitude [m]
%   V_pattern - Pattern speed [m/s]
%   h_field   - Field elevation [m]
%   geom      - Geometry structure
%
% Outputs:
%   landing   - Structure with landing results

g = 9.81;

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║         X-29 ELEVATED LANDING SIMULATION (ODE45)               ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Field conditions
[T_field, a_field, P_field, rho_field] = atmosisa(h_field);

fprintf('FIELD CONDITIONS:\n');
fprintf('  Elevation:     %.0f m (%.0f ft)\n', h_field, h_field/0.3048);
fprintf('  Density:       %.4f kg/m³\n', rho_field);
fprintf('  Pattern alt:   %.0f m (%.0f ft AGL)\n', h_pattern, (h_pattern - h_field)/0.3048);
fprintf('  Pattern speed: %.0f m/s (%.0f kts)\n', V_pattern, V_pattern * 1.944);
fprintf('\n');

% Landing configuration
delta_CD_gear = 0.015;   % Gear down drag increment

% Flaperon configuration for landing
delta_flaperon_landing = 20;  % Landing flaperon setting [deg]

% Get CL_max from Raymer model at alpha=12 (clean) + flaperon increment
[CL_max_clean, ~, ~, ~, ~, ~] = get_X29_aero(0.2, 12, h_field, geom);
[delta_CL_flap, delta_CD_flap, ~] = get_flaperon_effects(delta_flaperon_landing, 0.2, geom);
CL_max_land = CL_max_clean + delta_CL_flap;

% Critical speeds (at field density)
V_stall = sqrt(2 * W_initial * g / (rho_field * S_ref * CL_max_land));
V_approach = 1.2 * V_stall;
V_touchdown = 1.1 * V_stall;
V_min = 1.05 * V_stall;  % Minimum safe speed

fprintf('CRITICAL SPEEDS (at field elevation):\n');
fprintf('  V_stall:      %.1f m/s (%.0f kts)\n', V_stall, V_stall * 1.944);
fprintf('  V_approach:   %.1f m/s (%.0f kts)\n', V_approach, V_approach * 1.944);
fprintf('  V_touchdown:  %.1f m/s (%.0f kts)\n', V_touchdown, V_touchdown * 1.944);
fprintf('\n');

% Approach parameters
gamma_approach = -3;  % Glideslope angle [deg]
h_flare = 15;         % Flare initiation height AGL [m]
h_touchdown = h_field;

% Simulation
dt = 0.5;
t_max = 600;  % Max 10 minutes

% Initialize - ensure starting conditions are safe
t = 0;
h = h_pattern;
V = max(V_pattern, V_approach);  % Don't start slower than approach speed
x = 0;
W = W_initial;

fprintf('Starting approach at V=%.1f m/s\n', V);

% Storage arrays
t_arr = [t];
h_arr = [h];
V_arr = [V];
x_arr = [x];
W_arr = [W];

% PHASE 1: Final Approach
fprintf('\nPHASE 1: FINAL APPROACH\n');
h_agl = h - h_field;

while h_agl > h_flare && t < t_max
    [~, a_curr, ~, rho_curr] = atmosisa(h);
    M = V / a_curr;
    q = 0.5 * rho_curr * V^2;

    % Idle thrust
    prop = get_F404_thrust_calibrated(M, h, 35);
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Aero on glideslope
    CL = W * g * cosd(gamma_approach) / (q * S_ref);
    CL = min(CL, CL_max_land);
    % Get CL_alpha from Raymer model and estimate alpha
    [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, geom);
    alpha_est = CL / CL_alpha;  % Use proper CL_alpha from model
    alpha_est = max(0, min(alpha_est, 15));

    % Get drag from Raymer model + gear + flaperon increments
    [~, CD_clean, ~, ~, ~, ~] = get_X29_aero(M, alpha_est, h, geom);
    [~, dCD_flap, ~] = get_flaperon_effects(delta_flaperon_landing, M, geom);
    CD = CD_clean + delta_CD_gear + dCD_flap;
    D = q * S_ref * CD;

    % Flight path - descend at constant glideslope
    dh_dt = V * sind(gamma_approach);  % Negative (descending)
    dx_dt = V * cosd(gamma_approach);  % Positive (forward)

    % Speed control - try to maintain approach speed
    % Adjust pitch/drag to maintain target speed
    V_target = V_approach;
    dV_dt = -0.5 * (V - V_target);  % Proportional control to target speed
    dV_dt = max(-2, min(dV_dt, 2));  % Limit rate

    % Ensure we don't go below minimum safe speed
    if V + dV_dt * dt < V_min
        dV_dt = 0;
    end

    % Update
    V = V + dV_dt * dt;
    V = max(V, V_min);  % Enforce minimum speed
    h = h + dh_dt * dt;
    x = x + dx_dt * dt;
    W = W - fuel_flow * dt;
    t = t + dt;

    h_agl = h - h_field;

    % Store
    t_arr(end+1) = t;
    h_arr(end+1) = h;
    V_arr(end+1) = V;
    x_arr(end+1) = x;
    W_arr(end+1) = W;
end

fprintf('  Flare altitude reached at t=%.1fs, x=%.0f m\n', t, x);
x_approach = x;

% PHASE 2: Flare
fprintf('\nPHASE 2: FLARE\n');
V_flare = V;

while h > h_touchdown + 0.5 && t < t_max
    [~, a_curr, ~, rho_curr] = atmosisa(h);
    q = 0.5 * rho_curr * V^2;

    prop = get_F404_thrust_calibrated(V/a_curr, h, 35);
    fuel_flow = prop.fuel_flow_kg_s;

    % Pitch up during flare - reduce sink rate
    CL = min(CL_max_land * 0.9, W * g / (q * S_ref) * 1.1);
    M_flare = V / a_curr;
    % Get CL_alpha from Raymer model and estimate alpha
    [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M_flare, 5, h, geom);
    alpha_est = CL / CL_alpha;  % Use proper CL_alpha from model
    alpha_est = max(0, min(alpha_est, 15));

    % Get drag from Raymer model + gear + flaperon increments
    [~, CD_clean, ~, ~, ~, ~] = get_X29_aero(M_flare, alpha_est, h, geom);
    [~, dCD_flap, ~] = get_flaperon_effects(delta_flaperon_landing, M_flare, geom);
    CD = CD_clean + delta_CD_gear + dCD_flap;
    L = q * S_ref * CL;
    D = q * S_ref * CD;

    % Flare arc - gradually reduce descent rate
    sink_rate = -2;  % Target sink rate at touchdown [m/s]
    gamma_flare = asind(sink_rate / max(V, V_min));

    dh_dt = V * sind(gamma_flare);
    dx_dt = V * cosd(gamma_flare);

    % Slight deceleration during flare
    dV_dt = -0.5;
    if V + dV_dt * dt < V_touchdown
        dV_dt = 0;  % Don't go below touchdown speed
    end

    V = V + dV_dt * dt;
    V = max(V, V_touchdown * 0.95);  % Don't go below touchdown speed
    h = h + dh_dt * dt;
    x = x + dx_dt * dt;
    W = W - fuel_flow * dt;
    t = t + dt;

    t_arr(end+1) = t;
    h_arr(end+1) = h;
    V_arr(end+1) = V;
    x_arr(end+1) = x;
    W_arr(end+1) = W;
end

V_td = V;
fprintf('  Touchdown at t=%.1fs, V=%.1f m/s (%.0f kts)\n', t, V_td, V_td * 1.944);
x_flare = x - x_approach;

% PHASE 3: Ground Roll
fprintf('\nPHASE 3: GROUND ROLL\n');

mu_brake = 0.3;  % Braking friction coefficient
t_free_roll = 2; % Seconds before brakes applied

t_brake_start = t + t_free_roll;

while V > 2 && t < t_max
    [~, a_curr, ~, rho_curr] = atmosisa(h_field);
    q = 0.5 * rho_curr * V^2;

    % Nose down, spoilers deployed
    CL = 0.1;
    CD = 0.08;  % Spoilers + aero braking
    L = q * S_ref * CL;
    D = q * S_ref * CD;

    % Braking
    if t > t_brake_start
        N = W * g - L;
        F_brake = mu_brake * N;
    else
        F_brake = 0;
    end

    % Idle thrust
    prop = get_F404_thrust_calibrated(V/a_curr, h_field, 35);
    fuel_flow = prop.fuel_flow_kg_s;

    % Deceleration
    a_decel = -(D + F_brake) / W;
    V = V + a_decel * dt;
    V = max(V, 0);
    x = x + V * dt;
    W = W - fuel_flow * dt;
    t = t + dt;

    t_arr(end+1) = t;
    h_arr(end+1) = h_field;
    V_arr(end+1) = V;
    x_arr(end+1) = x;
    W_arr(end+1) = W;
end

x_rollout = x - x_approach - x_flare;
fuel_total = W_initial - W;

fprintf('  Ground roll: %.0f m (%.0f ft)\n', x_rollout, x_rollout/0.3048);

% Package output
landing = struct();
landing.t = t_arr(:);
landing.h = h_arr(:);
landing.V = V_arr(:);
landing.x = x_arr(:);
landing.W = W_arr(:);
landing.total_time = t;
landing.total_distance = x;
landing.total_fuel = fuel_total;
landing.W_initial = W_initial;
landing.W_final = W;
landing.V_touchdown = V_td;
landing.approach_distance = x_approach;
landing.flare_distance = x_flare;
landing.rollout_distance = x_rollout;
landing.delta_flaperon = delta_flaperon_landing;
landing.CL_max = CL_max_land;

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('LANDING SUMMARY:\n');
fprintf('  Approach:       %.0f m\n', x_approach);
fprintf('  Flare:          %.0f m\n', x_flare);
fprintf('  Ground roll:    %.0f m (%.0f ft)\n', x_rollout, x_rollout/0.3048);
fprintf('  Total distance: %.0f m (%.0f ft)\n', x, x/0.3048);
fprintf('  Total time:     %.1f s\n', t);
fprintf('  Fuel burned:    %.1f kg\n', fuel_total);
fprintf('  Final weight:   %.0f kg\n', W);
fprintf('  Flaperon:       %.1f deg\n', delta_flaperon_landing);
fprintf('  CL_max:         %.2f\n', CL_max_land);
fprintf('════════════════════════════════════════════════════════════════\n');

end
