function descent = sim_descent_to_pattern(W_initial, S_ref, h_start, V_start, h_target, geom)
%SIM_DESCENT_TO_PATTERN Simulate descent from cruise to pattern altitude
%
% Descends from subsonic cruise to specified pattern altitude.
% Uses part power with controlled descent rate.
% Maintains reasonable approach speed for landing.
%
% Inputs:
%   W_initial - Initial weight [kg]
%   S_ref     - Wing reference area [m²]
%   h_start   - Starting altitude [m]
%   V_start   - Starting velocity [m/s]
%   h_target  - Target pattern altitude [m] (should be AGL + field elevation)
%   geom      - Geometry structure
%
% Outputs:
%   descent   - Structure with descent results

g = 9.81;
PLA_descent = 50;   % Part power for controlled descent
V_approach = 90;    % Target approach speed [m/s] (~175 kts)

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║         DESCENT TO PATTERN (ODE45 Integration)                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

[~, a_start, ~, ~] = atmosisa(h_start);
M_start = V_start / a_start;

fprintf('DESCENT PROFILE:\n');
fprintf('  Start:  h=%.0f m (%.0f ft), M=%.2f, V=%.0f m/s\n', ...
    h_start, h_start/0.3048, M_start, V_start);
fprintf('  Target: h=%.0f m (%.0f ft)\n', h_target, h_target/0.3048);
fprintf('  Target V_approach: %.0f m/s (%.0f kts)\n', V_approach, V_approach * 1.944);
fprintf('\n');

% Simulation parameters
dt = 1.0;
t_max = 1200;  % Max 20 minutes

% Initialize
t = 0;
h = h_start;
V = V_start;
W = W_initial;
x = 0;

% Storage
t_arr = [0];
h_arr = [h_start];
V_arr = [V_start];
M_arr = [M_start];
W_arr = [W_initial];
x_arr = [0];

% Descent angle
gamma = -4;  % degrees (negative = descending)

fprintf('SIMULATING DESCENT...\n');

while h > h_target && t < t_max
    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);
    M = V / a;
    q = 0.5 * rho * V^2;

    % Adjust power based on speed - want to decelerate toward V_approach
    % As we get lower, slow down to approach speed
    h_frac = (h - h_target) / (h_start - h_target);  % 1 at start, 0 at target
    V_target = V_approach + h_frac * (V_start - V_approach);  % Linear blend

    if V > V_target * 1.1
        PLA_curr = 35;  % Idle to slow down
    elseif V < V_target * 0.9
        PLA_curr = 70;  % More power to speed up
    else
        PLA_curr = PLA_descent;
    end

    prop = get_F404_thrust_calibrated(M, h, PLA_curr);
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % CL for descending flight
    CL_req = W * g * cosd(gamma) / (q * S_ref);
    CL_req = min(CL_req, 1.0);

    % Get CL_alpha from Raymer model and estimate alpha
    [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, geom);
    alpha_est = CL_req / CL_alpha;  % Use proper CL_alpha from model
    alpha_est = max(0, min(alpha_est, 15));

    % Get drag from Raymer aero model
    [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_est, h, geom);
    D = q * S_ref * CD;

    % Flight path
    dh_dt = V * sind(gamma);
    dx_dt = V * cosd(gamma);

    % Speed control
    dV_dt = (T - D - W * g * sind(gamma)) / W;

    % Apply speed target tracking
    V_error = V - V_target;
    dV_dt = dV_dt - 0.1 * V_error;  % Proportional control

    dV_dt = max(-3, min(dV_dt, 3));  % Limit speed change

    % Update
    t = t + dt;
    V = V + dV_dt * dt;
    V = max(V, V_approach * 0.9);  % Don't go below 90% approach speed
    h = h + dh_dt * dt;
    x = x + dx_dt * dt;
    W = W - fuel_flow * dt;

    % Limit
    if h < h_target
        h = h_target;
    end

    % Store
    t_arr(end+1) = t;
    h_arr(end+1) = h;
    V_arr(end+1) = V;
    [~, a_new, ~, ~] = atmosisa(h);
    M_arr(end+1) = V / a_new;
    W_arr(end+1) = W;
    x_arr(end+1) = x;

    % Progress
    if mod(t, 120) < dt
        fprintf('  t=%4.0fs: h=%5.0f m, M=%.2f, V=%4.0f m/s, V_target=%.0f m/s\n', ...
            t, h, M_arr(end), V, V_target);
    end
end

fuel_burned = W_initial - W;

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('DESCENT SUMMARY:\n');
fprintf('  Duration:     %.1f s (%.1f min)\n', t, t/60);
fprintf('  Distance:     %.1f km (%.1f NM)\n', x/1000, x/1852);
fprintf('  Alt change:   %.0f m to %.0f m\n', h_start, h);
fprintf('  Mach change:  %.2f to %.2f\n', M_start, M_arr(end));
fprintf('  Final speed:  %.0f m/s (%.0f kts) - target was %.0f m/s\n', V, V*1.944, V_approach);
fprintf('  Fuel burned:  %.1f kg\n', fuel_burned);
fprintf('  Final weight: %.0f kg\n', W);
fprintf('════════════════════════════════════════════════════════════════\n');

% Package output
descent = struct();
descent.t = t_arr(:);
descent.h = h_arr(:);
descent.V = V_arr(:);
descent.M = M_arr(:);
descent.W = W_arr(:);
descent.x = x_arr(:);
descent.total_time = t;
descent.total_fuel = fuel_burned;
descent.W_initial = W_initial;
descent.W_final = W;
descent.distance = x;
descent.h_final = h;
descent.V_final = V;

end
