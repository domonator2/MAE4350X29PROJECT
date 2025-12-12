function descent = sim_descent_to_subsonic(W_initial, S_ref, h_start, V_start, h_target, M_target, geom)
%SIM_DESCENT_TO_SUBSONIC Descend and decelerate from supersonic to subsonic cruise
%
% Two-phase descent:
%   Phase 1: Descending deceleration from M1.4 to M0.9 (slight descent)
%   Phase 2: Continue descent to target altitude while decelerating to M_target
%
% Inputs:
%   W_initial - Initial weight [kg]
%   S_ref     - Wing reference area [m²]
%   h_start   - Starting altitude [m]
%   V_start   - Starting velocity [m/s]
%   h_target  - Target altitude [m]
%   M_target  - Target Mach number
%   geom      - Geometry structure
%
% Outputs:
%   descent   - Structure with descent results

g = 9.81;
PLA_idle = 35;  % Idle power for descent

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║      DESCENT TO SUBSONIC CRUISE (ODE45 Integration)            ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

[~, a_start, ~, ~] = atmosisa(h_start);
M_start = V_start / a_start;

[~, a_target, ~, ~] = atmosisa(h_target);
V_target = M_target * a_target;

fprintf('DESCENT PROFILE:\n');
fprintf('  Start:  h=%.0f m, M=%.2f, V=%.0f m/s\n', h_start, M_start, V_start);
fprintf('  Target: h=%.0f m, M=%.2f, V=%.0f m/s\n', h_target, M_target, V_target);
fprintf('\n');

% Simulation parameters
dt = 1.0;
t_max = 600;  % 10 min max

% Initialize
t = 0;
h = h_start;
V = V_start;
W = W_initial;
x = 0;  % Horizontal distance

% Storage
t_arr = [0];
h_arr = [h_start];
V_arr = [V_start];
M_arr = [M_start];
W_arr = [W_initial];
x_arr = [0];

% Descent path angle - gradual descent
gamma = -3;  % degrees (negative = descending)

fprintf('SIMULATING DESCENT...\n');

while h > h_target && t < t_max
    % Atmosphere
    [~, a, ~, rho] = atmosisa(h);
    M = V / a;
    q = 0.5 * rho * V^2;

    % Propulsion at idle
    prop = get_F404_thrust_calibrated(M, h, PLA_idle);
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Aerodynamics - estimate CL for descending flight
    CL_req = W * g * cosd(gamma) / (q * S_ref);
    CL_req = min(CL_req, 1.0);

    % Get CL_alpha from Raymer model and estimate alpha
    [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, geom);
    alpha_est = CL_req / CL_alpha;  % Use proper CL_alpha from model
    alpha_est = max(0, min(alpha_est, 15));

    % Get drag from Raymer aero model
    [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_est, h, geom);

    D = q * S_ref * CD;
    L = q * S_ref * CL_req;

    % Equations of motion for descending flight
    % dV/dt = (T*cos(alpha) - D - W*g*sin(gamma)) / (W/g)
    % For descent, we want to decelerate, so we adjust gamma

    % If still supersonic, decelerate faster with steeper descent
    if M > 0.95
        gamma = -2;  % Shallow descent, more drag
    elseif M > M_target + 0.1
        gamma = -3;  % Moderate descent
    else
        gamma = -4;  % Steeper descent to lose altitude faster
    end

    % Flight path
    dV_dt = (T - D - W * g * sind(gamma)) / W;
    dh_dt = V * sind(gamma);
    dx_dt = V * cosd(gamma);

    % Limit deceleration rate
    dV_dt = max(dV_dt, -5);  % Don't decelerate too fast

    % Update state
    V = V + dV_dt * dt;
    h = h + dh_dt * dt;
    x = x + dx_dt * dt;
    W = W - fuel_flow * dt;
    t = t + dt;

    % Limit V to target
    if V < V_target * 0.95
        V = V_target * 0.95;
    end

    % Store
    t_arr(end+1) = t;
    h_arr(end+1) = h;
    V_arr(end+1) = V;
    [~, a_curr, ~, ~] = atmosisa(h);
    M_arr(end+1) = V / a_curr;
    W_arr(end+1) = W;
    x_arr(end+1) = x;

    % Progress
    if mod(t, 60) < dt
        fprintf('  t=%.0fs: h=%.0f m, M=%.2f, V=%.0f m/s\n', t, h, M_arr(end), V);
    end
end

fuel_burned = W_initial - W;

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('DESCENT TO SUBSONIC SUMMARY:\n');
fprintf('  Duration:     %.1f s (%.1f min)\n', t, t/60);
fprintf('  Distance:     %.1f km (%.1f NM)\n', x/1000, x/1852);
fprintf('  Alt change:   %.0f m to %.0f m\n', h_start, h);
fprintf('  Mach change:  %.2f to %.2f\n', M_start, M_arr(end));
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
descent.distance = x;
descent.W_initial = W_initial;
descent.W_final = W;
descent.h_final = h;
descent.V_final = V;

end
