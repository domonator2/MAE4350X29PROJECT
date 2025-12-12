function accel = sim_accel(W_initial, S_ref, h, M_start, M_target, geom)
%SIM_ACCEL Simulate transonic acceleration at constant altitude
%
% Accelerates from M_start to M_target at constant altitude h.
% If already at or above M_target, returns minimal fuel burn passthrough.
%
% Inputs:
%   W_initial - Initial weight [kg]
%   S_ref     - Wing reference area [m²]
%   h         - Altitude [m]
%   M_start   - Starting Mach number
%   M_target  - Target Mach number
%   geom      - Geometry structure
%
% Outputs:
%   accel     - Structure with acceleration results

g = 9.81;

fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║      TRANSONIC ACCELERATION (Constant Altitude)               ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

% Atmosphere
[~, a, ~, rho] = atmosisa(h);
V_start = M_start * a;
V_target = M_target * a;

fprintf('ACCELERATION PROFILE:\n');
fprintf('  Altitude: %.0f m (%.0f ft)\n', h, h/0.3048);
fprintf('  Start:    M=%.2f, V=%.0f m/s\n', M_start, V_start);
fprintf('  Target:   M=%.2f, V=%.0f m/s\n', M_target, V_target);
fprintf('\n');

% Check if already at target
if M_start >= M_target * 0.99
    fprintf('  Already at target Mach - minimal acceleration needed.\n');

    % Just burn a tiny bit of fuel for ~10 seconds
    prop = get_F404_thrust_calibrated(M_start, h, 100);  % Cruise power
    fuel_burned = prop.fuel_flow_kg_s * 10;

    accel = struct();
    accel.t = [0; 10];
    accel.h = [h; h];
    accel.V = [V_start; V_start];
    accel.M = [M_start; M_start];
    accel.W = [W_initial; W_initial - fuel_burned];
    accel.x = [0; V_start * 10];
    accel.total_time = 10;
    accel.total_fuel = fuel_burned;
    accel.W_initial = W_initial;
    accel.W_final = W_initial - fuel_burned;
    accel.distance = V_start * 10;

    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('ACCELERATION COMPLETE (already at target)\n');
    fprintf('════════════════════════════════════════════════════════════════\n');
    return;
end

% Simulation parameters
dt = 0.5;
t_max = 600;  % Max 10 minutes for transonic accel (let it take as long as needed)

% Initialize
t = 0;
V = V_start;
W = W_initial;
x = 0;

% Storage
t_arr = [0];
V_arr = [V_start];
M_arr = [M_start];
W_arr = [W_initial];
x_arr = [0];

fprintf('SIMULATING ACCELERATION...\n');

while V < V_target && t < t_max
    % Current Mach
    M = V / a;
    q = 0.5 * rho * V^2;

    % Max AB thrust
    PLA = 130;
    prop = get_F404_thrust_calibrated(M, h, PLA);
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Level flight CL
    CL = W * g / (q * S_ref);
    CL = min(CL, 1.0);

    % Get CL_alpha from Raymer model and estimate alpha
    [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, geom);
    alpha_est = CL / CL_alpha;  % Use proper CL_alpha from model
    alpha_est = max(0, min(alpha_est, 15));

    % Get drag from Raymer aero model
    [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_est, h, geom);
    D = q * S_ref * CD;

    % Acceleration (level flight)
    dV_dt = (T - D) / W;
    dV_dt = max(0, min(dV_dt, 10));  % Limit

    % Update
    t = t + dt;
    V = V + dV_dt * dt;
    x = x + V * dt;
    W = W - fuel_flow * dt;

    % Limit to target
    if V > V_target
        V = V_target;
    end

    % Store
    t_arr(end+1) = t;
    V_arr(end+1) = V;
    M_arr(end+1) = V / a;
    W_arr(end+1) = W;
    x_arr(end+1) = x;

    % Progress
    if mod(t, 30) < dt
        fprintf('  t=%3.0fs: M=%.2f, V=%4.0f m/s, dV/dt=%.1f m/s²\n', t, V/a, V, dV_dt);
    end
end

fuel_burned = W_initial - W;

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('ACCELERATION SUMMARY:\n');
fprintf('  Duration:     %.1f s\n', t);
fprintf('  Mach:         %.2f → %.2f\n', M_start, V/a);
fprintf('  Distance:     %.1f km (%.1f NM)\n', x/1000, x/1852);
fprintf('  Fuel burned:  %.1f kg\n', fuel_burned);
fprintf('  Final weight: %.0f kg\n', W);
fprintf('════════════════════════════════════════════════════════════════\n');

% Package output
accel = struct();
accel.t = t_arr(:);
accel.h = h * ones(size(t_arr(:)));
accel.V = V_arr(:);
accel.M = M_arr(:);
accel.W = W_arr(:);
accel.x = x_arr(:);
accel.total_time = t;
accel.total_fuel = fuel_burned;
accel.W_initial = W_initial;
accel.W_final = W;
accel.distance = x;

end
