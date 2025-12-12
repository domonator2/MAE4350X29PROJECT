function cruise = simulate_subsonic_cruise_efficient(M_cruise, h_cruise, W_initial, duration, S_ref, geom)
%SIMULATE_SUBSONIC_CRUISE_EFFICIENT Fuel-efficient subsonic cruise with L=W trim
%
% Solves for angle of attack (alpha) that satisfies L = W at each timestep.
% Uses get_X29_aero() Raymer model for realistic aerodynamics.
%
% Inputs:
%   M_cruise  - Cruise Mach number
%   h_cruise  - Cruise altitude [m]
%   W_initial - Initial weight [kg]
%   duration  - Cruise duration [s]
%   S_ref     - Wing reference area [m²]
%   geom      - Geometry structure
%
% Outputs:
%   cruise    - Structure with cruise results

g = 9.81;

fprintf('========================================\n');
fprintf('X-29 SUBSONIC CRUISE - FUEL EFFICIENT\n');
fprintf('========================================\n');
fprintf('Mach: %.2f, Altitude: %.0f m (%.0f ft)\n', M_cruise, h_cruise, h_cruise/0.3048);
fprintf('Initial weight: %.0f kg\n', W_initial);
fprintf('Duration: %.0f s (%.1f min)\n', duration, duration/60);
fprintf('----------------------------------------\n');

% Atmosphere
[T_atm, a, P, rho] = atmosisa(h_cruise);
V = M_cruise * a;
q = 0.5 * rho * V^2;

fprintf('Flight conditions:\n');
fprintf('  Temperature:   %.1f K\n', T_atm);
fprintf('  Density:       %.4f kg/m³\n', rho);
fprintf('  Velocity:      %.1f m/s (%.0f kts)\n', V, V * 1.944);
fprintf('  Dynamic press: %.0f Pa\n', q);
fprintf('----------------------------------------\n');

fprintf('Aerodynamics: Using get_X29_aero() Raymer model\n');

% Get initial aero estimate to display
alpha_init = 5.0;  % Initial guess
[CL_init, CD_init, ~, ~, ~] = get_X29_aero(M_cruise, alpha_init, h_cruise, geom);
LD_init = CL_init / CD_init;
fprintf('  Initial estimate (α=%.1f°): CL=%.3f, CD=%.4f, L/D=%.2f\n', ...
    alpha_init, CL_init, CD_init, LD_init);
fprintf('----------------------------------------\n');

% Simulation
dt = 1.0;  % 1 second timesteps
n_steps = ceil(duration / dt);

% Preallocate
t = zeros(n_steps, 1);
W = zeros(n_steps, 1);
CL = zeros(n_steps, 1);
CD = zeros(n_steps, 1);
LD = zeros(n_steps, 1);
alpha = zeros(n_steps, 1);
D = zeros(n_steps, 1);
T = zeros(n_steps, 1);
PLA = zeros(n_steps, 1);
fuel_flow = zeros(n_steps, 1);
fuel_burned = zeros(n_steps, 1);
TSFC = zeros(n_steps, 1);

% Initial conditions
W(1) = W_initial;
fuel_burned(1) = 0;

fprintf('Simulating cruise (solving L=W with Raymer aero at each step)...\n');
fprintf('t=   0s: W=%6.0fkg, ', W_initial);

for i = 1:n_steps
    t(i) = (i-1) * dt;

    % Required CL for L = W (trim condition)
    CL_required = W(i) * g / (q * S_ref);

    % Solve for alpha that gives CL_required using get_X29_aero
    [alpha(i), CL(i), CD(i)] = solve_alpha_for_CL(CL_required, M_cruise, h_cruise, geom);

    % Drag
    D(i) = q * S_ref * CD(i);

    % L/D
    LD(i) = CL(i) / CD(i);

    % Thrust required = Drag (steady level flight)
    T_required = D(i);
    T(i) = T_required;

    % Find PLA for required thrust using physics model
    [PLA(i), fuel_flow(i), TSFC(i)] = solve_PLA_for_thrust(T_required, M_cruise, h_cruise);

    % Fuel burn
    if i < n_steps
        fuel_burned(i+1) = fuel_burned(i) + fuel_flow(i) * dt;
        W(i+1) = W_initial - fuel_burned(i+1);
    end

    % Print progress
    if mod(i-1, 60) == 0 && i > 1
        fprintf('t=%4.0fs: W=%6.0fkg, α=%.2f°, CL=%.4f, CD=%.4f, L/D=%.2f, D=%.1fkN, PLA=%.1f°, FF=%.3fkg/s\n', ...
            t(i), W(i), alpha(i), CL(i), CD(i), LD(i), D(i)/1000, PLA(i), fuel_flow(i));
    end
end

% Final entry
i = n_steps;
fprintf('t=%4.0fs: W=%6.0fkg, α=%.2f°, CL=%.4f, CD=%.4f, L/D=%.2f, D=%.1fkN, PLA=%.1f°, FF=%.3fkg/s\n', ...
    t(i), W(i), alpha(i), CL(i), CD(i), LD(i), D(i)/1000, PLA(i), fuel_flow(i));

fprintf('----------------------------------------\n');
fprintf('SUBSONIC CRUISE COMPLETE\n');
fprintf('----------------------------------------\n');
fprintf('Final weight:     %.0f kg\n', W(end));
fprintf('Total fuel burned: %.1f kg\n', fuel_burned(end));
fprintf('Average fuel flow: %.3f kg/s (%.1f kg/min)\n', mean(fuel_flow), mean(fuel_flow)*60);
fprintf('Average L/D:       %.2f\n', mean(LD));
fprintf('Alpha range:       %.2f° to %.2f°\n', min(alpha), max(alpha));
fprintf('Specific range:    %.3f km/kg\n', V * duration / fuel_burned(end) / 1000);
fprintf('========================================\n');

% Package output
cruise = struct();
cruise.t = t;
cruise.W = W;
cruise.CL = CL;
cruise.CD = CD;
cruise.LD = LD;
cruise.alpha = alpha;
cruise.D = D;
cruise.T = T;
cruise.PLA = PLA;
cruise.fuel_flow = fuel_flow;
cruise.fuel_burned = fuel_burned;
cruise.TSFC = TSFC;
cruise.conditions.M = M_cruise;
cruise.conditions.h = h_cruise;
cruise.conditions.V = V;
cruise.conditions.q = q;

end

function [alpha_sol, CL_sol, CD_sol] = solve_alpha_for_CL(CL_target, M, h, geom)
%SOLVE_ALPHA_FOR_CL Find alpha that gives target CL using get_X29_aero
%
% Uses bisection method to find angle of attack

    % Search bounds
    alpha_min = 0;
    alpha_max = 15;

    % Limit CL_target to achievable range
    CL_target = max(0.1, min(CL_target, 1.2));

    % Bisection search
    for iter = 1:20
        alpha_mid = (alpha_min + alpha_max) / 2;
        [CL_mid, CD_mid, ~, ~, ~] = get_X29_aero(M, alpha_mid, h, geom);

        if CL_mid < CL_target
            alpha_min = alpha_mid;
        else
            alpha_max = alpha_mid;
        end

        if abs(CL_mid - CL_target) < 0.001
            break;
        end
    end

    alpha_sol = alpha_mid;
    [CL_sol, CD_sol, ~, ~, ~] = get_X29_aero(M, alpha_sol, h, geom);
end

function [PLA_sol, fuel_flow, TSFC_sol] = solve_PLA_for_thrust(T_required, M, h)
%SOLVE_PLA_FOR_THRUST Binary search to find PLA that gives required thrust

    % Search bounds
    PLA_min = 50;   % Minimum reasonable cruise PLA
    PLA_max = 100;  % Max dry power

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
    prop = get_F404_thrust_calibrated(M, h, PLA_sol);
    fuel_flow = prop.fuel_flow_kg_s;

    % TSFC = fuel_flow / thrust (convert to lb/hr/lbf)
    T_lbf = prop.thrust_N / 4.448;
    Wf_lb_hr = fuel_flow * 2.205 * 3600;
    TSFC_sol = Wf_lb_hr / max(T_lbf, 1);
end
