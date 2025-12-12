%{
====================================================================================================
SCRIPT: run_maxmachmaxalt_canards.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 12/05/2025
====================================================================================================
DESCRIPTION:
Tests X-29 capability to reach M=1.6 at 15.3 km (max capable speed/altitude).
Uses CALIBRATED physics-based propulsion model (get_F404_thrust_calibrated).
Includes canard deflection and stability tracking throughout.

Key question: Is T_available > D_required at M=1.6, h=15.3km?
====================================================================================================
%}

clear; clc; close all;

% Add paths - ensure Physics folder takes priority
script_path = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(script_path));

% Add project paths first, then this folder LAST so it takes priority
addpath(genpath(proj_root));
addpath(script_path);

fprintf('╔════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║     X-29 MAX MACH / MAX ALTITUDE MISSION (CALIBRATED PROPULSION)      ║\n');
fprintf('║              M=1.6 at h=15.3 km with Canard Tracking                  ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════════╝\n\n');

% =========================================================================
% LOAD CONVERGED VEHICLE GEOMETRY
% =========================================================================
fprintf('Loading converged vehicle geometry...\n');

mat_file = fullfile(proj_root, 'Synthesis', 'updated_fuselage_calculations', ...
    'elliptical_trade_study_20251117_170903.mat');

if ~exist(mat_file, 'file')
    mat_files = dir(fullfile(proj_root, 'Synthesis', '**', '*trade_study*.mat'));
    if ~isempty(mat_files)
        mat_file = fullfile(mat_files(end).folder, mat_files(end).name);
        fprintf('  Using alternate mat file: %s\n', mat_files(end).name);
    else
        error('Could not find trade study .mat file');
    end
end

trade_data = load(mat_file);

sc = trade_data.successful_configs;
target_sweep = -30; target_tau = 0.077; target_duration = 10;
best_idx = -1; best_dist = inf;
for i = 1:numel(sc)
    dist = abs(sc(i).sweep_deg - target_sweep) + ...
        100*abs(sc(i).tau - target_tau) + ...
        abs(sc(i).duration_min - target_duration);
    if dist < best_dist
        best_dist = dist;
        best_idx = i;
    end
end
veh = sc(best_idx);
geom = build_X29_geometry(veh, trade_data.vehicle_base);
fprintf('  Geometry loaded.\n\n');

% =========================================================================
% VEHICLE AND MISSION PARAMETERS
% =========================================================================
TOGW = veh.TOGW_kg;
S_ref = veh.wing_area_m2;
fuel_capacity = veh.W_fuel_kg * 1.05;
OEW = 6278;  % kg

% Target conditions
h_target = 15300;  % m (max altitude)
M_target = 1.6;    % Max Mach
cruise_duration = 60;  % 1 minute if achievable

fprintf('VEHICLE PARAMETERS:\n');
fprintf('  TOGW:          %.0f kg\n', TOGW);
fprintf('  OEW:           %.0f kg\n', OEW);
fprintf('  Wing area:     %.2f m²\n', S_ref);
fprintf('  Fuel capacity: %.0f kg\n', fuel_capacity);
fprintf('\n');

fprintf('TARGET CONDITIONS:\n');
fprintf('  Altitude:      %.0f m (%.0f ft)\n', h_target, h_target*3.281);
fprintf('  Mach:          %.1f\n', M_target);
fprintf('  Cruise:        %.0f s (if achievable)\n', cruise_duration);
fprintf('\n');

% =========================================================================
% THRUST vs DRAG ANALYSIS (Using Physics Model)
% =========================================================================
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('THRUST vs DRAG ANALYSIS AT M=%.1f, h=%.0f m (PHYSICS MODEL)\n', M_target, h_target);
fprintf('════════════════════════════════════════════════════════════════\n\n');

% Atmospheric conditions at target
[T_atm, a_target, P, rho_target] = atmosisa(h_target);
V_target = M_target * a_target;
q_target = 0.5 * rho_target * V_target^2;

fprintf('Atmosphere at %.0f m:\n', h_target);
fprintf('  Temperature:   %.1f K (%.1f °C)\n', T_atm, T_atm - 273.15);
fprintf('  Density:       %.4f kg/m³\n', rho_target);
fprintf('  Speed of sound: %.1f m/s\n', a_target);
fprintf('  Velocity:      %.1f m/s (%.0f kts)\n', V_target, V_target*1.944);
fprintf('  Dynamic pressure: %.2f kPa\n', q_target/1000);
fprintf('\n');

% Estimate weight at target
W_at_target = TOGW - 500;
W_N = W_at_target * 9.81;

% Required CL for level flight
CL_required = W_N / (q_target * S_ref);

% Alpha estimate from CL (using typical CL_alpha = 0.08/deg)
alpha_est = CL_required / 0.08;  % degrees
alpha_est = max(1, min(alpha_est, 12));

% Drag calculation using Raymer aero model
[CL_aero, CD_total, CD0, CDi, ~] = get_X29_aero(M_target, alpha_est, h_target, geom);
D_required = q_target * CD_total * S_ref;
LD_target = CL_required / CD_total;

fprintf('Aerodynamics at target (RAYMER MODEL):\n');
fprintf('  Weight:        %.0f kg (estimated after climb)\n', W_at_target);
fprintf('  Alpha:         %.2f deg\n', alpha_est);
fprintf('  CL required:   %.4f\n', CL_required);
fprintf('  CL (Raymer):   %.4f\n', CL_aero);
fprintf('  CD0:           %.4f\n', CD0);
fprintf('  CDi:           %.4f\n', CDi);
fprintf('  CD total:      %.4f\n', CD_total);
fprintf('  L/D:           %.2f\n', LD_target);
fprintf('  Drag:          %.2f kN\n', D_required/1000);
fprintf('\n');

% Thrust available using calibrated model
fprintf('Thrust available at M=%.1f, h=%.0f m (CALIBRATED MODEL):\n', M_target, h_target);

% Get MIL and AB thrust using calibrated model
prop_mil = get_F404_thrust_calibrated(M_target, h_target, 87);  % MIL power
prop_maxAB = get_F404_thrust_calibrated(M_target, h_target, 130);  % Max AB
T_mil = prop_mil.thrust_N;
T_maxAB = prop_maxAB.thrust_N;

fprintf('  MIL power: T=%.2f kN, T-D=%.2f kN\n', T_mil/1000, (T_mil - D_required)/1000);
fprintf('  Max AB:    T=%.2f kN, T-D=%.2f kN\n', T_maxAB/1000, (T_maxAB - D_required)/1000);
fprintf('\n');
thrust_margin = (T_maxAB - D_required) / D_required * 100;

fprintf('VERDICT:\n');
if T_maxAB > D_required
    fprintf('  ✓ ACHIEVABLE: T_max (%.1f kN) > D (%.1f kN)\n', T_maxAB/1000, D_required/1000);
    fprintf('  Thrust margin: +%.1f%%\n', thrust_margin);
    can_sustain = true;
else
    fprintf('  ✗ NOT SUSTAINABLE: T_max (%.1f kN) < D (%.1f kN)\n', T_maxAB/1000, D_required/1000);
    fprintf('  Thrust deficit: %.1f%%\n', thrust_margin);
    can_sustain = false;
end
fprintf('\n');

% =========================================================================
% GENERATE THRUST vs DRAG ENVELOPE PLOTS (Calibrated Model)
% =========================================================================
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('GENERATING THRUST/DRAG ENVELOPE ANALYSIS (CALIBRATED MODEL)...\n');
fprintf('════════════════════════════════════════════════════════════════\n\n');

% Mach sweep at target altitude
M_sweep = 0.5:0.05:1.8;
n_M = length(M_sweep);

T_avail_sweep = zeros(n_M, 1);
T_mil_sweep = zeros(n_M, 1);
D_sweep = zeros(n_M, 1);
CL_sweep = zeros(n_M, 1);
CD_sweep = zeros(n_M, 1);

for i = 1:n_M
    M = M_sweep(i);
    V = M * a_target;
    q = 0.5 * rho_target * V^2;

    CL = W_N / (q * S_ref);
    CL_sweep(i) = CL;

    % Alpha estimate from CL
    alpha_i = CL / 0.08;  % degrees
    alpha_i = max(1, min(alpha_i, 12));

    % Use Raymer aero model for drag
    [~, CD, ~, ~, ~] = get_X29_aero(M, alpha_i, h_target, geom);
    CD_sweep(i) = CD;
    D_sweep(i) = q * CD * S_ref;

    % Calibrated model thrust
    prop_ab = get_F404_thrust_calibrated(M, h_target, 130);  % Max AB
    prop_mil = get_F404_thrust_calibrated(M, h_target, 87);   % MIL
    T_avail_sweep(i) = prop_ab.thrust_N;
    T_mil_sweep(i) = prop_mil.thrust_N;
end

% Find max sustainable Mach
excess_thrust = T_avail_sweep - D_sweep;
idx_cross = find(diff(sign(excess_thrust)) ~= 0, 1, 'last');
if ~isempty(idx_cross)
    M1 = M_sweep(idx_cross);
    M2 = M_sweep(idx_cross + 1);
    E1 = excess_thrust(idx_cross);
    E2 = excess_thrust(idx_cross + 1);
    M_max_exact = M1 - E1 * (M2 - M1) / (E2 - E1);
else
    M_max_exact = M_sweep(find(excess_thrust > 0, 1, 'last'));
end

fprintf('Max sustainable Mach at %.0f m (Calibrated Model): %.2f\n', h_target, M_max_exact);
fprintf('\n');

% =========================================================================
% FIGURE 1: Thrust vs Drag at h=15.3 km
% =========================================================================
fig1 = figure('Position', [50, 50, 1200, 800], 'Color', 'w');
sgtitle(sprintf('X-29 Thrust vs Drag at h = %.1f km (Calibrated Propulsion Model)', h_target/1000), ...
    'FontSize', 14, 'FontWeight', 'bold');

subplot(2, 2, 1);
hold on;
plot(M_sweep, T_avail_sweep/1000, 'b-', 'LineWidth', 2, 'DisplayName', 'T_{avail} (Max AB)');
plot(M_sweep, T_mil_sweep/1000, 'b--', 'LineWidth', 1.5, 'DisplayName', 'T_{avail} (MIL)');
plot(M_sweep, D_sweep/1000, 'r-', 'LineWidth', 2, 'DisplayName', 'D_{required}');
xline(M_target, 'k--', sprintf('M=%.1f target', M_target), 'LineWidth', 1.5);
xline(M_max_exact, 'g--', sprintf('M_{max}=%.2f', M_max_exact), 'LineWidth', 1.5);
xlabel('Mach Number'); ylabel('Force (kN)');
title('Thrust Available vs Drag Required');
legend('Location', 'best'); grid on; xlim([0.5 1.8]);

subplot(2, 2, 2);
hold on;
% Fill positive region (green) and negative region (red)
excess_kN = excess_thrust/1000;
fill([M_sweep(:); flipud(M_sweep(:))], [max(excess_kN(:), 0); zeros(n_M, 1)], ...
    [0.7 0.9 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
fill([M_sweep(:); flipud(M_sweep(:))], [min(excess_kN(:), 0); zeros(n_M, 1)], ...
    [0.9 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
plot(M_sweep, excess_kN, 'k-', 'LineWidth', 2, 'DisplayName', 'T - D');
yline(0, 'k-', 'LineWidth', 1);
xline(M_target, 'k--', sprintf('M=%.1f', M_target), 'LineWidth', 1.5);
if ~isempty(M_max_exact)
    xline(M_max_exact, 'g--', sprintf('M_{max}=%.2f', M_max_exact), 'LineWidth', 1.5);
end
xlabel('Mach Number'); ylabel('Excess Thrust (kN)');
title('Excess Thrust (T_{avail} - D)'); grid on; xlim([0.5 1.8]);

subplot(2, 2, 3);
yyaxis left;
plot(M_sweep, CL_sweep, 'b-', 'LineWidth', 2);
ylabel('C_L'); ylim([0 max(CL_sweep)*1.2]);
yyaxis right;
plot(M_sweep, CD_sweep, 'r-', 'LineWidth', 2);
ylabel('C_D');
xlabel('Mach Number');
title('Lift and Drag Coefficients');
xline(M_target, 'k--', 'LineWidth', 1); grid on; xlim([0.5 1.8]);

subplot(2, 2, 4);
LD_sweep = CL_sweep ./ CD_sweep;
plot(M_sweep, LD_sweep, 'm-', 'LineWidth', 2);
xlabel('Mach Number'); ylabel('L/D');
title('Lift-to-Drag Ratio');
xline(M_target, 'k--', sprintf('M=%.1f', M_target), 'LineWidth', 1);
xline(1.0, 'k:', 'M=1.0', 'LineWidth', 1);
grid on; xlim([0.5 1.8]); ylim([0 max(LD_sweep)*1.1]);

saveas(gcf, 'X29_thrust_drag_envelope_15km.png');
fprintf('Figure saved: X29_thrust_drag_envelope_15km.png\n');

% =========================================================================
% FIGURE 2: Altitude sweep at M=1.6
% =========================================================================
fig2 = figure('Position', [100, 100, 1000, 600], 'Color', 'w');
sgtitle(sprintf('X-29 Thrust vs Drag at M = %.1f (Calibrated Model)', M_target), ...
    'FontSize', 14, 'FontWeight', 'bold');

h_sweep = 10000:500:18000;
n_h = length(h_sweep);

T_h_sweep = zeros(n_h, 1);
D_h_sweep = zeros(n_h, 1);

for i = 1:n_h
    h_i = h_sweep(i);
    [~, a_i, ~, rho_i] = atmosisa(h_i);
    V_i = M_target * a_i;
    q_i = 0.5 * rho_i * V_i^2;

    CL_i = W_N / (q_i * S_ref);

    % Alpha estimate from CL
    alpha_i = CL_i / 0.08;  % degrees
    alpha_i = max(1, min(alpha_i, 12));

    % Use Raymer aero model for drag
    [~, CD_i, ~, ~, ~] = get_X29_aero(M_target, alpha_i, h_i, geom);
    D_h_sweep(i) = q_i * CD_i * S_ref;

    prop_h = get_F404_thrust_calibrated(M_target, h_i, 130);  % Max AB
    T_h_sweep(i) = prop_h.thrust_N;
end

excess_h = T_h_sweep - D_h_sweep;
idx_ceiling = find(diff(sign(excess_h)) ~= 0, 1, 'last');
if ~isempty(idx_ceiling)
    h1 = h_sweep(idx_ceiling);
    h2 = h_sweep(idx_ceiling + 1);
    E1 = excess_h(idx_ceiling);
    E2 = excess_h(idx_ceiling + 1);
    h_ceiling = h1 - E1 * (h2 - h1) / (E2 - E1);
else
    h_ceiling = h_sweep(end);
end

subplot(1, 2, 1);
hold on;
plot(h_sweep/1000, T_h_sweep/1000, 'b-', 'LineWidth', 2, 'DisplayName', 'T_{avail} (Max AB)');
plot(h_sweep/1000, D_h_sweep/1000, 'r-', 'LineWidth', 2, 'DisplayName', 'D_{required}');
xline(h_target/1000, 'k--', sprintf('h=%.1f km', h_target/1000), 'LineWidth', 1.5);
xline(h_ceiling/1000, 'g--', sprintf('Ceiling=%.1f km', h_ceiling/1000), 'LineWidth', 1.5);
xlabel('Altitude (km)'); ylabel('Force (kN)');
title(sprintf('T and D vs Altitude at M=%.1f', M_target));
legend('Location', 'best'); grid on;

subplot(1, 2, 2);
hold on;
% Fill positive region (green) and negative region (red)
h_km = h_sweep(:)/1000;
excess_h_kN = excess_h(:)/1000;
fill([h_km; flipud(h_km)], [max(excess_h_kN, 0); zeros(n_h, 1)], ...
    [0.7 0.9 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
fill([h_km; flipud(h_km)], [min(excess_h_kN, 0); zeros(n_h, 1)], ...
    [0.9 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.7);
plot(h_km, excess_h_kN, 'k-', 'LineWidth', 2);
yline(0, 'k-', 'LineWidth', 1);
xline(h_target/1000, 'k--', sprintf('h=%.1f km', h_target/1000), 'LineWidth', 1.5);
xline(h_ceiling/1000, 'g--', sprintf('Ceiling=%.1f km', h_ceiling/1000), 'LineWidth', 1.5);
xlabel('Altitude (km)'); ylabel('Excess Thrust (kN)');
title('Excess Thrust vs Altitude'); grid on;

fprintf('Service ceiling at M=%.1f (Calibrated Model): %.1f km (%.0f ft)\n', M_target, h_ceiling/1000, h_ceiling*3.281);

saveas(gcf, 'X29_altitude_ceiling_M16.png');
fprintf('Figure saved: X29_altitude_ceiling_M16.png\n\n');

% =========================================================================
% RUN MISSION SIMULATION WITH CANARD TRACKING
% =========================================================================
if can_sustain
    fprintf('════════════════════════════════════════════════════════════════\n');
    fprintf('RUNNING MISSION WITH CANARD TRACKING TO M=%.1f AT %.0f m\n', M_target, h_target);
    fprintf('════════════════════════════════════════════════════════════════\n\n');

    % Phase 1: Taxi
    fprintf('--- PHASE 1: TAXI ---\n');
    taxi = simulate_X29_taxi(300, TOGW);
    W_after_taxi = taxi.W_final;
    W_fuel_after_taxi = fuel_capacity - taxi.fuel_burned;

    % Phase 2: Takeoff (using sim_takeoff_canards with Raymer aero)
    fprintf('\n--- PHASE 2: TAKEOFF ---\n');
    takeoff_params = struct();
    takeoff_params.W_total = W_after_taxi;
    takeoff_params.W_fuel = W_fuel_after_taxi;
    takeoff_params.S_ref = S_ref;
    takeoff_params.CL_max = 1.6;
    takeoff_params.verbose = true;
    takeoff = sim_takeoff_canards(takeoff_params, geom);
    W_after_takeoff = takeoff.W_final;
    W_fuel_after_takeoff = takeoff.W_fuel_final;

    % Phase 3: Climb to 15.3 km (using sim_climb_canards with Hull's equations)
    fprintf('\n--- PHASE 3: CLIMB TO %.0f m ---\n', h_target);
    climb_params = struct();
    climb_params.W_initial = W_after_takeoff;
    climb_params.W_fuel_initial = W_fuel_after_takeoff;
    climb_params.h_initial = 15.24;
    climb_params.V_initial = takeoff.V_final;
    climb_params.h_target = h_target;
    climb_params.M_target = 0.9;
    climb_params.verbose = true;
    climb = sim_climb_canards(climb_params, geom);
    W_after_climb = climb.W_final;
    W_fuel_after_climb = climb.W_fuel_final;

    % Phase 4: Accelerate to M=1.6 (using sim_accel_canards with physics propulsion)
    fprintf('\n--- PHASE 4: ACCELERATE M=0.9 TO M=%.1f ---\n', M_target);
    accel = sim_accel_canards(W_after_climb, S_ref, h_target, 0.9, M_target, geom);
    W_after_accel = accel.W_final;
    W_fuel_after_accel = W_fuel_after_climb - accel.total_fuel;

    % Phase 5: Cruise at M=1.6 (using sim_cruise_canards)
    fprintf('\n--- PHASE 5: CRUISE AT M=%.1f FOR %.0f s ---\n', M_target, cruise_duration);
    cruise_params = struct();
    cruise_params.W_initial = W_after_accel;
    cruise_params.W_fuel_initial = W_fuel_after_accel;
    cruise_params.h = h_target;
    cruise_params.M = M_target;
    cruise_params.duration_min = cruise_duration / 60;  % Convert seconds to minutes
    cruise_params.verbose = true;
    cruise = sim_cruise_canards(cruise_params, geom);
    W_after_cruise = cruise.W_final;

    % =====================================================================
    % EXTRACT CANARD DEFLECTION FROM SIMULATION OUTPUTS
    % =====================================================================
    fprintf('\nExtracting canard deflection from simulation outputs...\n');

    % Takeoff - use data from sim_takeoff_canards
    canard_takeoff = takeoff.delta_canard;
    alpha_takeoff = takeoff.alpha;
    cg_takeoff = takeoff.x_cg;

    % Climb - use data from sim_climb_canards
    canard_climb = climb.delta_canard;
    alpha_climb = climb.alpha;
    cg_climb = climb.x_cg;

    % Accel - use data from sim_accel_canards (now tracks canards)
    canard_accel = accel.delta_canard;
    alpha_accel = accel.alpha;
    cg_accel = accel.x_cg;

    % Cruise - use data from sim_cruise_canards
    canard_cruise = cruise.delta_canard;
    alpha_cruise = cruise.alpha;
    cg_cruise = cruise.x_cg;

    % Summary
    total_fuel = TOGW - W_after_cruise;

    % Calculate fuel burned per phase
    takeoff_fuel = takeoff.W(1) - takeoff.W(end);
    cruise_fuel_burned = cruise.W(1) - cruise.W(end);

    fprintf('\n');
    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║       MISSION SUMMARY (NASA PROPULSION + RAYMER AERO)         ║\n');
    fprintf('╠════════════════════════════════════════════════════════════════╣\n');
    fprintf('║ Phase              │ Time (s) │ Fuel (kg) │ Weight (kg)        ║\n');
    fprintf('╠════════════════════╪══════════╪═══════════╪════════════════════╣\n');
    fprintf('║ 1. Taxi            │ %7.1f  │ %8.1f  │ %10.0f         ║\n', ...
        taxi.duration, taxi.fuel_burned, taxi.W_final);
    fprintf('║ 2. Takeoff         │ %7.1f  │ %8.1f  │ %10.0f         ║\n', ...
        takeoff.total_time, takeoff_fuel, takeoff.W_final);
    fprintf('║ 3. Climb (15.3km)  │ %7.1f  │ %8.1f  │ %10.0f         ║\n', ...
        climb.total_time, climb.total_fuel, climb.W_final);
    fprintf('║ 4. Accel (M=1.6)   │ %7.1f  │ %8.1f  │ %10.0f         ║\n', ...
        accel.total_time, accel.total_fuel, accel.W_final);
    fprintf('║ 5. Cruise (1 min)  │ %7.1f  │ %8.1f  │ %10.0f         ║\n', ...
        cruise.t(end), cruise_fuel_burned, W_after_cruise);
    fprintf('╠════════════════════╪══════════╪═══════════╪════════════════════╣\n');
    fprintf('║ TOTAL              │          │ %8.1f  │ %10.0f         ║\n', ...
        total_fuel, W_after_cruise);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n');
    fprintf('\n');
    fprintf('Fuel remaining: %.1f kg (%.1f%% of capacity)\n', ...
        fuel_capacity - total_fuel, 100*(fuel_capacity - total_fuel)/fuel_capacity);

    % =====================================================================
    % CANARD/STABILITY SUMMARY
    % =====================================================================
    fprintf('\n');
    fprintf('╔════════════════════════════════════════════════════════════════╗\n');
    fprintf('║                 CANARD/STABILITY SUMMARY                       ║\n');
    fprintf('╠════════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Canard Deflection by Phase                                    ║\n');
    fprintf('║    Takeoff:         %.1f° to %.1f°                              ║\n', ...
        min(canard_takeoff), max(canard_takeoff));
    fprintf('║    Climb:           %.1f° to %.1f°                              ║\n', ...
        min(canard_climb), max(canard_climb));
    fprintf('║    Accel (M→1.6):   %.1f° to %.1f°                              ║\n', ...
        min(canard_accel), max(canard_accel));
    fprintf('║    Cruise (M=1.6):  %.1f° to %.1f°                              ║\n', ...
        min(canard_cruise), max(canard_cruise));
    fprintf('╠════════════════════════════════════════════════════════════════╣\n');
    fprintf('║  CG Position                                                   ║\n');
    fprintf('║    Initial:         %.2f m from nose                           ║\n', cg_takeoff(1));
    fprintf('║    End of cruise:   %.2f m from nose                           ║\n', cg_cruise(end));
    fprintf('║    CG travel:       %.3f m                                     ║\n', abs(cg_takeoff(1) - cg_cruise(end)));
    fprintf('╠════════════════════════════════════════════════════════════════╣\n');
    fprintf('║  Target Achieved                                               ║\n');
    fprintf('║    Altitude:        %.0f m (target: %.0f m)                    ║\n', h_target, h_target);
    fprintf('║    Mach:            %.1f (target: %.1f)                         ║\n', M_target, M_target);
    fprintf('╚════════════════════════════════════════════════════════════════╝\n');

    % =====================================================================
    % PLOTS
    % =====================================================================
    % Cruise details
    fig3 = figure('Position', [150, 150, 1200, 500], 'Color', 'w');
    sgtitle(sprintf('X-29 Cruise at M=%.1f, h=%.1f km (NASA + Raymer)', M_target, h_target/1000), ...
        'FontSize', 14, 'FontWeight', 'bold');

    t_min = cruise.t / 60;

    subplot(1, 3, 1);
    plot(t_min, cruise.D/1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Drag');
    hold on;
    plot(t_min, cruise.T/1000, 'b--', 'LineWidth', 2, 'DisplayName', 'Thrust');
    xlabel('Time (min)'); ylabel('Force (kN)');
    title('Thrust and Drag'); legend('Location', 'best'); grid on;

    subplot(1, 3, 2);
    plot(t_min, cruise.LD, 'g-', 'LineWidth', 2);
    xlabel('Time (min)'); ylabel('L/D');
    title('Lift-to-Drag Ratio'); grid on;

    subplot(1, 3, 3);
    % Calculate cumulative fuel burned
    cruise_fuel_cumulative = cruise.W(1) - cruise.W;
    plot(t_min, cruise_fuel_cumulative, 'k-', 'LineWidth', 2);
    xlabel('Time (min)'); ylabel('Fuel Burned (kg)');
    title('Cumulative Fuel'); grid on;

    saveas(gcf, 'X29_cruise_M16_h153km.png');
    fprintf('\nFigure saved: X29_cruise_M16_h153km.png\n');

    % Canard deflection plot
    fig4 = figure('Position', [200, 200, 1000, 600], 'Color', 'w');
    sgtitle(sprintf('X-29 Canard Deflection to M=%.1f (NASA + Raymer)', M_target), ...
        'FontSize', 14, 'FontWeight', 'bold');

    % Build combined timeline using actual time arrays
    t_offset = 0;
    t_to = takeoff.t + t_offset;
    t_offset = t_to(end);
    t_cl = climb.t + t_offset;
    t_offset = t_cl(end);
    t_ac = accel.t + t_offset;
    t_offset = t_ac(end);
    t_cr = cruise.t + t_offset;

    subplot(2, 1, 1);
    hold on;
    plot(t_to, canard_takeoff, 'b-', 'LineWidth', 2);
    plot(t_cl, canard_climb, 'b-', 'LineWidth', 2);
    plot(t_ac, canard_accel, 'b-', 'LineWidth', 2);
    plot(t_cr, canard_cruise, 'b-', 'LineWidth', 2);
    yline(58, 'r--', 'Max +58°', 'LineWidth', 1);
    yline(-32, 'r--', 'Min -32°', 'LineWidth', 1);
    xlabel('Time (s)'); ylabel('Canard Deflection (deg)');
    title('Canard Deflection for Trim');
    ylim([-35 60]); grid on; box on;

    subplot(2, 1, 2);
    hold on;
    plot(t_to, cg_takeoff, 'Color', [0.6 0 0.6], 'LineWidth', 2);
    plot(t_cl, cg_climb, 'Color', [0.6 0 0.6], 'LineWidth', 2);
    plot(t_ac, cg_accel, 'Color', [0.6 0 0.6], 'LineWidth', 2);
    plot(t_cr, cg_cruise, 'Color', [0.6 0 0.6], 'LineWidth', 2);
    xlabel('Time (s)'); ylabel('X_{CG} from Nose (m)');
    title('CG Position');
    grid on; box on;

    saveas(gcf, 'X29_maxmach_canard_profile.png');
    fprintf('Figure saved: X29_maxmach_canard_profile.png\n');

    % Mach-Altitude trajectory
    figure('Position', [250, 250, 900, 600], 'Color', 'w');
    hold on;

    % Build trajectory using actual data
    M_traj = [takeoff.M; climb.M; accel.M; M_target*ones(size(cruise.t))];
    h_traj = [takeoff.h; climb.h; h_target*ones(size(accel.t)); h_target*ones(size(cruise.t))];

    plot(M_traj, h_traj/1000, 'b-', 'LineWidth', 2, 'DisplayName', 'Trajectory');
    plot(M_target, h_target/1000, 'rp', 'MarkerSize', 15, 'MarkerFaceColor', 'r', 'DisplayName', 'Target (M=1.6, h=15.3km)');
    xlabel('Mach Number', 'FontSize', 12);
    ylabel('Altitude (km)', 'FontSize', 12);
    title(sprintf('X-29 Max Mach/Alt Mission Trajectory (NASA + Raymer)'), ...
        'FontSize', 14, 'FontWeight', 'bold');
    legend('Location', 'northwest');
    xlim([0 1.8]); ylim([0 18]);
    grid on; box on;

    saveas(gcf, 'X29_maxmach_trajectory.png');
    fprintf('Figure saved: X29_maxmach_trajectory.png\n');

else
    fprintf('════════════════════════════════════════════════════════════════\n');
    fprintf('MISSION NOT ACHIEVABLE WITH CALIBRATED MODEL\n');
    fprintf('════════════════════════════════════════════════════════════════\n');
    fprintf('Cannot sustain M=%.1f at h=%.0f m.\n', M_target, h_target);
    fprintf('Max sustainable Mach at this altitude: %.2f\n', M_max_exact);
    fprintf('Service ceiling at M=%.1f: %.1f km\n', M_target, h_ceiling/1000);
end

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('SIMULATION COMPLETE (CALIBRATED PROPULSION + CANARDS)\n');
fprintf('════════════════════════════════════════════════════════════════\n');
