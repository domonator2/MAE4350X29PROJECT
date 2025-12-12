%{
====================================================================================================
SCRIPT: run_Edwards_to_Luke_mission.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 12/05/2025
====================================================================================================
DESCRIPTION:
Realistic X-29 mission simulation from Edwards AFB to Luke AFB.

Route Information:
  - Start: Edwards AFB, CA (elevation 2,303 ft / 702 m)
  - End:   Luke AFB, AZ (elevation 1,085 ft / 331 m)
  - Distance: 286.5 NM (530.6 km)

Mission Profile:
  1. Taxi + Warmup (5 min at Edwards AFB)
  2. Takeoff (from 702 m field elevation)
  3. Climb (to 12.2 km cruise altitude)
  4. Transonic Acceleration (M=0.9 to M=1.4)
  5. Supersonic Cruise (M=1.4, h=12.2 km)
  6. Descent to Subsonic Cruise (M=1.4->0.6, h=12.2->9.12 km)
  7. Subsonic Cruise (M=0.6, h=9.12 km) - fuel efficient, L=W trim
  8. Descent to Pattern (h=9.12 km to pattern altitude)
  9. Landing (at Luke AFB, 331 m elevation)

The simulation calculates phase distances and adjusts cruise durations to match
the total route distance of 286.5 NM.
====================================================================================================
%}

clear; clc; close all;

% Add paths - get project root
% Script is in: Trajectory_Performance/X29_Simulation_Canards_CycleAnalysis_RealisticTrajectory/
% Project root is 2 levels up: MAE4350-X29-PROJECT/
script_path = fileparts(mfilename('fullpath'));
traj_perf_path = fileparts(script_path);  % Trajectory_Performance
proj_root = fileparts(traj_perf_path);     % MAE4350-X29-PROJECT
addpath(genpath(proj_root));

fprintf('╔════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║      X-29 REALISTIC TRAJECTORY: EDWARDS AFB → LUKE AFB                 ║\n');
fprintf('║              Physics-Based Propulsion + Cycle Analysis                 ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════════╝\n\n');

% =========================================================================
% MISSION PARAMETERS
% =========================================================================
% Airfield data
h_edwards_ft = 2303;      % Edwards AFB elevation [ft]
h_luke_ft = 1085;         % Luke AFB elevation [ft]
h_edwards = h_edwards_ft * 0.3048;  % [m] = 702 m
h_luke = h_luke_ft * 0.3048;        % [m] = 331 m

% Route distance
distance_nm = 286.5;      % Total route distance [NM]
distance_km = distance_nm * 1.852;  % [km] = 530.6 km
distance_m = distance_km * 1000;    % [m]

fprintf('ROUTE INFORMATION:\n');
fprintf('  Edwards AFB:    %.0f ft (%.0f m) elevation\n', h_edwards_ft, h_edwards);
fprintf('  Luke AFB:       %.0f ft (%.0f m) elevation\n', h_luke_ft, h_luke);
fprintf('  Distance:       %.1f NM (%.1f km)\n', distance_nm, distance_km);
fprintf('\n');

% Cruise parameters
M_supersonic = 1.4;       % Supersonic cruise Mach
h_supersonic = 12200;     % Supersonic cruise altitude [m] (40,000 ft)
M_subsonic = 0.6;         % Subsonic cruise Mach
h_subsonic = 9120;        % Subsonic cruise altitude [m] (29,921 ft)

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

% Find configuration: sweep=-30°, tau=0.077, duration=10 min
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

fprintf('  Configuration: sweep=%.0f°, tau=%.4f, duration=%.0f min\n', ...
    veh.sweep_deg, veh.tau, veh.duration_min);

geom = build_X29_geometry(veh, trade_data.vehicle_base);
fprintf('  Geometry loaded successfully.\n\n');

% =========================================================================
% VEHICLE PARAMETERS
% =========================================================================
TOGW = veh.TOGW_kg;
S_ref = veh.wing_area_m2;
fuel_capacity = veh.W_fuel_kg * 1.05;
OEW = 6278;

fprintf('VEHICLE PARAMETERS:\n');
fprintf('  TOGW:          %.0f kg\n', TOGW);
fprintf('  OEW:           %.0f kg\n', OEW);
fprintf('  Wing area:     %.2f m²\n', S_ref);
fprintf('  Fuel capacity: %.0f kg\n', fuel_capacity);
fprintf('\n');

% =========================================================================
% PHASE 1: TAXI + WARMUP (at Edwards AFB)
% =========================================================================
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 1: TAXI + WARMUP (Edwards AFB, %.0f m)\n', h_edwards);
fprintf('════════════════════════════════════════════════════════════════\n');

taxi = simulate_X29_taxi(300, TOGW);  % 5 minutes
W_after_taxi = taxi.W_final;
dist_taxi = 0;  % No horizontal distance during taxi

% =========================================================================
% PHASE 2: TAKEOFF (from Edwards AFB elevation)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 2: TAKEOFF (from %.0f m field elevation)\n', h_edwards);
fprintf('════════════════════════════════════════════════════════════════\n');

% Modify takeoff for Edwards elevation
takeoff = sim_takeoff_elevated(W_after_taxi, S_ref, h_edwards, geom);
W_after_takeoff = takeoff.W_final;
dist_takeoff = takeoff.total_distance;

% =========================================================================
% PHASE 3: CLIMB (to supersonic cruise altitude)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 3: CLIMB (%.0f m to %.0f m)\n', h_edwards + 15, h_supersonic);
fprintf('════════════════════════════════════════════════════════════════\n');

% Use sim_climb_elevated for climb from elevated airfield
% Climb to h_supersonic at M=0.9 (acceleration phase will take it to M_supersonic)
climb = sim_climb_elevated(W_after_takeoff, S_ref, h_edwards + 15.24, takeoff.V_LOF, h_supersonic, 0.9, geom);
W_after_climb = climb.W_final;
dist_climb = climb.x(end);  % Horizontal distance from Hull's equations integration

% =========================================================================
% PHASE 4: TRANSONIC ACCELERATION
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 4: TRANSONIC ACCELERATION (M=0.9 to M=%.1f)\n', M_supersonic);
fprintf('════════════════════════════════════════════════════════════════\n');

% Check if climb already reached target Mach
if climb.M_final >= M_supersonic * 0.99
    M_accel_start = climb.M_final;
else
    M_accel_start = 0.9;
end
accel = sim_accel(W_after_climb, S_ref, h_supersonic, M_accel_start, M_supersonic, geom);
W_after_accel = accel.W_final;
dist_accel = accel.distance;

% =========================================================================
% CALCULATE REMAINING DISTANCE FOR CRUISES
% =========================================================================
% Distance used so far (non-cruise phases)
dist_used = dist_takeoff + dist_climb + dist_accel;

% Estimate descent distances
[~, a_supersonic, ~, ~] = atmosisa(h_supersonic);
[~, a_subsonic, ~, ~] = atmosisa(h_subsonic);
V_supersonic = M_supersonic * a_supersonic;
V_subsonic = M_subsonic * a_subsonic;

% Descent 1: Supersonic to subsonic cruise (12.2 km -> 9.12 km)
% Approximate distance based on descent angle
gamma_desc1 = -3;  % degrees
h_drop1 = h_supersonic - h_subsonic;
dist_descent1_est = h_drop1 / tand(abs(gamma_desc1));

% Descent 2: Subsonic cruise to pattern (9.12 km -> 300 m)
gamma_desc2 = -5;  % degrees
h_drop2 = h_subsonic - 300;
dist_descent2_est = h_drop2 / tand(abs(gamma_desc2));

% Landing distance estimate
dist_landing_est = 10000;  % ~10 km for approach + rollout

% Total non-cruise distance
dist_non_cruise = dist_used + dist_descent1_est + dist_descent2_est + dist_landing_est;

% Remaining distance for both cruise phases
dist_cruise_total = distance_m - dist_non_cruise;

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('DISTANCE BUDGET:\n');
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('  Takeoff:              %.1f km\n', dist_takeoff/1000);
fprintf('  Climb:                %.1f km\n', dist_climb/1000);
fprintf('  Accel:                %.1f km\n', dist_accel/1000);
fprintf('  Descent 1 (est):      %.1f km\n', dist_descent1_est/1000);
fprintf('  Descent 2 (est):      %.1f km\n', dist_descent2_est/1000);
fprintf('  Landing (est):        %.1f km\n', dist_landing_est/1000);
fprintf('  ─────────────────────────────\n');
fprintf('  Non-cruise total:     %.1f km\n', dist_non_cruise/1000);
fprintf('  Available for cruise: %.1f km\n', dist_cruise_total/1000);
fprintf('  Route total:          %.1f km (%.1f NM)\n', distance_km, distance_nm);

% FIXED REQUIREMENT: 10 minutes supersonic cruise (PRIORITY)
t_supersonic_cruise = 600;  % 10 minutes [s]
dist_supersonic_cruise = V_supersonic * t_supersonic_cruise;

% Subsonic cruise covers remaining distance (may be zero or small)
% Correction: Actual descent distances exceed estimates by ~17 km
% At M=0.6, 30kft: V ≈ 182 m/s, so 17 km = 93s of cruise time
dist_correction = 17000;  % m - empirical correction for descent distance underestimate
dist_subsonic_cruise = dist_cruise_total - dist_supersonic_cruise - dist_correction;
if dist_subsonic_cruise < 0
    fprintf('\n  NOTE: Route too short for full 10-min supersonic + subsonic cruise.\n');
    fprintf('  Supersonic cruise will exceed planned cruise distance.\n');
    fprintf('  Total distance will be ~%.1f km (%.1f NM) instead of %.1f NM.\n', ...
        (dist_non_cruise + dist_supersonic_cruise)/1000, ...
        (dist_non_cruise + dist_supersonic_cruise)/1852, distance_nm);
    % Keep 10-min supersonic, skip subsonic cruise
    dist_subsonic_cruise = 0;
    t_subsonic_cruise = 0;
else
    t_subsonic_cruise = dist_subsonic_cruise / V_subsonic;
end

fprintf('\n  Supersonic cruise:    %.1f km (%.1f s = %.1f min)\n', ...
    dist_supersonic_cruise/1000, t_supersonic_cruise, t_supersonic_cruise/60);
fprintf('  Subsonic cruise:      %.1f km (%.1f s = %.1f min)\n', ...
    dist_subsonic_cruise/1000, t_subsonic_cruise, t_subsonic_cruise/60);

% =========================================================================
% PHASE 5: SUPERSONIC CRUISE
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 5: SUPERSONIC CRUISE (M=%.1f, h=%.0f m, %.1f min)\n', ...
    M_supersonic, h_supersonic, t_supersonic_cruise/60);
fprintf('════════════════════════════════════════════════════════════════\n');

cruise_ss = simulate_X29_cruise_quasi_steady(M_supersonic, h_supersonic, ...
    W_after_accel, t_supersonic_cruise, 1, geom);
W_after_ss_cruise = cruise_ss.W(end);
dist_ss_cruise_actual = V_supersonic * cruise_ss.t(end);

% Pattern altitude
h_pattern = h_luke + 300;  % Pattern altitude is 300m AGL above Luke

if t_subsonic_cruise > 0
    % =========================================================================
    % PHASE 6: DESCENT TO SUBSONIC CRUISE
    % =========================================================================
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('PHASE 6: DESCENT TO SUBSONIC CRUISE (M=%.1f->%.1f, h=%.0f->%.0f m)\n', ...
        M_supersonic, M_subsonic, h_supersonic, h_subsonic);
    fprintf('════════════════════════════════════════════════════════════════\n');

    descent1 = sim_descent_to_subsonic(W_after_ss_cruise, S_ref, ...
        h_supersonic, V_supersonic, h_subsonic, M_subsonic, geom);
    W_after_descent1 = descent1.W_final;
    dist_descent1_actual = descent1.distance;

    % =========================================================================
    % PHASE 7: SUBSONIC CRUISE (Fuel Efficient, L=W Trim)
    % =========================================================================
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('PHASE 7: SUBSONIC CRUISE (M=%.1f, h=%.0f m, %.1f min)\n', ...
        M_subsonic, h_subsonic, t_subsonic_cruise/60);
    fprintf('════════════════════════════════════════════════════════════════\n');

    cruise_sub = simulate_subsonic_cruise_efficient(M_subsonic, h_subsonic, ...
        W_after_descent1, t_subsonic_cruise, S_ref, geom);
    W_after_sub_cruise = cruise_sub.W(end);
    dist_sub_cruise_actual = V_subsonic * cruise_sub.t(end);

    % =========================================================================
    % PHASE 8: DESCENT TO PATTERN
    % =========================================================================
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('PHASE 8: DESCENT TO PATTERN (h=%.0f m to %.0f m)\n', h_subsonic, h_pattern);
    fprintf('════════════════════════════════════════════════════════════════\n');

    V_sub_end = M_subsonic * a_subsonic;
    descent2 = sim_descent_to_pattern(W_after_sub_cruise, S_ref, h_subsonic, V_sub_end, h_pattern, geom);
else
    % SKIP subsonic cruise - descend directly from supersonic cruise to pattern
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('PHASES 6-7: SKIPPED (No subsonic cruise - route too short)\n');
    fprintf('════════════════════════════════════════════════════════════════\n');

    dist_descent1_actual = 0;
    dist_sub_cruise_actual = 0;

    % Create dummy structures for skipped phases
    descent1 = struct('total_time', 0, 'total_fuel', 0, 'W_final', W_after_ss_cruise, 'distance', 0);
    cruise_sub = struct('t', 0, 'fuel_burned', 0, 'W', W_after_ss_cruise, 'LD', 0, 'TSFC', 0);
    W_after_sub_cruise = W_after_ss_cruise;

    % =========================================================================
    % PHASE 8: DESCENT DIRECTLY TO PATTERN (from supersonic cruise)
    % =========================================================================
    fprintf('\n════════════════════════════════════════════════════════════════\n');
    fprintf('PHASE 8: DESCENT TO PATTERN (h=%.0f m to %.0f m)\n', h_supersonic, h_pattern);
    fprintf('════════════════════════════════════════════════════════════════\n');

    descent2 = sim_descent_to_pattern(W_after_ss_cruise, S_ref, h_supersonic, V_supersonic, h_pattern, geom);
end
W_after_descent2 = descent2.W_final;
dist_descent2_actual = descent2.distance;

% =========================================================================
% PHASE 9: LANDING (at Luke AFB)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 9: LANDING (at Luke AFB, %.0f m elevation)\n', h_luke);
fprintf('════════════════════════════════════════════════════════════════\n');

landing = sim_landing_elevated(W_after_descent2, S_ref, descent2.h_final, ...
    descent2.V_final, h_luke, geom);
W_final = landing.W_final;
dist_landing_actual = landing.total_distance;

% =========================================================================
% MISSION SUMMARY
% =========================================================================
total_distance = dist_takeoff + dist_climb + dist_accel + ...
    dist_ss_cruise_actual + dist_descent1_actual + ...
    dist_sub_cruise_actual + dist_descent2_actual + dist_landing_actual;
total_fuel = TOGW - W_final;
fuel_remaining = fuel_capacity - total_fuel;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    EDWARDS → LUKE MISSION SUMMARY                      ║\n');
fprintf('╠════════════════════════════════════════════════════════════════════════╣\n');
fprintf('║ Phase                    │ Dist (km) │ Time (s) │ Fuel (kg) │ Weight  ║\n');
fprintf('╠══════════════════════════╪═══════════╪══════════╪═══════════╪═════════╣\n');
fprintf('║ 1. Taxi (Edwards)        │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    0, taxi.duration, taxi.fuel_burned, taxi.W_final);
fprintf('║ 2. Takeoff               │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_takeoff/1000, takeoff.total_time, takeoff.fuel_burned, takeoff.W_final);
fprintf('║ 3. Climb                 │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_climb/1000, climb.total_time, climb.total_fuel, climb.W_final);
fprintf('║ 4. Transonic Accel       │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_accel/1000, accel.total_time, accel.total_fuel, accel.W_final);
fprintf('║ 5. Supersonic Cruise     │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_ss_cruise_actual/1000, cruise_ss.t(end), cruise_ss.fuel_burned(end), W_after_ss_cruise);
fprintf('║ 6. Descent to Subsonic   │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_descent1_actual/1000, descent1.total_time, descent1.total_fuel, descent1.W_final);
fprintf('║ 7. Subsonic Cruise       │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_sub_cruise_actual/1000, cruise_sub.t(end), cruise_sub.fuel_burned(end), W_after_sub_cruise);
fprintf('║ 8. Descent to Pattern    │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_descent2_actual/1000, descent2.total_time, descent2.total_fuel, descent2.W_final);
fprintf('║ 9. Landing (Luke)        │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    dist_landing_actual/1000, landing.total_time, landing.total_fuel, landing.W_final);
fprintf('╠══════════════════════════╪═══════════╪══════════╪═══════════╪═════════╣\n');
fprintf('║ TOTAL                    │ %8.1f  │ %7.1f  │ %8.1f  │ %7.0f  ║\n', ...
    total_distance/1000, ...
    taxi.duration + takeoff.total_time + climb.total_time + accel.total_time + ...
    cruise_ss.t(end) + descent1.total_time + cruise_sub.t(end) + ...
    descent2.total_time + landing.total_time, ...
    total_fuel, W_final);
fprintf('╚════════════════════════════════════════════════════════════════════════╝\n');

fprintf('\n');
fprintf('MISSION METRICS:\n');
fprintf('  Route distance:      %.1f NM (%.1f km) - Target: %.1f NM\n', ...
    total_distance/1852, total_distance/1000, distance_nm);
fprintf('  Total mission time:  %.1f min\n', ...
    (taxi.duration + takeoff.total_time + climb.total_time + accel.total_time + ...
    cruise_ss.t(end) + descent1.total_time + cruise_sub.t(end) + ...
    descent2.total_time + landing.total_time) / 60);
fprintf('  Total fuel burned:   %.1f kg\n', total_fuel);
fprintf('  Fuel remaining:      %.1f kg (%.1f%%)\n', fuel_remaining, 100*fuel_remaining/fuel_capacity);
fprintf('\n');

fprintf('CRUISE COMPARISON:\n');
fprintf('  Supersonic (M=%.1f, h=%.0f km):\n', M_supersonic, h_supersonic/1000);
fprintf('    Duration: %.1f min, Distance: %.1f km, Fuel: %.1f kg\n', ...
    cruise_ss.t(end)/60, dist_ss_cruise_actual/1000, cruise_ss.fuel_burned(end));
fprintf('    Avg L/D: %.2f, Avg TSFC: %.3f lb/hr/lbf\n', mean(cruise_ss.LD), mean(cruise_ss.TSFC));
fprintf('  Subsonic (M=%.1f, h=%.1f km):\n', M_subsonic, h_subsonic/1000);
fprintf('    Duration: %.1f min, Distance: %.1f km, Fuel: %.1f kg\n', ...
    cruise_sub.t(end)/60, dist_sub_cruise_actual/1000, cruise_sub.fuel_burned(end));
fprintf('    Avg L/D: %.2f, Avg TSFC: %.3f lb/hr/lbf\n', mean(cruise_sub.LD), mean(cruise_sub.TSFC));
fprintf('\n');

% =========================================================================
% GENERATE PLOTS
% =========================================================================
fprintf('Generating plots...\n');

% Build combined trajectory arrays
% Time offsets for each phase
t_offset = 0;
t_taxi_end = taxi.duration;
t_takeoff_end = t_taxi_end + takeoff.total_time;
t_climb_end = t_takeoff_end + climb.total_time;
t_accel_end = t_climb_end + accel.total_time;
t_ss_cruise_end = t_accel_end + cruise_ss.t(end);
t_descent1_end = t_ss_cruise_end + descent1.total_time;
t_sub_cruise_end = t_descent1_end + cruise_sub.t(end);
t_descent2_end = t_sub_cruise_end + descent2.total_time;

% Combine altitude data
t_combined = [];
h_combined = [];
M_combined = [];
W_combined = [];
x_combined = [];

% Helper to append phase data
x_offset = 0;

% Taxi (ground)
t_combined = [t_combined; linspace(0, taxi.duration, 10)'];
h_combined = [h_combined; h_edwards * ones(10, 1)];
M_combined = [M_combined; zeros(10, 1)];
W_combined = [W_combined; linspace(TOGW, taxi.W_final, 10)'];
x_combined = [x_combined; zeros(10, 1)];

% Takeoff (compute M from V and h, W from linear interpolation)
t_combined = [t_combined; t_taxi_end + takeoff.t(:)];
h_combined = [h_combined; takeoff.h(:)];
% Compute Mach from V and h
M_takeoff = zeros(size(takeoff.V));
for ii = 1:length(takeoff.V)
    [~, a_to, ~, ~] = atmosisa(takeoff.h(ii));
    M_takeoff(ii) = takeoff.V(ii) / a_to;
end
M_combined = [M_combined; M_takeoff(:)];
% Linear interpolate weight during takeoff
W_takeoff = linspace(takeoff.W_initial, takeoff.W_final, length(takeoff.t))';
W_combined = [W_combined; W_takeoff];
x_combined = [x_combined; takeoff.x(:)];
x_offset = takeoff.x(end);

% Climb
t_combined = [t_combined; t_takeoff_end + climb.t(:)];
h_combined = [h_combined; climb.h(:)];
M_combined = [M_combined; climb.M(:)];
W_combined = [W_combined; climb.W(:)];
x_combined = [x_combined; x_offset + climb.x(:)];
x_offset = x_offset + climb.x(end);

% Acceleration (accel.h may be scalar - expand if needed)
t_combined = [t_combined; t_climb_end + accel.t(:)];
if isscalar(accel.h)
    h_accel = accel.h * ones(length(accel.t), 1);
else
    h_accel = accel.h(:);
end
h_combined = [h_combined; h_accel];
M_combined = [M_combined; accel.M(:)];
W_combined = [W_combined; accel.W(:)];
x_combined = [x_combined; x_offset + accel.x(:)];
x_offset = x_offset + accel.x(end);

% Supersonic cruise
t_combined = [t_combined; t_accel_end + cruise_ss.t(:)];
h_combined = [h_combined; h_supersonic * ones(length(cruise_ss.t), 1)];
[~, a_ss, ~, ~] = atmosisa(h_supersonic);
M_combined = [M_combined; M_supersonic * ones(length(cruise_ss.t), 1)];
W_combined = [W_combined; cruise_ss.W(:)];
V_ss = M_supersonic * a_ss;
x_ss_phase = cumsum([0; diff(cruise_ss.t(:))] * V_ss);
x_combined = [x_combined; x_offset + x_ss_phase];
x_offset = x_offset + x_ss_phase(end);

% Descent 1 (to subsonic)
if descent1.total_time > 0
    t_combined = [t_combined; t_ss_cruise_end + descent1.t(:)];
    h_combined = [h_combined; descent1.h(:)];
    M_combined = [M_combined; descent1.M(:)];
    W_combined = [W_combined; descent1.W(:)];
    x_combined = [x_combined; x_offset + descent1.x(:)];
    x_offset = x_offset + descent1.x(end);
end

% Subsonic cruise
if cruise_sub.t(end) > 0
    t_combined = [t_combined; t_descent1_end + cruise_sub.t(:)];
    h_combined = [h_combined; h_subsonic * ones(length(cruise_sub.t), 1)];
    M_combined = [M_combined; M_subsonic * ones(length(cruise_sub.t), 1)];
    W_combined = [W_combined; cruise_sub.W(:)];
    [~, a_sub, ~, ~] = atmosisa(h_subsonic);
    V_sub = M_subsonic * a_sub;
    x_sub_phase = cumsum([0; diff(cruise_sub.t(:))] * V_sub);
    x_combined = [x_combined; x_offset + x_sub_phase];
    x_offset = x_offset + x_sub_phase(end);
end

% Descent 2 (to pattern)
t_combined = [t_combined; t_sub_cruise_end + descent2.t(:)];
h_combined = [h_combined; descent2.h(:)];
M_combined = [M_combined; descent2.M(:)];
W_combined = [W_combined; descent2.W(:)];
x_combined = [x_combined; x_offset + descent2.x(:)];
x_offset = x_offset + descent2.x(end);

% Landing
t_combined = [t_combined; t_descent2_end + landing.t(:)];
h_combined = [h_combined; landing.h(:)];
% Estimate Mach for landing
[~, a_land, ~, ~] = atmosisa(h_luke);
M_land = landing.V(:) / a_land;
M_combined = [M_combined; M_land];
W_combined = [W_combined; landing.W(:)];
x_combined = [x_combined; x_offset + landing.x(:)];

% Verify all arrays have the same length
n_t = length(t_combined);
n_h = length(h_combined);
n_M = length(M_combined);
n_W = length(W_combined);
n_x = length(x_combined);

fprintf('  Array sizes: t=%d, h=%d, M=%d, W=%d, x=%d\n', n_t, n_h, n_M, n_W, n_x);

% Truncate to minimum length if there's a mismatch
n_min = min([n_t, n_h, n_M, n_W, n_x]);
if n_min < n_t
    fprintf('  WARNING: Truncating arrays to %d points due to size mismatch\n', n_min);
    t_combined = t_combined(1:n_min);
    h_combined = h_combined(1:n_min);
    M_combined = M_combined(1:n_min);
    W_combined = W_combined(1:n_min);
    x_combined = x_combined(1:n_min);
end

% Create figure with subplots
figure('Name', 'Edwards to Luke Mission Profile', 'Position', [100, 100, 1200, 800]);

% Altitude vs Distance
subplot(2, 2, 1);
plot(x_combined/1000, h_combined/1000, 'b-', 'LineWidth', 1.5);
hold on;
% Mark key points
plot(0, h_edwards/1000, 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(x_combined(end)/1000, h_luke/1000, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
xlabel('Distance (km)');
ylabel('Altitude (km)');
title('Mission Profile: Altitude vs Distance');
grid on;
legend('Flight Path', 'Edwards AFB', 'Luke AFB', 'Location', 'best');

% Mach vs Time
subplot(2, 2, 2);
plot(t_combined/60, M_combined, 'r-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Mach Number');
title('Mach Number vs Time');
grid on;

% Weight vs Time
subplot(2, 2, 3);
plot(t_combined/60, W_combined, 'g-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Weight (kg)');
title('Aircraft Weight vs Time');
grid on;

% Altitude vs Time
subplot(2, 2, 4);
plot(t_combined/60, h_combined/1000, 'b-', 'LineWidth', 1.5);
xlabel('Time (min)');
ylabel('Altitude (km)');
title('Altitude vs Time');
grid on;

% Save figure
saveas(gcf, fullfile(script_path, 'Edwards_to_Luke_mission_profile.png'));
fprintf('  Saved: Edwards_to_Luke_mission_profile.png\n');

% Second figure: Cruise performance comparison
figure('Name', 'Cruise Performance Comparison', 'Position', [150, 150, 1000, 400]);

subplot(1, 2, 1);
bar([1, 2], [mean(cruise_ss.LD), mean(cruise_sub.LD)]);
set(gca, 'XTickLabel', {'Supersonic (M=1.4)', 'Subsonic (M=0.6)'});
ylabel('L/D Ratio');
title('Lift-to-Drag Ratio Comparison');
grid on;

subplot(1, 2, 2);
bar([1, 2], [cruise_ss.fuel_burned(end), cruise_sub.fuel_burned(end)]);
set(gca, 'XTickLabel', {'Supersonic', 'Subsonic'});
ylabel('Fuel Burned (kg)');
title('Fuel Consumption by Cruise Phase');
grid on;

saveas(gcf, fullfile(script_path, 'Edwards_to_Luke_cruise_comparison.png'));
fprintf('  Saved: Edwards_to_Luke_cruise_comparison.png\n');

% =========================================================================
% GEOGRAPHIC VISUALIZATION (Geoplot 3D)
% =========================================================================
fprintf('\nGenerating geographic flight path visualization...\n');

% Airfield coordinates
% Edwards AFB, CA
lat_edwards = 34.9054;
lon_edwards = -117.8839;

% Luke AFB, AZ
lat_luke = 33.5386;
lon_luke = -112.3830;

% Calculate flight path (great circle interpolation)
n_geo_points = length(x_combined);

% Interpolate lat/lon based on distance traveled
x_total_m = x_combined(end);

% Linear interpolation along great circle path
lat_path = zeros(n_geo_points, 1);
lon_path = zeros(n_geo_points, 1);

for i = 1:n_geo_points
    frac = x_combined(i) / x_total_m;
    frac = max(0, min(1, frac));
    lat_path(i) = lat_edwards + frac * (lat_luke - lat_edwards);
    lon_path(i) = lon_edwards + frac * (lon_luke - lon_edwards);
end

% Altitude in meters for 3D plot
alt_path = h_combined;

% =========================================================================
% FLIGHT PROFILE PLOT - Distance vs Altitude (Side View)
% =========================================================================
figure('Name', 'Edwards to Luke Flight Profile', 'Position', [200, 200, 1200, 500]);

% Downsample for faster plotting (every 10th point)
ds = 10;  % Downsample factor
idx_ds = 1:ds:n_geo_points;
if idx_ds(end) ~= n_geo_points
    idx_ds = [idx_ds, n_geo_points];  % Include last point
end

lon_ds = lon_path(idx_ds);
lat_ds = lat_path(idx_ds);
alt_ds = alt_path(idx_ds);
M_ds = M_combined(idx_ds);

% Convert to distance-based coordinates for visualization
x_km = x_combined(idx_ds) / 1000;  % Distance in km
alt_km = alt_ds / 1000;            % Altitude in km

% Categorize by Mach regime for batch plotting
idx_ground = M_ds < 0.3;
idx_subsonic = M_ds >= 0.3 & M_ds < 0.9;
idx_transonic = M_ds >= 0.9 & M_ds < 1.0;
idx_supersonic = M_ds >= 1.0;

hold on;

% Ensure column vectors
x_km = x_km(:);
alt_km = alt_km(:);

% Plot flight path as thick line
plot(x_km, alt_km, '-', 'Color', [0.2 0.2 0.2], 'LineWidth', 2.5);

% Overlay colored markers by regime
if any(idx_ground)
    plot(x_km(idx_ground), alt_km(idx_ground), 'o', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6);
end
if any(idx_subsonic)
    plot(x_km(idx_subsonic), alt_km(idx_subsonic), 'o', 'Color', [0 0.6 1], 'MarkerFaceColor', [0 0.6 1], 'MarkerSize', 6);
end
if any(idx_transonic)
    plot(x_km(idx_transonic), alt_km(idx_transonic), 'o', 'Color', [1 0.6 0], 'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 7);
end
if any(idx_supersonic)
    plot(x_km(idx_supersonic), alt_km(idx_supersonic), 'o', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 6);
end

% Mark airfields
plot(0, h_edwards/1000, 'g^', 'MarkerSize', 14, 'MarkerFaceColor', 'g', 'LineWidth', 2);
plot(x_km(end), h_luke/1000, 'rv', 'MarkerSize', 14, 'MarkerFaceColor', 'r', 'LineWidth', 2);

% Add ground line
plot([0, max(x_km)], [0, 0], 'k-', 'LineWidth', 2);

% Labels and formatting
xlabel('Distance from Edwards AFB (km)', 'FontSize', 12);
ylabel('Altitude (km)', 'FontSize', 12);
title(sprintf('X-29 Flight Profile: Edwards AFB → Luke AFB (%.0f NM)', distance_nm), 'FontSize', 14);

grid on;
box on;

% Set axis limits
xlim([-10, max(x_km)+10]);
ylim([-0.5, max(alt_km)*1.2]);

% Add text annotations
text(15, h_edwards/1000 + 1.2, sprintf('Edwards AFB\n(%.0f m)', h_edwards), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'left');
text(x_km(end)-15, h_luke/1000 + 1.2, sprintf('Luke AFB\n(%.0f m)', h_luke), 'FontSize', 9, 'FontWeight', 'bold', 'HorizontalAlignment', 'right');

% Add altitude reference lines
yline(12.2, '--', 'Cruise: 12.2 km (40,000 ft)', 'Color', [0.4 0.4 0.4], 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'FontSize', 9);
yline(9.12, ':', 'Subsonic: 9.1 km (30,000 ft)', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'LabelHorizontalAlignment', 'left', 'FontSize', 9);

% Create legend
h1 = plot(NaN, NaN, 's', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 10);
h2 = plot(NaN, NaN, 's', 'Color', [0 0.6 1], 'MarkerFaceColor', [0 0.6 1], 'MarkerSize', 10);
h3 = plot(NaN, NaN, 's', 'Color', [1 0.6 0], 'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 10);
h4 = plot(NaN, NaN, 's', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 10);
legend([h1, h2, h3, h4], {'Ground/Taxi (M<0.3)', 'Subsonic (M<0.9)', 'Transonic (0.9-1.0)', 'Supersonic (M>1.0)'}, ...
    'Location', 'east');

set(gca, 'FontSize', 11);

% Add note about vertical exaggeration
text(x_km(end)/2, -0.3, 'Note: Vertical scale exaggerated ~40x for visibility', ...
    'FontSize', 8, 'HorizontalAlignment', 'center', 'FontAngle', 'italic', 'Color', [0.5 0.5 0.5]);

saveas(gcf, fullfile(script_path, 'Edwards_to_Luke_flight_profile.png'));
fprintf('  Saved: Edwards_to_Luke_flight_profile.png\n');

% =========================================================================
% 2D GEOGRAPHIC PLOT (overhead view with Mach coloring)
% =========================================================================
figure('Name', 'Edwards to Luke - Overhead View', 'Position', [250, 250, 1000, 700]);

% Use downsampled data (already computed above)
% lon_ds, lat_ds, M_ds, idx_ground, idx_subsonic, idx_transonic, idx_supersonic

% Try geoplot, fallback to regular plot
try
    gx = geoaxes;
    hold(gx, 'on');

    % Plot base path
    geoplot(gx, lat_ds, lon_ds, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);

    % Overlay colored scatter points by regime
    if any(idx_ground)
        geoscatter(gx, lat_ds(idx_ground), lon_ds(idx_ground), 30, [0.5 0.5 0.5], 'filled');
    end
    if any(idx_subsonic)
        geoscatter(gx, lat_ds(idx_subsonic), lon_ds(idx_subsonic), 30, [0 0.6 1], 'filled');
    end
    if any(idx_transonic)
        geoscatter(gx, lat_ds(idx_transonic), lon_ds(idx_transonic), 40, [1 0.6 0], 'filled');
    end
    if any(idx_supersonic)
        geoscatter(gx, lat_ds(idx_supersonic), lon_ds(idx_supersonic), 30, [1 0 0], 'filled');
    end

    % Mark airfields
    geoplot(gx, lat_edwards, lon_edwards, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    geoplot(gx, lat_luke, lon_luke, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

    geobasemap(gx, 'streets-light');
    title(gx, 'X-29 Flight Path Colored by Mach Number', 'FontSize', 14);

    % Create legend
    h1 = plot(NaN, NaN, 's', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 10);
    h2 = plot(NaN, NaN, 's', 'Color', [0 0.6 1], 'MarkerFaceColor', [0 0.6 1], 'MarkerSize', 10);
    h3 = plot(NaN, NaN, 's', 'Color', [1 0.6 0], 'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 10);
    h4 = plot(NaN, NaN, 's', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 10);
    legend([h1, h2, h3, h4], {'Ground/Taxi', 'Subsonic (M<0.9)', 'Transonic', 'Supersonic (M≥1.0)'}, ...
        'Location', 'southwest');

catch
    % Fallback to regular 2D plot
    fprintf('  geoplot not available, using standard 2D plot...\n');

    hold on;

    % Plot base path
    plot(lon_ds, lat_ds, '-', 'Color', [0.8 0.8 0.8], 'LineWidth', 1);

    % Overlay colored scatter points by regime
    if any(idx_ground)
        scatter(lon_ds(idx_ground), lat_ds(idx_ground), 30, [0.5 0.5 0.5], 'filled');
    end
    if any(idx_subsonic)
        scatter(lon_ds(idx_subsonic), lat_ds(idx_subsonic), 30, [0 0.6 1], 'filled');
    end
    if any(idx_transonic)
        scatter(lon_ds(idx_transonic), lat_ds(idx_transonic), 40, [1 0.6 0], 'filled');
    end
    if any(idx_supersonic)
        scatter(lon_ds(idx_supersonic), lat_ds(idx_supersonic), 30, [1 0 0], 'filled');
    end

    plot(lon_edwards, lat_edwards, 'ko', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
    plot(lon_luke, lat_luke, 'ks', 'MarkerSize', 12, 'MarkerFaceColor', 'r');

    xlabel('Longitude (°)');
    ylabel('Latitude (°)');
    title('X-29 Flight Path Colored by Mach Number', 'FontSize', 14);
    grid on;
    axis equal;

    h1 = plot(NaN, NaN, 's', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 10);
    h2 = plot(NaN, NaN, 's', 'Color', [0 0.6 1], 'MarkerFaceColor', [0 0.6 1], 'MarkerSize', 10);
    h3 = plot(NaN, NaN, 's', 'Color', [1 0.6 0], 'MarkerFaceColor', [1 0.6 0], 'MarkerSize', 10);
    h4 = plot(NaN, NaN, 's', 'Color', [1 0 0], 'MarkerFaceColor', [1 0 0], 'MarkerSize', 10);
    legend([h1, h2, h3, h4], {'Ground/Taxi', 'Subsonic (M<0.9)', 'Transonic', 'Supersonic (M≥1.0)'}, ...
        'Location', 'southwest');
end

saveas(gcf, fullfile(script_path, 'Edwards_to_Luke_geoplot_mach.png'));
fprintf('  Saved: Edwards_to_Luke_geoplot_mach.png\n');

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('EDWARDS → LUKE MISSION SIMULATION COMPLETE\n');
fprintf('════════════════════════════════════════════════════════════════\n');
