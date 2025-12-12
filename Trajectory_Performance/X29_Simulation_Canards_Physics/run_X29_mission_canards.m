%{
====================================================================================================
SCRIPT: run_X29_mission_canards.m
AUTHOR: Dominic Larin (Chief Engineer), Damian Do (Trajectory Lead), Carlos Gonzalez (Trajectory Secondary)
INITIATED: 12/05/2025
LAST MODIFIED: 12/06/2025
====================================================================================================
DESCRIPTION:
Runs complete X-29 mission simulation with canard/stability tracking.
Uses sim_*_canards functions with RAYMER AERODYNAMICS MODEL (get_X29_aero) for accurate
drag calculations based on Mach, altitude, and angle of attack.

Mission Profile:
  1. Taxi + Warmup (5 min at idle)
  2. Takeoff (sim_takeoff_canards - Raymer aero + gear drag)
  3. Climb (sim_climb_canards - Raymer aero)
  4. Transonic Acceleration (sim_accel - Raymer aero)
  5. Cruise (sim_cruise_canards - Raymer aero, 10 min at M=1.4, h=12.2km)
  6. Descent (sim_descent_canards - Raymer aero)
  7. Landing (sim_landing_canards - Raymer aero + gear/flap drag)

Tracking Throughout Mission:
  - CG position (fuel burn effect)
  - Canard deflection for trim
  - Cm_alpha (pitch stability derivative)
  - Alpha (angle of attack)
  - Fuel consumption based on physics-based F404 model
====================================================================================================
%}

clear; clc; close all;

% Add paths - get project root (2 levels up from this script)
% X29_Simulation_Canards_Physics -> Trajectory_Performance -> MAE4350-X29-PROJECT
script_path = fileparts(mfilename('fullpath'));
proj_root = fileparts(fileparts(script_path));

% Add the rest of the project first
addpath(genpath(proj_root));

% Then add this folder LAST so its functions take priority over other folders
% (addpath adds to the BEGINNING of the path, so the last call has highest priority)
addpath(script_path);

fprintf('╔════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║       X-29 FULL MISSION SIMULATION WITH CANARD/STABILITY TRACKING      ║\n');
fprintf('║              Using Raymer Aero Model (get_X29_aero)                    ║\n');
fprintf('╚════════════════════════════════════════════════════════════════════════╝\n\n');

% =========================================================================
% LOAD CONVERGED VEHICLE GEOMETRY FROM TRADE STUDY
% =========================================================================
fprintf('Loading converged vehicle geometry...\n');

% Load trade study data
mat_file = fullfile(proj_root, 'Synthesis', 'updated_fuselage_calculations', ...
    'elliptical_trade_study_20251117_170903.mat');

if ~exist(mat_file, 'file')
    % Try alternate location
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

% Build geometry structure for Raymer aero model
geom = build_X29_geometry(veh, trade_data.vehicle_base);

fprintf('  Geometry loaded successfully.\n\n');

% =========================================================================
% VEHICLE PARAMETERS (from converged trade study)
% =========================================================================
TOGW = veh.TOGW_kg;
S_ref = veh.wing_area_m2;
fuel_capacity = 1900;
OEW = 6278;  % kg (from Calculations.xlsx)

fprintf('VEHICLE (from converged trade study):\n');
fprintf('  TOGW:          %.0f kg\n', TOGW);
fprintf('  OEW:           %.0f kg\n', OEW);
fprintf('  Wing area:     %.2f m²\n', S_ref);
fprintf('  Fuel capacity: %.0f kg\n', fuel_capacity);
fprintf('  L/D cruise:    %.2f (synthesis)\n', veh.LD_cruise);
fprintf('\n');

% =========================================================================
% PHASE 1: TAXI + WARMUP
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 1: TAXI + WARMUP\n');
fprintf('════════════════════════════════════════════════════════════════\n');

taxi = simulate_X29_taxi(300, TOGW);  % 5 minutes (no canards version needed)
W_after_taxi = taxi.W_final;
W_fuel_after_taxi = fuel_capacity - taxi.fuel_burned;

% =========================================================================
% PHASE 2: TAKEOFF (using sim_takeoff_canards with Raymer aero)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 2: TAKEOFF\n');
fprintf('════════════════════════════════════════════════════════════════\n');

takeoff_params = struct();
takeoff_params.W_total = W_after_taxi;
takeoff_params.W_fuel = W_fuel_after_taxi;
takeoff_params.S_ref = S_ref;
takeoff_params.CL_max = 1.6;
takeoff_params.verbose = true;

takeoff = sim_takeoff_canards(takeoff_params, geom);
W_after_takeoff = takeoff.W_final;
W_fuel_after_takeoff = takeoff.W_fuel_final;

% Add compatibility fields
takeoff.fuel_burned = takeoff.W(1) - takeoff.W(end);

% =========================================================================
% PHASE 3: CLIMB (using sim_climb_canards with ODE45 Hull equations)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 3: CLIMB\n');
fprintf('════════════════════════════════════════════════════════════════\n');

% Use sim_climb_canards with ODE45 Hull equations + canard tracking
climb_params = struct();
climb_params.W_initial = W_after_takeoff;
climb_params.W_fuel_initial = W_fuel_after_takeoff;
climb_params.h_initial = 15.24;  % 50 ft after takeoff
climb_params.V_initial = takeoff.V_final;
climb_params.h_target = 12200;
climb_params.M_target = 0.9;
climb_params.verbose = true;

climb = sim_climb_canards(climb_params, geom);
W_after_climb = climb.W_final;
W_fuel_after_climb = climb.W_fuel_final;

% =========================================================================
% PHASE 4: TRANSONIC ACCELERATION (with canards tracking)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 4: TRANSONIC ACCELERATION\n');
fprintf('════════════════════════════════════════════════════════════════\n');

accel = sim_accel_canards(W_after_climb, S_ref, 12200, 0.9, 1.4, geom);
W_after_accel = accel.W_final;
W_fuel_after_accel = W_fuel_after_climb - accel.total_fuel;

% =========================================================================
% PHASE 5: CRUISE (using sim_cruise_canards with Raymer aero)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 5: CRUISE (10 minutes)\n');
fprintf('════════════════════════════════════════════════════════════════\n');

cruise_params = struct();
cruise_params.W_initial = W_after_accel;
cruise_params.W_fuel_initial = W_fuel_after_accel;
cruise_params.h = 12200;
cruise_params.M = 1.4;
cruise_params.duration_min = 10;
cruise_params.verbose = true;

cruise = sim_cruise_canards(cruise_params, geom);
W_after_cruise = cruise.W_final;
W_fuel_after_cruise = cruise.W_fuel_final;

% Add compatibility fields for downstream code
cruise.fuel_burned = cruise.W(1) - cruise.W;  % Array of cumulative fuel burned
cruise.conditions.M = 1.4;
cruise.aero.CD0 = mean(cruise.CD) * 0.6;  % Approximate parasitic portion
cruise.aero.K = 0.15;  % Approximate induced drag factor

% =========================================================================
% PHASE 6: DESCENT (using sim_descent_canards with Raymer aero)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 6: DESCENT\n');
fprintf('════════════════════════════════════════════════════════════════\n');

descent_params = struct();
descent_params.W_initial = W_after_cruise;
descent_params.W_fuel_initial = W_fuel_after_cruise;
descent_params.h_initial = 12200;
descent_params.h_final = 610;  % Pattern altitude (2000 ft)
descent_params.verbose = true;

descent = sim_descent_canards(descent_params, geom);
W_after_descent = descent.W_final;
W_fuel_after_descent = descent.W_fuel_final;

% Add compatibility fields for downstream code
descent.total_time = descent.t(end);
descent.total_fuel = descent.W(1) - descent.W(end);
descent.h_final = descent.h(end);
descent.V_final = descent.V(end);

% =========================================================================
% PHASE 7: LANDING (using sim_landing_canards with Raymer aero)
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('PHASE 7: LANDING\n');
fprintf('════════════════════════════════════════════════════════════════\n');

landing_params = struct();
landing_params.W_total = W_after_descent;
landing_params.W_fuel = W_fuel_after_descent;
landing_params.h_pattern = descent.h_final;
landing_params.CL_max = 1.6;
landing_params.verbose = true;

landing = sim_landing_canards(landing_params, geom);
W_final = landing.W_final;

% Add compatibility fields for downstream code
if isfield(landing, 'fuel_consumed_total')
    landing.total_fuel = landing.fuel_consumed_total;
else
    landing.total_fuel = landing.W(1) - landing.W(end);
end

% =========================================================================
% STITCH TIMELINE TOGETHER
% =========================================================================
fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('STITCHING MISSION TIMELINE...\n');
fprintf('════════════════════════════════════════════════════════════════\n');

% Time offsets
t_offset_taxi = 0;
t_offset_takeoff = taxi.duration;
t_offset_climb = t_offset_takeoff + takeoff.total_time;
t_offset_accel = t_offset_climb + climb.total_time;
t_offset_cruise = t_offset_accel + accel.total_time;
t_offset_descent = t_offset_cruise + cruise.t(end);
t_offset_landing = t_offset_descent + descent.total_time;

% Create unified timeline arrays
% Taxi (simplified - just start and end points)
t_taxi = [0; taxi.duration];
h_taxi = [0; 0];
V_taxi = [0; 0];
M_taxi = [0; 0];
W_taxi = [TOGW; taxi.W_final];

% Takeoff (using full arrays from sim_takeoff_canards)
t_takeoff = takeoff.t + t_offset_takeoff;
h_takeoff = takeoff.h;
V_takeoff = takeoff.V;
% Compute Mach from V and h
M_takeoff = zeros(size(takeoff.V));
for ii = 1:length(takeoff.V)
    [~, a_to, ~, ~] = atmosisa(max(0, takeoff.h(ii)));
    M_takeoff(ii) = takeoff.V(ii) / a_to;
end
W_takeoff = takeoff.W;

% Climb (using full arrays from sim_climb_canards)
t_climb = climb.t + t_offset_climb;
h_climb = climb.h;
V_climb = climb.V;
M_climb = climb.M;
W_climb = climb.W;

% Acceleration
t_accel = accel.t + t_offset_accel;
% Handle scalar vs array for h
if isscalar(accel.h)
    h_accel = accel.h * ones(size(accel.t));
else
    h_accel = accel.h;
end
V_accel = accel.V;
M_accel = accel.M;
W_accel = accel.W;

% Cruise (using full arrays from sim_cruise_canards)
t_cruise = cruise.t + t_offset_cruise;
h_cruise = cruise.h;
V_cruise_vec = cruise.V;
M_cruise_vec = cruise.M;
W_cruise = cruise.W;

% Descent
t_descent = descent.t + t_offset_descent;
h_descent = descent.h;
V_descent = descent.V;
M_descent = descent.M;
W_descent = descent.W;

% Landing
t_landing = landing.t + t_offset_landing;
h_landing = landing.h;
V_landing = landing.V;
M_landing = zeros(size(landing.V));
for i = 1:length(landing.V)
    [~, a_land, ~, ~] = atmosisa(max(0, landing.h(i)));
    M_landing(i) = landing.V(i) / a_land;
end
W_landing = landing.W;

% =========================================================================
% CG AND CANARD TRACKING WITH SMOOTHING
% =========================================================================
fprintf('Calculating CG and canard deflection throughout mission...\n');

% Transition period parameters
T_blend = 10.0;  % Transition period duration [s] for blending between phases

% Helper function to calculate canard deflection for trim
calculate_canard_trim = @(M, alpha_deg, Cm_alpha, Cm_delta_c) ...
    max(-32, min(58, -Cm_alpha * (alpha_deg * pi/180) / Cm_delta_c * 180/pi));

% Helper function for linear transition blending
% t_in_phase: time since phase started
% alpha_start: alpha at end of previous phase
% alpha_target: calculated alpha for current conditions
% T_blend: transition period duration
blend_alpha = @(t_in_phase, alpha_start, alpha_target, T_blend) ...
    alpha_start + (alpha_target - alpha_start) * min(1, t_in_phase / T_blend);

% =========================================================================
% TAXI - Ground, alpha = 0
% =========================================================================
cg_taxi = zeros(size(W_taxi));
canard_taxi = zeros(size(W_taxi));
alpha_taxi = zeros(size(W_taxi));
Cm_alpha_taxi = zeros(size(W_taxi));
for i = 1:length(W_taxi)
    fuel_rem = W_taxi(i) - OEW;
    [cg_taxi(i), ~] = calculate_X29_CG(W_taxi(i), fuel_rem, fuel_capacity);
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(0);
    Cm_alpha_taxi(i) = Cm_alpha;
    alpha_taxi(i) = 0;  % On ground
    canard_taxi(i) = 0;
end
alpha_phase_start = 0;  % Alpha at end of taxi (start of takeoff)

% =========================================================================
% TAKEOFF - Rotation to liftoff (with transition blending)
% =========================================================================
cg_takeoff = zeros(size(W_takeoff));
canard_takeoff = zeros(size(W_takeoff));
alpha_takeoff = zeros(size(W_takeoff));
Cm_alpha_takeoff = zeros(size(W_takeoff));
t_phase_start = t_takeoff(1);  % Time at start of this phase

for i = 1:length(W_takeoff)
    fuel_rem = W_takeoff(i) - OEW;
    [cg_takeoff(i), ~] = calculate_X29_CG(W_takeoff(i), fuel_rem, fuel_capacity);
    M_val = max(0.1, M_takeoff(i));
    V_val = max(10, V_takeoff(i));
    h_val = max(0, h_takeoff(i));

    % Calculate target alpha for takeoff (rotation to 8-10 deg)
    if i == 1
        alpha_target = 0;  % Start of ground roll
    else
        % During rotation/liftoff, calculate trim alpha
        [alpha_target, ~, ~] = get_trim_alpha(W_takeoff(i), V_val, h_val, S_ref, 5, 1.0);
        alpha_target = min(alpha_target, 10);  % Cap for takeoff
    end

    % Apply transition blending
    t_in_phase = t_takeoff(i) - t_phase_start;
    alpha_takeoff(i) = blend_alpha(t_in_phase, alpha_phase_start, alpha_target, T_blend);

    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(max(0.3, M_val));
    Cm_alpha_takeoff(i) = Cm_alpha;
    if abs(Cm_delta_c) > 0.01
        canard_takeoff(i) = calculate_canard_trim(M_val, alpha_takeoff(i), Cm_alpha, Cm_delta_c);
    else
        canard_takeoff(i) = 10;
    end
end
alpha_phase_start = alpha_takeoff(end);  % Carry forward to next phase

% =========================================================================
% CLIMB - Solved alpha with transition blending
% =========================================================================
cg_climb = zeros(size(W_climb));
canard_climb = zeros(size(W_climb));
alpha_climb = zeros(size(W_climb));
Cm_alpha_climb = zeros(size(W_climb));
t_phase_start = t_climb(1);  % Time at start of this phase

for i = 1:length(W_climb)
    fuel_rem = W_climb(i) - OEW;
    [cg_climb(i), ~] = calculate_X29_CG(W_climb(i), fuel_rem, fuel_capacity);
    M_val = M_climb(i);
    V_val = V_climb(i);
    h_val = h_climb(i);

    % Estimate climb angle from altitude change rate
    if i > 1
        dh = h_climb(i) - h_climb(i-1);
        dt_c = t_climb(i) - t_climb(i-1);
        if dt_c > 0 && V_val > 0
            gamma_climb = asind(dh / (dt_c * V_val));
            gamma_climb = max(-10, min(30, gamma_climb));
        else
            gamma_climb = 10;
        end
    else
        gamma_climb = 10;
    end

    % Calculate target alpha for trim during climb
    [alpha_target, ~, ~] = get_trim_alpha(W_climb(i), V_val, h_val, S_ref, gamma_climb, 1.0);

    % Apply transition blending at phase start, then follow target
    t_in_phase = t_climb(i) - t_phase_start;
    alpha_climb(i) = blend_alpha(t_in_phase, alpha_phase_start, alpha_target, T_blend);

    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(max(0.3, M_val));
    Cm_alpha_climb(i) = Cm_alpha;
    if abs(Cm_delta_c) > 0.01
        canard_climb(i) = calculate_canard_trim(M_val, alpha_climb(i), Cm_alpha, Cm_delta_c);
    else
        canard_climb(i) = 8;
    end
end
alpha_phase_start = alpha_climb(end);  % Carry forward to next phase

% =========================================================================
% ACCEL - Level acceleration (with transition blending)
% =========================================================================
cg_accel = zeros(size(W_accel));
canard_accel = zeros(size(W_accel));
alpha_accel = zeros(size(W_accel));
Cm_alpha_accel = zeros(size(W_accel));
t_phase_start = t_accel(1);  % Time at start of this phase

for i = 1:length(W_accel)
    fuel_rem = W_accel(i) - OEW;
    [cg_accel(i), ~] = calculate_X29_CG(W_accel(i), fuel_rem, fuel_capacity);
    M_val = M_accel(i);
    V_val = V_accel(i);
    h_val = h_accel(i);

    % Level acceleration - calculate trim alpha
    [alpha_target, ~, ~] = get_trim_alpha(W_accel(i), V_val, h_val, S_ref, 0, 1.0);

    % Apply transition blending at phase start, then follow target
    t_in_phase = t_accel(i) - t_phase_start;
    alpha_accel(i) = blend_alpha(t_in_phase, alpha_phase_start, alpha_target, T_blend);

    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_val);
    Cm_alpha_accel(i) = Cm_alpha;
    if abs(Cm_delta_c) > 0.01
        canard_accel(i) = calculate_canard_trim(M_val, alpha_accel(i), Cm_alpha, Cm_delta_c);
    else
        canard_accel(i) = 5;
    end
end
alpha_phase_start = alpha_accel(end);  % Carry forward to next phase

% =========================================================================
% CRUISE - Use simulation alpha with transition blending
% =========================================================================
cg_cruise = zeros(size(W_cruise));
canard_cruise = zeros(size(W_cruise));
alpha_cruise = zeros(size(W_cruise));
Cm_alpha_cruise = zeros(size(W_cruise));
[Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(1.4);
t_phase_start = t_cruise(1);  % Time at start of this phase

for i = 1:length(W_cruise)
    fuel_rem = W_cruise(i) - OEW;
    [cg_cruise(i), ~] = calculate_X29_CG(W_cruise(i), fuel_rem, fuel_capacity);

    % Use cruise simulation's calculated alpha as target
    alpha_target = cruise.alpha(i);

    % Apply transition blending at phase start, then follow target
    t_in_phase = t_cruise(i) - t_phase_start;
    alpha_cruise(i) = blend_alpha(t_in_phase, alpha_phase_start, alpha_target, T_blend);

    Cm_alpha_cruise(i) = Cm_alpha;
    if abs(Cm_delta_c) > 0.01
        canard_cruise(i) = calculate_canard_trim(1.4, alpha_cruise(i), Cm_alpha, Cm_delta_c);
    else
        canard_cruise(i) = 5;
    end
end
alpha_phase_start = alpha_cruise(end);  % Carry forward to next phase

% =========================================================================
% DESCENT - Solved alpha with transition blending
% =========================================================================
cg_descent = zeros(size(W_descent));
canard_descent = zeros(size(W_descent));
alpha_descent = zeros(size(W_descent));
Cm_alpha_descent = zeros(size(W_descent));
t_phase_start = t_descent(1);  % Time at start of this phase

for i = 1:length(W_descent)
    fuel_rem = W_descent(i) - OEW;
    [cg_descent(i), ~] = calculate_X29_CG(W_descent(i), fuel_rem, fuel_capacity);
    M_val = M_descent(i);
    V_val = V_descent(i);
    h_val = h_descent(i);

    % Estimate descent angle
    if i > 1
        dh = h_descent(i) - h_descent(i-1);
        dt_d = t_descent(i) - t_descent(i-1);
        if dt_d > 0 && V_val > 0
            gamma_desc = asind(dh / (dt_d * V_val));
            gamma_desc = max(-30, min(5, gamma_desc));
        else
            gamma_desc = -3;
        end
    else
        gamma_desc = -3;
    end

    % Calculate target alpha
    [alpha_target, ~, ~] = get_trim_alpha(W_descent(i), V_val, h_val, S_ref, gamma_desc, 1.0);

    % Apply transition blending at phase start, then follow target
    t_in_phase = t_descent(i) - t_phase_start;
    alpha_descent(i) = blend_alpha(t_in_phase, alpha_phase_start, alpha_target, T_blend);

    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(max(0.3, M_val));
    Cm_alpha_descent(i) = Cm_alpha;
    if abs(Cm_delta_c) > 0.01
        canard_descent(i) = calculate_canard_trim(M_val, alpha_descent(i), Cm_alpha, Cm_delta_c);
    else
        canard_descent(i) = 3;
    end
end
alpha_phase_start = alpha_descent(end);  % Carry forward to next phase

% =========================================================================
% LANDING - Solved alpha with transition blending
% =========================================================================
cg_landing = zeros(size(W_landing));
canard_landing = zeros(size(W_landing));
alpha_landing = zeros(size(W_landing));
Cm_alpha_landing = zeros(size(W_landing));
t_phase_start = t_landing(1);  % Time at start of this phase

for i = 1:length(W_landing)
    fuel_rem = W_landing(i) - OEW;
    [cg_landing(i), ~] = calculate_X29_CG(W_landing(i), fuel_rem, fuel_capacity);
    M_val = max(0.15, M_landing(i));
    V_val = max(10, V_landing(i));
    h_val = max(0, h_landing(i));

    % Landing phases: approach (higher alpha), flare, touchdown, rollout
    % Estimate from flight path
    if i > 1
        dh = h_landing(i) - h_landing(i-1);
        dt_l = t_landing(i) - t_landing(i-1);
        if dt_l > 0 && V_val > 0
            gamma_land = asind(dh / (dt_l * V_val));
            gamma_land = max(-10, min(5, gamma_land));
        else
            gamma_land = -3;
        end
    else
        gamma_land = -3;  % Approach
    end

    % Calculate target alpha (higher during approach/flare)
    if h_val > 5  % In the air
        [alpha_target, ~, ~] = get_trim_alpha(W_landing(i), V_val, h_val, S_ref, gamma_land, 1.0);
        alpha_target = min(alpha_target, 12);  % Cap for approach
    else  % On ground
        alpha_target = 0;  % Nosedown after touchdown
    end

    % Apply transition blending at phase start, then follow target
    t_in_phase = t_landing(i) - t_phase_start;
    alpha_landing(i) = blend_alpha(t_in_phase, alpha_phase_start, alpha_target, T_blend);

    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_val);
    Cm_alpha_landing(i) = Cm_alpha;
    if abs(Cm_delta_c) > 0.01
        canard_landing(i) = calculate_canard_trim(M_val, alpha_landing(i), Cm_alpha, Cm_delta_c);
    else
        canard_landing(i) = 15;
    end
end

% Get initial and final CG data for summary
[~, cg_data_init] = calculate_X29_CG(TOGW, fuel_capacity, fuel_capacity);
[~, cg_data_final] = calculate_X29_CG(W_final, W_final - OEW, fuel_capacity);

fprintf('  Initial CG: %.2f m (%.2f ft) from nose cone - %.1f%% MAC\n', ...
    cg_data_init.x_cg_m, cg_data_init.x_cg_ft, cg_data_init.x_cg_MAC_percent);
fprintf('  Final CG:   %.2f m (%.2f ft) from nose cone - %.1f%% MAC\n', ...
    cg_data_final.x_cg_m, cg_data_final.x_cg_ft, cg_data_final.x_cg_MAC_percent);
fprintf('  CG travel:  %.3f m (%.2f ft) - %.1f%% MAC\n', ...
    abs(cg_data_init.x_cg_m - cg_data_final.x_cg_m), ...
    abs(cg_data_init.x_cg_ft - cg_data_final.x_cg_ft), ...
    abs(cg_data_init.x_cg_MAC_percent - cg_data_final.x_cg_MAC_percent));

% =========================================================================
% MISSION SUMMARY
% =========================================================================
total_time = t_landing(end);
total_fuel = TOGW - W_final;
fuel_remaining = fuel_capacity - total_fuel;

fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                         MISSION SUMMARY                                ║\n');
fprintf('╠════════════════════════════════════════════════════════════════════════╣\n');
fprintf('║ Phase              │ Time (s) │ Fuel (kg) │ Weight (kg) │ Alt (m)     ║\n');
fprintf('╠════════════════════╪══════════╪═══════════╪═════════════╪═════════════╣\n');
fprintf('║ 1. Taxi + Warmup   │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    taxi.duration, taxi.fuel_burned, taxi.W_final, 0);
fprintf('║ 2. Takeoff         │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    takeoff.total_time, takeoff.fuel_burned, takeoff.W_final, 15);
fprintf('║ 3. Climb           │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    climb.total_time, climb.total_fuel, climb.W_final, 12200);
fprintf('║ 4. Accel (M→1.4)   │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    accel.total_time, accel.total_fuel, accel.W_final, 12200);
fprintf('║ 5. Cruise (10 min) │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    cruise.t(end), cruise.fuel_burned(end), W_after_cruise, 12200);
fprintf('║ 6. Descent         │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    descent.total_time, descent.total_fuel, descent.W_final, descent.h_final);
fprintf('║ 7. Landing         │ %7.1f  │ %8.1f  │ %10.0f  │ %10.0f  ║\n', ...
    landing.total_time, landing.total_fuel, landing.W_final, 0);
fprintf('╠════════════════════╪══════════╪═══════════╪═════════════╪═════════════╣\n');
fprintf('║ TOTAL              │ %7.1f  │ %8.1f  │ %10.0f  │             ║\n', ...
    total_time, total_fuel, W_final);
fprintf('╚════════════════════════════════════════════════════════════════════════╝\n');
fprintf('\n');
fprintf('  Total mission time:  %.1f s (%.2f min)\n', total_time, total_time/60);
fprintf('  Total fuel burned:   %.1f kg\n', total_fuel);
fprintf('  Fuel remaining:      %.1f kg (%.1f%% of capacity)\n', fuel_remaining, 100*fuel_remaining/fuel_capacity);
fprintf('  Final weight:        %.0f kg\n', W_final);
fprintf('\n');

% =========================================================================
% CANARD/STABILITY SUMMARY
% =========================================================================
fprintf('╔════════════════════════════════════════════════════════════════════════╗\n');
fprintf('║                    CANARD/STABILITY SUMMARY                            ║\n');
fprintf('╠════════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  CG Position                                                           ║\n');
fprintf('║    Initial:         %.2f m from nose (%.1f%% MAC)                      ║\n', ...
    cg_data_init.x_cg_m, cg_data_init.x_cg_MAC_percent);
fprintf('║    Final:           %.2f m from nose (%.1f%% MAC)                      ║\n', ...
    cg_data_final.x_cg_m, cg_data_final.x_cg_MAC_percent);
fprintf('║    Total travel:    %.3f m (%.1f%% MAC)                                ║\n', ...
    abs(cg_data_init.x_cg_m - cg_data_final.x_cg_m), ...
    abs(cg_data_init.x_cg_MAC_percent - cg_data_final.x_cg_MAC_percent));
fprintf('╠════════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Canard Deflection by Phase                                            ║\n');
fprintf('║    Takeoff:         %.1f° to %.1f°                                      ║\n', ...
    min(canard_takeoff), max(canard_takeoff));
fprintf('║    Climb:           %.1f° to %.1f°                                      ║\n', ...
    min(canard_climb), max(canard_climb));
fprintf('║    Accel:           %.1f° to %.1f°                                      ║\n', ...
    min(canard_accel), max(canard_accel));
fprintf('║    Cruise:          %.1f° to %.1f°                                      ║\n', ...
    min(canard_cruise), max(canard_cruise));
fprintf('║    Descent:         %.1f° to %.1f°                                      ║\n', ...
    min(canard_descent), max(canard_descent));
fprintf('║    Landing:         %.1f° to %.1f°                                      ║\n', ...
    min(canard_landing), max(canard_landing));
fprintf('╠════════════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Cm_alpha (Pitch Stability)                                            ║\n');
fprintf('║    Subsonic (M<0.9):   %.3f /rad (relaxed stable)                     ║\n', ...
    mean(Cm_alpha_climb(M_climb < 0.9)));
fprintf('║    Transonic:          %.3f /rad                                       ║\n', ...
    mean(Cm_alpha_accel));
fprintf('║    Supersonic (M=1.4): %.3f /rad                                       ║\n', ...
    mean(Cm_alpha_cruise));
fprintf('╚════════════════════════════════════════════════════════════════════════╝\n');

% =========================================================================
% PLOTTING
% =========================================================================
fprintf('\nGenerating plots...\n');

total_time = t_landing(end);

% Phase boundaries for shading
phase_times = [0, t_offset_takeoff, t_offset_climb, t_offset_accel, t_offset_cruise, ...
    t_offset_descent, t_offset_landing, total_time];
phase_names = {'Taxi', 'TO', 'Climb', 'Accel', 'Cruise', 'Desc', 'Land'};
phase_colors = [0.9 0.9 0.9; 0.8 0.9 0.8; 0.8 0.8 0.95; 0.95 0.85 0.8; ...
    0.85 0.95 0.85; 0.95 0.95 0.8; 0.9 0.8 0.8];

% ----- Figure 1: CG vs Time -----
figure('Position', [50, 50, 900, 500], 'Color', 'w');
hold on;

all_cg = [cg_taxi; cg_takeoff; cg_climb; cg_accel; cg_cruise; cg_descent; cg_landing];
cg_min = min(all_cg) - 0.05;
cg_max = max(all_cg) + 0.05;

for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [cg_min cg_min cg_max cg_max], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

plot(t_taxi, cg_taxi, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_takeoff, cg_takeoff, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_climb, cg_climb, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_accel, cg_accel, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_cruise, cg_cruise, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_descent, cg_descent, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_landing, cg_landing, 'Color', [0.6 0 0.6], 'LineWidth', 2);

plot(0, cg_taxi(1), 'go', 'MarkerSize', 10, 'MarkerFaceColor', 'g');
plot(t_landing(end), cg_landing(end), 'rs', 'MarkerSize', 10, 'MarkerFaceColor', 'r');

xlabel('Time (s)');
ylabel('X_{CG} from Nose (m)');
title(sprintf('X-29 CG Position vs Time (TOGW=%.0f kg)', TOGW));
ylim([cg_min cg_max]);
xlim([0 total_time]);
grid on; box on;

text(0.95, 0.95, sprintf('CG Travel: %.2f m\n(%.1f%% MAC)', ...
    abs(cg_taxi(1) - cg_landing(end)), ...
    abs(cg_data_init.x_cg_MAC_percent - cg_data_final.x_cg_MAC_percent)), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'BackgroundColor', 'w');
saveas(gcf, 'X29_canard_cg_vs_time.png');
fprintf('  Figure saved: X29_canard_cg_vs_time.png\n');

% ----- Figure 2: Canard Deflection vs Time -----
figure('Position', [100, 100, 900, 500], 'Color', 'w');
hold on;

for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [-35 -35 60 60], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

plot(t_taxi, canard_taxi, 'b-', 'LineWidth', 2);
plot(t_takeoff, canard_takeoff, 'b-', 'LineWidth', 2);
plot(t_climb, canard_climb, 'b-', 'LineWidth', 2);
plot(t_accel, canard_accel, 'b-', 'LineWidth', 2);
plot(t_cruise, canard_cruise, 'b-', 'LineWidth', 2);
plot(t_descent, canard_descent, 'b-', 'LineWidth', 2);
plot(t_landing, canard_landing, 'b-', 'LineWidth', 2);

yline(58, 'r--', 'Max deflection (+58°)', 'LineWidth', 1);
yline(-32, 'r--', 'Min deflection (-32°)', 'LineWidth', 1);

xlabel('Time (s)');
ylabel('Canard Deflection (deg)');
title(sprintf('X-29 Canard Deflection for Trim vs Time (TOGW=%.0f kg)', TOGW));
ylim([-35 60]);
xlim([0 total_time]);
grid on; box on;
saveas(gcf, 'X29_canard_deflection_vs_time.png');
fprintf('  Figure saved: X29_canard_deflection_vs_time.png\n');

% ----- Figure 3: Cm_alpha vs Time -----
figure('Position', [150, 150, 900, 500], 'Color', 'w');
hold on;

all_Cm = [Cm_alpha_taxi; Cm_alpha_takeoff; Cm_alpha_climb; Cm_alpha_accel; ...
    Cm_alpha_cruise; Cm_alpha_descent; Cm_alpha_landing];
Cm_min = min(all_Cm) - 0.01;
Cm_max = max(all_Cm) + 0.01;

for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [Cm_min Cm_min Cm_max Cm_max], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

plot(t_taxi, Cm_alpha_taxi, 'r-', 'LineWidth', 2);
plot(t_takeoff, Cm_alpha_takeoff, 'r-', 'LineWidth', 2);
plot(t_climb, Cm_alpha_climb, 'r-', 'LineWidth', 2);
plot(t_accel, Cm_alpha_accel, 'r-', 'LineWidth', 2);
plot(t_cruise, Cm_alpha_cruise, 'r-', 'LineWidth', 2);
plot(t_descent, Cm_alpha_descent, 'r-', 'LineWidth', 2);
plot(t_landing, Cm_alpha_landing, 'r-', 'LineWidth', 2);

yline(0, 'k--', 'Neutral stability', 'LineWidth', 1);

xlabel('Time (s)');
ylabel('C_{m_\alpha} (1/rad)');
title(sprintf('X-29 Pitch Stability Derivative vs Time (TOGW=%.0f kg)', TOGW));
ylim([Cm_min Cm_max]);
xlim([0 total_time]);
grid on; box on;

text(0.95, 0.95, 'Cm_\alpha > 0: Unstable (RSS)', ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 10, 'BackgroundColor', 'w');
saveas(gcf, 'X29_canard_Cm_alpha_vs_time.png');
fprintf('  Figure saved: X29_canard_Cm_alpha_vs_time.png\n');

% ----- Figure 4: Alpha vs Time -----
figure('Position', [200, 200, 900, 500], 'Color', 'w');
hold on;

all_alpha = [alpha_taxi; alpha_takeoff; alpha_climb; alpha_accel; ...
    alpha_cruise; alpha_descent; alpha_landing];
alpha_min = min(all_alpha) - 1;
alpha_max = max(all_alpha) + 1;

for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [alpha_min alpha_min alpha_max alpha_max], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end

plot(t_taxi, alpha_taxi, 'g-', 'LineWidth', 2);
plot(t_takeoff, alpha_takeoff, 'g-', 'LineWidth', 2);
plot(t_climb, alpha_climb, 'g-', 'LineWidth', 2);
plot(t_accel, alpha_accel, 'g-', 'LineWidth', 2);
plot(t_cruise, alpha_cruise, 'g-', 'LineWidth', 2);
plot(t_descent, alpha_descent, 'g-', 'LineWidth', 2);
plot(t_landing, alpha_landing, 'g-', 'LineWidth', 2);

xlabel('Time (s)');
ylabel('\alpha (deg)');
title(sprintf('X-29 Angle of Attack vs Time (TOGW=%.0f kg)', TOGW));
ylim([alpha_min alpha_max]);
xlim([0 total_time]);
grid on; box on;
saveas(gcf, 'X29_canard_alpha_vs_time.png');
fprintf('  Figure saved: X29_canard_alpha_vs_time.png\n');

% ----- Figure 5: Combined Mission Profile -----
figure('Position', [50, 50, 1200, 800], 'Color', 'w');

% Subplot 1: Altitude
subplot(3,2,1);
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [0 0 15000 15000], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, h_taxi, 'b-', 'LineWidth', 2);
plot(t_takeoff, h_takeoff, 'b-', 'LineWidth', 2);
plot(t_climb, h_climb, 'b-', 'LineWidth', 2);
plot(t_accel, h_accel, 'b-', 'LineWidth', 2);
plot(t_cruise, h_cruise, 'b-', 'LineWidth', 2);
plot(t_descent, h_descent, 'b-', 'LineWidth', 2);
plot(t_landing, h_landing, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Altitude (m)');
title('Altitude Profile');
ylim([0 14000]); xlim([0 total_time]);
grid on; box on;

% Subplot 2: Mach
subplot(3,2,2);
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [0 0 2 2], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, M_taxi, 'r-', 'LineWidth', 2);
plot(t_takeoff, M_takeoff, 'r-', 'LineWidth', 2);
plot(t_climb, M_climb, 'r-', 'LineWidth', 2);
plot(t_accel, M_accel, 'r-', 'LineWidth', 2);
plot(t_cruise, M_cruise_vec, 'r-', 'LineWidth', 2);
plot(t_descent, M_descent, 'r-', 'LineWidth', 2);
plot(t_landing, M_landing, 'r-', 'LineWidth', 2);
yline(1.0, 'k--', 'M=1.0', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Mach');
title('Mach Profile');
ylim([0 1.6]); xlim([0 total_time]);
grid on; box on;

% Subplot 3: Weight
subplot(3,2,3);
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [5500 5500 9000 9000], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, W_taxi, 'g-', 'LineWidth', 2);
plot(t_takeoff, W_takeoff, 'g-', 'LineWidth', 2);
plot(t_climb, W_climb, 'g-', 'LineWidth', 2);
plot(t_accel, W_accel, 'g-', 'LineWidth', 2);
plot(t_cruise, W_cruise, 'g-', 'LineWidth', 2);
plot(t_descent, W_descent, 'g-', 'LineWidth', 2);
plot(t_landing, W_landing, 'g-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Weight (kg)');
title(sprintf('Weight Profile (Fuel burned: %.0f kg)', total_fuel));
ylim([6000 8500]); xlim([0 total_time]);
grid on; box on;

% Subplot 4: CG
subplot(3,2,4);
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [cg_min cg_min cg_max cg_max], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, cg_taxi, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_takeoff, cg_takeoff, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_climb, cg_climb, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_accel, cg_accel, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_cruise, cg_cruise, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_descent, cg_descent, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_landing, cg_landing, 'Color', [0.6 0 0.6], 'LineWidth', 2);
xlabel('Time (s)'); ylabel('X_{CG} (m)');
title('CG Position');
ylim([cg_min cg_max]); xlim([0 total_time]);
grid on; box on;

% Subplot 5: Canard Deflection
subplot(3,2,5);
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [-35 -35 60 60], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, canard_taxi, 'b-', 'LineWidth', 2);
plot(t_takeoff, canard_takeoff, 'b-', 'LineWidth', 2);
plot(t_climb, canard_climb, 'b-', 'LineWidth', 2);
plot(t_accel, canard_accel, 'b-', 'LineWidth', 2);
plot(t_cruise, canard_cruise, 'b-', 'LineWidth', 2);
plot(t_descent, canard_descent, 'b-', 'LineWidth', 2);
plot(t_landing, canard_landing, 'b-', 'LineWidth', 2);
yline(58, 'r--', 'LineWidth', 1);
yline(-32, 'r--', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('\delta_c (deg)');
title('Canard Deflection for Trim');
ylim([-35 60]); xlim([0 total_time]);
grid on; box on;

% Subplot 6: Alpha
subplot(3,2,6);
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [alpha_min alpha_min alpha_max alpha_max], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, alpha_taxi, 'g-', 'LineWidth', 2);
plot(t_takeoff, alpha_takeoff, 'g-', 'LineWidth', 2);
plot(t_climb, alpha_climb, 'g-', 'LineWidth', 2);
plot(t_accel, alpha_accel, 'g-', 'LineWidth', 2);
plot(t_cruise, alpha_cruise, 'g-', 'LineWidth', 2);
plot(t_descent, alpha_descent, 'g-', 'LineWidth', 2);
plot(t_landing, alpha_landing, 'g-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('\alpha (deg)');
title('Angle of Attack');
ylim([alpha_min alpha_max]); xlim([0 total_time]);
grid on; box on;

sgtitle(sprintf('X-29 Mission Profile with Canard/Stability Tracking (TOGW=%.0f kg)', TOGW));
saveas(gcf, 'X29_canard_mission_profile.png');
fprintf('  Figure saved: X29_canard_mission_profile.png\n');

% =========================================================================
% STANDARD MISSION PLOTS (matching run_full_mission.m)
% =========================================================================
fprintf('\nGenerating standard mission plots...\n');

% Calculate cumulative fuel burned for each phase (matching time array lengths)
fuel_taxi = linspace(0, taxi.fuel_burned, length(t_taxi))';
fuel_takeoff = taxi.fuel_burned + (takeoff.W(1) - takeoff.W);
fuel_climb = taxi.fuel_burned + takeoff.fuel_burned + (climb.W(1) - climb.W);
fuel_accel = taxi.fuel_burned + takeoff.fuel_burned + climb.total_fuel + (accel.W(1) - accel.W);
fuel_cruise_vec = taxi.fuel_burned + takeoff.fuel_burned + climb.total_fuel + accel.total_fuel + cruise.fuel_burned;
fuel_descent = taxi.fuel_burned + takeoff.fuel_burned + climb.total_fuel + accel.total_fuel + ...
    cruise.fuel_burned(end) + (descent.W(1) - descent.W);
fuel_landing = taxi.fuel_burned + takeoff.fuel_burned + climb.total_fuel + accel.total_fuel + ...
    cruise.fuel_burned(end) + descent.total_fuel + (landing.W(1) - landing.W);

% ----- Mission Plot 1: Altitude vs Time -----
figure('Position', [50, 50, 800, 500], 'Color', 'w');
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [0 0 15000 15000], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, h_taxi, 'b-', 'LineWidth', 2);
plot(t_takeoff, h_takeoff, 'b-', 'LineWidth', 2);
plot(t_climb, h_climb, 'b-', 'LineWidth', 2);
plot(t_accel, h_accel, 'b-', 'LineWidth', 2);
plot(t_cruise, h_cruise, 'b-', 'LineWidth', 2);
plot(t_descent, h_descent, 'b-', 'LineWidth', 2);
plot(t_landing, h_landing, 'b-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Altitude (m)');
title(sprintf('X-29 Altitude Profile (TOGW=%.0f kg) - Physics Propulsion', TOGW));
ylim([0 14000]); xlim([0 total_time]);
grid on; box on;
saveas(gcf, 'X29_mission_altitude.png');
fprintf('  Figure saved: X29_mission_altitude.png\n');

% ----- Mission Plot 2: Mach vs Time -----
figure('Position', [100, 100, 800, 500], 'Color', 'w');
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [0 0 2 2], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, M_taxi, 'r-', 'LineWidth', 2);
plot(t_takeoff, M_takeoff, 'r-', 'LineWidth', 2);
plot(t_climb, M_climb, 'r-', 'LineWidth', 2);
plot(t_accel, M_accel, 'r-', 'LineWidth', 2);
plot(t_cruise, M_cruise_vec, 'r-', 'LineWidth', 2);
plot(t_descent, M_descent, 'r-', 'LineWidth', 2);
plot(t_landing, M_landing, 'r-', 'LineWidth', 2);
yline(1.0, 'k--', 'M = 1.0', 'LineWidth', 1);
yline(1.4, 'k:', 'M = 1.4', 'LineWidth', 1);
xlabel('Time (s)'); ylabel('Mach Number');
title(sprintf('X-29 Mach Profile (TOGW=%.0f kg) - Physics Propulsion', TOGW));
ylim([0 1.6]); xlim([0 total_time]);
grid on; box on;
saveas(gcf, 'X29_mission_mach.png');
fprintf('  Figure saved: X29_mission_mach.png\n');

% ----- Mission Plot 3: Weight vs Time -----
figure('Position', [150, 150, 800, 500], 'Color', 'w');
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [5500 5500 9000 9000], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, W_taxi, 'g-', 'LineWidth', 2);
plot(t_takeoff, W_takeoff, 'g-', 'LineWidth', 2);
plot(t_climb, W_climb, 'g-', 'LineWidth', 2);
plot(t_accel, W_accel, 'g-', 'LineWidth', 2);
plot(t_cruise, W_cruise, 'g-', 'LineWidth', 2);
plot(t_descent, W_descent, 'g-', 'LineWidth', 2);
plot(t_landing, W_landing, 'g-', 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Weight (kg)');
title(sprintf('X-29 Weight History (TOGW=%.0f kg) - Physics Propulsion', TOGW));
ylim([6000 8500]); xlim([0 total_time]);
grid on; box on;
saveas(gcf, 'X29_mission_weight.png');
fprintf('  Figure saved: X29_mission_weight.png\n');

% ----- Mission Plot 4: Fuel Burned vs Time -----
figure('Position', [200, 200, 800, 500], 'Color', 'w');
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [0 0 2000 2000], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, fuel_taxi, 'm-', 'LineWidth', 2);
plot(t_takeoff, fuel_takeoff, 'm-', 'LineWidth', 2);
plot(t_climb, fuel_climb, 'm-', 'LineWidth', 2);
plot(t_accel, fuel_accel, 'm-', 'LineWidth', 2);
plot(t_cruise, fuel_cruise_vec, 'm-', 'LineWidth', 2);
plot(t_descent, fuel_descent, 'm-', 'LineWidth', 2);
plot(t_landing, fuel_landing, 'm-', 'LineWidth', 2);
yline(fuel_capacity, 'r--', 'Fuel Capacity', 'LineWidth', 1.5);
xlabel('Time (s)'); ylabel('Fuel Burned (kg)');
title(sprintf('X-29 Cumulative Fuel Consumption (TOGW=%.0f kg) - Physics Propulsion', TOGW));
ylim([0 2000]); xlim([0 total_time]);
grid on; box on;
saveas(gcf, 'X29_mission_fuel.png');
fprintf('  Figure saved: X29_mission_fuel.png\n');

% ----- Mission Plot 5: Flight Envelope (Mach vs Altitude) -----
figure('Position', [250, 250, 800, 500], 'Color', 'w');
hold on;
plot(M_taxi, h_taxi/1000, 'k-', 'LineWidth', 2, 'DisplayName', 'Taxi');
plot(M_takeoff, h_takeoff/1000, 'g-', 'LineWidth', 2, 'DisplayName', 'Takeoff');
plot(M_climb, h_climb/1000, 'b-', 'LineWidth', 2, 'DisplayName', 'Climb');
plot(M_accel, h_accel/1000, 'r-', 'LineWidth', 2, 'DisplayName', 'Accel');
plot(M_cruise_vec, h_cruise/1000, 'm-', 'LineWidth', 2, 'DisplayName', 'Cruise');
plot(M_descent, h_descent/1000, 'c-', 'LineWidth', 2, 'DisplayName', 'Descent');
plot(M_landing, h_landing/1000, 'y-', 'LineWidth', 2, 'DisplayName', 'Landing');
plot(0, 0, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'g', 'DisplayName', 'Start');
plot(0, 0, 'ks', 'MarkerSize', 10, 'MarkerFaceColor', 'r', 'DisplayName', 'End');
xlabel('Mach Number'); ylabel('Altitude (km)');
title(sprintf('X-29 Flight Envelope (TOGW=%.0f kg) - Physics Propulsion', TOGW));
legend('Location', 'northwest');
xlim([0 1.6]); ylim([0 14]);
grid on; box on;
saveas(gcf, 'X29_mission_envelope.png');
fprintf('  Figure saved: X29_mission_envelope.png\n');

% ----- Mission Plot 6: CG vs Time -----
figure('Position', [300, 300, 800, 500], 'Color', 'w');
hold on;
for i = 1:7
    patch([phase_times(i) phase_times(i+1) phase_times(i+1) phase_times(i)], ...
        [cg_min cg_min cg_max cg_max], phase_colors(i,:), 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(t_taxi, cg_taxi, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_takeoff, cg_takeoff, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_climb, cg_climb, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_accel, cg_accel, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_cruise, cg_cruise, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_descent, cg_descent, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(t_landing, cg_landing, 'Color', [0.6 0 0.6], 'LineWidth', 2);
plot(0, cg_taxi(1), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');
plot(t_landing(end), cg_landing(end), 'rs', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('Time (s)'); ylabel('X_{CG} from Nose Cone (m)');
title(sprintf('X-29 Center of Gravity vs Time (TOGW=%.0f kg) - Physics Propulsion', TOGW));
ylim([cg_min cg_max]); xlim([0 total_time]);
grid on; box on;
text(0.95, 0.95, sprintf('CG Travel: %.2f m', abs(cg_taxi(1) - cg_landing(end))), ...
    'Units', 'normalized', 'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', ...
    'FontSize', 9, 'BackgroundColor', 'w');
saveas(gcf, 'X29_mission_cg.png');
fprintf('  Figure saved: X29_mission_cg.png\n');

% =========================================================================
% CRUISE DETAILED ANALYSIS PLOTS
% =========================================================================
fprintf('\nGenerating cruise detailed analysis plots...\n');

% Extract cruise data
t_cruise_plot = cruise.t;
t_cruise_min = t_cruise_plot / 60;

% Get atmospheric conditions at cruise
[~, a_cr, ~, rho_cr] = atmosisa(12200);
V_cr = 1.4 * a_cr;
q_cr = 0.5 * rho_cr * V_cr^2;

% Drag breakdown
CD0_cruise = cruise.aero.CD0;
K_cruise = cruise.aero.K;
CDi_cruise = K_cruise * cruise.CL.^2;

% ----- Cruise Plot 1: L/D vs Time -----
figure('Position', [50, 50, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.LD, 'b-', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('L/D');
title(sprintf('X-29 Cruise L/D vs Time (M=1.4, h=12.2 km) - Physics Propulsion'));
grid on;
ylim([0 max(cruise.LD)*1.2]);
text(0.05, 0.95, sprintf('L/D_{avg} = %.2f', mean(cruise.LD)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_LD.png');
fprintf('  Figure saved: X29_cruise_LD.png\n');

% ----- Cruise Plot 2: Alpha vs Time -----
figure('Position', [100, 100, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.alpha, 'r-', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('\alpha (deg)');
title('X-29 Cruise Angle of Attack vs Time (M=1.4) - Physics Propulsion');
grid on;
text(0.05, 0.95, sprintf('\\alpha_{range} = %.2f° to %.2f°', min(cruise.alpha), max(cruise.alpha)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_alpha.png');
fprintf('  Figure saved: X29_cruise_alpha.png\n');

% ----- Cruise Plot 3: Thrust Required vs Time -----
figure('Position', [150, 150, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.D/1000, 'g-', 'LineWidth', 2, 'DisplayName', 'Drag = T_{req}');
hold on;
plot(t_cruise_min, cruise.T/1000, 'k--', 'LineWidth', 1.5, 'DisplayName', 'Thrust');
xlabel('Time (min)'); ylabel('Force (kN)');
title('X-29 Cruise Thrust Required vs Time (M=1.4) - Physics Propulsion');
legend('Location', 'best');
grid on;
text(0.05, 0.95, sprintf('T_{req,avg} = %.1f kN', mean(cruise.D)/1000), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_thrust.png');
fprintf('  Figure saved: X29_cruise_thrust.png\n');

% ----- Cruise Plot 4: Weight vs Time (Cruise) -----
figure('Position', [200, 200, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.W, 'm-', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('Weight (kg)');
title('X-29 Cruise Weight vs Time (M=1.4) - Physics Propulsion');
grid on;
text(0.05, 0.95, sprintf('\\Delta W = %.1f kg', cruise.W(1) - cruise.W(end)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_weight.png');
fprintf('  Figure saved: X29_cruise_weight.png\n');

% ----- Cruise Plot 5: CL and CD vs Time -----
figure('Position', [250, 250, 800, 500], 'Color', 'w');
yyaxis left;
plot(t_cruise_min, cruise.CL, 'b-', 'LineWidth', 2);
ylabel('C_L');
yyaxis right;
plot(t_cruise_min, cruise.CD, 'r-', 'LineWidth', 2);
ylabel('C_D');
xlabel('Time (min)');
title('X-29 Cruise CL and CD vs Time (M=1.4) - Physics Propulsion');
grid on;
legend('C_L', 'C_D', 'Location', 'best');
saveas(gcf, 'X29_cruise_CL_CD.png');
fprintf('  Figure saved: X29_cruise_CL_CD.png\n');

% ----- Cruise Plot 6: Fuel Flow vs Time -----
figure('Position', [300, 300, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.fuel_flow, 'Color', [0.8 0.4 0], 'LineWidth', 2);
xlabel('Time (min)'); ylabel('Fuel Flow (kg/s)');
title('X-29 Cruise Fuel Flow Rate vs Time (M=1.4) - Physics Propulsion');
grid on;
text(0.05, 0.95, sprintf('FF_{avg} = %.3f kg/s', mean(cruise.fuel_flow)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_fuel_flow.png');
fprintf('  Figure saved: X29_cruise_fuel_flow.png\n');

% ----- Cruise Plot 7: PLA vs Time -----
figure('Position', [350, 350, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.PLA, 'Color', [0.5 0 0.5], 'LineWidth', 2);
xlabel('Time (min)'); ylabel('PLA (deg)');
title('X-29 Cruise Power Lever Angle vs Time (M=1.4) - Physics Propulsion');
grid on;
ylim([80 140]);
yline(87, 'k--', 'MIL Power', 'LineWidth', 1);
yline(130, 'r--', 'Max AB', 'LineWidth', 1);
text(0.05, 0.85, sprintf('PLA_{avg} = %.1f°', mean(cruise.PLA)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_PLA.png');
fprintf('  Figure saved: X29_cruise_PLA.png\n');

% ----- Cruise Plot 8: Cumulative Fuel Burned -----
figure('Position', [400, 400, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise.fuel_burned, 'k-', 'LineWidth', 2);
xlabel('Time (min)'); ylabel('Fuel Burned (kg)');
title('X-29 Cruise Cumulative Fuel Burned vs Time (M=1.4) - Physics Propulsion');
grid on;
text(0.05, 0.95, sprintf('Total = %.1f kg', cruise.fuel_burned(end)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_fuel_burned.png');
fprintf('  Figure saved: X29_cruise_fuel_burned.png\n');

% ----- Cruise Plot 9: Drag Breakdown -----
figure('Position', [450, 450, 800, 500], 'Color', 'w');
area(t_cruise_min, [CD0_cruise*ones(size(t_cruise_min)), CDi_cruise], 'LineWidth', 1);
xlabel('Time (min)'); ylabel('C_D');
title('X-29 Cruise Drag Breakdown (M=1.4) - Physics Propulsion');
legend('C_{D0} (Parasitic)', 'C_{Di} (Induced)', 'Location', 'best');
grid on;
text(0.5, 0.95, sprintf('C_{D0}/C_D = %.1f%%', 100*CD0_cruise/mean(cruise.CD)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center');
saveas(gcf, 'X29_cruise_drag_breakdown.png');
fprintf('  Figure saved: X29_cruise_drag_breakdown.png\n');

% ----- Cruise Plot 10: Drag Polar -----
figure('Position', [50, 50, 800, 500], 'Color', 'w');
CL_range = linspace(0, 0.5, 100);
CD_polar = CD0_cruise + K_cruise * CL_range.^2;
plot(CD_polar, CL_range, 'b-', 'LineWidth', 2);
hold on;
plot(cruise.CD, cruise.CL, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('C_D'); ylabel('C_L');
title('X-29 Drag Polar (Supersonic M=1.4) - Physics Propulsion');
grid on;
legend('Polar', 'Operating Points', 'Location', 'best');
saveas(gcf, 'X29_cruise_drag_polar.png');
fprintf('  Figure saved: X29_cruise_drag_polar.png\n');

% ----- Cruise Plot 11: L/D vs CL -----
figure('Position', [100, 100, 800, 500], 'Color', 'w');
LD_polar = CL_range ./ CD_polar;
plot(CL_range, LD_polar, 'b-', 'LineWidth', 2);
hold on;
plot(cruise.CL, cruise.LD, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
xlabel('C_L'); ylabel('L/D');
title('X-29 L/D vs C_L (Supersonic M=1.4) - Physics Propulsion');
grid on;
[LD_max, idx_max] = max(LD_polar);
CL_opt = CL_range(idx_max);
plot(CL_opt, LD_max, 'g^', 'MarkerSize', 12, 'MarkerFaceColor', 'g');
text(CL_opt+0.02, LD_max, sprintf('(L/D)_{max}=%.2f\nC_L=%.3f', LD_max, CL_opt), 'FontSize', 9);
legend('Theoretical', 'Operating', 'Optimum', 'Location', 'best');
saveas(gcf, 'X29_cruise_LD_vs_CL.png');
fprintf('  Figure saved: X29_cruise_LD_vs_CL.png\n');

% ----- Cruise Plot 12: Specific Range -----
figure('Position', [150, 150, 800, 500], 'Color', 'w');
SR_cruise = V_cr ./ (cruise.fuel_flow * 1000);
plot(t_cruise_min, SR_cruise, 'Color', [0 0.6 0], 'LineWidth', 2);
xlabel('Time (min)'); ylabel('Specific Range (km/kg)');
title('X-29 Cruise Specific Range vs Time (M=1.4) - Physics Propulsion');
grid on;
text(0.05, 0.95, sprintf('SR_{avg} = %.3f km/kg', mean(SR_cruise)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_specific_range.png');
fprintf('  Figure saved: X29_cruise_specific_range.png\n');

% ----- Cruise Plot 13: TSFC -----
% Calculate TSFC from fuel flow and thrust: TSFC = fuel_flow [lb/hr] / thrust [lbf]
cruise_fuel_flow_lb_hr = cruise.fuel_flow * 7936.64;  % kg/s to lb/hr
cruise_thrust_lbf = cruise.T / 4.44822;  % N to lbf
cruise_TSFC = cruise_fuel_flow_lb_hr ./ cruise_thrust_lbf;
figure('Position', [200, 200, 800, 500], 'Color', 'w');
plot(t_cruise_min, cruise_TSFC, 'Color', [0.8 0.2 0.2], 'LineWidth', 2);
xlabel('Time (min)'); ylabel('TSFC (lb/hr/lbf)');
title('X-29 Cruise TSFC vs Time (M=1.4) - Calibrated Propulsion');
grid on;
text(0.05, 0.95, sprintf('TSFC_{avg} = %.3f', mean(cruise_TSFC)), ...
    'Units', 'normalized', 'FontSize', 10, 'VerticalAlignment', 'top');
saveas(gcf, 'X29_cruise_TSFC.png');
fprintf('  Figure saved: X29_cruise_TSFC.png\n');

% =========================================================================
% CRUISE SUMMARY TABLE
% =========================================================================
fprintf('\n');
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║               SUPERSONIC CRUISE SUMMARY                        ║\n');
fprintf('║              (Physics-Based Propulsion Model)                  ║\n');
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Flight Conditions                                             ║\n');
fprintf('║    Mach number:           %.2f                                 ║\n', 1.4);
fprintf('║    Altitude:              %.0f m (%.0f ft)                    ║\n', 12200, 12200*3.281);
fprintf('║    True airspeed:         %.1f m/s (%.0f kts)                 ║\n', V_cr, V_cr*1.944);
fprintf('║    Dynamic pressure:      %.1f kPa                             ║\n', q_cr/1000);
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Aerodynamics                                                  ║\n');
fprintf('║    C_D0:                  %.4f                                ║\n', CD0_cruise);
fprintf('║    K:                     %.4f                                ║\n', K_cruise);
fprintf('║    C_L range:             %.4f to %.4f                       ║\n', min(cruise.CL), max(cruise.CL));
fprintf('║    C_D range:             %.4f to %.4f                       ║\n', min(cruise.CD), max(cruise.CD));
fprintf('║    L/D range:             %.2f to %.2f                         ║\n', min(cruise.LD), max(cruise.LD));
fprintf('║    Alpha range:           %.2f° to %.2f°                       ║\n', min(cruise.alpha), max(cruise.alpha));
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Propulsion (0-D Cycle Analysis)                               ║\n');
fprintf('║    Thrust required:       %.1f to %.1f kN                     ║\n', min(cruise.D)/1000, max(cruise.D)/1000);
fprintf('║    PLA range:             %.1f° to %.1f°                       ║\n', min(cruise.PLA), max(cruise.PLA));
fprintf('║    Fuel flow:             %.3f to %.3f kg/s                   ║\n', min(cruise.fuel_flow), max(cruise.fuel_flow));
fprintf('║    TSFC:                  %.3f lb/hr/lbf                       ║\n', mean(cruise_TSFC));
fprintf('╠════════════════════════════════════════════════════════════════╣\n');
fprintf('║  Mission Performance                                           ║\n');
fprintf('║    Duration:              %.0f s (%.1f min)                    ║\n', cruise.t(end), cruise.t(end)/60);
fprintf('║    Initial weight:        %.0f kg                              ║\n', cruise.W(1));
fprintf('║    Final weight:          %.0f kg                              ║\n', cruise.W(end));
fprintf('║    Fuel burned:           %.1f kg                              ║\n', cruise.fuel_burned(end));
fprintf('║    Specific range:        %.3f km/kg                           ║\n', mean(SR_cruise));
fprintf('║    Distance covered:      %.1f km (%.1f NM)                   ║\n', V_cr*cruise.t(end)/1000, V_cr*cruise.t(end)/1852);
fprintf('╚════════════════════════════════════════════════════════════════╝\n');

% =========================================================================
% ALTITUDE VS DISTANCE PLOT
% =========================================================================
fprintf('\nGenerating Altitude vs Distance plot...\n');

% Use the already-extracted altitude and time arrays:
% h_takeoff, h_climb, h_accel, h_cruise, h_descent, h_landing
% t_takeoff, t_climb, t_accel, t_cruise, t_descent, t_landing

% Calculate horizontal distance for each phase

% Takeoff (short, just ground roll + initial climb)
n_to = length(t_takeoff);
x_to_phase = linspace(0, 500, n_to)';  % ~500m takeoff roll

% Climb - account for climb angle
n_cl = length(t_climb);
dt_cl = [0; diff(t_climb)];
dh_cl = [0; diff(h_climb)];
V_cl = climb.V;
gamma_cl = zeros(n_cl, 1);
for i = 2:n_cl
    denom = dt_cl(i) * V_cl(i);
    if denom > 0
        sin_gamma = dh_cl(i) / denom;
        sin_gamma = max(-1, min(1, sin_gamma));  % Clamp to valid range
        gamma_cl(i) = asind(sin_gamma);
    end
end
gamma_cl(isnan(gamma_cl) | abs(gamma_cl) > 45) = 10;
x_cl_phase = cumsum(dt_cl .* V_cl .* cosd(gamma_cl));

% Accel (constant altitude at 12.2 km)
n_ac = length(t_accel);
dt_ac = [0; diff(t_accel)];
V_ac = accel.V;
x_ac_phase = cumsum(dt_ac .* V_ac);

% Cruise (constant velocity at M=1.4, h=12.2km)
n_cr = length(t_cruise);
[~, a_cr_temp, ~, ~] = atmosisa(12200);
V_cruise_val = 1.4 * a_cr_temp;
x_cr_phase = (t_cruise - t_cruise(1)) * V_cruise_val;

% Descent - account for descent angle
n_de = length(t_descent);
dt_de = [0; diff(t_descent)];
dh_de = [0; diff(h_descent)];
V_de = descent.V;
gamma_de = zeros(n_de, 1);
for i = 2:n_de
    denom = dt_de(i) * V_de(i);
    if denom > 0
        sin_gamma = dh_de(i) / denom;
        sin_gamma = max(-1, min(1, sin_gamma));
        gamma_de(i) = asind(sin_gamma);
    end
end
gamma_de(isnan(gamma_de) | abs(gamma_de) > 30) = -3;
x_de_phase = cumsum(dt_de .* V_de .* cosd(gamma_de));

% Landing
n_la = length(t_landing);
dt_la = [0; diff(t_landing)];
V_la = landing.V;
x_la_phase = cumsum(dt_la .* V_la);

% Cumulative distances (offsets for each phase)
x0_to = 0;
x0_cl = x_to_phase(end);
x0_ac = x0_cl + x_cl_phase(end);
x0_cr = x0_ac + x_ac_phase(end);
x0_de = x0_cr + x_cr_phase(end);
x0_la = x0_de + x_de_phase(end);

x_total_to = x0_to + x_to_phase;
x_total_cl = x0_cl + x_cl_phase;
x_total_ac = x0_ac + x_ac_phase;
x_total_cr = x0_cr + x_cr_phase;
x_total_de = x0_de + x_de_phase;
x_total_la = x0_la + x_la_phase;

total_distance_km = (x0_la + x_la_phase(end)) / 1000;

% Plot Altitude vs Distance
figure('Position', [100, 100, 1000, 500], 'Color', 'w');
hold on;
plot(x_total_to/1000, h_takeoff/1000, 'g-', 'LineWidth', 2.5, 'DisplayName', 'Takeoff');
plot(x_total_cl/1000, h_climb/1000, 'b-', 'LineWidth', 2.5, 'DisplayName', 'Climb');
plot(x_total_ac/1000, h_accel/1000, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Accel');
plot(x_total_cr/1000, h_cruise/1000, 'm-', 'LineWidth', 2.5, 'DisplayName', 'Cruise');
plot(x_total_de/1000, h_descent/1000, 'c-', 'LineWidth', 2.5, 'DisplayName', 'Descent');
plot(x_total_la/1000, h_landing/1000, 'y-', 'LineWidth', 2.5, 'DisplayName', 'Landing');

% Mark waypoints
scatter([x0_cl, x0_ac, x0_cr, x0_de, x0_la]/1000, ...
    [h_takeoff(end), h_climb(end), h_accel(end), h_cruise(end), h_descent(end)]/1000, ...
    80, 'ko', 'filled', 'HandleVisibility', 'off');

xlabel('Horizontal Distance (km)', 'FontSize', 12);
ylabel('Altitude (km)', 'FontSize', 12);
title(sprintf('X-29 Mission Profile - Total Distance: %.0f km (%.0f nm)', ...
    total_distance_km, total_distance_km/1.852), 'FontSize', 14, 'FontWeight', 'bold');
legend('Location', 'northeast', 'FontSize', 10);
grid on; box on;
xlim([0, total_distance_km * 1.02]);
ylim([0, 14]);

saveas(gcf, 'X29_mission_altitude_vs_distance.png');
fprintf('  Figure saved: X29_mission_altitude_vs_distance.png\n');
fprintf('  Total mission distance: %.1f km (%.1f nm)\n', total_distance_km, total_distance_km/1.852);

fprintf('\n════════════════════════════════════════════════════════════════\n');
fprintf('MISSION SIMULATION WITH CANARD TRACKING COMPLETE\n');
fprintf('════════════════════════════════════════════════════════════════\n');
