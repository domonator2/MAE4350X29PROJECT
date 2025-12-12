%% AERO_MODEL_FOR_DAVID - Simple X-29 Aerodynamics Calculator
% This script demonstrates how to calculate CL, CD, and L/D for the X-29
% using the Raymer aerodynamics model with optional flaperon/strake effects.
%
% Author: Dominic Larin (Chief Engineer)
% For: Andrea (Aero Lead)
% Date: 12/08/2025

clear; clc; close all;

% Add project paths
addpath(genpath(fileparts(pwd)));

%% =========================================================================
% SECTION 1: BASIC INPUTS
% =========================================================================
fprintf('=== X-29 AERODYNAMICS CALCULATOR ===\n\n');

% Flight conditions
M = 1.4;              % Mach number [-]
alpha_deg = 4;        % Angle of attack [deg]
h = 12200;            % Altitude [m] (40,000 ft)

fprintf('Flight Conditions:\n');
fprintf('  Mach:     %.2f\n', M);
fprintf('  Alpha:    %.1f deg\n', alpha_deg);
fprintf('  Altitude: %.0f m (%.0f ft)\n\n', h, h*3.281);

%% =========================================================================
% SECTION 2: BUILD X-29 GEOMETRY
% =========================================================================
% Load trade study data to get converged vehicle parameters
script_path = fileparts(mfilename('fullpath'));
proj_root = fileparts(script_path);

mat_file = fullfile(proj_root, 'Synthesis', 'updated_fuselage_calculations', ...
    'elliptical_trade_study_20251117_170903.mat');

if ~exist(mat_file, 'file')
    mat_files = dir(fullfile(proj_root, 'Synthesis', '**', '*trade_study*.mat'));
    if ~isempty(mat_files)
        mat_file = fullfile(mat_files(end).folder, mat_files(end).name);
    else
        error('Could not find trade study .mat file');
    end
end

trade_data = load(mat_file);

% Find X-29 baseline config: sweep=-30 deg, tau=0.077
sc = trade_data.successful_configs;
target_sweep = -30; target_tau = 0.077; target_duration = 10;
best_idx = 1; best_dist = inf;
for i = 1:numel(sc)
    dist = abs(sc(i).sweep_deg - target_sweep) + 100*abs(sc(i).tau - target_tau);
    if dist < best_dist
        best_dist = dist;
        best_idx = i;
    end
end
veh = sc(best_idx);

% Build geometry structure
geom = build_X29_geometry(veh, trade_data.vehicle_base);

fprintf('Vehicle Geometry Loaded:\n');
fprintf('  Wing Area:  %.2f m^2\n', geom.wing.S);
fprintf('  Wing Span:  %.2f m\n', geom.wing.b);
fprintf('  Sweep:      %.0f deg\n', veh.sweep_deg);
fprintf('  Slenderness (tau): %.3f\n\n', veh.tau);

%% =========================================================================
% SECTION 3: CALCULATE CLEAN AERODYNAMICS (NO FLAPS)
% =========================================================================
fprintf('--- CLEAN CONFIGURATION (No Flaps) ---\n');

[CL, CD, CD0, K, L_over_D, CL_alpha] = get_X29_aero(M, alpha_deg, h, geom);

fprintf('  CL:       %.4f\n', CL);
fprintf('  CD:       %.4f\n', CD);
fprintf('  CD0:      %.4f (zero-lift drag)\n', CD0);
fprintf('  CDi:      %.4f (induced drag)\n', CD - CD0);
fprintf('  K:        %.4f (induced drag factor)\n', K);
fprintf('  L/D:      %.2f\n', L_over_D);
fprintf('  CL_alpha: %.4f /deg\n\n', CL_alpha);

%% =========================================================================
% SECTION 4: CALCULATE WITH FLAPERON/STRAKE EFFECTS
% =========================================================================
fprintf('--- WITH FLAPERON + STRAKE (Takeoff Config) ---\n');

% Typical takeoff deflections
delta_flaperon = 15;  % deg (trailing edge down)
delta_strake = 6;     % deg

% Get clean aero as a structure for apply_flap_effects
aero_clean = raymer_aero_model_elliptical(M, alpha_deg, h, geom);

% Apply flap effects
aero_flapped = apply_flap_effects(aero_clean, delta_flaperon, geom, delta_strake);

fprintf('  Flaperon: %.1f deg, Strake: %.1f deg\n', delta_flaperon, delta_strake);
fprintf('  CL:       %.4f (was %.4f, delta = +%.4f)\n', aero_flapped.CL, CL, aero_flapped.CL - CL);
fprintf('  CD:       %.4f (was %.4f, delta = +%.4f)\n', aero_flapped.CD, CD, aero_flapped.CD - CD);
fprintf('  L/D:      %.2f (was %.2f)\n\n', aero_flapped.CL/aero_flapped.CD, L_over_D);

%% =========================================================================
% SECTION 5: LANDING CONFIGURATION
% =========================================================================
fprintf('--- WITH FLAPERON + STRAKE (Landing Config) ---\n');

% Typical landing deflections (higher)
delta_flaperon_land = 20;  % deg
delta_strake_land = 10;    % deg

aero_landing = apply_flap_effects(aero_clean, delta_flaperon_land, geom, delta_strake_land);

fprintf('  Flaperon: %.1f deg, Strake: %.1f deg\n', delta_flaperon_land, delta_strake_land);
fprintf('  CL:       %.4f (delta = +%.4f from clean)\n', aero_landing.CL, aero_landing.CL - CL);
fprintf('  CD:       %.4f (delta = +%.4f from clean)\n', aero_landing.CD, aero_landing.CD - CD);
fprintf('  L/D:      %.2f\n\n', aero_landing.CL/aero_landing.CD);

%% =========================================================================
% SECTION 6: MACH SWEEP (Example)
% =========================================================================
fprintf('--- MACH SWEEP (alpha = 4 deg, h = 12.2 km) ---\n');

Mach_range = [0.4, 0.6, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6];
fprintf('  Mach    CL       CD       L/D\n');
fprintf('  ----   ------   ------   -----\n');

for i = 1:length(Mach_range)
    [CL_i, CD_i, ~, ~, LD_i, ~] = get_X29_aero(Mach_range(i), alpha_deg, h, geom);
    fprintf('  %.1f    %.4f   %.4f   %.2f\n', Mach_range(i), CL_i, CD_i, LD_i);
end

%% =========================================================================
% SECTION 7: ALPHA SWEEP (Example)
% =========================================================================
fprintf('\n--- ALPHA SWEEP (M = 1.4, h = 12.2 km) ---\n');

alpha_range = 0:2:12;
fprintf('  Alpha   CL       CD       L/D\n');
fprintf('  -----  ------   ------   -----\n');

for i = 1:length(alpha_range)
    [CL_i, CD_i, ~, ~, LD_i, ~] = get_X29_aero(M, alpha_range(i), h, geom);
    fprintf('  %4.0f    %.4f   %.4f   %.2f\n', alpha_range(i), CL_i, CD_i, LD_i);
end

fprintf('\n=== DONE ===\n');
