function extracted_data = extract_converged_vehicle_data(duration_min, tau, sweep_deg, varargin)
%EXTRACT_CONVERGED_VEHICLE_DATA Extract all parameters from a converged vehicle design
%
% Usage:
%   data = extract_converged_vehicle_data(duration_min, tau, sweep_deg)
%   data = extract_converged_vehicle_data(duration_min, tau, sweep_deg, 'save', true)
%   data = extract_converged_vehicle_data(duration_min, tau, sweep_deg, 'save', true, 'filename', 'my_vehicle.mat')
%
% Inputs:
%   duration_min - Mission duration in minutes (e.g., 10)
%   tau          - Slenderness parameter (0.01 - 0.10, e.g., 0.0724)
%   sweep_deg    - Leading edge sweep angle in degrees (e.g., -30)
%
% Optional Parameters:
%   'save'     - Save to file (default: false)
%   'filename' - Custom filename (default: auto-generated with timestamp)
%   'verbose'  - Print summary to console (default: true)
%
% Outputs:
%   extracted_data - Structure containing all vehicle parameters
%
% Example:
%   % Extract X-29 baseline configuration
%   data = extract_converged_vehicle_data(10, 0.0724, -30);
%
%   % Extract and save to file
%   data = extract_converged_vehicle_data(10, 0.07, -30, 'save', true);
%
%   % Extract with custom filename
%   data = extract_converged_vehicle_data(5, 0.06, -30, 'save', true, 'filename', 'low_tau_vehicle.mat');

%% Parse Optional Arguments
p = inputParser;
addParameter(p, 'save', false, @islogical);
addParameter(p, 'filename', '', @(x) ischar(x) || isstring(x));
addParameter(p, 'verbose', true, @islogical);
parse(p, varargin{:});

save_file = p.Results.save;
custom_filename = p.Results.filename;
verbose = p.Results.verbose;

%% Add paths
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(fileparts(script_dir));
addpath(genpath(project_root));

%% Setup Vehicle Configuration (Calibrated Parameters)
vehicle = struct();
vehicle.S_ref_m2 = 17.544;
vehicle.E_TW = 7.323;
vehicle.mu_a = 0.0;
vehicle.f_sys = 0.2;
vehicle.C_sys_kg = 1250;
vehicle.I_str_kg_m2 = 18.7;
vehicle.rho_pay_kg_m3 = 150;
vehicle.rho_fuel_kg_m3 = 807.5;
vehicle.k_vv = 0.1;
vehicle.k_vs = 0.03;
vehicle.k_ve = 0.2641;
vehicle.v_fix_m3 = 11.0;
vehicle.n_ult = 12;
vehicle.Kw = [];
vehicle.lambda_wing = 0.4;
vehicle.tc_max = 0.103;
vehicle.tau = tau;
vehicle.sweep_deg = sweep_deg;

%% Setup Mission
mission = struct();
mission.M_cruise = 1.4;
mission.h_cruise_km = 12.2;
mission.h_cruise_m = 12200;
mission.payload_kg = 0;
mission.crew_count = 1;
mission.range_km = mission.M_cruise * sqrt(1.4 * 287.05 * (288.15 * (1 - 0.0065*mission.h_cruise_m/288.15))) * duration_min * 60 / 1000;

%% Initial Guess
initialGuess = struct();
initialGuess.TOGW_kg = 8000;
initialGuess.S_pln_m2 = 38;

%% Convergence Options
opts = struct();
opts.tolerance = 1e-4;
opts.maxIterations = 100;
opts.units = 'SI';
opts.verbose = false;

%% Run Convergence
if verbose
    fprintf('=================================================================\n');
    fprintf('  Running Convergence for Vehicle Configuration\n');
    fprintf('=================================================================\n');
    fprintf('  Tau:       %.4f\n', tau);
    fprintf('  Sweep:     %d°\n', sweep_deg);
    fprintf('  Duration:  %d min\n', duration_min);
    fprintf('  Range:     %.1f km\n\n', mission.range_km);
end

[design, status] = hypersonic_convergence_elliptical(mission, vehicle, initialGuess, opts);

if ~strcmp(status, 'Converged')
    error('Convergence failed! Status: %s', status);
end

if verbose
    fprintf('✓ Convergence successful!\n\n');
end

%% Get Geometry
geom = test_elliptical_fuselage(design.S_pln_m2, tau, sweep_deg);

%% Extract All Parameters
extracted_data = struct();

%=== DESIGN VARIABLES ===
extracted_data.Tau = vehicle.tau;
extracted_data.Sweep_deg = vehicle.sweep_deg;
extracted_data.Mission_Time_min = mission.range_km / (mission.M_cruise * sqrt(1.4 * 287.05 * (288.15 * (1 - 0.0065*mission.h_cruise_m/288.15))) / 1000 * 60);

%=== VEHICLE PARAMETERS ===
extracted_data.Crew = mission.crew_count;
extracted_data.TOGW_kg = design.TOGW_kg;
extracted_data.OEW_kg = design.OEW_kg;
extracted_data.Fuel_kg = design.TOGW_kg - design.OEW_kg;
extracted_data.S_pln_m2 = design.S_pln_m2;
extracted_data.S_wing_m2 = geom.wing.S;
extracted_data.WS_kg_m2 = design.TOGW_kg / vehicle.S_ref_m2;

%=== WEIGHT BREAKDOWN ===
extracted_data.W_structural_kg = design.W_structural_kg;
extracted_data.W_engine_kg = design.W_engine_kg;
extracted_data.W_systems_kg = design.W_systems_kg;

% Optional fields
if isfield(design, 'W_fixed_kg')
    extracted_data.W_fixed_kg = design.W_fixed_kg;
else
    extracted_data.W_fixed_kg = 0;
end

if isfield(design, 'W_crew_kg')
    extracted_data.W_crew_kg = design.W_crew_kg;
else
    extracted_data.W_crew_kg = mission.crew_count * 91;
end

% Detailed structural components using Roskam weight fractions
% These are stored in the design structure from convergence
if isfield(design, 'W_wing_kg')
    extracted_data.W_wing_kg = design.W_wing_kg;
    extracted_data.W_fuselage_kg = design.W_fuselage_kg;
    extracted_data.W_empennage_kg = design.W_empennage_kg;
    extracted_data.W_gear_kg = design.W_gear_kg;
    extracted_data.W_nacelles_kg = design.W_nacelles_kg;
    if isfield(design, 'W_vtail_kg')
        extracted_data.W_vtail_kg = design.W_vtail_kg;
    end
    if isfield(design, 'W_canard_kg')
        extracted_data.W_canard_kg = design.W_canard_kg;
    end
else
    % Fallback: use Roskam fractions of TOGW
    extracted_data.W_wing_kg = 0.0982 * extracted_data.TOGW_kg;
    extracted_data.W_fuselage_kg = 0.10716 * extracted_data.TOGW_kg;
    extracted_data.W_empennage_kg = 0.030556 * extracted_data.TOGW_kg;
    extracted_data.W_gear_kg = 0.039498 * extracted_data.TOGW_kg;
    extracted_data.W_nacelles_kg = 0.049438 * extracted_data.TOGW_kg;
    extracted_data.W_vtail_kg = 0.014591 * extracted_data.TOGW_kg;
    extracted_data.W_canard_kg = 0.015965 * extracted_data.TOGW_kg;
end

%=== PERFORMANCE METRICS ===
extracted_data.LD_cruise = design.LD_cruise;
extracted_data.Mission_Range_km = mission.range_km;
extracted_data.Fuel_Fraction = extracted_data.Fuel_kg / extracted_data.TOGW_kg;

% Thrust-to-weight ratios
if isfield(design, 'T_W_required')
    extracted_data.Total_TW_cruise = design.T_W_required;
else
    extracted_data.Total_TW_cruise = NaN;
end
extracted_data.Engine_TW_cruise = vehicle.E_TW;

%=== GEOMETRY PARAMETERS ===
% Fuselage
extracted_data.Fuse_Length_m = geom.fuse.L;
extracted_data.Fuse_Width_m = geom.fuse.W;
extracted_data.Fuse_Height_m = geom.fuse.H;
extracted_data.Fuse_Volume_m3 = geom.fuse.V;
extracted_data.Fuse_Swet_m2 = geom.fuse.S_wet;

% Wing
extracted_data.Wing_Span_m = geom.wing.span;
extracted_data.Wing_AR = geom.wing.AR;
extracted_data.Wing_CR_m = geom.wing.c_r;
extracted_data.Wing_CT_m = geom.wing.c_t;
extracted_data.Wing_MAC_m = geom.wing.c_mac;
extracted_data.Wing_Lambda_LE_deg = geom.wing.Lambda_LE_deg;
extracted_data.Wing_Swet_m2 = geom.wing.S_wet;

% Canard
extracted_data.Canard_S_m2 = geom.canard.S;
extracted_data.Canard_Span_m = geom.canard.span;
extracted_data.Canard_Swet_m2 = geom.canard.S_wet;

% Vertical tail
extracted_data.Vtail_S_m2 = geom.vtail.S;
extracted_data.Vtail_Span_m = geom.vtail.span;
extracted_data.Vtail_Swet_m2 = geom.vtail.S_wet;

% Overall
extracted_data.Total_Length_m = geom.fuse.L;
extracted_data.Total_Span_m = geom.wing.span;
extracted_data.V_total_m3 = geom.V_tot;
extracted_data.Wetted_Area_Total_m2 = geom.S_wet_total;
extracted_data.Kw = geom.Kw;

%=== STRUCTURAL INDEX ===
if isfield(design, 'I_str_kg_m2')
    extracted_data.Istr_kg_m2 = design.I_str_kg_m2;
else
    extracted_data.Istr_kg_m2 = extracted_data.W_structural_kg / extracted_data.S_pln_m2;
end

%=== CONVERGENCE METRICS ===
if isfield(design, 'iterations')
    extracted_data.Iterations = design.iterations;
else
    extracted_data.Iterations = NaN;
end

if isfield(design, 'error_final')
    extracted_data.OWE_Error = design.error_final;
else
    extracted_data.OWE_Error = NaN;
end

%=== AERODYNAMIC DETAILS ===
% Get cruise aero at design point
h_m = mission.h_cruise_m;
T0 = 288.15;
T_ratio = 1 - (0.0065 * h_m / T0);
T = T0 * T_ratio;
p = 101325 * T_ratio^(9.81/(287.05*0.0065));
rho = p / (287.05 * T);
a = sqrt(1.4 * 287.05 * T);

atm = struct();
atm.T = T;
atm.p = p;
atm.rho = rho;
atm.a = a;
atm.mu = 1.458e-6 * T^1.5 / (T + 110.4);

% Use raymer_aero_model with elliptical fuselage geometry
opts_aero = struct();
opts_aero.units = 'SI';
opts_aero.Ewave = 2.0;  % Wave drag efficiency (Raymer default for conventional aircraft)
opts_aero.disable_empirical_drags = true;  % Elliptical fuselage (not spherocylinder)

aero = raymer_aero_model(mission.M_cruise, 2.0, geom, atm, opts_aero);
extracted_data.CL_cruise = aero.CL;
extracted_data.CD_cruise = aero.CD;
extracted_data.CD0_cruise = aero.CD0;
extracted_data.K_induced = aero.K;

%=== VOLUME BUDGET ===
extracted_data.V_crew_m3 = 1.2;
extracted_data.V_void_m3 = vehicle.k_vv * geom.V_tot;
extracted_data.V_systems_m3 = vehicle.k_vs * geom.V_tot;
% k_ve is specific volume per thrust (m³/tf), not volume fraction
T_engine_kN = vehicle.E_TW * extracted_data.W_engine_kg * 9.81 / 1000;
T_engine_tf = T_engine_kN / 9.81;
extracted_data.V_engine_m3 = vehicle.k_ve * T_engine_tf;
extracted_data.V_fixed_m3 = vehicle.v_fix_m3;
extracted_data.V_fuel_available_m3 = geom.V_tot - extracted_data.V_crew_m3 - extracted_data.V_void_m3 - ...
    extracted_data.V_systems_m3 - extracted_data.V_engine_m3 - extracted_data.V_fixed_m3;
extracted_data.V_fuel_used_m3 = extracted_data.Fuel_kg / vehicle.rho_fuel_kg_m3;
extracted_data.V_fuel_fit_percent = 100 * extracted_data.V_fuel_used_m3 / extracted_data.V_fuel_available_m3;

%=== AERODYNAMIC MODEL CALIBRATION ===
extracted_data.Ewave = 2.0;  % Wave drag efficiency factor (Raymer default)
extracted_data.Amax_m2 = geom.Amax;
extracted_data.Fineness_Ratio = geom.fuse.L / sqrt(geom.fuse.W * geom.fuse.H);

%% Display Summary (if verbose)
if verbose
    fprintf('=================================================================\n');
    fprintf('CONVERGED VEHICLE PARAMETERS\n');
    fprintf('=================================================================\n\n');

    fprintf('DESIGN VARIABLES:\n');
    fprintf('  Tau:                %.4f\n', extracted_data.Tau);
    fprintf('  Sweep (LE):         %d°\n', extracted_data.Sweep_deg);
    fprintf('  Mission duration:   %.1f min\n', extracted_data.Mission_Time_min);
    fprintf('  Range:              %.1f km\n\n', extracted_data.Mission_Range_km);

    fprintf('WEIGHT SUMMARY:\n');
    fprintf('  TOGW:               %d kg\n', round(extracted_data.TOGW_kg));
    fprintf('  OEW:                %d kg\n', round(extracted_data.OEW_kg));
    fprintf('  Fuel:               %d kg (%.1f%%)\n', round(extracted_data.Fuel_kg), extracted_data.Fuel_Fraction*100);
    fprintf('  - Structural:       %d kg\n', round(extracted_data.W_structural_kg));
    fprintf('  - Engine:           %d kg\n', round(extracted_data.W_engine_kg));
    fprintf('  - Systems:          %d kg\n', round(extracted_data.W_systems_kg));
    fprintf('  - Crew:             %d kg\n\n', round(extracted_data.W_crew_kg));

    fprintf('GEOMETRY:\n');
    fprintf('  S_pln:              %.1f m²\n', extracted_data.S_pln_m2);
    fprintf('  S_wing:             %.2f m²\n', extracted_data.S_wing_m2);
    fprintf('  Length:             %.2f m\n', extracted_data.Total_Length_m);
    fprintf('  Span:               %.2f m\n', extracted_data.Total_Span_m);
    fprintf('  V_total:            %.2f m³\n', extracted_data.V_total_m3);
    fprintf('  Wing loading:       %.1f kg/m²\n\n', extracted_data.WS_kg_m2);

    fprintf('PERFORMANCE:\n');
    fprintf('  L/D cruise:         %.2f\n', extracted_data.LD_cruise);
    fprintf('  CL cruise:          %.4f\n', extracted_data.CL_cruise);
    fprintf('  CD cruise:          %.5f (CD0=%.5f, K=%.4f)\n', extracted_data.CD_cruise, extracted_data.CD0_cruise, extracted_data.K_induced);
    fprintf('  T/W cruise:         %.3f\n\n', extracted_data.Total_TW_cruise);

    fprintf('VOLUME BUDGET:\n');
    fprintf('  Total volume:       %.2f m³\n', extracted_data.V_total_m3);
    fprintf('  - Fuel available:   %.2f m³\n', extracted_data.V_fuel_available_m3);
    fprintf('  - Fuel used:        %.2f m³ (%.1f%% filled)\n', extracted_data.V_fuel_used_m3, extracted_data.V_fuel_fit_percent);
    fprintf('  - Engine:           %.2f m³\n', extracted_data.V_engine_m3);
    fprintf('  - Fixed:            %.2f m³\n', extracted_data.V_fixed_m3);
    fprintf('  - Systems:          %.2f m³\n', extracted_data.V_systems_m3);
    fprintf('  - Crew:             %.2f m³\n', extracted_data.V_crew_m3);
    fprintf('  - Void:             %.2f m³\n\n', extracted_data.V_void_m3);
end

%% Save Data (if requested)
if save_file
    % Determine filename
    if isempty(custom_filename)
        timestamp = datestr(now, 'yyyymmdd_HHMMSS');
        filename = sprintf('vehicle_tau%.4f_sweep%d_dur%.0fmin_%s.mat', ...
            extracted_data.Tau, extracted_data.Sweep_deg, extracted_data.Mission_Time_min, timestamp);
    else
        filename = custom_filename;
        if ~endsWith(filename, '.mat')
            filename = [filename, '.mat'];
        end
    end

    % Find data_files directory
    script_dir = fileparts(mfilename('fullpath'));
    project_root = fileparts(fileparts(script_dir));
    data_files_dir = fullfile(project_root, 'data_files');

    if ~exist(data_files_dir, 'dir')
        mkdir(data_files_dir);
    end

    % Save as MAT file
    mat_filename = fullfile(data_files_dir, filename);
    save(mat_filename, 'extracted_data', 'design', 'geom', 'vehicle', 'mission');

    if verbose
        fprintf('=================================================================\n');
        fprintf('✓ Saved data to:\n  %s\n', mat_filename);
        fprintf('=================================================================\n');
    end
end

end
