% RUN_ELLIPTICAL_LOAD_CASES
% Driver script to compute elliptical spanwise lift distributions
% for many Mach/altitude combinations using raymer_aero_model_elliptical
% and the converged X-29 geometry.

clear; clc; close all;

%% 0. Optional: add project root to path (adjust as needed)
this_file    = mfilename('fullpath');
project_root = fileparts(this_file);
addpath(genpath(project_root));

fprintf('=== Elliptical Load Solver: Multiple Mach/Altitude Cases ===\n\n');

%% 1. Build converged X-29 geometry (same as run_mach_sweep_elliptical)

% X-29 baseline configuration parameters
tau          = 0.077;    % slenderness parameter
LambdaLE_deg = -30;      % leading-edge forward sweep [deg]

fprintf('Extracting converged geometry (tau = %.3f, Lambda_{LE} = %.1f deg)...\n', ...
    tau, LambdaLE_deg);

converged = extract_converged_vehicle_data(10, tau, LambdaLE_deg);

geom = struct();

% Wing
geom.wing = struct();
geom.wing.S             = converged.S_wing_m2;
geom.wing.b             = converged.Wing_Span_m;
geom.wing.span          = converged.Wing_Span_m;
geom.wing.AR            = converged.Wing_AR;
geom.wing.Lambda_LE_deg = converged.Wing_Lambda_LE_deg;

% Quarter-chord sweep (for completeness; not required directly here)
lam    = 0.4;  % taper ratio (assumption used elsewhere in project)
Kshift = (4/geom.wing.AR) * (1 - lam) / (1 + lam);
tanLE  = tand(geom.wing.Lambda_LE_deg);
geom.wing.Lambda_c25_deg = atand(tanLE - 0.25 * Kshift);

% Fuselage
geom.fuse = struct();
geom.fuse.L = converged.Fuse_Length_m;
geom.fuse.W = converged.Fuse_Width_m;
geom.fuse.H = converged.Fuse_Height_m;
geom.fuse.D = converged.Fuse_Width_m;  % width as equivalent diameter

% Wetted area and exposure (required by raymer_aero_model_elliptical)
geom.S_wet_total        = converged.Wetted_Area_Total_m2;
geom.Sexposed_over_Sref = 0.85;  % exposed wing area fraction (project constant)
geom.fuselage_d_over_b  = converged.Fuse_Width_m / converged.Wing_Span_m;

fprintf('  Wing Area   S = %.3f m^2\n', geom.wing.S);
fprintf('  Wing Span   b = %.3f m\n',   geom.wing.b);
fprintf('  Aspect Ratio AR = %.3f\n',   geom.wing.AR);
fprintf('  Fuselage L = %.3f m, W = %.3f m, H = %.3f m\n', ...
    geom.fuse.L, geom.fuse.W, geom.fuse.H);
fprintf('  d/b = %.4f\n\n', geom.fuse.D / geom.wing.b);

%% 2. Define Mach/altitude cases

% EDIT THESE VECTORS to match the DCFCs
Mach_list = [0.6, 0.9, 1.2, 1.6];         % example Mach numbers
h_list    = [10e3, 10e3, 12.2e3, 15.3e3]; % [m] example altitudes

if numel(Mach_list) ~= numel(h_list)
    error('Mach_list and h_list must have the same length.');
end

nCases = numel(Mach_list);

% Wing lift fraction and AoA schedule
lambda_w  = 0.8;              % fraction of total lift carried by wing
alpha_vec = (0:2:10).';       % [deg] AoA grid (edit if needed)

%% 3. Loop over cases and compute loads

Cases = struct([]);

for k = 1:nCases
    M_k = Mach_list(k);
    h_k = h_list(k);

    fprintf('Case %d/%d: M = %.2f, h = %.1f km\n', ...
        k, nCases, M_k, h_k/1000);

    [y, eta, L_prime, L_total_vec, alpha_deg_vec, CL_vec] = ...
        compute_elliptical_load_case(M_k, h_k, geom, lambda_w, alpha_vec);

    % Store results in struct for later plotting or export
    Cases(k).Mach          = M_k;
    Cases(k).h_m           = h_k;
    Cases(k).y             = y;
    Cases(k).eta           = eta;
    Cases(k).L_prime       = L_prime;
    Cases(k).L_total_vec   = L_total_vec;
    Cases(k).alpha_deg_vec = alpha_deg_vec;
    Cases(k).CL_vec        = CL_vec;

    % Identify the AoA with maximum total wing lift (useful for max-load case)
    [L_max, idx_max]             = max(L_total_vec);
    Cases(k).L_max               = L_max;
    Cases(k).alpha_deg_max_load  = alpha_deg_vec(idx_max);
    Cases(k).L_prime_max_load    = L_prime(:, idx_max);

    fprintf('  Max lift in this case: L = %.3f MN at alpha = %.2f deg\n\n', ...
        L_max / 1e6, alpha_deg_vec(idx_max));
end

%% 4. Plot spanwise load for the maximum-load AoA of each case

figure;
hold on; grid on; box on;

colors = lines(nCases);

for k = 1:nCases
    y_plot   = Cases(k).y;
    L_plot   = Cases(k).L_prime_max_load;
    alpha_ml = Cases(k).alpha_deg_max_load;

    plot(y_plot, L_plot, ...
        'LineWidth', 1.6, 'Color', colors(k,:), ...
        'DisplayName', sprintf('M = %.2f, h = %.1f km, \\alpha_{max} = %.1f^\\circ', ...
                               Cases(k).Mach, Cases(k).h_m/1000, alpha_ml));
end

xlabel('Spanwise coordinate y [m]');
ylabel('Lift per unit span L''(y) [N/m]');
title('Elliptical spanwise lift at max-load \alpha for each Mach/altitude case');
legend('Location', 'best');

%% 5. Example: total wing lift vs angle of attack for one case

% Choose a case index to visualize (for example, last case)
k_plot = nCases;

figure;
plot(Cases(k_plot).alpha_deg_vec, Cases(k_plot).L_total_vec, 'o-', 'LineWidth', 1.6);
grid on; box on;
xlabel('Angle of attack \alpha [deg]');
ylabel('Total wing lift L_{total} [N]');
title(sprintf('Total Wing Lift vs \\alpha, M = %.2f, h = %.1f km', ...
      Cases(k_plot).Mach, Cases(k_plot).h_m/1000));

%% 6. Example: write CSV files for structural handoff (one per case)

write_csv = true;  % set false if CSV outputs are not needed

if write_csv
    out_dir = fullfile(project_root, 'elliptical_load_outputs');
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end

    for k = 1:nCases
        y_case      = Cases(k).y;
        Lprime_case = Cases(k).L_prime_max_load;
        alpha_ml    = Cases(k).alpha_deg_max_load;

        T = table(y_case, Lprime_case, ...
            'VariableNames', {'y_m', 'Lprime_N_per_m'});

        fname = sprintf('elliptical_load_M%.2f_h%.1fkm_alpha%.1fdeg.csv', ...
            Cases(k).Mach, Cases(k).h_m/1000, alpha_ml);
        fname = strrep(fname, '.', 'p'); % avoid multiple dots in filename

        writetable(T, fullfile(out_dir, fname));
    end

    fprintf('CSV load tables written to:\n  %s\n', out_dir);
end
