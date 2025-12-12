%{
====================================================================================================
SCRIPT NAME: run_mach_sweep_elliptical.m
AUTHOR: Dominic Larin
INITIATED: 11/21/2025
====================================================================================================
SCRIPT DESCRIPTION:
Performs Mach sweep analysis using raymer_aero_model_elliptical at constant altitude and angle
of attack. Generates plots showing variation of aerodynamic coefficients and performance across
the flight envelope.

Test Conditions:
- Altitude: 12.2 km (40,000 ft) - X-29 cruise altitude
- Angle of Attack: 2 degrees (typical cruise condition)
- Mach Range: 0.3 to 1.8 (subsonic through supersonic)
- Focus plots at: M = 0.6, 0.9, 1.2 (subsonic, transonic, supersonic)

Outputs:
- Plots of CL, CD, L/D vs Mach number
- Drag polar plots at selected Mach numbers
- Regime transition identification
====================================================================================================
%}

clear; close all; clc;

% Add all project paths
project_root = fileparts(pwd);  % Go up one level from Aerodynamics to project root
addpath(genpath(project_root));

%% Test Configuration
h = 12200;              % [m] Altitude (12.2 km cruise altitude)
alpha_deg = 6;        % [deg] Angle of attack (typical cruise)

% Use X-29 baseline configuration from test_elliptical_fuselage
Spln = 39.1;          % [m^2] Total planform area (X-29 baseline)
tau = 0.077;             % [-] Slenderness parameter (X-29 baseline)
LambdaLE_deg = -30;   % [deg] Leading edge sweep (X-29 forward sweep)

fprintf('=== Mach Sweep Analysis for Raymer Elliptical Aero Model ===\n');
fprintf('Altitude: %.1f km (%.0f ft)\n', h/1000, h*3.28084);
fprintf('Angle of Attack: %.1f deg\n', alpha_deg);
fprintf('Configuration: Spln=%.2f m², tau=%.3f, ΛLE=%.1f°\n\n', Spln, tau, LambdaLE_deg);

%% Generate Geometry
fprintf('Extracting converged vehicle geometry (10 min, tau=%.3f, sweep=%.0f)...\n', tau, LambdaLE_deg);

% Extract converged vehicle data
converged = extract_converged_vehicle_data(10, tau, LambdaLE_deg);

% Build geometry structure
geom = struct();

% Wing parameters
geom.wing = struct();
geom.wing.S = converged.S_wing_m2;
geom.wing.b = converged.Wing_Span_m;
geom.wing.AR = converged.Wing_AR;
geom.wing.Lambda_LE_deg = converged.Wing_Lambda_LE_deg;

% Calculate quarter-chord sweep using Kshift method (from get_geometry.m)
lam = 0.4;                          % [-] Wing taper ratio
Kshift = (4/geom.wing.AR) * (1-lam)/(1+lam);  % Sweep shift factor
tanLE = tand(geom.wing.Lambda_LE_deg);
geom.wing.Lambda_c25_deg = atand(tanLE - 0.25*Kshift);  % [deg] Quarter-chord sweep

% Fuselage parameters
geom.fuse = struct();
geom.fuse.L = converged.Fuse_Length_m;
geom.fuse.W = converged.Fuse_Width_m;
geom.fuse.H = converged.Fuse_Height_m;
geom.fuse.D = converged.Fuse_Width_m; % Use width as equivalent diameter

% Total wetted area and exposure
geom.S_wet_total = converged.Wetted_Area_Total_m2;
geom.Sexposed_over_Sref = 0.85;     % [-] Exposed wing area fraction
geom.fuselage_d_over_b = converged.Fuse_Width_m / converged.Wing_Span_m;

fprintf('  Wing Area: %.2f m²\n', geom.wing.S);
fprintf('  Wing Span: %.2f m\n', geom.wing.b);
fprintf('  Fuselage: L=%.2f m, W=%.2f m, H=%.2f m\n', ...
    geom.fuse.L, geom.fuse.W, geom.fuse.H);
fprintf('  d/b = %.4f\n\n', geom.fuse.D/geom.wing.b);

%% Mach Sweep
M_array = linspace(0.3, 1.8, 61);  % 61 points from M=0.3 to M=1.8
nM = length(M_array);

% Preallocate output arrays
CL = zeros(nM, 1);
CD = zeros(nM, 1);
CD0 = zeros(nM, 1);
CDi = zeros(nM, 1);
K = zeros(nM, 1);
L_over_D = zeros(nM, 1);
LD_max = zeros(nM, 1);
regime = cell(nM, 1);

fprintf('Running Mach sweep from M=%.1f to M=%.1f...\n', M_array(1), M_array(end));
tic;

for i = 1:nM
    M = M_array(i);

    % Call elliptical aero model
    aero = raymer_aero_model_elliptical(M, alpha_deg, h, geom);

    % Store results
    CL(i) = aero.CL;
    CD(i) = aero.CD;
    CD0(i) = aero.CD0;
    CDi(i) = aero.CDi;
    K(i) = aero.K;
    L_over_D(i) = aero.L_over_D;
    LD_max(i) = aero.LD_max;
    regime{i} = aero.regime;

    % Progress indicator
    if mod(i, 10) == 0
        fprintf('  M=%.2f: CL=%.4f, CD=%.4f, L/D=%.2f (%s)\n', ...
            M, CL(i), CD(i), L_over_D(i), regime{i});
    end
end

elapsed = toc;
fprintf('Sweep completed in %.2f seconds.\n\n', elapsed);

%% Summary Statistics
fprintf('=== Performance Summary ===\n');
[max_LD, idx_max_LD] = max(L_over_D);
fprintf('Maximum L/D: %.2f at M=%.2f (%s)\n', max_LD, M_array(idx_max_LD), regime{idx_max_LD});

% Find cruise condition (M=1.4)
[~, idx_cruise] = min(abs(M_array - 1.4));
fprintf('Cruise (M=1.4): CL=%.4f, CD=%.4f, L/D=%.2f\n', ...
    CL(idx_cruise), CD(idx_cruise), L_over_D(idx_cruise));

% Regime transitions
transitions = find(diff(strcmp(regime(1:end-1), regime(2:end))) ~= 0);
fprintf('\nRegime Transitions:\n');
for i = 1:length(transitions)
    idx = transitions(i);
    fprintf('  M=%.3f: %s → %s\n', M_array(idx), regime{idx}, regime{idx+1});
end

%% Plot 1: Aerodynamic Coefficients vs Mach
figure('Position', [100, 100, 1200, 800], 'Name', 'Mach Sweep - Aero Coefficients');

% CL vs Mach
subplot(2,2,1)
plot(M_array, CL, 'b-', 'LineWidth', 2);
hold on;
plot(1.4, CL(idx_cruise), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('Mach Number', 'FontSize', 12);
ylabel('C_L', 'FontSize', 12);
title(sprintf('Lift Coefficient (\\alpha = %.1f°)', alpha_deg), 'FontSize', 14);
legend('C_L', 'Cruise (M=1.4)', 'Location', 'best');

% CD vs Mach
subplot(2,2,2)
plot(M_array, CD, 'r-', 'LineWidth', 2);
hold on;
plot(M_array, CD0, 'g--', 'LineWidth', 1.5);
plot(M_array, CDi, 'm--', 'LineWidth', 1.5);
plot(1.4, CD(idx_cruise), 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'k');
grid on;
xlabel('Mach Number', 'FontSize', 12);
ylabel('C_D', 'FontSize', 12);
title('Drag Coefficients', 'FontSize', 14);
legend('C_D Total', 'C_{D0} (Parasitic)', 'C_{Di} (Induced)', 'Cruise', 'Location', 'best');

% L/D vs Mach
subplot(2,2,3)
plot(M_array, L_over_D, 'b-', 'LineWidth', 2);
hold on;
plot(M_array, LD_max, 'g--', 'LineWidth', 1.5);
plot(1.4, L_over_D(idx_cruise), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('Mach Number', 'FontSize', 12);
ylabel('L/D', 'FontSize', 12);
title('Lift-to-Drag Ratio', 'FontSize', 14);
legend('L/D (actual)', '(L/D)_{max}', 'Cruise (M=1.4)', 'Location', 'best');

% K factor vs Mach
subplot(2,2,4)
plot(M_array, K, 'k-', 'LineWidth', 2);
hold on;
plot(1.4, K(idx_cruise), 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
grid on;
xlabel('Mach Number', 'FontSize', 12);
ylabel('K Factor', 'FontSize', 12);
title('Induced Drag Factor', 'FontSize', 14);
legend('K', 'Cruise (M=1.4)', 'Location', 'best');

sgtitle(sprintf('X-29 Elliptical Model - Mach Sweep (h=%.1f km, \\alpha=%.1f°)', h/1000, alpha_deg), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Plot 2: Drag Polars at Selected Mach Numbers
M_selected = [0.6, 0.9, 1.2];  % Subsonic, transonic, supersonic
alpha_array = linspace(-2, 8, 21);  % Alpha sweep for drag polar
nAlpha = length(alpha_array);

figure('Position', [150, 150, 1400, 500], 'Name', 'Drag Polars at Selected Mach');

for iM = 1:length(M_selected)
    M_target = M_selected(iM);

    % Preallocate
    CL_polar = zeros(nAlpha, 1);
    CD_polar = zeros(nAlpha, 1);
    LD_polar = zeros(nAlpha, 1);

    % Alpha sweep at this Mach
    for iAlpha = 1:nAlpha
        alpha = alpha_array(iAlpha);
        aero = raymer_aero_model_elliptical(M_target, alpha, h, geom);
        CL_polar(iAlpha) = aero.CL;
        CD_polar(iAlpha) = aero.CD;
        LD_polar(iAlpha) = aero.L_over_D;
    end

    % Plot drag polar (CL vs CD)
    subplot(1,3,iM)
    plot(CD_polar, CL_polar, 'b-', 'LineWidth', 2);
    hold on;

    % Mark alpha = 2 deg (cruise)
    [~, idx_cruise_alpha] = min(abs(alpha_array - alpha_deg));
    plot(CD_polar(idx_cruise_alpha), CL_polar(idx_cruise_alpha), 'ro', ...
        'MarkerSize', 10, 'MarkerFaceColor', 'r');

    grid on;
    xlabel('C_D', 'FontSize', 12);
    ylabel('C_L', 'FontSize', 12);
    title(sprintf('M = %.1f', M_target), 'FontSize', 14);
    legend('Drag Polar', sprintf('\\alpha=%.0f°', alpha_deg), 'Location', 'best');
    axis([0 max(CD_polar)*1.1 min(CL_polar) max(CL_polar)*1.1]);
end

sgtitle(sprintf('X-29 Elliptical Model - Drag Polars (h=%.1f km)', h/1000), ...
    'FontSize', 16, 'FontWeight', 'bold');

%% Plot 3: L/D vs Alpha at Selected Mach Numbers
figure('Position', [200, 200, 900, 600], 'Name', 'L/D vs Alpha');

colors = {'b', 'r', 'g'};
markers = {'o', 's', '^'};

for iM = 1:length(M_selected)
    M_target = M_selected(iM);

    % Preallocate
    LD_alpha = zeros(nAlpha, 1);

    % Alpha sweep at this Mach
    for iAlpha = 1:nAlpha
        alpha = alpha_array(iAlpha);
        aero = raymer_aero_model_elliptical(M_target, alpha, h, geom);
        LD_alpha(iAlpha) = aero.L_over_D;
    end

    % Plot
    plot(alpha_array, LD_alpha, [colors{iM} '-'], 'LineWidth', 2, 'DisplayName', sprintf('M=%.1f', M_target));
    hold on;

    % Mark cruise alpha
    [~, idx_cruise_alpha] = min(abs(alpha_array - alpha_deg));
    plot(alpha_array(idx_cruise_alpha), LD_alpha(idx_cruise_alpha), ...
        [colors{iM} markers{iM}], 'MarkerSize', 10, 'MarkerFaceColor', colors{iM}, ...
        'HandleVisibility', 'off');
end

grid on;
xlabel('Angle of Attack [deg]', 'FontSize', 12);
ylabel('L/D', 'FontSize', 12);
title(sprintf('Lift-to-Drag Ratio vs Angle of Attack (h=%.1f km)', h/1000), 'FontSize', 14);
legend('Location', 'best', 'FontSize', 11);

%% Save Results
fprintf('\n=== Saving Results ===\n');

% Save figure
saveas(gcf, 'X29_Elliptical_Mach_Sweep_LD_vs_Alpha.png');
fprintf('Saved: X29_Elliptical_Mach_Sweep_LD_vs_Alpha.png\n');

% Save workspace
save('mach_sweep_elliptical_results.mat', 'M_array', 'CL', 'CD', 'CD0', 'CDi', ...
    'K', 'L_over_D', 'LD_max', 'regime', 'geom', 'h', 'alpha_deg');
fprintf('Saved: mach_sweep_elliptical_results.mat\n');

% Export data to CSV
results_table = table(M_array', CL, CD, CD0, CDi, K, L_over_D, LD_max, regime, ...
    'VariableNames', {'Mach', 'CL', 'CD', 'CD0', 'CDi', 'K', 'L_over_D', 'LD_max', 'Regime'});
writetable(results_table, 'mach_sweep_elliptical_results.csv');
fprintf('Saved: mach_sweep_elliptical_results.csv\n');

fprintf('\nAnalysis complete!\n');


%% Supersonic percent error analysis vs literature (M = 1.2 to 1.6)
% Compare test configuration at alpha = 6 deg to:
%  - Budd Mach-based CL, CD  baseline  AoA not accounted for
%  - Budd Mach–AoA lookup CL, CD at alpha = 6 deg
%  - Saltzman and Hicks ACC L/D  true in flight

fprintf('\n=== Supersonic Percent Error Analysis vs Literature (M = 1.2 to 1.6) ===\n');

% Supersonic Mach grid for comparison
M_sup = 1.2:0.05:1.6;
nSup  = numel(M_sup);

% Model results for test configuration  alpha_deg = 6
CL_model_sup = interp1(M_array, CL,       M_sup, 'linear');
CD_model_sup = interp1(M_array, CD,       M_sup, 'linear');
LD_model_sup = interp1(M_array, L_over_D, M_sup, 'linear');

% Preallocate literature arrays
CL_LLL_sup    = nan(size(M_sup));   % Budd Mach-based  baseline
CD_LLL_sup    = nan(size(M_sup));
LD_LLL_sup    = nan(size(M_sup));   % used for plots

CL_lookup_sup = nan(size(M_sup));   % Budd lookup Mach,alpha
CD_lookup_sup = nan(size(M_sup));
LD_lookup_sup = nan(size(M_sup));   % used for plots

LD_ACC_sup    = nan(size(M_sup));   % Saltzman and Hicks ACC L/D

for i = 1:nSup
    M = M_sup(i);

    if M >= 1.2 && M <= 1.5
        % Within shared Mach range where lookup data are valid
        [LLL, Lookup, ACC] = x29a_aero_all_sources(M, alpha_deg, true);

        % Budd Mach-based  baseline  AoA not accounted for
        CL_LLL_sup(i) = LLL.CL;
        CD_LLL_sup(i) = LLL.CD;
        LD_LLL_sup(i) = LLL.L_over_D;

        % Budd lookup at alpha = 6 deg
        CL_lookup_sup(i) = Lookup.CL;
        CD_lookup_sup(i) = Lookup.CD;
        LD_lookup_sup(i) = Lookup.L_over_D;

        % Saltzman and Hicks ACC L/D
        LD_ACC_sup(i) = ACC.L_over_D;

    else
        % For M > 1.5 there is no lookup data
        [CD_LLL, CL_LLL, LD_LLL] = x29a_aero_coeffs_Lit_LLL(M, true);
        CL_LLL_sup(i) = CL_LLL;
        CD_LLL_sup(i) = CD_LLL;
        LD_LLL_sup(i) = LD_LLL;

        LD_ACC_sup(i) = X29A_LtoD_ACC_9p14km(M);
        % CL_lookup_sup and CD_lookup_sup remain NaN here
    end
end

% Percent error helper
percent_error = @(measured, actual) 100 * (measured - actual) ./ actual;

% Baseline Budd Mach-based errors
PE_CL_LLL = percent_error(CL_model_sup, CL_LLL_sup);
PE_CD_LLL = percent_error(CD_model_sup, CD_LLL_sup);

% Budd Mach–AoA lookup errors  alpha = 6 deg
PE_CL_lookup = percent_error(CL_model_sup, CL_lookup_sup);
PE_CD_lookup = percent_error(CD_model_sup, CD_lookup_sup);

% ACC L/D errors
PE_LD_ACC = percent_error(LD_model_sup, LD_ACC_sup);

% Mean absolute percent errors  ignore NaNs
MAPE_CL_LLL    = mean(abs(PE_CL_LLL(~isnan(PE_CL_LLL))));
MAPE_CD_LLL    = mean(abs(PE_CD_LLL(~isnan(PE_CD_LLL))));
MAPE_CL_lookup = mean(abs(PE_CL_lookup(~isnan(PE_CL_lookup))));
MAPE_CD_lookup = mean(abs(PE_CD_lookup(~isnan(PE_CD_lookup))));
MAPE_LD_ACC    = mean(abs(PE_LD_ACC(~isnan(PE_LD_ACC))));

% Columns with constant MAPE values for table display
MAPE_CL_LLL_col    = repmat(MAPE_CL_LLL,    nSup, 1);
MAPE_CD_LLL_col    = repmat(MAPE_CD_LLL,    nSup, 1);
MAPE_CL_lookup_col = repmat(MAPE_CL_lookup, nSup, 1);
MAPE_CD_lookup_col = repmat(MAPE_CD_lookup, nSup, 1);
MAPE_LD_ACC_col    = repmat(MAPE_LD_ACC,    nSup, 1);

%% Tables  2 CL  2 CD  1 L/D

% 1  CL vs Budd Mach-based  baseline
fprintf('\nCL comparison vs Budd Mach-based dataset  baseline  AoA not accounted for\n');
T_CL_LLL = table(M_sup.', CL_LLL_sup.', CL_model_sup.', PE_CL_LLL.', MAPE_CL_LLL_col, ...
    'VariableNames', {'Mach', 'CL_Lit_baseline', 'CL_model_measured', ...
                      'PercentError', 'MAPE_CL_baseline'});
disp(T_CL_LLL);

% 2  CD vs Budd Mach-based  baseline
fprintf('\nCD comparison vs Budd Mach-based dataset  baseline  AoA not accounted for\n');
T_CD_LLL = table(M_sup.', CD_LLL_sup.', CD_model_sup.', PE_CD_LLL.', MAPE_CD_LLL_col, ...
    'VariableNames', {'Mach', 'CD_Lit_baseline', 'CD_model_measured', ...
                      'PercentError', 'MAPE_CD_baseline'});
disp(T_CD_LLL);

% 3  CL vs Budd Mach–AoA lookup  alpha = 6 deg
fprintf('\nCL comparison vs Budd Mach–AoA lookup dataset  alpha = 6 deg\n');
T_CL_lookup = table(M_sup.', CL_lookup_sup.', CL_model_sup.', PE_CL_lookup.', MAPE_CL_lookup_col, ...
    'VariableNames', {'Mach', 'CL_Lit_lookup', 'CL_model_measured', ...
                      'PercentError', 'MAPE_CL_lookup'});
disp(T_CL_lookup);

% 4  CD vs Budd Mach–AoA lookup  alpha = 6 deg
fprintf('\nCD comparison vs Budd Mach–AoA lookup dataset  alpha = 6 deg\n');
T_CD_lookup = table(M_sup.', CD_lookup_sup.', CD_model_sup.', PE_CD_lookup.', MAPE_CD_lookup_col, ...
    'VariableNames', {'Mach', 'CD_Lit_lookup', 'CD_model_measured', ...
                      'PercentError', 'MAPE_CD_lookup'});
disp(T_CD_lookup);

% 5  L/D vs Saltzman and Hicks ACC dataset
fprintf('\nL/D comparison vs Saltzman and Hicks ACC dataset  true in-flight L/D\n');
T_LD_ACC = table(M_sup.', LD_ACC_sup.', LD_model_sup.', PE_LD_ACC.', MAPE_LD_ACC_col, ...
    'VariableNames', {'Mach', 'LD_ACC_actual', 'LD_model_measured', ...
                      'PercentError', 'MAPE_LD_ACC'});
disp(T_LD_ACC);

%% Supersonic CL, CD, and L/D vs Mach with literature overlay

figure('Position', [200, 100, 1100, 750], 'Name', 'Supersonic Aero Coefficients vs Mach');

% ----- CL vs Mach (top left) -----
subplot(2,2,1);
hold on; grid on;
plot(M_sup, CL_model_sup, 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Model  test config  \alpha = 6°');
plot(M_sup, CL_LLL_sup,   'bo-', 'LineWidth', 1.5, ...
    'DisplayName', 'Budd Mach-based  baseline');
plot(M_sup, CL_lookup_sup,'rs--','LineWidth', 1.5, ...
    'DisplayName', 'Budd lookup  Mach,\alpha = 6°');
ylabel('C_L');
xlim([1.2 1.6]);
title('Supersonic C_L vs Mach');
legend('Location', 'best');

% ----- CD vs Mach (top right) -----
subplot(2,2,2);
hold on; grid on;
plot(M_sup, CD_model_sup, 'k-', 'LineWidth', 2, ...
    'DisplayName', 'Model  test config  \alpha = 6°');
plot(M_sup, CD_LLL_sup,   'bo-', 'LineWidth', 1.5, ...
    'DisplayName', 'Budd Mach-based  baseline');
plot(M_sup, CD_lookup_sup,'rs--','LineWidth', 1.5, ...
    'DisplayName', 'Budd lookup  Mach,\alpha = 6°');
ylabel('C_D');
xlim([1.2 1.6]);
title('Supersonic C_D vs Mach');
legend('Location', 'best');

% ----- L/D vs Mach (bottom row, spans both columns) -----
subplot(2,2,[3 4]);
hold on; grid on;
plot(M_sup, LD_model_sup, 'k-',  'LineWidth', 2, ...
    'DisplayName', 'Model  test config  \alpha = 6°');
plot(M_sup, LD_ACC_sup,   'gs--','LineWidth', 1.5, ...
    'DisplayName', 'Saltzman and Hicks ACC');
xlabel('Mach');
ylabel('L/D');
xlim([1.2 1.6]);
title('Supersonic L/D vs Mach');
legend('Location', 'best');
