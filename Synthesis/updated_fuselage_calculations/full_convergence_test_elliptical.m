%{
====================================================================================================
SCRIPT NAME: full_convergence_test_elliptical.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 10/13/2025
====================================================================================================
SCRIPT DESCRIPTION:
Complete trade study convergence test using the new elliptical fuselage model.
This script replicates execute_design_trade_study.m exactly, but uses test_elliptical_fuselage
instead of get_geometry to validate that the new geometry model produces viable aircraft.

Outputs:
1. Wing loading vs Sref plots
2. Solution space colored by tau
3. Comparison with original spherocylinder results (if available)
4. Full convergence statistics and geometry validation

This script follows the EXACT workflow from execute_design_trade_study.m with all the same
parameters, constants, and disciplinary analysis to ensure apples-to-apples comparison.
====================================================================================================
%}

%% 1. SETUP
clear; clc; close all;
theme('light');
set(groot, 'defaultFigureColor', 'white');
set(groot, 'defaultAxesColor', 'white');
set(groot, 'defaultFigureInvertHardcopy', 'off');

% Add paths
current_dir = fileparts(mfilename('fullpath'));
parent_dir = fileparts(fileparts(current_dir));
addpath(genpath(parent_dir));

fprintf('=================================================================\n');
fprintf('    MAE 4350 X-29 Trade Study - ELLIPTICAL FUSELAGE TEST\n');
fprintf('=================================================================\n');

%% 2. DEFINE TRADE MATRIX PARAMETERS
% EXACT same as execute_design_trade_study.m
fprintf('Setting up trade matrix parameters...\n');

trade.sweep_LE_deg = [-30];  % Focus on X-29 sweep
trade.tau_min = 0.01;
trade.tau_max = 0.1;
trade.tau_steps = 19;
trade.tau_values = linspace(trade.tau_min, trade.tau_max, trade.tau_steps);
trade.cruise_duration_min = [1, 5, 10, 15, 20, 25, 30];

n_sweep = length(trade.sweep_LE_deg);
n_tau = length(trade.tau_values);
n_duration = length(trade.cruise_duration_min);
total_design_points = n_sweep * n_tau * n_duration;

fprintf('Trade Matrix Summary:\n');
fprintf('  Wing Sweep Angles: %d values from %.0f° to %.0f°\n', n_sweep, min(trade.sweep_LE_deg), max(trade.sweep_LE_deg));
fprintf('  Slenderness (τ): %d values from %.3f to %.3f\n', n_tau, trade.tau_min, trade.tau_max);
fprintf('  Cruise Durations: %d values from %d to %d minutes\n', n_duration, min(trade.cruise_duration_min), max(trade.cruise_duration_min));
fprintf('  Total Design Points: %d\n\n', total_design_points);

%% 3. MISSION AND AIRCRAFT PARAMETERS
% EXACT same as execute_design_trade_study.m
fprintf('Setting up mission and aircraft parameters...\n');

mission.M_cruise = 1.4;
mission.h_cruise_km = 12.2;
mission.h_cruise_m = mission.h_cruise_km * 1000;
mission.M_max = 1.6;
mission.h_max_km = 15.3;

aircraft.engine_type = 'GE F404';

S_pln_base = 35.0;
TOGW_base = 8000;

%% 4. INITIALIZE RESULTS STORAGE
fprintf('Initializing storage for all trade study results...\n');

successful_configs = [];
failed_configs = [];

%% 5. SETUP VEHICLE CONFIGURATION
% EXACT same parameters as execute_design_trade_study.m
fprintf('Setting up vehicle configuration...\n');

vehicle = struct();
vehicle.S_ref_m2 = 17.544;  % X-29 exposed wing reference area
vehicle.E_TW = 7.323;
vehicle.mu_a = 0.0;
vehicle.f_sys = 0.2;
vehicle.C_sys_kg = 1250;
vehicle.rho_pay_kg_m3 = 150;
vehicle.rho_fuel_kg_m3 = 807.5;
vehicle.k_vv = 0.1;
vehicle.k_vs = 0.03;
vehicle.k_ve = 0.2641;
vehicle.v_fix_m3 = 11.5;  % Increased from 9.0 to force larger S_pln (more volume needed)
vehicle.n_ult = 12;
vehicle.Kw = 1.00;
vehicle.lambda_wing = 0.4;
vehicle.tc_max = 0.103;

mission_enhanced = struct();
mission_enhanced.payload_kg = 0;
mission_enhanced.crew_count = 1;

initialGuess = struct();
initialGuess.TOGW_kg = 8000;
initialGuess.S_pln_m2 = 35;

opts = struct();
opts.tolerance = 1e-4;
opts.maxIterations = 100;
opts.units = 'SI';
opts.verbose = false;

fprintf('Vehicle configuration complete.\n\n');

%% 6. EXECUTE CONVERGENCE TRADE STUDY
fprintf('\nRunning Elliptical Fuselage Convergence Analysis...\n');
fprintf('This uses test_elliptical_fuselage() instead of get_geometry()\n');

design_counter = 0;
progress_bar_length = 50;
fprintf('Progress: [%s] 0%% (0/%d)', repmat(' ', 1, progress_bar_length), total_design_points);

for i_sweep = 1:n_sweep
    for i_tau = 1:n_tau
        for i_duration = 1:n_duration
            design_counter = design_counter + 1;

            % Update progress bar
            if mod(design_counter, max(1, floor(total_design_points/5))) == 0 || design_counter == total_design_points
                percent = floor(100 * design_counter / total_design_points);
                filled_length = floor(progress_bar_length * design_counter / total_design_points);
                progress_str = [repmat('=', 1, filled_length), repmat(' ', 1, progress_bar_length - filled_length)];
                fprintf('\rProgress: [%s] %3d%% (%d/%d)', progress_str, percent, design_counter, total_design_points);
            end

            try
                % Current design point parameters
                current_sweep = trade.sweep_LE_deg(i_sweep);
                current_tau = trade.tau_values(i_tau);
                current_duration = trade.cruise_duration_min(i_duration);

                % Update vehicle configuration
                vehicle_config = vehicle;
                vehicle_config.tau = current_tau;
                vehicle_config.sweep_deg = current_sweep;

                % Update mission
                mission_config = mission;
                mission_config.payload_kg = mission_enhanced.payload_kg;
                mission_config.crew_count = mission_enhanced.crew_count;
                mission_config.range_km = calculate_range_from_duration(current_duration, mission.M_cruise);

                % Run convergence with ELLIPTICAL FUSELAGE
                % This calls hypersonic_convergence_elliptical which uses test_elliptical_fuselage
                [design, status] = hypersonic_convergence_elliptical(mission_config, vehicle_config, initialGuess, opts);

                if strcmp(status, 'Converged')
                    % Store comprehensive result data
                    result = struct();

                    % Basic configuration parameters
                    result.sweep_deg = current_sweep;
                    result.tau = current_tau;
                    result.duration_min = current_duration;
                    result.range_km = mission_config.range_km;
                    result.payload_kg = mission_config.payload_kg;
                    result.crew_count = mission_config.crew_count;

                    % Converged design parameters
                    result.TOGW_kg = design.TOGW_kg;
                    result.OEW_kg = design.OEW_kg;
                    result.S_pln_m2 = design.S_pln_m2;
                    result.wing_loading = design.wing_loading;
                    result.fuel_fraction = design.fuel_fraction;
                    result.LD_cruise = design.LD_cruise;
                    result.T_W_required = design.T_W_required;
                    result.iterations = design.iterations;

                    % Derived parameters
                    result.fuel_weight_kg = result.TOGW_kg * result.fuel_fraction;
                    result.engine_TW = vehicle_config.E_TW;
                    result.vehicle_TW = result.T_W_required;

                    % Weight breakdown
                    result.W_structural_kg = design.W_structural_kg;
                    if isfield(design, 'W_wing_kg')
                        result.W_wing_kg = design.W_wing_kg;
                        result.W_fuselage_kg = design.W_fuselage_kg;
                        result.W_empennage_kg = design.W_empennage_kg;
                        result.W_gear_kg = design.W_gear_kg;
                        result.W_nacelles_kg = design.W_nacelles_kg;
                        result.structural_method = design.structural_method;
                    end
                    result.W_engine_kg = design.W_engine_kg;
                    result.W_systems_kg = design.W_systems_kg;
                    result.I_str_kg_m2 = design.I_str_kg_m2;

                    % Get ELLIPTICAL geometry data
                    try
                        geom = test_elliptical_fuselage(result.S_pln_m2, current_tau, current_sweep);
                        result.volume_m3 = geom.V_tot;
                        result.length_m = geom.Lref;
                        result.wingspan_m = geom.wing.span;
                        result.wetted_area_m2 = geom.S_wet_total;
                        result.Kw = geom.Kw;
                        result.wing_area_m2 = geom.Sref;

                        % NEW: Store elliptical fuselage geometry
                        result.fuse_width_m = geom.fuse.W;
                        result.fuse_height_m = geom.fuse.H;
                        result.fuse_aspect_ratio = geom.fuse.aspect_ratio;
                        result.fuse_fineness_ratio = geom.fuse.fineness_ratio;
                        result.fuse_diameter_equiv_m = geom.fuse.D;

                    catch ME_geom
                        warning('Geometry extraction failed for tau=%.3f, sweep=%d: %s', current_tau, current_sweep, ME_geom.message);
                    end

                    % Store in successful configs
                    successful_configs = [successful_configs; result];

                else
                    % Store failure
                    failure = struct();
                    failure.sweep_deg = current_sweep;
                    failure.tau = current_tau;
                    failure.duration_min = current_duration;
                    failure.status = status;
                    failed_configs = [failed_configs; failure];
                end

            catch ME
                % Store error
                error_record = struct();
                error_record.sweep_deg = current_sweep;
                error_record.tau = current_tau;
                error_record.duration_min = current_duration;
                error_record.error_message = ME.message;
                failed_configs = [failed_configs; error_record];
            end
        end
    end
end

fprintf('\rProgress: [%s] 100%% (%d/%d)\n', repmat('=', 1, progress_bar_length), total_design_points, total_design_points);
fprintf('Convergence sweep complete!\n\n');

%% 7. RESULTS ANALYSIS
fprintf('=================================================================\n');
fprintf('                    CONVERGENCE RESULTS\n');
fprintf('=================================================================\n');

n_successful = length(successful_configs);
n_failed = length(failed_configs);
success_rate = 100 * n_successful / total_design_points;

fprintf('Convergence Statistics:\n');
fprintf('  Successful: %d / %d (%.1f%%)\n', n_successful, total_design_points, success_rate);
fprintf('  Failed: %d / %d (%.1f%%)\n\n', n_failed, total_design_points, 100-success_rate);

if n_successful > 0
    fprintf('SUCCESSFUL DESIGNS SUMMARY:\n');
    fprintf('  TOGW range: %.0f to %.0f kg (mean: %.0f)\n', ...
            min([successful_configs.TOGW_kg]), max([successful_configs.TOGW_kg]), mean([successful_configs.TOGW_kg]));
    fprintf('  OEW range: %.0f to %.0f kg (mean: %.0f)\n', ...
            min([successful_configs.OEW_kg]), max([successful_configs.OEW_kg]), mean([successful_configs.OEW_kg]));
    fprintf('  Spln range: %.1f to %.1f m² (mean: %.1f)\n', ...
            min([successful_configs.S_pln_m2]), max([successful_configs.S_pln_m2]), mean([successful_configs.S_pln_m2]));
    fprintf('  Wing loading range: %.0f to %.0f kg/m² (mean: %.0f)\n', ...
            min([successful_configs.wing_loading]), max([successful_configs.wing_loading]), mean([successful_configs.wing_loading]));
    fprintf('  L/D range: %.2f to %.2f (mean: %.2f)\n', ...
            min([successful_configs.LD_cruise]), max([successful_configs.LD_cruise]), mean([successful_configs.LD_cruise]));

    fprintf('\nWEIGHT BREAKDOWN (at tau=0.07 baseline):\n');
    % Find X-29 baseline case (tau = 0.07, 10 min cruise)
    baseline_idx = find(abs([successful_configs.tau] - 0.07) < 0.001 & [successful_configs.duration_min] == 10, 1);
    if ~isempty(baseline_idx)
        baseline = successful_configs(baseline_idx);
        fprintf('  Structural: %.0f kg\n', baseline.W_structural_kg);
        fprintf('  Engine: %.0f kg\n', baseline.W_engine_kg);
        fprintf('  Systems: %.0f kg\n', baseline.W_systems_kg);
        fprintf('  ----\n');
        fprintf('  OEW: %.0f kg (actual X-29: 6260 kg)\n', baseline.OEW_kg);
        fprintf('  Fuel: %.0f kg\n', baseline.TOGW_kg * baseline.fuel_fraction);
        fprintf('  ----\n');
        fprintf('  TOGW: %.0f kg (actual X-29: 8074 kg)\n', baseline.TOGW_kg);
    end

    fprintf('\nELLIPTICAL FUSELAGE GEOMETRY CHECK:\n');
    fprintf('  Width range: %.2f to %.2f m (mean: %.2f)\n', ...
            min([successful_configs.fuse_width_m]), max([successful_configs.fuse_width_m]), mean([successful_configs.fuse_width_m]));
    fprintf('  Height range: %.2f to %.2f m (mean: %.2f)\n', ...
            min([successful_configs.fuse_height_m]), max([successful_configs.fuse_height_m]), mean([successful_configs.fuse_height_m]));
    fprintf('  H/W ratio range: %.2f to %.2f (mean: %.2f)\n', ...
            min([successful_configs.fuse_aspect_ratio]), max([successful_configs.fuse_aspect_ratio]), mean([successful_configs.fuse_aspect_ratio]));

    % Check for realistic geometry
    unrealistic = ([successful_configs.fuse_width_m] < 0.5) | ([successful_configs.fuse_height_m] < 0.3);
    if any(unrealistic)
        fprintf('  WARNING: %d cases have questionable geometry\n', sum(unrealistic));
    else
        fprintf('  ✓ All geometries are physically realistic!\n');
    end
end

%% 8. GENERATE PLOTS
fprintf('\nGenerating comprehensive plots...\n');

if n_successful > 0
    % Extract data for plotting
    tau_plot = [successful_configs.tau];
    TOGW_plot = [successful_configs.TOGW_kg];
    OEW_plot = [successful_configs.OEW_kg];
    Spln_plot = [successful_configs.S_pln_m2];
    wing_loading_plot = [successful_configs.wing_loading];
    Sref_plot = [successful_configs.wing_area_m2];
    LD_plot = [successful_configs.LD_cruise];
    duration_plot = [successful_configs.duration_min];
    fuse_W = [successful_configs.fuse_width_m];
    fuse_H = [successful_configs.fuse_height_m];
    fuse_aspect = [successful_configs.fuse_aspect_ratio];

    %% FIGURE 1: CARPET PLOT - TOGW vs Spln with cruise duration lines
    figure('Position', [100 100 1200 900]);

    % Plot all points colored by tau
    scatter(Spln_plot, TOGW_plot, 80, tau_plot, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    hold on;

    % Draw lines connecting points of the same cruise duration across tau values
    line_colors = lines(n_duration);
    for i_dur = 1:n_duration
        dur = trade.cruise_duration_min(i_dur);
        idx_dur = (duration_plot == dur);

        % Sort by tau for proper line drawing
        [tau_sorted, sort_idx] = sort(tau_plot(idx_dur));
        Spln_sorted = Spln_plot(idx_dur);
        Spln_sorted = Spln_sorted(sort_idx);
        TOGW_sorted = TOGW_plot(idx_dur);
        TOGW_sorted = TOGW_sorted(sort_idx);

        % Draw line through this cruise duration
        plot(Spln_sorted, TOGW_sorted, '-', 'Color', line_colors(i_dur,:), ...
             'LineWidth', 2, 'DisplayName', sprintf('%d min', dur));
    end

    % Add X-29 actual data point
    plot(38.8, 8074, 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
         'LineWidth', 3, 'DisplayName', 'X-29 Actual');

    % Formatting
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Slenderness (\tau)';
    cb.Label.FontSize = 12;
    caxis([trade.tau_min trade.tau_max]);

    xlabel('Planform Area S_{pln} [m²]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('TOGW [kg]', 'FontSize', 13, 'FontWeight', 'bold');
    title('Carpet Plot: TOGW vs Planform Area - Elliptical Fuselage', 'FontSize', 14, 'FontWeight', 'bold');

    legend('Location', 'northwest', 'FontSize', 10);
    grid on;
    box on;

    saveas(gcf, 'carpet_plot_TOGW_vs_Spln_elliptical.png');

    %% FIGURE 2: CARPET PLOT - TOGW vs Wing Area with cruise duration lines and wing loading constraint
    figure('Position', [100 100 1200 900]);

    % Define axis limits
    xlim_range = [10, 70];

    % Identify feasible and infeasible points based on wing loading constraint
    WL_limit = 550;  % kg/m²
    wing_loading_actual = TOGW_plot ./ Sref_plot;
    feasible_idx = wing_loading_actual <= WL_limit;
    infeasible_idx = wing_loading_actual > WL_limit;

    % Create wing loading constraint line data
    Sref_range = linspace(xlim_range(1), xlim_range(2), 100);
    TOGW_limit = WL_limit * Sref_range;

    % Set y-axis upper limit
    ylim_upper = 15000;  % kg

    % Fill the infeasible region above the constraint line with grey
    fill([Sref_range, fliplr(Sref_range)], [TOGW_limit, ones(size(TOGW_limit))*ylim_upper], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.3, 'DisplayName', 'Infeasible Region (WL > 550)');
    hold on;

    % Plot infeasible points in dark grey
    if any(infeasible_idx)
        scatter(Sref_plot(infeasible_idx), TOGW_plot(infeasible_idx), 80, [0.3 0.3 0.3], 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'DisplayName', 'Infeasible Points');
    end

    % Plot feasible points colored by tau
    scatter(Sref_plot(feasible_idx), TOGW_plot(feasible_idx), 80, tau_plot(feasible_idx), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);

    % Draw lines connecting points of the same cruise duration across tau values
    line_colors = lines(n_duration);
    for i_dur = 1:n_duration
        dur = trade.cruise_duration_min(i_dur);
        idx_dur = (duration_plot == dur);

        % Sort by tau for proper line drawing
        [tau_sorted, sort_idx] = sort(tau_plot(idx_dur));
        Sref_sorted = Sref_plot(idx_dur);
        Sref_sorted = Sref_sorted(sort_idx);
        TOGW_sorted = TOGW_plot(idx_dur);
        TOGW_sorted = TOGW_sorted(sort_idx);

        % Draw line through this cruise duration
        plot(Sref_sorted, TOGW_sorted, '-', 'Color', line_colors(i_dur,:), ...
             'LineWidth', 2, 'DisplayName', sprintf('%d min', dur));
    end

    % Draw wing loading constraint line: TOGW = 550 * Sref
    plot(Sref_range, TOGW_limit, 'r--', 'LineWidth', 2.5, 'DisplayName', 'WL = 550 kg/m² Limit');

    % Add X-29 actual data point
    plot(17.544, 8074, 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
         'LineWidth', 3, 'DisplayName', 'X-29 Actual');

    % Formatting
    xlim(xlim_range);
    ylim([5000, ylim_upper]);
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Slenderness (\tau)';
    cb.Label.FontSize = 12;
    caxis([trade.tau_min trade.tau_max]);

    xlabel('Wing Reference Area S_{ref} [m²]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('TOGW [kg]', 'FontSize', 13, 'FontWeight', 'bold');
    title('Carpet Plot: TOGW vs Wing Area - Elliptical Fuselage', 'FontSize', 14, 'FontWeight', 'bold');

    legend('Location', 'northwest', 'FontSize', 10);
    grid on;
    box on;

    saveas(gcf, 'carpet_plot_TOGW_vs_WingArea_elliptical.png');

    %% FIGURE 3: CARPET PLOT - OEW vs Spln with cruise duration lines
    figure('Position', [100 100 1200 900]);

    % Plot all points colored by tau
    scatter(Spln_plot, OEW_plot, 80, tau_plot, 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    hold on;

    % Draw lines connecting points of the same cruise duration across tau values
    line_colors = lines(n_duration);
    for i_dur = 1:n_duration
        dur = trade.cruise_duration_min(i_dur);
        idx_dur = (duration_plot == dur);

        % Sort by tau for proper line drawing
        [tau_sorted, sort_idx] = sort(tau_plot(idx_dur));
        Spln_sorted = Spln_plot(idx_dur);
        Spln_sorted = Spln_sorted(sort_idx);
        OEW_sorted = OEW_plot(idx_dur);
        OEW_sorted = OEW_sorted(sort_idx);

        % Draw line through this cruise duration
        plot(Spln_sorted, OEW_sorted, '-', 'Color', line_colors(i_dur,:), ...
             'LineWidth', 2, 'DisplayName', sprintf('%d min', dur));
    end

    % Add X-29 actual data point
    plot(38.8, 6260, 'p', 'MarkerSize', 20, 'MarkerFaceColor', 'y', 'MarkerEdgeColor', 'k', ...
         'LineWidth', 3, 'DisplayName', 'X-29 Actual');

    % Formatting
    colormap(jet);
    cb = colorbar;
    cb.Label.String = 'Slenderness (\tau)';
    cb.Label.FontSize = 12;
    caxis([trade.tau_min trade.tau_max]);

    xlabel('Planform Area S_{pln} [m²]', 'FontSize', 13, 'FontWeight', 'bold');
    ylabel('OEW [kg]', 'FontSize', 13, 'FontWeight', 'bold');
    title('Carpet Plot: Operating Empty Weight vs Planform Area - Elliptical Fuselage', 'FontSize', 14, 'FontWeight', 'bold');

    legend('Location', 'northwest', 'FontSize', 10);
    grid on;
    box on;

    saveas(gcf, 'carpet_plot_OEW_vs_Spln_elliptical.png');

    %% FIGURE 4: Fuselage Geometry Evolution
    figure('Position', [100 100 1400 1000]);

    % Focus on 10-minute duration
    idx_10min = (duration_plot == 10);

    subplot(2,2,1);
    scatter(tau_plot(idx_10min), fuse_W(idx_10min), 50, 'b', 'filled');
    xlabel('\tau (Slenderness)');
    ylabel('Fuselage Width [m]');
    title('Fuselage Width vs Tau');
    grid on;
    yline(0.75, 'r--', 'Min Constraint');
    yline(1.30, 'r--', 'Max Constraint');

    subplot(2,2,2);
    scatter(tau_plot(idx_10min), fuse_H(idx_10min), 50, 'g', 'filled');
    xlabel('\tau (Slenderness)');
    ylabel('Fuselage Height [m]');
    title('Fuselage Height vs Tau');
    grid on;

    subplot(2,2,3);
    scatter(tau_plot(idx_10min), fuse_aspect(idx_10min), 50, 'r', 'filled');
    xlabel('\tau (Slenderness)');
    ylabel('H/W Aspect Ratio');
    title('Fuselage Aspect Ratio vs Tau');
    grid on;
    yline(1.0, 'k--', 'Circular');

    subplot(2,2,4);
    scatter(tau_plot(idx_10min), LD_plot(idx_10min), 50, tau_plot(idx_10min), 'filled');
    colorbar;
    xlabel('\tau (Slenderness)');
    ylabel('L/D at Cruise');
    title('Cruise L/D vs Tau');
    grid on;

    sgtitle('Elliptical Fuselage Geometry and Performance (10-min cruise)', 'FontSize', 14, 'FontWeight', 'bold');
    saveas(gcf, 'fuselage_geometry_evolution_elliptical.png');

    fprintf('Plots generated successfully!\n');
end

%% 9. SAVE RESULTS
fprintf('\nSaving results...\n');
timestamp = datestr(now, 'yyyymmdd_HHMMSS');
filename = sprintf('elliptical_fuselage_trade_study_%s.mat', timestamp);
save(filename, 'successful_configs', 'failed_configs', 'trade', 'mission', 'vehicle');
fprintf('Results saved to: %s\n', filename);

fprintf('=================================================================\n');
fprintf('Elliptical fuselage convergence test complete!\n');
fprintf('=================================================================\n');

%% HELPER FUNCTION: Calculate range from duration
function range_km = calculate_range_from_duration(duration_min, M_cruise)
    % Get atmospheric properties at cruise altitude
    h_cruise_m = 12200;
    atm = standard_atmosphere_1976(h_cruise_m);
    V_cruise = M_cruise * atm.a;  % [m/s]
    range_m = V_cruise * duration_min * 60;  % [m]
    range_km = range_m / 1000;  % [km]
end
