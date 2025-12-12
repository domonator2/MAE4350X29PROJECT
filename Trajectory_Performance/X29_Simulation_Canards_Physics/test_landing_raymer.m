%% TEST LANDING SIMULATION - RAYMER METHODOLOGY
% Verifies landing simulation against report Table 31 values
clear; clc; close all;

% Add paths
addpath(genpath(pwd));
addpath(genpath('../../Propulsion/PythonScriptF404-400'));
addpath(genpath('../../Stability'));
addpath(genpath('../X29_Simulation'));  % For calculate_X29_CG

fprintf('=====================================================\n');
fprintf('  X-29 LANDING SIMULATION TEST (Raymer Methodology)\n');
fprintf('=====================================================\n\n');

%% Test Case: Match report assumptions (W = 6800 kg)
params = struct();
params.W_total = 6800;          % Landing weight [kg] (~15,000 lb)
params.W_fuel = 500;            % Remaining fuel [kg]
params.h_pattern = 610;         % Pattern altitude [m] (2000 ft)
params.CL_max = 1.6;            % With flaperons
params.verbose = true;

% Run landing simulation
landing = sim_landing_canards(params);

%% Compare with Report Table 31
fprintf('\n\n');
fprintf('=====================================================\n');
fprintf('  COMPARISON WITH REPORT TABLE 31\n');
fprintf('=====================================================\n');
fprintf('%25s | %12s | %12s | %10s\n', 'Parameter', 'Report', 'Simulation', 'Error (%)');
fprintf('%s\n', repmat('-', 1, 70));

% Report values (from Table 31)
report.V_stall_ms = 62.4;
report.V_approach_ms = 74.9;
report.V_touchdown_ms = 68.6;
report.R_flare_m = 2627;
report.h_flare_m = 3.6;
report.S_approach_m = 222;
report.S_flare_m = 137;
report.S_free_roll_m = 137;
report.S_brake_m = 800;
report.S_ground_roll_m = 937;
report.S_total_m = 1296;

% Compare
compare_val('Stall Speed (m/s)', report.V_stall_ms, landing.V_stall);
compare_val('Approach Speed (m/s)', report.V_approach_ms, landing.V_approach);
compare_val('Touchdown Speed (m/s)', report.V_touchdown_ms, landing.V_touchdown);
compare_val('Flare Radius (m)', report.R_flare_m, landing.R_flare);
compare_val('Flare Height (m)', report.h_flare_m, landing.h_flare_init);
compare_val('Approach Dist (m)', report.S_approach_m, landing.S_approach_from_obstacle);
compare_val('Flare Dist (m)', report.S_flare_m, landing.S_flare);
compare_val('Free Roll Dist (m)', report.S_free_roll_m, landing.S_free_roll);
compare_val('Braking Dist (m)', report.S_brake_m, landing.S_brake);
compare_val('Ground Roll (m)', report.S_ground_roll_m, landing.S_ground_roll);
compare_val('Total Landing (m)', report.S_total_m, landing.S_total_from_obstacle);

fprintf('%s\n', repmat('=', 1, 70));

%% Generate plots
figure('Position', [50, 50, 1400, 900], 'Color', 'w');
sgtitle('X-29 Landing Trajectory (Raymer Methodology)', 'FontSize', 14, 'FontWeight', 'bold');

% Plot 1: Altitude vs Distance
subplot(2,3,1);
plot(landing.x/1000, landing.h, 'b-', 'LineWidth', 2);
hold on;
xline(landing.S_approach_from_pattern/1000, 'r--', 'Flare', 'LabelVerticalAlignment', 'bottom');
xlabel('Horizontal Distance [km]');
ylabel('Altitude [m]');
title('Approach Profile');
grid on;

% Plot 2: Velocity vs Time
subplot(2,3,2);
plot(landing.t, landing.V * 1.94384, 'b-', 'LineWidth', 2);
hold on;
yline(landing.V_approach * 1.94384, 'g--', 'V_{approach}');
yline(landing.V_touchdown * 1.94384, 'r--', 'V_{TD}');
yline(landing.V_stall * 1.94384, 'k:', 'V_{stall}');
xlabel('Time [s]');
ylabel('Velocity [kts]');
title('Velocity Profile');
legend('Velocity', 'Location', 'northeast');
grid on;

% Plot 3: Ground Roll Detail
subplot(2,3,3);
idx_ground = landing.h <= 0.1;
t_ground = landing.t(idx_ground) - landing.t(find(idx_ground, 1));
V_ground = landing.V(idx_ground);
plot(t_ground, V_ground * 1.94384, 'b-', 'LineWidth', 2);
hold on;
xline(2, 'r--', 'Brakes Applied');
xlabel('Time after Touchdown [s]');
ylabel('Velocity [kts]');
title('Ground Roll Deceleration');
grid on;

% Plot 4: Canard Deflection
subplot(2,3,4);
plot(landing.t, landing.delta_canard, 'b-', 'LineWidth', 2);
xlabel('Time [s]');
ylabel('Canard Deflection [deg]');
title('Canard Schedule During Landing');
grid on;

% Plot 5: Weight and Fuel
subplot(2,3,5);
yyaxis left;
plot(landing.t, landing.W, 'b-', 'LineWidth', 2);
ylabel('Total Weight [kg]');
yyaxis right;
plot(landing.t, landing.W_fuel, 'r-', 'LineWidth', 2);
ylabel('Fuel Weight [kg]');
xlabel('Time [s]');
title('Weight During Landing');
legend('Total', 'Fuel', 'Location', 'east');
grid on;

% Plot 6: Landing Distance Breakdown (Bar chart)
subplot(2,3,6);
distances = [landing.S_approach_from_obstacle, landing.S_flare, ...
             landing.S_free_roll, landing.S_brake];
bar(distances/0.3048, 'FaceColor', [0.2 0.6 0.8]);
set(gca, 'XTickLabel', {'Approach', 'Flare', 'Free Roll', 'Braking'});
ylabel('Distance [ft]');
title('Landing Distance Breakdown');
grid on;

% Add total distance annotation
text(2.5, max(distances)/0.3048 * 1.1, ...
    sprintf('Total: %.0f ft', landing.S_total_from_obstacle/0.3048), ...
    'FontSize', 11, 'FontWeight', 'bold');

saveas(gcf, 'X29_landing_raymer_analysis.png');
fprintf('\nPlot saved: X29_landing_raymer_analysis.png\n');

%% Helper function
function compare_val(name, report_val, sim_val)
    err_pct = (sim_val - report_val) / report_val * 100;
    fprintf('%25s | %12.1f | %12.1f | %+9.1f%%\n', name, report_val, sim_val, err_pct);
end
