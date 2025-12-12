%% PLOT MISSION ALTITUDE VS DISTANCE
% Creates altitude vs horizontal distance plot for entire X-29 mission
clear; clc; close all;

% Add paths
addpath(genpath(pwd));
addpath(genpath('../../Propulsion/PythonScriptF404-400'));
addpath(genpath('../../Stability'));
addpath(genpath('../X29_Simulation'));

fprintf('==============================================\n');
fprintf('  X-29 FULL MISSION - ALTITUDE VS DISTANCE\n');
fprintf('==============================================\n\n');

%% Mission Parameters
% Initial conditions
W_TO = 8074;            % Takeoff weight [kg] (~17,800 lb)
W_fuel_initial = 2182;  % Initial fuel [kg]

% Mission profile
h_cruise = 12200;       % Cruise altitude [m] (40,000 ft)
M_cruise = 1.4;         % Cruise Mach
t_cruise = 5 * 60;      % Cruise duration [s] (5 minutes)
h_pattern = 610;        % Pattern altitude [m] (2000 ft)

%% Phase 1: Takeoff
fprintf('Phase 1: Takeoff...\n');
takeoff_params = struct();
takeoff_params.W_total = W_TO;
takeoff_params.W_fuel = W_fuel_initial;
takeoff_params.verbose = false;

takeoff = sim_takeoff_canards(takeoff_params);

% Takeoff distance (ground roll + rotation + initial climb)
dt_takeoff = diff(takeoff.t);
V_avg_takeoff = (takeoff.V(1:end-1) + takeoff.V(2:end)) / 2;
x_takeoff = [0; cumsum(dt_takeoff .* V_avg_takeoff)];
h_takeoff = takeoff.h;
W_after_takeoff = takeoff.W(end);
W_fuel_after_takeoff = takeoff.W_fuel(end);

fprintf('  Ground roll + climb to %.0f m\n', takeoff.h(end));
fprintf('  Distance: %.1f km\n', x_takeoff(end)/1000);

%% Phase 2: Climb to cruise altitude
fprintf('Phase 2: Climb to %.0f m...\n', h_cruise);
climb_params = struct();
climb_params.W_total = W_after_takeoff;
climb_params.W_fuel = W_fuel_after_takeoff;
climb_params.h_initial = takeoff.h(end);
climb_params.h_target = h_cruise;
climb_params.M_climb = 0.9;
climb_params.verbose = false;

climb = sim_climb_canards(climb_params);

% Calculate horizontal distance during climb
dt_climb = diff(climb.t);
dh_climb = diff(climb.h);
V_avg_climb = (climb.V(1:end-1) + climb.V(2:end)) / 2;

% Flight path angle
gamma_climb = asind(dh_climb ./ (dt_climb .* V_avg_climb + 1e-6));
gamma_climb(isnan(gamma_climb) | abs(gamma_climb) > 45) = 10;

% Horizontal distance
dx_climb = dt_climb .* V_avg_climb .* cosd(gamma_climb);
x_climb = [0; cumsum(dx_climb)];

h_climb = climb.h;
W_after_climb = climb.W(end);
W_fuel_after_climb = climb.W_fuel(end);

fprintf('  Climbed to %.0f m\n', climb.h(end));
fprintf('  Distance: %.1f km\n', x_climb(end)/1000);

%% Phase 3: Cruise
fprintf('Phase 3: Cruise at M%.1f for %.0f min...\n', M_cruise, t_cruise/60);
cruise_params = struct();
cruise_params.W_initial = W_after_climb;
cruise_params.W_fuel_initial = W_fuel_after_climb;
cruise_params.h = h_cruise;
cruise_params.M = M_cruise;
cruise_params.duration_min = t_cruise / 60;  % Convert to minutes
cruise_params.verbose = false;

cruise = sim_cruise_canards(cruise_params);

% Calculate cruise distance
[~, a_cruise, ~, ~] = atmosisa(h_cruise);
V_cruise = M_cruise * a_cruise;
x_cruise = cruise.t * V_cruise;  % Constant velocity cruise

h_cruise_arr = cruise.h;
W_after_cruise = cruise.W(end);
W_fuel_after_cruise = cruise.W_fuel(end);

fprintf('  Cruise distance: %.1f km\n', x_cruise(end)/1000);

%% Phase 4: Descent
fprintf('Phase 4: Descent to pattern altitude...\n');
descent_params = struct();
descent_params.W_total = W_after_cruise;
descent_params.W_fuel = W_fuel_after_cruise;
descent_params.h_initial = h_cruise;
descent_params.h_final = h_pattern;
descent_params.verbose = false;

descent = sim_descent_canards(descent_params);

% Calculate descent distance
dt_descent = diff(descent.t);
dh_descent = diff(descent.h);
V_avg_descent = (descent.V(1:end-1) + descent.V(2:end)) / 2;

% Flight path angle
gamma_descent = asind(dh_descent ./ (dt_descent .* V_avg_descent + 1e-6));
gamma_descent(isnan(gamma_descent) | abs(gamma_descent) > 30) = -3;

% Horizontal distance
dx_descent = dt_descent .* V_avg_descent .* cosd(gamma_descent);
x_descent = [0; cumsum(dx_descent)];

h_descent = descent.h;
W_after_descent = descent.W(end);
W_fuel_after_descent = descent.W_fuel(end);

fprintf('  Descended to %.0f m\n', descent.h(end));
fprintf('  Distance: %.1f km\n', x_descent(end)/1000);

%% Phase 5: Landing
fprintf('Phase 5: Landing...\n');
landing_params = struct();
landing_params.W_total = W_after_descent;
landing_params.W_fuel = W_fuel_after_descent;
landing_params.h_pattern = h_pattern;
landing_params.verbose = false;

landing = sim_landing_canards(landing_params);

x_landing = landing.x;
h_landing = landing.h;

fprintf('  Landing distance: %.1f km\n', x_landing(end)/1000);

%% Combine all phases with cumulative distance
% Starting distances for each phase
x0_takeoff = 0;
x0_climb = x_takeoff(end);
x0_cruise = x0_climb + x_climb(end);
x0_descent = x0_cruise + x_cruise(end);
x0_landing = x0_descent + x_descent(end);

% Cumulative distances
x_total_takeoff = x0_takeoff + x_takeoff;
x_total_climb = x0_climb + x_climb;
x_total_cruise = x0_cruise + x_cruise;
x_total_descent = x0_descent + x_descent;
x_total_landing = x0_landing + x_landing;

% Total mission distance
total_distance_m = x0_landing + x_landing(end);
total_distance_km = total_distance_m / 1000;
total_distance_nm = total_distance_km / 1.852;

%% Print Summary
fprintf('\n==============================================\n');
fprintf('  MISSION DISTANCE SUMMARY\n');
fprintf('==============================================\n');
fprintf('  Phase          Distance (km)   Distance (nm)\n');
fprintf('  ------------------------------------------------\n');
fprintf('  Takeoff:       %8.1f        %8.1f\n', x_takeoff(end)/1000, x_takeoff(end)/1852);
fprintf('  Climb:         %8.1f        %8.1f\n', x_climb(end)/1000, x_climb(end)/1852);
fprintf('  Cruise:        %8.1f        %8.1f\n', x_cruise(end)/1000, x_cruise(end)/1852);
fprintf('  Descent:       %8.1f        %8.1f\n', x_descent(end)/1000, x_descent(end)/1852);
fprintf('  Landing:       %8.1f        %8.1f\n', x_landing(end)/1000, x_landing(end)/1852);
fprintf('  ------------------------------------------------\n');
fprintf('  TOTAL:         %8.1f km     %8.1f nm\n', total_distance_km, total_distance_nm);
fprintf('==============================================\n');

%% Create Altitude vs Distance Plot
figure('Position', [100, 100, 1400, 600], 'Color', 'w');

% Main plot
subplot(1,1,1);
hold on;

% Plot each phase with different colors
p1 = plot(x_total_takeoff/1000, h_takeoff/1000, 'b-', 'LineWidth', 2.5);
p2 = plot(x_total_climb/1000, h_climb/1000, 'g-', 'LineWidth', 2.5);
p3 = plot(x_total_cruise/1000, h_cruise_arr/1000, 'r-', 'LineWidth', 2.5);
p4 = plot(x_total_descent/1000, h_descent/1000, 'm-', 'LineWidth', 2.5);
p5 = plot(x_total_landing/1000, h_landing/1000, 'c-', 'LineWidth', 2.5);

% Mark phase transitions
scatter([x0_climb, x0_cruise, x0_descent, x0_landing]/1000, ...
        [h_takeoff(end), h_climb(end), h_cruise_arr(end), h_descent(end)]/1000, ...
        100, 'ko', 'filled');

% Labels
xlabel('Horizontal Distance [km]', 'FontSize', 12);
ylabel('Altitude [km]', 'FontSize', 12);
title(sprintf('X-29 Mission Profile - Total Distance: %.0f km (%.0f nm)', ...
    total_distance_km, total_distance_nm), 'FontSize', 14, 'FontWeight', 'bold');

legend([p1, p2, p3, p4, p5], {'Takeoff', 'Climb', 'Cruise (M1.4)', 'Descent', 'Landing'}, ...
    'Location', 'northeast', 'FontSize', 11);

grid on;
xlim([0, total_distance_km * 1.02]);
ylim([0, max(h_cruise_arr)/1000 * 1.1]);

% Add distance annotations
text(x0_climb/1000, -0.5, sprintf('%.0f km', x_takeoff(end)/1000), ...
    'HorizontalAlignment', 'center', 'FontSize', 9);
text((x0_climb + x0_cruise)/2/1000, -0.5, sprintf('%.0f km', x_climb(end)/1000), ...
    'HorizontalAlignment', 'center', 'FontSize', 9);
text((x0_cruise + x0_descent)/2/1000, h_cruise/1000 + 0.5, sprintf('Cruise: %.0f km', x_cruise(end)/1000), ...
    'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
text((x0_descent + x0_landing)/2/1000, -0.5, sprintf('%.0f km', x_descent(end)/1000), ...
    'HorizontalAlignment', 'center', 'FontSize', 9);

% Save plot
saveas(gcf, 'X29_mission_altitude_vs_distance.png');
fprintf('\nPlot saved: X29_mission_altitude_vs_distance.png\n');

%% Also create a table with key waypoints
fprintf('\n==============================================\n');
fprintf('  MISSION WAYPOINTS\n');
fprintf('==============================================\n');
fprintf('  Point              Distance(km)  Altitude(m)\n');
fprintf('  ------------------------------------------------\n');
fprintf('  Takeoff Start      %8.1f      %8.0f\n', 0, 0);
fprintf('  End Takeoff        %8.1f      %8.0f\n', x0_climb/1000, h_takeoff(end));
fprintf('  Top of Climb       %8.1f      %8.0f\n', x0_cruise/1000, h_climb(end));
fprintf('  End Cruise         %8.1f      %8.0f\n', x0_descent/1000, h_cruise_arr(end));
fprintf('  End Descent        %8.1f      %8.0f\n', x0_landing/1000, h_descent(end));
fprintf('  Touchdown          %8.1f      %8.0f\n', total_distance_km, 0);
fprintf('==============================================\n');
