%% CHECK IDLE FUEL FLOW DURING APPROACH
% Diagnose why fuel burn is low during landing approach
clear; clc;

addpath(genpath(pwd));
addpath(genpath('../../Propulsion/PythonScriptF404-400'));

fprintf('==============================================\n');
fprintf('  IDLE FUEL FLOW ANALYSIS (PLA = 30)\n');
fprintf('==============================================\n\n');

% Approach conditions
V_approach = 74.8;  % m/s
h_pattern = 610;    % m
gamma_a = 3;        % degrees

% Calculate Mach at approach
[~, a_sl, ~, ~] = atmosisa(0);
[~, a_pattern, ~, ~] = atmosisa(h_pattern);

M_approach = V_approach / a_pattern;

fprintf('Approach Conditions:\n');
fprintf('  Velocity:    %.1f m/s (%.0f kts)\n', V_approach, V_approach * 1.94384);
fprintf('  Altitude:    %.0f m (%.0f ft)\n', h_pattern, h_pattern * 3.28084);
fprintf('  Mach:        %.3f\n', M_approach);
fprintf('  Glideslope:  %.0f degrees\n', gamma_a);

%% Test physics model at idle
fprintf('\n--- Physics Model at Idle (PLA=30) ---\n');

try
    prop_idle = get_F404_thrust_physics(M_approach, h_pattern, 30);
    fprintf('  Thrust:      %.0f N (%.0f lbf)\n', prop_idle.thrust_N, prop_idle.thrust_lbf);
    fprintf('  TSFC:        %.4f kg/N/s (%.3f lb/hr/lbf)\n', prop_idle.TSFC_SI, prop_idle.TSFC_imperial);
    fprintf('  Fuel Flow:   %.4f kg/s (%.1f kg/hr)\n', prop_idle.fuel_flow_kg_s, prop_idle.fuel_flow_kg_s * 3600);
catch ME
    fprintf('  ERROR: %s\n', ME.message);
end

%% Test at different PLA values
fprintf('\n--- Fuel Flow vs PLA at Approach Conditions ---\n');
fprintf('%6s | %10s | %12s | %12s\n', 'PLA', 'Thrust (N)', 'TSFC', 'Fuel (kg/s)');
fprintf('%s\n', repmat('-', 1, 50));

for PLA = [30, 40, 50, 60, 70, 80, 87]
    prop = get_F404_thrust_physics(M_approach, h_pattern, PLA);
    fprintf('%6d | %10.0f | %12.4f | %12.4f\n', PLA, prop.thrust_N, prop.TSFC_SI, prop.fuel_flow_kg_s);
end

%% Calculate expected fuel burn
fprintf('\n--- Expected Fuel Burn During Approach ---\n');

% Approach time (from pattern to flare)
h_flare = 3.6;  % m
t_approach = (h_pattern - h_flare) / (V_approach * sind(gamma_a));

fprintf('  Approach distance: %.0f m\n', (h_pattern - h_flare) / sind(gamma_a));
fprintf('  Approach time:     %.1f s (%.2f min)\n', t_approach, t_approach/60);

% At idle
prop_idle = get_F404_thrust_physics(M_approach, h_pattern, 30);
fuel_idle = prop_idle.fuel_flow_kg_s * t_approach;
fprintf('\n  At IDLE (PLA=30):\n');
fprintf('    Fuel flow:  %.4f kg/s\n', prop_idle.fuel_flow_kg_s);
fprintf('    Fuel burn:  %.1f kg (%.1f lb)\n', fuel_idle, fuel_idle * 2.205);

% At typical approach power (maybe PLA=50 for some drag)
prop_50 = get_F404_thrust_physics(M_approach, h_pattern, 50);
fuel_50 = prop_50.fuel_flow_kg_s * t_approach;
fprintf('\n  At PLA=50:\n');
fprintf('    Fuel flow:  %.4f kg/s\n', prop_50.fuel_flow_kg_s);
fprintf('    Fuel burn:  %.1f kg (%.1f lb)\n', fuel_50, fuel_50 * 2.205);

%% What does NASA say about idle fuel flow?
fprintf('\n==============================================\n');
fprintf('  REFERENCE: F404 IDLE FUEL FLOW\n');
fprintf('==============================================\n');
fprintf('  Typical turbofan idle fuel flow: 500-800 lb/hr\n');
fprintf('  That is: 0.063 - 0.101 kg/s\n');
fprintf('  Over %.0f s approach: %.0f - %.0f kg\n', t_approach, 0.063*t_approach, 0.101*t_approach);
fprintf('==============================================\n');
