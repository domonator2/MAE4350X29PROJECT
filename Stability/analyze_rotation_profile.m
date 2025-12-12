function analyze_rotation_profile()
%ANALYZE_ROTATION_PROFILE Show canard behavior through takeoff phases
%
% Compares canard deflection requirements across takeoff phases:
%   1. Ground roll - moments about main gear
%   2. Rotation/Liftoff - transition
%   3. Climb - accelerating, trimmed at higher V
%   4. Cruise - steady level flight
%
% Author: Dominic Larin
% Date: December 2025

%% Constants
g = 9.81;
rho = 1.225;        % Sea level
ft2m = 0.3048;

%% Aircraft Parameters
S_ref = 17.54;              % Wing reference area [m²]
S_canard = 3.4923;          % Canard area [m²]
c_bar = 8.0 * ft2m;         % MAC [m]

% Positions from nose [m]
x_wing_ac = 40.91 * ft2m;   % 12.47 m
x_canard_ac = 6.49;         % m
x_main_gear = 29.50 * ft2m; % 8.99 m

% Weight
W_total = 8000;             % [kg]
W_fuel = 1877;              % Full fuel [kg]
W = W_total * g;            % [N]

% CG
[x_cg, ~] = calculate_X29_CG(W_total, W_fuel, 1877);

%% Aerodynamic Parameters
CL_alpha_wing = 0.08;       % [1/deg]
CL_0_wing = 0.05;           % Zero-alpha lift
CL_alpha_canard = 0.0263;   % [1/deg]

fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('  X-29 TAKEOFF CANARD SCHEDULE\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('\n');
fprintf('  CG = %.2f m from nose\n', x_cg);
fprintf('  Main gear = %.2f m from nose\n', x_main_gear);
fprintf('\n');

%% Phase 1: Ground Roll - Moment about main gear for rotation
fprintf('PHASE 1: GROUND ROLL (Rotation initiation)\n');
fprintf('─────────────────────────────────────────\n');

V_rotate = 90;  % m/s
q_rotate = 0.5 * rho * V_rotate^2;
alpha_ground = 2;  % deg

% Moment arms about MAIN GEAR
arm_cg_mg = x_main_gear - x_cg;              % Positive if CG ahead of MG
arm_canard_mg = x_main_gear - x_canard_ac;   % Positive (canard far ahead)
arm_wing_mg = x_main_gear - x_wing_ac;       % Negative (wing behind MG)

fprintf('  Moment arms about main gear:\n');
fprintf('    CG:      %+.2f m\n', arm_cg_mg);
fprintf('    Canard:  %+.2f m\n', arm_canard_mg);
fprintf('    Wing:    %+.2f m\n', arm_wing_mg);

% Weight creates nose-down moment (CG ahead of main gear)
M_weight_mg = -W * arm_cg_mg;

% Wing lift at ground roll alpha
CL_wing_ground = CL_0_wing + CL_alpha_wing * alpha_ground;
L_wing_ground = q_rotate * S_ref * CL_wing_ground;
M_wing_mg = L_wing_ground * arm_wing_mg;

% Solve for canard δ needed to get M_net = 0 (just starting rotation)
% M_canard + M_aero = -(M_weight + M_wing)
M_to_overcome = -(M_weight_mg + M_wing_mg);

% Get Cm derivatives
a_sl = 340.3;
M_mach = V_rotate / a_sl;
[Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M_mach);

% Effectiveness: dM/dδ
alpha_rad = alpha_ground * pi/180;
M_aero_alpha = q_rotate * S_ref * c_bar * Cm_alpha * alpha_rad;
M_to_overcome = M_to_overcome - M_aero_alpha;

canard_moment_effectiveness = q_rotate * S_canard * CL_alpha_canard * arm_canard_mg;
aero_moment_effectiveness = q_rotate * S_ref * c_bar * Cm_delta_c * (pi/180);
total_effectiveness = canard_moment_effectiveness + aero_moment_effectiveness;

delta_rotate = M_to_overcome / total_effectiveness;

fprintf('\n  At V = %.0f m/s, α = %.0f°:\n', V_rotate, alpha_ground);
fprintf('    M_weight (about MG): %+.1f kN·m (nose down)\n', M_weight_mg/1000);
fprintf('    M_wing (about MG):   %+.1f kN·m\n', M_wing_mg/1000);
fprintf('    M_aero(α):           %+.1f kN·m (instability helps!)\n', M_aero_alpha/1000);
fprintf('    → δ_canard for rotation: %.1f°\n', delta_rotate);

%% Phase 2: Trimmed flight at various speeds
fprintf('\n');
fprintf('PHASE 2-4: TRIMMED FLIGHT (L=W, M=0)\n');
fprintf('─────────────────────────────────────────\n');

% Call calculate_trim at various speeds
speeds = [90, 100, 120, 150, 200];
fprintf('\n');
fprintf('  Speed      α_trim   δ_canard   CL_total   L/D\n');
fprintf('  [m/s]       [deg]     [deg]       [-]      [-]\n');
fprintf('  ─────────────────────────────────────────────────\n');

trim_data = struct();
for i = 1:length(speeds)
    V = speeds(i);
    % Suppress output from calculate_trim
    evalc('trim = calculate_trim(V, 0, W_total, W_fuel);');

    trim_data(i).V = V;
    trim_data(i).alpha = trim.alpha_deg;
    trim_data(i).delta = trim.delta_canard_deg;
    trim_data(i).CL = trim.CL_total;
    trim_data(i).LD = trim.LD;
    trim_data(i).valid = trim.valid;

    valid_str = '';
    if ~trim.valid
        valid_str = ' *';
    end

    fprintf('  %5.0f      %6.2f    %6.1f     %6.3f   %5.1f%s\n', ...
        V, trim.alpha_deg, trim.delta_canard_deg, trim.CL_total, trim.LD, valid_str);
end
fprintf('  ─────────────────────────────────────────────────\n');
fprintf('  * = outside limits (α>15° or δ>58°)\n');

%% Comparison Summary
fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('  CANARD DEFLECTION SCHEDULE SUMMARY\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('\n');
fprintf('  TAKEOFF SEQUENCE:\n');
fprintf('  ─────────────────\n');
fprintf('  1. Ground roll start:     δ ≈ 0° (building speed)\n');
fprintf('  2. Rotation (V=90 m/s):   δ ≈ %.0f° (initiating pitch-up)\n', delta_rotate);
fprintf('  3. Liftoff (V=90 m/s):    δ ≈ %.0f° (trimmed flight)\n', trim_data(1).delta);
fprintf('  4. Climb (V=120 m/s):     δ ≈ %.0f° (accelerating)\n', trim_data(3).delta);
fprintf('  5. Cruise (V=200 m/s):    δ ≈ %.0f° (steady flight)\n', trim_data(5).delta);
fprintf('\n');

% The key insight
fprintf('  KEY FINDING:\n');
fprintf('  ─────────────\n');
if trim_data(1).delta < delta_rotate
    fprintf('  Canard deflection DECREASES from %.0f° to %.0f° after liftoff!\n', ...
        delta_rotate, trim_data(1).delta);
    fprintf('  Why? On ground, canard must overcome weight moment about main gear.\n');
    fprintf('        In flight, canard only needs to trim (L=W, M=0).\n');
else
    fprintf('  Canard deflection stays similar (%.0f° → %.0f°)\n', ...
        delta_rotate, trim_data(1).delta);
    fprintf('  Both require significant canard authority at low speed.\n');
end
fprintf('\n');
fprintf('  As speed increases, canard deflection DECREASES:\n');
fprintf('  • Higher q → less α needed for lift\n');
fprintf('  • Less α → less nose-down wing moment\n');
fprintf('  • Less canard deflection needed for trim\n');
fprintf('\n');
fprintf('  At 200 m/s: δ = %.1f° (vs %.1f° at 90 m/s)\n', ...
    trim_data(5).delta, trim_data(1).delta);
fprintf('═══════════════════════════════════════════════════════════════════\n');

%% Plot
figure('Name', 'X-29 Canard Schedule', 'Position', [100 100 800 500]);

speeds_plot = [speeds(1), speeds];
delta_plot = [delta_rotate, [trim_data.delta]];
labels = {'Rotation', 'Liftoff', '100 m/s', '120 m/s', '150 m/s', '200 m/s'};

bar(delta_plot);
hold on;
yline(58, 'r--', 'LineWidth', 2);
text(1, 58+2, 'Max TEU (58°)', 'Color', 'r', 'FontSize', 10);

set(gca, 'XTickLabel', labels);
ylabel('Canard Deflection [deg]');
xlabel('Flight Phase');
title('X-29 Canard Deflection Through Takeoff');
grid on;

% Add value labels on bars
for i = 1:length(delta_plot)
    text(i, delta_plot(i)+1, sprintf('%.0f°', delta_plot(i)), ...
        'HorizontalAlignment', 'center', 'FontWeight', 'bold');
end

saveas(gcf, 'X29_canard_schedule.png');
fprintf('\nFigure saved as X29_canard_schedule.png\n');

end

%% Helper function
function out = ternary(cond, true_val, false_val)
    if cond
        out = true_val;
    else
        out = false_val;
    end
end
