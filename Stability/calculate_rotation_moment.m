function [results] = calculate_rotation_moment(V, W_total, W_fuel, delta_canard, alpha_wing)
%CALCULATE_ROTATION_MOMENT Calculate moments about main gear for takeoff rotation
%
% Performs moment balance about main gear pivot to determine if canard
% can rotate the aircraft nose-up during takeoff ground roll.
%
% Sign convention: Positive moment = nose UP
%
% Inputs:
%   V            - Airspeed [m/s]
%   W_total      - Total aircraft weight [kg]
%   W_fuel       - Current fuel weight [kg] (for CG calculation)
%   delta_canard - Canard deflection [deg] (positive = TEU = nose up)
%   alpha_wing   - Wing angle of attack [deg] (ground roll ~ 2 deg)
%
% Outputs:
%   results - Structure with moment breakdown and rotation feasibility
%
% Author: Dominic Larin
% Date: December 2025

%% Default inputs for testing
if nargin < 1, V = 70; end              % Typical rotation speed [m/s]
if nargin < 2, W_total = 8000; end      % TOGW [kg]
if nargin < 3, W_fuel = 1877; end       % Full fuel [kg]
if nargin < 4, delta_canard = 20; end   % Canard deflection [deg]
if nargin < 5, alpha_wing = 2; end      % Ground roll alpha [deg]

%% Constants
g = 9.81;           % Gravity [m/s²]
rho = 1.225;        % Sea level density [kg/m³]
ft2m = 0.3048;

%% Aircraft Geometry (from calculate_X29_CG.m - Updated Dec 2025)
% Positions from nose tip [m]
x_nose_gear = 11.50 * ft2m;   % 3.51 m (CORRECTED)
x_main_gear = 29.50 * ft2m;   % 8.99 m (CORRECTED)
x_canard_ac = 6.49;           % Canard AC (from canard_properties)
x_wing_ac = 40.91 * ft2m;     % 12.47 m - Wing AC (~quarter chord of wing)

% Reference areas and dimensions
S_ref = 17.54;          % Wing reference area [m²]
S_canard = 3.4923;      % Canard area [m²]
c_bar = 8.0 * ft2m;     % Mean aerodynamic chord [m] (2.44 m)

%% Aerodynamic Parameters
% Canard (from canard_properties.m)
CL_alpha_canard = 0.0263;   % [1/deg]
canard_delta_max_up = 58;   % [deg]
canard_delta_max_down = 32; % [deg]

% Wing (ground roll, with ground effect)
CL_alpha_wing = 0.08;       % [1/deg] - approximate
CL_0_wing = 0.1;            % Zero-alpha lift (camber + ground effect)

% Get Cm derivatives at current Mach
a_sl = 340.3;  % Speed of sound at sea level [m/s]
M = V / a_sl;
[Cm_alpha, Cm_delta_c, Cm_delta_s, Cm_delta_f] = get_Cm_derivatives(M);

%% Calculate CG Position
% Use the CG calculation function
[x_cg, cg_data] = calculate_X29_CG(W_total, W_fuel, 1877);

%% Dynamic Pressure
q = 0.5 * rho * V^2;

%% Forces

% Weight
W = W_total * g;    % [N]

% Canard lift
CL_canard = CL_alpha_canard * delta_canard;
% Check for stall (rough limit ~1.2)
if abs(CL_canard) > 1.2
    CL_canard = sign(CL_canard) * 1.2;
    canard_stalled = true;
else
    canard_stalled = false;
end
L_canard = q * S_canard * CL_canard;

% Wing lift (small during ground roll due to low alpha)
CL_wing = CL_0_wing + CL_alpha_wing * alpha_wing;
L_wing = q * S_ref * CL_wing;

% Total lift
L_total = L_canard + L_wing;

% Ground reaction (weight minus lift, distributed between gears)
R_total = W - L_total;
if R_total < 0
    R_total = 0;  % Airborne
end

%% Moment Arms (about main gear pivot)
% Positive = ahead of main gear (creates nose-down moment when pushed down)
arm_cg = x_main_gear - x_cg;           % CG arm (negative if CG ahead of main gear)
arm_canard = x_main_gear - x_canard_ac; % Canard arm (positive, canard ahead)
arm_wing = x_main_gear - x_wing_ac;     % Wing arm (negative, wing behind)
arm_nose_gear = x_main_gear - x_nose_gear; % Nose gear arm (positive)

%% Moment Balance About Main Gear
% M = sum of (Force × Arm)
% Positive moment = nose UP

% Weight moment (acts at CG)
% If CG is ahead of main gear: arm_cg > 0, weight pushes nose down → negative moment
M_weight = -W * arm_cg;  % Nose down if CG ahead of main gear

% Canard lift moment
% Canard ahead of main gear: arm_canard > 0, lift up → nose UP → positive
M_canard = L_canard * arm_canard;

% Wing lift moment
% Wing behind main gear: arm_wing < 0, lift up → nose DOWN → but sign works out
M_wing = L_wing * arm_wing;

% Aerodynamic pitching moment (about CG, then transfer to main gear)
% Cm = Cm_alpha * alpha + Cm_delta_c * delta_canard
% Note: alpha and delta in radians for Cm calculation
alpha_rad = alpha_wing * pi/180;
delta_rad = delta_canard * pi/180;
Cm_total = Cm_alpha * alpha_rad + Cm_delta_c * delta_rad;
M_aero_cg = q * S_ref * c_bar * Cm_total;  % Moment about CG

% Transfer to main gear pivot: M_aero already acts about CG, no transfer needed
% This is the "free" pitching moment from aircraft instability
M_aero = M_aero_cg;  % Positive = nose up

% Ground reaction at nose gear (before rotation)
% Estimate nose gear load from static equilibrium
% Sum moments about main gear = 0 for static case
% R_nose * arm_nose_gear = W * arm_cg (approximately, ignoring aero)
if arm_nose_gear > 0
    R_nose_static = W * arm_cg / arm_nose_gear;
    R_nose_static = max(0, R_nose_static);  % Can't be negative (tension)
else
    R_nose_static = 0;
end

% During rotation attempt, nose gear reaction reduces as nose lifts
% At rotation point, we want R_nose → 0
% The moment from nose gear reaction (resisting rotation):
M_nose_gear = 0;  % At rotation initiation, we're trying to make this zero

%% Net Moment
M_net = M_canard + M_wing + M_weight + M_aero;

% For rotation: M_net must be positive (nose up)
can_rotate = M_net > 0;

%% Calculate Required Canard Deflection for Rotation
% Set M_net = 0 and solve for delta_canard
% M_canard + M_aero(delta) = -(M_wing + M_weight + M_aero(alpha))
% This is iterative because M_aero depends on delta, but we can linearize:
% M_canard(delta) = L_canard(delta) * arm_canard = q*S_c*CL_alpha*delta * arm_canard
% M_aero(delta) = q*S*c_bar*Cm_delta_c*delta
% Combined: delta * (q*S_c*CL_alpha*arm_canard + q*S*c_bar*Cm_delta_c) = -(M_wing + M_weight + M_aero_alpha)
M_aero_alpha_only = q * S_ref * c_bar * Cm_alpha * alpha_rad;
M_to_overcome = -(M_wing + M_weight + M_aero_alpha_only);
canard_effectiveness = q * S_canard * CL_alpha_canard * (pi/180) * arm_canard + q * S_ref * c_bar * Cm_delta_c * (pi/180);
delta_required = M_to_overcome / canard_effectiveness;  % [deg]

% Check if within limits
delta_feasible = (delta_required <= canard_delta_max_up) && (delta_required >= -canard_delta_max_down);

%% Angular Acceleration (if rotating)
% I_yy for X-29 (estimate based on fighter aircraft)
I_yy = 30000;  % kg·m² (rough estimate)

if M_net > 0
    alpha_dot_dot = M_net / I_yy;  % rad/s²
    alpha_dot_dot_deg = alpha_dot_dot * 180/pi;  % deg/s²
else
    alpha_dot_dot = 0;
    alpha_dot_dot_deg = 0;
end

%% Package Results
results = struct();

% Input conditions
results.V = V;
results.V_kts = V * 1.944;
results.W_total = W_total;
results.W_fuel = W_fuel;
results.delta_canard = delta_canard;
results.alpha_wing = alpha_wing;

% Positions
results.x_cg = x_cg;
results.x_main_gear = x_main_gear;
results.x_canard_ac = x_canard_ac;
results.x_wing_ac = x_wing_ac;

% Moment arms (about main gear)
results.arm_cg = arm_cg;
results.arm_canard = arm_canard;
results.arm_wing = arm_wing;

% Aerodynamics
results.q = q;
results.CL_canard = CL_canard;
results.CL_wing = CL_wing;
results.L_canard = L_canard;
results.L_wing = L_wing;
results.L_total = L_total;
results.canard_stalled = canard_stalled;

% Moments (about main gear)
results.M_canard = M_canard;
results.M_wing = M_wing;
results.M_weight = M_weight;
results.M_aero = M_aero;
results.M_net = M_net;

% Stability derivatives
results.Cm_alpha = Cm_alpha;
results.Cm_delta_c = Cm_delta_c;
results.Cm_total = Cm_total;

% Rotation analysis
results.can_rotate = can_rotate;
results.delta_required = delta_required;
results.delta_feasible = delta_feasible;
results.alpha_dot_dot_deg_s2 = alpha_dot_dot_deg;

% Weight on gears
results.R_nose_static = R_nose_static;
results.nose_gear_load_percent = R_nose_static / W * 100;

%% Print Summary
fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('  X-29 ROTATION MOMENT ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('\n');
fprintf('CONDITIONS:\n');
fprintf('  Airspeed:        %.1f m/s (%.0f kts)\n', V, V*1.944);
fprintf('  Weight:          %.0f kg\n', W_total);
fprintf('  Canard δ:        %.1f deg\n', delta_canard);
fprintf('  Wing α:          %.1f deg\n', alpha_wing);
fprintf('  Dynamic pressure: %.0f Pa (%.2f kPa)\n', q, q/1000);
fprintf('\n');

fprintf('POSITIONS (from nose):\n');
fprintf('  CG:              %.2f m\n', x_cg);
fprintf('  Main gear:       %.2f m\n', x_main_gear);
fprintf('  Canard AC:       %.2f m\n', x_canard_ac);
fprintf('  Wing AC:         %.2f m\n', x_wing_ac);
fprintf('\n');

fprintf('MOMENT ARMS (about main gear, + = ahead):\n');
fprintf('  CG:              %+.2f m\n', arm_cg);
fprintf('  Canard:          %+.2f m\n', arm_canard);
fprintf('  Wing:            %+.2f m\n', arm_wing);
fprintf('\n');

fprintf('FORCES:\n');
fprintf('  Weight:          %.1f kN\n', W/1000);
fprintf('  Canard lift:     %.2f kN (CL = %.3f%s)\n', L_canard/1000, CL_canard, ...
    ternary(canard_stalled, ' STALLED', ''));
fprintf('  Wing lift:       %.2f kN (CL = %.3f)\n', L_wing/1000, CL_wing);
fprintf('  Total lift:      %.2f kN (%.1f%% of weight)\n', L_total/1000, L_total/W*100);
fprintf('\n');

fprintf('STABILITY DERIVATIVES (M=%.2f):\n', M);
fprintf('  Cm_alpha:        %.3f /rad (unstable > 0)\n', Cm_alpha);
fprintf('  Cm_delta_c:      %.3f /rad\n', Cm_delta_c);
fprintf('\n');

fprintf('MOMENTS ABOUT MAIN GEAR (+ = nose UP):\n');
fprintf('  M_weight:        %+.1f kN·m\n', M_weight/1000);
fprintf('  M_canard:        %+.1f kN·m\n', M_canard/1000);
fprintf('  M_wing:          %+.1f kN·m\n', M_wing/1000);
fprintf('  M_aero (Cm):     %+.1f kN·m  (instability + canard)\n', M_aero/1000);
fprintf('  ─────────────────────────────\n');
fprintf('  M_net:           %+.1f kN·m\n', M_net/1000);
fprintf('\n');

fprintf('ROTATION ANALYSIS:\n');
if can_rotate
    fprintf('  ✓ CAN ROTATE (M_net > 0)\n');
    fprintf('  Pitch acceleration: %.1f deg/s²\n', alpha_dot_dot_deg);
else
    fprintf('  ✗ CANNOT ROTATE (M_net < 0)\n');
end
fprintf('  Required δ for rotation: %.1f deg', delta_required);
if delta_feasible
    fprintf(' (within limits)\n');
else
    fprintf(' (EXCEEDS LIMITS)\n');
end
fprintf('  Static nose gear load: %.1f%% of weight\n', results.nose_gear_load_percent);
fprintf('═══════════════════════════════════════════════════════════════════\n');

end

%% Helper function
function out = ternary(cond, true_val, false_val)
    if cond
        out = true_val;
    else
        out = false_val;
    end
end
