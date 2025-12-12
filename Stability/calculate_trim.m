function [trim] = calculate_trim(V, h, W_total, W_fuel)
%CALCULATE_TRIM Find trimmed flight condition (L=W, M=0)
%
% Solves for alpha and delta_canard such that:
%   1. Total lift = Weight
%   2. Total pitching moment = 0
%
% Inputs:
%   V        - Airspeed [m/s]
%   h        - Altitude [m]
%   W_total  - Total aircraft weight [kg]
%   W_fuel   - Current fuel weight [kg]
%
% Outputs:
%   trim - Structure with trim solution
%
% Author: Dominic Larin
% Date: December 2025

%% Default inputs
if nargin < 1, V = 200; end           % Cruise speed [m/s]
if nargin < 2, h = 12200; end         % Cruise altitude [m]
if nargin < 3, W_total = 7500; end    % Weight [kg]
if nargin < 4, W_fuel = 1000; end     % Fuel [kg]

%% Constants
g = 9.81;
ft2m = 0.3048;

%% Atmosphere
[T_atm, a, P, rho] = atmosisa(h);
M = V / a;
q = 0.5 * rho * V^2;

%% Aircraft Parameters
S_ref = 17.54;              % Wing reference area [m²]
S_canard = 3.4923;          % Canard area [m²]
c_bar = 8.0 * ft2m;         % MAC [m]

% Positions from nose [m]
x_wing_ac = 40.91 * ft2m;   % 12.47 m
x_canard_ac = 6.49;         % m

% Get CG
[x_cg, ~] = calculate_X29_CG(W_total, W_fuel, 1877);

% Moment arms about CG (positive = ahead of CG)
arm_wing = x_cg - x_wing_ac;      % Negative (wing behind CG)
arm_canard = x_cg - x_canard_ac;  % Positive (canard ahead of CG)

%% Aerodynamic Parameters
% Wing
CL_alpha_wing = 0.08;       % [1/deg]
CL_0_wing = 0.05;           % Zero-alpha lift

% Canard
CL_alpha_canard = 0.0263;   % [1/deg]

% Stability derivatives
[Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M);

% Limits
delta_max_up = 58;
delta_max_down = 32;
alpha_max = 15;  % Stall limit

%% Weight
W = W_total * g;  % [N]

%% Solve trim equations
% Variables: x = [alpha_deg, delta_canard_deg]
%
% Equation 1 (Lift = Weight):
%   q*S_ref*(CL_0 + CL_alpha_wing*alpha) + q*S_canard*CL_alpha_canard*delta = W
%
% Equation 2 (Moment = 0 about CG):
%   L_wing*arm_wing + L_canard*arm_canard + q*S_ref*c_bar*(Cm_alpha*alpha_rad + Cm_delta_c*delta_rad) = 0

% Initial guess
x0 = [5; 10];  % [alpha_deg, delta_deg]

% Solve
options = optimoptions('fsolve', 'Display', 'off', 'TolFun', 1e-10);
[x_sol, fval, exitflag] = fsolve(@(x) trim_equations(x, q, S_ref, S_canard, c_bar, ...
    CL_0_wing, CL_alpha_wing, CL_alpha_canard, arm_wing, arm_canard, ...
    Cm_alpha, Cm_delta_c, W), x0, options);

alpha_trim = x_sol(1);
delta_trim = x_sol(2);

%% Check solution validity
valid = true;
warnings = {};

if alpha_trim > alpha_max
    valid = false;
    warnings{end+1} = sprintf('Alpha (%.1f°) exceeds stall limit (%.1f°)', alpha_trim, alpha_max);
end
if alpha_trim < -5
    valid = false;
    warnings{end+1} = sprintf('Alpha (%.1f°) is negative/unrealistic', alpha_trim);
end
if delta_trim > delta_max_up
    valid = false;
    warnings{end+1} = sprintf('Canard (%.1f°) exceeds TEU limit (%.1f°)', delta_trim, delta_max_up);
end
if delta_trim < -delta_max_down
    valid = false;
    warnings{end+1} = sprintf('Canard (%.1f°) exceeds TED limit (%.1f°)', delta_trim, -delta_max_down);
end

%% Calculate trim aerodynamics
CL_wing = CL_0_wing + CL_alpha_wing * alpha_trim;
CL_canard = CL_alpha_canard * delta_trim;
L_wing = q * S_ref * CL_wing;
L_canard = q * S_canard * CL_canard;
L_total = L_wing + L_canard;

% Drag (simple model)
CD0 = 0.025;
K = 0.10;
if M > 1.0
    CD0 = 0.045;
    K = 0.15;
end
CL_total = L_total / (q * S_ref);
CD = CD0 + K * CL_total^2;
D = q * S_ref * CD;
LD = L_total / D;

% Moments about CG
M_wing = L_wing * arm_wing;
M_canard = L_canard * arm_canard;
alpha_rad = alpha_trim * pi/180;
delta_rad = delta_trim * pi/180;
M_aero = q * S_ref * c_bar * (Cm_alpha * alpha_rad + Cm_delta_c * delta_rad);
M_total = M_wing + M_canard + M_aero;

%% Package results
trim = struct();

% Flight condition
trim.V = V;
trim.V_kts = V * 1.944;
trim.h = h;
trim.h_ft = h / ft2m;
trim.M = M;
trim.q = q;

% Weight & CG
trim.W_total = W_total;
trim.W_fuel = W_fuel;
trim.x_cg = x_cg;

% Trim solution
trim.alpha_deg = alpha_trim;
trim.delta_canard_deg = delta_trim;
trim.valid = valid;
trim.warnings = warnings;
trim.converged = (exitflag > 0);
trim.residual = norm(fval);

% Aerodynamics
trim.CL_wing = CL_wing;
trim.CL_canard = CL_canard;
trim.CL_total = CL_total;
trim.CD = CD;
trim.LD = LD;
trim.L_total = L_total;
trim.D = D;

% Moments (should be ~0 at trim)
trim.M_wing = M_wing;
trim.M_canard = M_canard;
trim.M_aero = M_aero;
trim.M_total = M_total;

% Stability
trim.Cm_alpha = Cm_alpha;
trim.Cm_delta_c = Cm_delta_c;

%% Print results
fprintf('\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('  X-29 TRIM ANALYSIS\n');
fprintf('═══════════════════════════════════════════════════════════════════\n');
fprintf('\n');

fprintf('FLIGHT CONDITION:\n');
fprintf('  Altitude:        %.0f m (%.0f ft)\n', h, h/ft2m);
fprintf('  Airspeed:        %.1f m/s (%.0f kts, M=%.2f)\n', V, V*1.944, M);
fprintf('  Weight:          %.0f kg\n', W_total);
fprintf('  Dynamic pressure: %.2f kPa\n', q/1000);
fprintf('  CG:              %.2f m from nose\n', x_cg);
fprintf('\n');

fprintf('TRIM SOLUTION:\n');
if trim.converged
    fprintf('  ✓ Converged (residual = %.2e)\n', trim.residual);
else
    fprintf('  ✗ Did not converge!\n');
end
fprintf('  Alpha:           %.2f deg\n', alpha_trim);
fprintf('  Canard δ:        %.2f deg\n', delta_trim);

if ~valid
    fprintf('\n  ⚠ WARNINGS:\n');
    for i = 1:length(warnings)
        fprintf('    - %s\n', warnings{i});
    end
end
fprintf('\n');

fprintf('AERODYNAMICS:\n');
fprintf('  CL_wing:         %.4f\n', CL_wing);
fprintf('  CL_canard:       %.4f\n', CL_canard);
fprintf('  CL_total:        %.4f\n', CL_total);
fprintf('  CD:              %.4f\n', CD);
fprintf('  L/D:             %.2f\n', LD);
fprintf('\n');

fprintf('FORCE BALANCE:\n');
fprintf('  Weight:          %.2f kN\n', W/1000);
fprintf('  Total Lift:      %.2f kN\n', L_total/1000);
fprintf('  Lift - Weight:   %.2f N (should be ~0)\n', L_total - W);
fprintf('\n');

fprintf('MOMENT BALANCE (about CG):\n');
fprintf('  M_wing:          %+.2f kN·m\n', M_wing/1000);
fprintf('  M_canard:        %+.2f kN·m\n', M_canard/1000);
fprintf('  M_aero:          %+.2f kN·m\n', M_aero/1000);
fprintf('  Total:           %+.2f N·m (should be ~0)\n', M_total);
fprintf('═══════════════════════════════════════════════════════════════════\n');

end

%% Trim equations to solve
function F = trim_equations(x, q, S_ref, S_canard, c_bar, ...
    CL_0_wing, CL_alpha_wing, CL_alpha_canard, arm_wing, arm_canard, ...
    Cm_alpha, Cm_delta_c, W)

    alpha_deg = x(1);
    delta_deg = x(2);

    alpha_rad = alpha_deg * pi/180;
    delta_rad = delta_deg * pi/180;

    % Lift coefficients
    CL_wing = CL_0_wing + CL_alpha_wing * alpha_deg;
    CL_canard = CL_alpha_canard * delta_deg;

    % Lift forces
    L_wing = q * S_ref * CL_wing;
    L_canard = q * S_canard * CL_canard;
    L_total = L_wing + L_canard;

    % Moments about CG
    M_wing = L_wing * arm_wing;
    M_canard = L_canard * arm_canard;
    M_aero = q * S_ref * c_bar * (Cm_alpha * alpha_rad + Cm_delta_c * delta_rad);
    M_total = M_wing + M_canard + M_aero;

    % Equations: want both = 0
    F(1) = (L_total - W) / W;           % Normalized lift equation
    F(2) = M_total / (q * S_ref * c_bar); % Normalized moment equation
end
