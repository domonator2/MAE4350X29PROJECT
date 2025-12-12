%{
====================================================================================================
FUNCTION NAME: sim_accel_canards.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 12/06/2025
====================================================================================================
FUNCTION DESCRIPTION:
Transonic/supersonic acceleration simulation for X-29 using ODE45 integration.
Level acceleration at constant altitude with canard tracking.

Uses CALIBRATED physics-based propulsion model via get_F404_thrust_calibrated.
Uses RAYMER aerodynamics model via get_X29_aero.

Equations of Motion (level flight, constant altitude):
  ẋ = V
  V̇ = (T - D) / m
  Ẇ = -fuel_flow * g

Models transonic drag rise (wave drag) which peaks around M=1.0-1.2.

====================================================================================================
INPUTS:
W0              - [kg] Weight at start of acceleration
S_ref           - [m²] Wing reference area
h               - [m] Altitude (constant during acceleration)
M0              - [-] Initial Mach number
M_target        - [-] Target Mach number
geom            - Geometry structure for Raymer aero model

OUTPUTS:
results         - [struct] Simulation results with canard tracking
====================================================================================================
%}

function results = sim_accel_canards(W0, S_ref, h, M0, M_target, geom)

% =========================================================================
% DEFAULTS
% =========================================================================
if nargin < 1 || isempty(W0)
    W0 = 7616;  % kg (after climb)
end
if nargin < 2 || isempty(S_ref)
    S_ref = 17.54;  % m²
end
if nargin < 3 || isempty(h)
    h = 12200;  % m (cruise altitude)
end
if nargin < 4 || isempty(M0)
    M0 = 0.9;  % Start Mach
end
if nargin < 5 || isempty(M_target)
    M_target = 1.4;  % Target cruise Mach
end
if nargin < 6 || isempty(geom)
    geom = [];
end

% =========================================================================
% CONSTANTS
% =========================================================================
g = 9.81;

% =========================================================================
% ATMOSPHERE AT ALTITUDE
% =========================================================================
[T_atm, a, P, rho] = atmosisa(h);
h_ft = h / 0.3048;

% Store in struct for EOM function
atm.a = a;
atm.rho = rho;
atm.h = h;
atm.h_ft = h_ft;

% Aero parameters struct
aero.S_ref = S_ref;
aero.AR = 4.0;
aero.geom = geom;
aero.h = h;

% =========================================================================
% HEADER
% =========================================================================
fprintf('╔════════════════════════════════════════════════════════════════╗\n');
fprintf('║   X-29 TRANSONIC ACCELERATION (Physics Propulsion + Raymer)   ║\n');
fprintf('╚════════════════════════════════════════════════════════════════╝\n\n');

fprintf('ACCELERATION PROFILE:\n');
fprintf('  Altitude:      %.0f m (%.0f ft) - constant\n', h, h_ft);
fprintf('  Initial Mach:  %.2f\n', M0);
fprintf('  Target Mach:   %.2f\n', M_target);
fprintf('  Propulsion:    Physics-based 0-D cycle analysis\n');
fprintf('  Aerodynamics:  Raymer model (get_X29_aero)\n');
fprintf('\n');

% Initial velocity
V0 = M0 * a;
V_target = M_target * a;
fprintf('  Initial V:     %.1f m/s (%.0f kts)\n', V0, V0*1.944);
fprintf('  Target V:      %.1f m/s (%.0f kts)\n', V_target, V_target*1.944);
fprintf('\n');

% =========================================================================
% ODE45 SETUP
% =========================================================================
fprintf('SIMULATING ACCELERATION (ODE45)...\n');

% Initial state: [x, V, W]
y0 = [0; V0; W0 * g];

% ODE options with event function for target Mach
ode_opts = odeset('RelTol', 1e-6, 'AbsTol', 1e-8, 'MaxStep', 1.0, ...
                  'Events', @(t,y) event_mach_target(t, y, M_target, atm));

tspan = [0, 600];  % Max 10 minutes

% Integrate
[t_ode, Y_ode, te, ye, ie] = ode45(@(t,y) accel_EOM_physics(t, y, aero, atm, g), tspan, y0, ode_opts);

% =========================================================================
% POST-PROCESS RESULTS
% =========================================================================
n_points = length(t_ode);

% Extract state vectors
x_hist = Y_ode(:, 1);
V_hist = Y_ode(:, 2);
W_hist = Y_ode(:, 3) / g;  % Convert to kg

% Compute derived quantities
M_hist = V_hist / a;
q_hist = 0.5 * rho * V_hist.^2;
h_hist = h * ones(n_points, 1);  % Constant altitude

% Compute thrust, drag, accel at each point for output
T_hist = zeros(n_points, 1);
D_hist = zeros(n_points, 1);
accel_hist = zeros(n_points, 1);
PLA_hist = zeros(n_points, 1);
fuel_rate_hist = zeros(n_points, 1);
alpha_hist = zeros(n_points, 1);
CL_hist = zeros(n_points, 1);
CD_hist = zeros(n_points, 1);
delta_canard_hist = zeros(n_points, 1);
x_cg_hist = zeros(n_points, 1);

for i = 1:n_points
    V = V_hist(i);
    W = W_hist(i);
    M = M_hist(i);
    q = q_hist(i);

    % Level flight: L = W, solve for CL
    CL = (W * g) / (q * S_ref);
    CL_hist(i) = CL;

    % Estimate alpha from CL using proper CL_alpha from Raymer model
    if ~isempty(geom)
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, geom);
        alpha_deg = CL / CL_alpha;
    else
        alpha_deg = CL / 0.08;  % Fallback
    end
    alpha_deg = max(0, min(15, alpha_deg));
    alpha_hist(i) = alpha_deg;

    % Aerodynamics using Raymer model
    if ~isempty(geom)
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, geom);
    else
        CD0 = get_CD0_accel(M);
        K = get_K_accel(M);
        CD = CD0 + K * CL^2;
    end
    CD_hist(i) = CD;
    D = q * S_ref * CD;
    D_hist(i) = D;

    % Thrust using PHYSICS model (max AB for acceleration)
    prop = get_F404_thrust_calibrated(M, h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;
    T_hist(i) = T;
    fuel_rate_hist(i) = fuel_flow;
    PLA_hist(i) = 130;  % Max AB

    % Acceleration
    accel_hist(i) = (T - D) / W;

    % Canard deflection for trim
    [Cm_alpha, Cm_delta_c, ~, ~] = get_Cm_derivatives(M);
    if abs(Cm_delta_c) > 0.01
        delta_rad = -Cm_alpha * (alpha_deg * pi/180) / Cm_delta_c;
        delta_deg = delta_rad * 180/pi;
    else
        delta_deg = 5;
    end
    delta_canard_hist(i) = max(-32, min(58, delta_deg));

    % CG position
    W_fuel_est = max(0, W - 6278);  % Estimate fuel from weight
    [x_cg_hist(i), ~] = calculate_X29_CG(W, W_fuel_est, 1877);
end

% Compute cumulative fuel burned from weight change
fuel_burned_hist = W_hist(1) - W_hist;

% Final values
t_final = t_ode(end);
V_final = V_hist(end);
M_final = M_hist(end);
W_final = W_hist(end);
total_fuel = W_hist(1) - W_final;

% =========================================================================
% SUMMARY
% =========================================================================
fprintf('\nACCELERATION SUMMARY:\n');
fprintf('  Mach:          %.2f → %.2f\n', M0, M_final);
fprintf('  Velocity:      %.0f → %.0f m/s\n', V0, V_final);
fprintf('  Time:          %.1f s\n', t_final);
fprintf('  Fuel burned:   %.1f kg\n', total_fuel);
fprintf('  Final weight:  %.0f kg\n', W_final);
fprintf('  Canard range:  %.1f° to %.1f°\n', min(delta_canard_hist), max(delta_canard_hist));

% =========================================================================
% PACKAGE RESULTS
% =========================================================================
results = struct();
results.t = t_ode;
results.x = x_hist;
results.h = h_hist;
results.V = V_hist;
results.M = M_hist;
results.q = q_hist;
results.W = W_hist;
results.T = T_hist;
results.D = D_hist;
results.CD = CD_hist;
results.CL = CL_hist;
results.accel = accel_hist;
results.fuel_burned = fuel_burned_hist;
results.PLA = PLA_hist;
results.alpha = alpha_hist;
results.delta_canard = delta_canard_hist;
results.x_cg = x_cg_hist;

results.total_time = t_final;
results.total_fuel = total_fuel;
results.W_final = W_final;
results.M_final = M_final;
results.V_final = V_final;

results.config.W0 = W0;
results.config.h = h;
results.config.M0 = M0;
results.config.M_target = M_target;

% Total horizontal distance
results.distance = x_hist(end);

end

% =========================================================================
% EVENT FUNCTION: Stop when Mach reaches target
% =========================================================================
function [value, isterminal, direction] = event_mach_target(t, y, M_target, atm)
    V = y(2);
    M = V / atm.a;
    value = M - M_target;
    isterminal = 1;
    direction = 1;  % Increasing
end

% =========================================================================
% EQUATIONS OF MOTION: Level acceleration (PHYSICS PROPULSION)
% =========================================================================
function dydt = accel_EOM_physics(t, y, aero, atm, g)
    % State: [x, V, W]
    x = y(1);
    V = max(10, y(2));
    W = y(3);  % Weight in Newtons

    % Atmosphere (constant altitude)
    a = atm.a;
    rho = atm.rho;
    h = atm.h;

    M = V / a;
    q = 0.5 * rho * V^2;

    % Level flight: L = W, solve for CL
    CL = W / (q * aero.S_ref);

    % Aerodynamics - use Raymer model
    if ~isempty(aero.geom)
        % Get CL_alpha from Raymer model for proper alpha estimation
        [~, ~, ~, ~, ~, CL_alpha] = get_X29_aero(M, 5, h, aero.geom);
        alpha_deg = CL / CL_alpha;
        alpha_deg = max(0, min(15, alpha_deg));
        [~, CD, ~, ~, ~, ~] = get_X29_aero(M, alpha_deg, h, aero.geom);
    else
        % Fallback to simple Mach-dependent model
        CD0 = get_CD0_accel(M);
        K = get_K_accel(M);
        CD = CD0 + K * CL^2;
    end
    D = q * aero.S_ref * CD;

    % Thrust - PHYSICS MODEL with max AB for acceleration
    prop = get_F404_thrust_calibrated(max(0.3, M), h, 130);  % PLA=130 for max AB
    T = prop.thrust_N;
    fuel_flow = prop.fuel_flow_kg_s;

    % Equations of motion (level flight)
    dxdt = V;
    dVdt = (T - D) * g / W;  % a = F/m = F*g/W
    dWdt = -fuel_flow * g;

    dydt = [dxdt; dVdt; dWdt];
end

% =========================================================================
% CD0 MODEL WITH TRANSONIC DRAG RISE (fallback if no geometry)
% =========================================================================
function CD0 = get_CD0_accel(M)
    if M < 0.85
        CD0 = 0.025;
    elseif M < 1.05
        x = (M - 0.85) / (1.05 - 0.85);
        CD0 = 0.025 + 0.035 * (3*x^2 - 2*x^3);
    elseif M < 1.2
        CD0 = 0.060;
    else
        CD0 = 0.060 - 0.005 * (M - 1.2);
        CD0 = max(0.050, CD0);
    end
end

% =========================================================================
% INDUCED DRAG FACTOR (K) MODEL (fallback if no geometry)
% =========================================================================
function K = get_K_accel(M)
    AR = 4.0;
    if M < 0.85
        e = 0.80;
        K = 1 / (pi * AR * e);
    elseif M < 1.2
        e = 0.70;
        K = 1 / (pi * AR * e);
    else
        K = 0.17;
    end
end
