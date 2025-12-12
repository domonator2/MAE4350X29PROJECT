function canard = canard_properties()
%CANARD_PROPERTIES X-29 canard geometry and aerodynamics for 3DOF pitch analysis
%
% Returns a structure containing canard parameters needed for
% longitudinal (pitch) stability and control analysis.
%
% Usage:
%   canard = canard_properties();
%
% Output:
%   canard - Structure with fields:
%       .geom       - Geometry parameters
%       .aero       - Aerodynamic coefficients
%       .control    - Control surface effectiveness
%       .limits     - Deflection limits
%
% Reference: NASA TM-4140, Grumman X-29 data
%
% Author: Dominic Larin
% Date: December 2025

canard = struct();

%% =========================================================================
%  GEOMETRY (Pitch-relevant parameters only)
%  =========================================================================
% All positions measured from nose tip (positive aft)

% Planform
canard.geom.S = 3.4923;              % [m²] Canard planform area
canard.geom.b = 4.393;               % [m] Canard span (full span, both sides)
canard.geom.AR = 1.47;               % [-] Aspect ratio (b²/S)
canard.geom.lambda = 0.32;           % [-] Taper ratio (c_tip/c_root)
canard.geom.c_root = 1.193;          % [m] Root chord
canard.geom.c_tip = 0.382;           % [m] Tip chord (= lambda * c_root)
canard.geom.c_mac = 0.864;           % [m] Mean aerodynamic chord

% Sweep angles
canard.geom.Lambda_LE = 42;          % [deg] Leading edge sweep

% Calculate other sweep angles using Kshift method
AR = canard.geom.AR;
lam = canard.geom.lambda;
Kshift = (4/AR) * (1-lam)/(1+lam);
tanLE = tand(canard.geom.Lambda_LE);
canard.geom.Lambda_c4 = atand(tanLE - 0.25*Kshift);   % [deg] Quarter-chord sweep
canard.geom.Lambda_TE = atand(tanLE - 1.00*Kshift);   % [deg] Trailing edge sweep

% Position (from nose tip, positive aft)
canard.geom.x_ac = 6.49;             % [m] X-position of canard AC from nose (at hinge line)

% Incidence
canard.geom.i_canard = 0;            % [deg] Fixed incidence angle (relative to FRL)

%% =========================================================================
%  AERODYNAMIC COEFFICIENTS (Pitch-relevant)
%  =========================================================================
% All coefficients are for the isolated canard unless noted

% Lift curve slope
canard.aero.CL_alpha_deg = 0.0263;   % [1/deg] Lift curve slope
canard.aero.CL_alpha = canard.aero.CL_alpha_deg * (180/pi);  % [1/rad] Lift curve slope

% Zero-lift properties
canard.aero.CL_0 = 0;                % [-] Zero-alpha lift coefficient (symmetric airfoil)
canard.aero.alpha_0L = 0;            % [deg] Zero-lift angle of attack

% Stall
canard.aero.CL_max = NaN;            % [-] Maximum lift coefficient
canard.aero.alpha_stall = NaN;       % [deg] Stall angle of attack

% Pitching moment (about canard AC)
canard.aero.Cm_ac = 0;               % [-] Moment coefficient about AC (symmetric airfoil ~ 0)

% Downwash effects on wing (canard-induced)
canard.aero.d_epsilon_d_alpha = NaN; % [-] Downwash gradient (dε/dα) at wing
canard.aero.epsilon_0 = 0;           % [deg] Downwash at zero canard lift

%% =========================================================================
%  CONTROL SURFACE EFFECTIVENESS (All-moving canard)
%  =========================================================================
% X-29 canard is all-moving (variable incidence for pitch control)

% For all-moving surface: CL_delta ≈ CL_alpha
canard.control.CL_delta_deg = canard.aero.CL_alpha_deg;  % [1/deg] ∂CL_canard/∂δ_canard
canard.control.CL_delta = canard.aero.CL_alpha;          % [1/rad] ∂CL_canard/∂δ_canard

% Pitching moment effectiveness (about aircraft CG) - calculated by canard_init
canard.control.Cm_delta = NaN;       % [1/rad] ∂Cm/∂δ_canard (about aircraft CG)
canard.control.Cm_delta_deg = NaN;   % [1/deg] ∂Cm/∂δ_canard

%% =========================================================================
%  DEFLECTION LIMITS
%  =========================================================================
canard.limits.delta_max_up = 58;     % [deg] Maximum trailing-edge up (nose up command)
canard.limits.delta_max_down = 32;   % [deg] Maximum trailing-edge down (nose down command)
canard.limits.delta_rate_max = NaN;  % [deg/s] Maximum deflection rate (if known)

%% =========================================================================
%  AIRCRAFT REFERENCE (for moment arm calculations)
%  =========================================================================
% These will be set by canard_init()
canard.aircraft.x_cg = NaN;          % [m] CG position from nose (to be set)
canard.aircraft.c_bar = NaN;         % [m] Aircraft reference MAC (to be set)
canard.aircraft.S_ref = NaN;         % [m²] Aircraft reference area (to be set)

%% =========================================================================
%  CALCULATED PARAMETERS (populated by canard_init)
%  =========================================================================
canard.calc.l_canard = NaN;          % [m] Moment arm: x_cg - x_ac_canard (positive = canard ahead of CG)
canard.calc.V_H = NaN;               % [-] Canard volume coefficient: (l_c * S_c) / (c_bar * S_ref)
canard.calc.eta_canard = 1.0;        % [-] Dynamic pressure ratio at canard (q_canard/q_inf)

%% =========================================================================
%  PRINT SUMMARY
%  =========================================================================
fprintf('\n=== X-29 CANARD PROPERTIES (3DOF Pitch) ===\n');
fprintf('Geometry:\n');
fprintf('  S = %.4f m²,  b = %.3f m,  AR = %.2f\n', canard.geom.S, canard.geom.b, canard.geom.AR);
fprintf('  c_root = %.3f m,  c_tip = %.3f m,  MAC = %.3f m\n', canard.geom.c_root, canard.geom.c_tip, canard.geom.c_mac);
fprintf('  Lambda_LE = %.1f°,  Lambda_c/4 = %.1f°,  Lambda_TE = %.1f°\n', ...
    canard.geom.Lambda_LE, canard.geom.Lambda_c4, canard.geom.Lambda_TE);
fprintf('Aerodynamics:\n');
fprintf('  CL_alpha = %.4f /deg (%.3f /rad)\n', canard.aero.CL_alpha_deg, canard.aero.CL_alpha);
fprintf('\nCall canard_init(canard, x_cg, c_bar, S_ref) to calculate moment arm and Cm_delta.\n');

end

%% =========================================================================
%  INITIALIZATION FUNCTION (call after setting aircraft parameters)
%  =========================================================================
function canard = canard_init(canard, x_cg, c_bar, S_ref)
%CANARD_INIT Calculate derived canard parameters for pitch analysis
%
% Inputs:
%   canard - Canard structure from canard_properties()
%   x_cg   - [m] Aircraft CG position from nose
%   c_bar  - [m] Aircraft reference MAC
%   S_ref  - [m²] Aircraft reference wing area
%
% Output:
%   canard - Updated structure with calculated parameters

% Store aircraft reference
canard.aircraft.x_cg = x_cg;
canard.aircraft.c_bar = c_bar;
canard.aircraft.S_ref = S_ref;

% Moment arm (positive = canard ahead of CG, which is typical)
canard.calc.l_canard = x_cg - canard.geom.x_ac;

% Canard volume coefficient
canard.calc.V_H = (canard.calc.l_canard * canard.geom.S) / (c_bar * S_ref);

% Calculate Cm_delta for pitch control effectiveness
% Cm_delta = -CL_delta * (l_canard / c_bar) * (S_canard / S_ref) * eta
canard.control.Cm_delta = -canard.control.CL_delta * ...
    (canard.calc.l_canard / c_bar) * (canard.geom.S / S_ref) * canard.calc.eta_canard;
canard.control.Cm_delta_deg = canard.control.Cm_delta * pi/180;

fprintf('\n=== CANARD INITIALIZED FOR PITCH ANALYSIS ===\n');
fprintf('  Moment arm (l_c):     %.3f m\n', canard.calc.l_canard);
fprintf('  Volume coefficient:   %.4f\n', canard.calc.V_H);
fprintf('  Cm_delta:             %.4f /rad (%.6f /deg)\n', ...
    canard.control.Cm_delta, canard.control.Cm_delta_deg);

end
