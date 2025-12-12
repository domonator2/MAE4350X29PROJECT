%{
====================================================================================================
FUNCTION NAME: hypersonic_convergence_elliptical.m
AUTHOR: Dominic Larin
INITIATED: 9/26/2025
LAST REVISION: 10/13/2025 - Modified to use elliptical fuselage model
====================================================================================================
FUNCTION DESCRIPTION:
This function implements an advanced physics-based convergence methodology for supersonic aircraft
sizing and design optimization using the ELLIPTICAL FUSELAGE MODEL. It integrates real aerodynamics,
propulsion, geometry, and trajectory analysis with simultaneous weight and volume budget convergence
to determine realistic vehicle configurations. The function uses a dual-budget approach with enhanced
structural modeling and disciplinary integration for comprehensive X-29 aircraft design analysis.

MODIFICATION: Uses test_elliptical_fuselage() instead of get_geometry() for all geometry calculations.
This validates that the new elliptical fuselage model produces viable aircraft designs.

Key Features:
- Dual-budget convergence (weight budget + volume budget)
- Integrated disciplinary analysis (aerodynamics, propulsion, geometry)
- Unified aerodynamics model with automatic regime selection
- Enhanced structural weight modeling using component buildup
- Realistic fuel fraction calculation using Breguet range equation
- Robust convergence algorithm with multiple retry strategies
- Comprehensive error handling and fallback calculations

Technical Methodology:
- Weight Budget: Balances structural, systems, engine, and inert weights with margins
- Volume Budget: Uses Kuechemann's slenderness parameter for fuel and systems packaging
- Convergence: Employs fsolve with trust-region-dogleg algorithm for equation solving
- Aerodynamics: Unified raymer_aero_model with subsonic/transonic/supersonic capability
- Propulsion: F404 engine performance data with realistic TSFC interpolation
====================================================================================================
INPUTS:
mission         - [struct] Mission parameters structure
  .payload_kg   - [kg] Payload weight
  .crew_count   - [-] Number of crew members
  .range_km     - [km] Mission range requirement

vehicle         - [struct] Vehicle configuration structure
  .S_ref_m2     - [m^2] Wing reference area
  .tau          - [-] Slenderness parameter (Kuechemann)
  .sweep_deg    - [deg] Leading edge sweep angle
  .E_TW         - [-] Engine thrust-to-weight ratio
  .mu_a         - [-] Margin factor
  .f_sys        - [-] Systems weight fraction
  .C_sys_kg     - [kg] Constant systems weight
  .rho_pay_kg_m3- [kg/m^3] Payload density
  .rho_fuel_kg_m3- [kg/m^3] Fuel density
  .k_vv         - [-] Void volume fraction
  .k_vs         - [-] Systems volume fraction
  .k_ve         - [-] Engine volume fraction
  .v_fix_m3     - [m^3] Fixed volume
  .n_ult        - [-] Ultimate load factor
  .Kw           - [-] Wing factor
  .lambda_wing  - [-] Wing taper ratio
  .tc_max       - [-] Max thickness-to-chord ratio

initialGuess    - [struct] Initial estimates structure
  .TOGW_kg      - [kg] Initial takeoff gross weight estimate
  .S_pln_m2     - [m^2] Initial planform area estimate

opts            - [struct] Options structure (optional)
  .tolerance    - [-] Convergence tolerance (default: 1e-3)
  .maxIterations- [-] Maximum iterations (default: 150)
  .units        - [string] Unit system 'SI' or 'Imperial' (default: 'SI')
  .verbose      - [logical] Enable detailed output (default: true)

OUTPUTS:
convergedDesign - [struct] Converged design parameters structure
  .TOGW_kg      - [kg] Converged takeoff gross weight
  .OEW_kg       - [kg] Operating empty weight
  .S_pln_m2     - [m^2] Converged planform area
  .fuel_fraction- [-] Mission fuel fraction
  .LD_cruise    - [-] Cruise lift-to-drag ratio
  .wing_loading - [kg/m^2] Wing loading
  .T_W_required - [-] Required thrust-to-weight ratio
  .iterations   - [-] Number of convergence iterations
  .error_final  - [-] Final convergence error
  .engine_performance - [struct] Engine performance data (if available)
  .geometry     - [struct] Vehicle geometry data
  .CL_cruise    - [-] Cruise lift coefficient
  .CD_cruise    - [-] Cruise drag coefficient

status          - [string] Convergence status
                  'Converged' - Successful convergence
                  'Failed: Maximum iterations reached' - Iteration limit exceeded
                  'Failed: Solver convergence error' - fsolve failed to converge

DEPENDENCIES:
External Functions:
- test_elliptical_fuselage.m: Vehicle geometry generation with elliptical fuselage (Geometry discipline)
- standard_atmosphere_1976.m: Atmospheric properties (Common utilities)
- raymer_aero_model.m: Unified aerodynamics analysis (Aerodynamics discipline)
- get_tsfc.m: Thrust specific fuel consumption (Propulsion discipline)
- get_F404_performance.m: Engine performance data (Propulsion discipline)

MATLAB Toolboxes Required:
- Optimization Toolbox (for fsolve nonlinear equation solver)
====================================================================================================
REVISION HISTORY:
9/26/2025 - D. Larin - Initial implementation with basic convergence framework
9/27/2025 - D. Larin - Enhanced budget equations and structural weight modeling
9/28/2025 - D. Larin - Integrated unified aerodynamics model and comprehensive documentation
====================================================================================================
%}

function [convergedDesign, status] = hypersonic_convergence_elliptical(mission, vehicle, initialGuess, opts)

%% Input Validation and Setup
% Check if options structure is provided
if nargin < 4
    % [-] Create empty options structure if not provided
    opts = struct();
end

% Default options
% [-] Set default convergence tolerance (relaxed for stability)
if ~isfield(opts, 'tolerance'), opts.tolerance = 1e-3; end
% [-] Set default maximum iterations (increased for complex problems)
if ~isfield(opts, 'maxIterations'), opts.maxIterations = 150; end
% [string] Set default unit system to SI
if ~isfield(opts, 'units'), opts.units = 'SI'; end
% [logical] Set default to enable verbose output
if ~isfield(opts, 'verbose'), opts.verbose = true; end

% Mission parameters
% [m] Cruise altitude for X-29 mission profile
h_cruise = 12200;
% [-] Cruise Mach number for supersonic operation
M_cruise = 1.4;
% [deg] Cruise angle of attack for optimal L/D
alpha_cruise = 2.0;

% Extract sweep from vehicle configuration
% Check if sweep angle is specified in vehicle structure
if isfield(vehicle, 'sweep_deg')
    % [deg] Use specified sweep angle
    sweep_deg = vehicle.sweep_deg;
else
    % [deg] Default X-29 forward sweep configuration
    sweep_deg = -30;
end

%% Initialize Convergence Variables
% [kg] Extract initial takeoff gross weight estimate
TOGW_guess = initialGuess.TOGW_kg;
% [m^2] Extract initial planform area estimate
S_pln_guess = initialGuess.S_pln_m2;
% [-] Initialize iteration counter
iterationCount = 0;
% [-] Initialize convergence error (start with large value)
error_TOGW = 1.0;
% [-] Extract convergence tolerance from options
tolerance = opts.tolerance;
% [-] Extract maximum iterations from options
maxIterations = opts.maxIterations;

% Print convergence header
% Check if verbose output is enabled
if opts.verbose
    % [-] Print convergence analysis header
    fprintf('\n=== Enhanced X-29 Hypersonic Convergence ===\n');
    % [-] Print initial guess values
    fprintf('Initial Guess: TOGW=%.0f kg, S_pln=%.1f m²\n', TOGW_guess, S_pln_guess);
    % [-] Print convergence target description
    fprintf('Target: L/D~4.5, Realistic wave drag, Integrated aero/geometry\n');
    % [-] Print convergence table header
    fprintf('%-5s | %-10s | %-10s | %-8s | %-10s | %-8s\n', ...
            'Iter', 'TOGW(kg)', 'S_pln(m²)', 'L/D', 'WingLoad', 'Error');
    % [-] Print table separator line
    fprintf('------|------------|------------|----------|------------|----------\n');
end

%% Main Convergence Loop
% Continue while error exceeds tolerance and under iteration limit
while (error_TOGW > tolerance) && (iterationCount < maxIterations)
    % [-] Increment iteration counter
    iterationCount = iterationCount + 1;

    % [struct] Generate vehicle geometry from current sizing parameters
    % NOTE: Using ELLIPTICAL FUSELAGE model instead of spherocylinder
    % NOTE: Moved before wing loading calculation so S_ref scales with vehicle size
    geom = test_elliptical_fuselage(S_pln_guess, vehicle.tau, sweep_deg);

    % Update wing reference area from geometry (scales with S_pln)
    % [m^2] Wing reference area calculated as 85% of planform area (for aerodynamics)
    S_ref_aero = geom.Sref;  % This is S_w = 0.85 * S_pln (used for aero calculations)

    % Calculate current wing loading using CONSTANT X-29 reference area
    % [kg/m^2] Wing loading based on CONSTANT S_ref for comparable metrics across all vehicles
    wing_loading = TOGW_guess / vehicle.S_ref_m2;  % Always use 17.544 m² for wing loading

    % Calculate required thrust-to-weight ratio from constraints
    % [-] Required T/W from takeoff, cruise, and maneuver constraints
    T_W_required = calculate_thrust_constraints(wing_loading, M_cruise);

    % Calculate L/D using unified physics-based aerodynamics model
    % Always use raymer_aero_model for consistency across all configurations

    % Get atmospheric conditions at cruise
    % [struct] Standard atmosphere properties at cruise altitude
    atm_cruise = standard_atmosphere_1976(h_cruise, struct('units', 'SI'));

    % Set up options for unified aerodynamic model
    % [struct] Aerodynamics options structure
    opts_aero = struct();
    % [string] Unit system for aerodynamic calculations
    opts_aero.units = 'SI';
    % [struct] Component breakdown for drag buildup analysis
    opts_aero.drag_buildup = geom.components;
    % [float] Wave drag efficiency - calibrated to match Oct 28 baseline convergence
    opts_aero.Ewave = 1.2;  % Increased from 0.5 to add drag (spherocylinder default 2.0)
    % [bool] Disable trim and misc drags - these were calibrated for spherocylinder, not elliptical
    opts_aero.disable_empirical_drags = true;

    % Call unified aerodynamic model - USE STANDARD VERSION (elliptical version has missing dependencies)
    % [struct] Complete aerodynamic analysis - works with elliptical fuselage geometry
    aero_result = raymer_aero_model(M_cruise, alpha_cruise, geom, atm_cruise, opts_aero);

    % Extract aerodynamic coefficients
    % [-] Lift coefficient from unified analysis
    CL = aero_result.CL;
    % [-] Drag coefficient from unified analysis
    CD = aero_result.CD;
    % [-] Lift-to-drag ratio from unified analysis (actual L/D at current condition)
    LD_actual = aero_result.L_over_D;

    % REMOVED: Artificial L/D bounds - let physics determine actual values
    % LD_actual = max(2.0, min(8.0, LD_actual));

    % Calculate real engine performance at cruise conditions using updated interpolated F404 model
    % Use the new get_tsfc function with proper thrust interpolation and realistic TSFC values
    % [kg/s/N] Thrust specific fuel consumption from interpolated F404 data
    TSFC_actual = get_tsfc(M_cruise, h_cruise);

    % Calculate thrust using F404 data table for consistency
    % [ft] Convert cruise altitude from meters to feet for F404 model
    h_cruise_ft = h_cruise * 3.28084;
    % [deg] Military power setting for afterburner cruise operation
    PLA_cruise = 109;
    % [struct] Get F404 engine performance data at cruise conditions
    engine_data = get_F404_performance(M_cruise, h_cruise_ft, PLA_cruise, 'WT');
    % [N] Available thrust from F404 engine at cruise conditions
    thrust_available = engine_data.thrust_N;

    % Calculate fuel fraction using real TSFC and L/D
    % [-] Mission fuel fraction from time-based TSFC calculation
    % Pass cruise time if available, otherwise calculate from range
    if isfield(mission, 't_cruise')
        t_cruise_min = mission.t_cruise;
    else
        % Calculate cruise time from range if not provided
        atm = get_atmosphere(mission.h_cruise_m, 'SI');
        V_cruise = M_cruise * atm.a;
        t_cruise_min = (mission.range_km * 1000) / (V_cruise * 60);
    end
    fuel_fraction = calculate_fuel_fraction(LD_actual, TSFC_actual, t_cruise_min, M_cruise, TOGW_guess, opts);

    % Ensure fuel fraction is reasonable to prevent numerical issues
    % [-] Check if fuel fraction is within reasonable bounds
    if isnan(fuel_fraction) || fuel_fraction <= 0 || fuel_fraction >= 0.8
        % Calculate a more reasonable fuel fraction based on mission
        % [-] Check mission range for short missions
        if mission.range_km <= 100
            % [-] Set fuel fraction for short range missions
            fuel_fraction = 0.12;
        % [-] Check mission range for medium missions
        elseif mission.range_km <= 500
            % [-] Set fuel fraction for medium range missions
            fuel_fraction = 0.25;
        % [-] Handle long range missions
        else
            % [-] Set fuel fraction for longer range missions
            fuel_fraction = 0.4;
        % [-] End mission range check
        end
        % [-] Check if verbose output is enabled
        if opts.verbose
            % [-] Print fuel fraction adjustment warning
            fprintf('  Warning: Fuel fraction adjusted to %.3f for numerical stability\n', fuel_fraction);
        % [-] End verbose check
        end
    % [-] End fuel fraction bounds check
    end

    % Calculate weight ratio
    % [-] Weight ratio from fuel fraction for mission sizing
    WR = 1 / (1 - fuel_fraction);

    % Get structural parameters (mu_a is already set as constant in vehicle struct)
    % [kg/m^2, -] Enhanced structural analysis for weight calculations
    [I_str, K_w] = analyze_structure_enhanced(vehicle.tau, S_pln_guess, sweep_deg, TOGW_guess, opts);

    % Setup enhanced budget equations
    % [function handle] Create budget equations function for fsolve
    budget_function = @(x) enhanced_budget_equations(x, mission, vehicle, T_W_required, WR, I_str, K_w, LD_actual, TOGW_guess, opts);

    % Initial guess for solver [OEW, S_pln] - improved scaling
    % [kg] Crew weight estimate
    W_crew = mission.crew_count * 90;
    % [kg] Operating empty weight initial guess from weight balance
    OEW_guess = TOGW_guess * (1 - fuel_fraction) - mission.payload_kg - W_crew;

    % Improve initial guess with better bounds and validation
    % [-] Check if OEW guess is valid
    if isnan(OEW_guess) || OEW_guess <= 0
        % [kg] Default OEW to 60% of TOGW if invalid
        OEW_guess = 0.6 * TOGW_guess;
    % [-] End OEW validation check
    end
    % REMOVED: Artificial OEW bounds - let physics determine actual values
    % OEW_guess = max(3000, min(15000, OEW_guess));
    % REMOVED: Artificial S_pln bounds - let physics determine actual values
    % S_pln_guess = max(10, min(120, S_pln_guess));

    % Scale initial guess for better numerical conditioning
    % [-] Scale variables for improved numerical conditioning
    x0 = [OEW_guess/1000, S_pln_guess/10];

    % Solve budget equations with robust settings
    % [struct] Optimization options for robust convergence
    solver_options = optimoptions('fsolve', 'Display', 'off', ...
                                 'MaxFunctionEvaluations', 2000, ...
                                 'FunctionTolerance', 1e-4, ...
                                 'OptimalityTolerance', 1e-4, ...
                                 'StepTolerance', 1e-6, ...
                                 'Algorithm', 'trust-region-dogleg');

    % [array, array, int] Solve budget equations using fsolve
    [x_sol, fval, exitflag] = fsolve(budget_function, x0, solver_options);

    % Enhanced retry logic with multiple strategies
    % [-] Initialize retry counter
    retry_count = 0;
    % [-] Maximum number of retry attempts
    max_retries = 4;

    % [-] Continue retrying while solver fails and within retry limit
    while exitflag < 1 && retry_count < max_retries
        % [-] Increment retry counter
        retry_count = retry_count + 1;

        % [-] Check if first retry attempt
        if retry_count == 1
            % Strategy 1: Perturb initial guess
            % [-] Apply first perturbation strategy
            x0_retry = x0 .* [1.3, 0.7];
        % [-] Check if second retry attempt
        elseif retry_count == 2
            % Strategy 2: Different perturbation
            % [-] Apply second perturbation strategy
            x0_retry = x0 .* [0.8, 1.2];
        % [-] Check if third retry attempt
        elseif retry_count == 3
            % Strategy 3: Use conservative baseline
            % [-] Use conservative baseline values (8000kg OEW, 35m² S_pln)
            x0_retry = [8.0, 3.5];
        % [-] Final retry attempt
        else
            % Strategy 4: Try levenberg-marquardt algorithm
            % [string] Switch to Levenberg-Marquardt algorithm
            solver_options.Algorithm = 'levenberg-marquardt';
            % [-] Use original initial guess
            x0_retry = x0;
        % [-] End retry strategy selection
        end

        % [array, array, int] Retry solving with modified initial guess
        [x_sol, fval, exitflag] = fsolve(budget_function, x0_retry, solver_options);

        % If still failing, try relaxed tolerance
        % [-] Check if solver still failing after second retry
        if exitflag < 1 && retry_count >= 2
            % [-] Relax function tolerance for better convergence
            solver_options.FunctionTolerance = 1e-3;
            % [-] Relax optimality tolerance for better convergence
            solver_options.OptimalityTolerance = 1e-3;
        % [-] End tolerance relaxation check
        end
    % [-] End retry loop
    end

    % Final check - if still failed, return gracefully
    % [-] Check if solver ultimately failed
    if exitflag < 1
        % [string] Create failure status message with iteration and retry info
        status = sprintf('Failed: Solver convergence error at iteration %d after %d retries', iterationCount, retry_count);
        % [struct] Return empty converged design structure
        convergedDesign = struct();
        % [-] Exit function with failure status
        return;
    % [-] End final solver check
    end

    % Extract solution (scale back up)
    % [kg] Scale operating empty weight back to original units
    OEW_new = x_sol(1) * 1000;
    % [m^2] Scale planform area back to original units
    S_pln_new = x_sol(2) * 10;

    % REMOVED: Bounds checking and clamping - let physics determine actual values
    % bounds_violated = false;
    % if S_pln_new < 15 || S_pln_new > 80
    %     bounds_violated = true;
    % end
    % if OEW_new < 2000 || OEW_new > 15000
    %     bounds_violated = true;
    % end
    % S_pln_new = max(15, min(80, S_pln_new));
    % OEW_new = max(2000, min(15000, OEW_new));

    % Calculate new TOGW
    % [kg] Crew weight for TOGW calculation
    W_crew = mission.crew_count * 90;

    % Ensure fuel fraction is real to prevent complex TOGW
    fuel_fraction = real(fuel_fraction);
    % REMOVED: Artificial fuel_fraction bounds - let physics determine actual values
    % fuel_fraction = max(0.01, min(0.8, fuel_fraction));

    % [kg] New takeoff gross weight from weight balance
    TOGW_new = real((OEW_new + mission.payload_kg + W_crew) / (1 - fuel_fraction));

    % Calculate convergence error
    % [-] Relative error between new and guess TOGW
    error_TOGW = abs(TOGW_new - TOGW_guess) / TOGW_guess;

    % Print iteration status
    % [-] Check if verbose output is enabled
    if opts.verbose
        % [-] Print iteration results in tabular format
        fprintf('%-5d | %-10.0f | %-10.2f | %-8.2f | %-10.1f | %-8.4f\n', ...
                iterationCount, TOGW_new, S_pln_new, LD_actual, ...
                wing_loading, error_TOGW);
    % [-] End verbose output check
    end

    % Update guesses with damping for stability
    % [-] Under-relaxation factor for convergence stability
    damping_factor = 0.7;
    % [kg] Update TOGW guess with under-relaxation
    TOGW_guess = damping_factor * TOGW_new + (1 - damping_factor) * TOGW_guess;
    % [m^2] Update planform area guess with under-relaxation
    S_pln_guess = damping_factor * S_pln_new + (1 - damping_factor) * S_pln_guess;
% [-] End main convergence loop
end

%% Check Convergence Status
% [-] Check if maximum iterations exceeded
if iterationCount >= maxIterations
    % [string] Create failure status message for iteration limit
    status = sprintf('Failed: Maximum iterations (%d) reached', maxIterations);
    % [struct] Return empty converged design structure
    convergedDesign = struct();
    % [-] Exit function with failure status
    return;
% [-] Convergence successful
else
    % [string] Set successful convergence status
    status = 'Converged';
% [-] End convergence status check
end

%% Package Converged Design
% [struct] Initialize converged design results structure
convergedDesign = struct();
% [kg] Store converged takeoff gross weight
convergedDesign.TOGW_kg = TOGW_guess;
% [kg] Store converged operating empty weight
convergedDesign.OEW_kg = OEW_new;
% [m^2] Store converged planform area
convergedDesign.S_pln_m2 = S_pln_guess;
% [m^2] Store CONSTANT wing reference area (17.544 m² for all vehicles)
convergedDesign.S_ref_m2 = vehicle.S_ref_m2;
% [-] Store mission fuel fraction
convergedDesign.fuel_fraction = fuel_fraction;
% [-] Store cruise lift-to-drag ratio
convergedDesign.LD_cruise = LD_actual;
% [kg/m^2] Calculate and store wing loading using CONSTANT reference area
convergedDesign.wing_loading = TOGW_guess / vehicle.S_ref_m2;
% [-] Store required thrust-to-weight ratio
convergedDesign.T_W_required = T_W_required;
% [-] Store convergence iteration count
convergedDesign.iterations = iterationCount;
% [-] Store final convergence error
convergedDesign.error_final = error_TOGW;
% REMOVED: bounds_violated tracking (no longer using artificial bounds)
% convergedDesign.bounds_violated = bounds_violated;
% convergedDesign.OEW_unbounded = x_sol(1) * 1000;
% convergedDesign.S_pln_unbounded = x_sol(2) * 10;

% Add structural weight data with detailed component breakdown
% [kg] Calculate final structural weight using converged parameters
S_wet_final = K_w * S_pln_guess;
W_structural_final = I_str * S_wet_final;

% FORCE USE OF I_STR METHOD with Roskam component weight fractions
% Total structural weight from I_str calibration
convergedDesign.W_structural_kg = W_structural_final;
convergedDesign.structural_method = 'I_str_with_Roskam_fractions';

% Apply Roskam's weight fractions to TOGW for detailed component breakdown
% These fractions are from Roskam's methods and represent % of TOGW:
% Wing:         9.82% of TOGW
% Fuselage:    10.716% of TOGW
% Vertical tail: 1.4591% of TOGW
% Canard:       1.5965% of TOGW
% Empennage:    3.0556% of TOGW (total tail surfaces)
% Landing gear: 3.9498% of TOGW
% Nacelle:      4.9438% of TOGW

convergedDesign.W_wing_kg = 0.0982 * TOGW_guess;
convergedDesign.W_fuselage_kg = 0.10716 * TOGW_guess;
convergedDesign.W_vtail_kg = 0.014591 * TOGW_guess;
convergedDesign.W_canard_kg = 0.015965 * TOGW_guess;
convergedDesign.W_empennage_kg = 0.030556 * TOGW_guess;  % Total empennage
convergedDesign.W_gear_kg = 0.039498 * TOGW_guess;
convergedDesign.W_nacelles_kg = 0.049438 * TOGW_guess;

% [kg/m^2] Store structural index used
convergedDesign.I_str_kg_m2 = I_str;

% Add engine performance data
% [-] Check if engine performance data is available
if exist('engine_data', 'var')
    % [struct] Store engine performance data
    convergedDesign.engine_performance = engine_data;
    % [N] Store available thrust
    convergedDesign.thrust_available = thrust_available;
    % [kg/s/N] Store actual TSFC
    convergedDesign.TSFC_actual = TSFC_actual;
    % [kg] Store engine weight
    convergedDesign.W_engine_kg = 1036; % GE F404 fixed weight
% [-] End engine data check
end

% Add systems weight
% [kg] Calculate systems weight from final OEW
W_systems_final = vehicle.C_sys_kg + vehicle.f_sys * OEW_new;
% [kg] Store systems weight
convergedDesign.W_systems_kg = W_systems_final;

% Add geometry and aerodynamic details
% [struct] Store geometry data
convergedDesign.geometry = geom;
% [-] Store cruise lift coefficient
convergedDesign.CL_cruise = CL;
% [-] Store cruise drag coefficient
convergedDesign.CD_cruise = CD;

% Add volume constraint diagnostics
% Calculate final volume budget status
% [m^3] Total vehicle volume
V_total_final = vehicle.tau * S_pln_guess^1.5;
% [m^3] Fuel volume
V_fuel_final = convergedDesign.TOGW_kg * fuel_fraction / vehicle.rho_fuel_kg_m3;
% [m^3] Payload volume
V_payload_final = mission.payload_kg / vehicle.rho_pay_kg_m3;
% [m^3] Crew volume
V_crew_final = mission.crew_count * 2.0;
% [m^3] Systems volume
V_systems_final = vehicle.k_vs * V_total_final;
% [m^3] Engine volume
V_engine_final = vehicle.k_ve * 1036 / 1000;
% [m^3] Void volume
V_void_final = vehicle.k_vv * V_total_final;
% [m^3] Fixed volume
V_fixed_final = vehicle.v_fix_m3;
% [m^3] Total volume used
V_used_final = V_fuel_final + V_payload_final + V_crew_final + V_systems_final + V_engine_final + V_void_final + V_fixed_final;
% [m^3] Available volume remaining
V_available_final = V_total_final - V_used_final;

% Store volume diagnostics
convergedDesign.V_total_m3 = V_total_final;
convergedDesign.V_used_m3 = V_used_final;
convergedDesign.V_available_m3 = V_available_final;
convergedDesign.V_fuel_m3 = V_fuel_final;
convergedDesign.volume_constraint_active = (V_available_final < 0.1);  % Active if <0.1 m³ margin

% [-] Check if verbose output is enabled
if opts.verbose
    % [-] Print convergence success message
    fprintf('\n✓ Convergence achieved in %d iterations\n', iterationCount);
    % [-] Print final design summary
    fprintf('Final Design: TOGW=%.0f kg, S_pln=%.1f m², L/D=%.2f\n', ...
            convergedDesign.TOGW_kg, convergedDesign.S_pln_m2, convergedDesign.LD_cruise);
% [-] End verbose output check
end


% [-] End main function
end

%% Enhanced Budget Equations Function
function F = enhanced_budget_equations(x, mission, vehicle, T_W, WR, I_str, K_w, LD_actual, TOGW_guess, opts)
    % Extract variables (scaled inputs)
    % [kg] Operating empty weight scaled back up from solver
    OEW = x(1) * 1000;
    % [m^2] Planform area scaled back up from solver
    S_pln = x(2) * 10;

    % Soft bounds checking with gradual penalties
    % [-] Initialize penalty accumulator
    penalty = 0;

    % OEW bounds with gradual penalty
    % [-] Check if OEW below minimum bound
    if OEW < 2000
        % [-] Add penalty for OEW below minimum
        penalty = penalty + (2000 - OEW)^2 / 1000;
    % [-] Check if OEW above maximum bound
    elseif OEW > 15000
        % [-] Add penalty for OEW above maximum
        penalty = penalty + (OEW - 15000)^2 / 1000;
    % [-] End OEW bounds check
    end

    % S_pln bounds with gradual penalty
    % [-] Check if planform area below minimum bound
    if S_pln < 15
        % [-] Add penalty for planform area below minimum
        penalty = penalty + (15 - S_pln)^2 / 10;
    % [-] Check if planform area above maximum bound
    elseif S_pln > 80
        % [-] Add penalty for planform area above maximum
        penalty = penalty + (S_pln - 80)^2 / 10;
    % [-] End planform area bounds check
    end

    % REMOVED: Penalty factor soft constraints - let fsolve find natural solutions
    % Always set penalty_factor to 0 (no artificial penalties)
    penalty_factor = 0;

    % Calculate total weight components
    % [kg] Crew weight estimate
    W_crew = mission.crew_count * 90;
    % [-] Fuel fraction from weight ratio
    fuel_fraction = (WR - 1) / WR;
    % [kg] Current TOGW from weight balance equation
    TOGW_current = (OEW + mission.payload_kg + W_crew) / (1 - fuel_fraction);

    % Woew - Weight from Empty Weight Budget
    % Structural weight
    % [m^2] Wetted area from wetted area ratio and planform area
    S_wet = K_w * S_pln;
    % [kg] Structural weight from index and wetted area
    W_structural = I_str * S_wet;

    % Engine weight - fixed for GE F404
    % [kg] GE F404 engine weight (fixed value)
    W_engine = 1036;

    % Systems weight
    % [kg] Systems weight from constant and variable components
    W_systems = vehicle.C_sys_kg + vehicle.f_sys * OEW;

    % Apply safety margin
    % [-] Safety margin factor from vehicle configuration
    mu_a = vehicle.mu_a;

    % Weight from empty weight budget
    % [kg] Total weight from empty weight budget with margin
    Woew = (W_structural + W_engine + W_systems) * (1 + mu_a);

    % Wowe - Weight from Volume Budget
    % Get total vehicle volume with consistent calculation
    % [m^3] Total vehicle volume using Kuechemann's slenderness parameter
    V_total = vehicle.tau * S_pln^1.5;

    % Calculate required volumes
    % [m^3] Fuel volume from weight and density
    V_fuel = TOGW_current * fuel_fraction / vehicle.rho_fuel_kg_m3;
    % [m^3] Payload volume from weight and density
    V_payload = mission.payload_kg / vehicle.rho_pay_kg_m3;
    % [m^3] Cockpit volume per crew member
    V_crew = mission.crew_count * 2.0;
    % [m^3] Systems volume as fraction of total volume
    V_systems = vehicle.k_vs * V_total;
    % [m^3] Engine volume from weight-based calculation
    V_engine = vehicle.k_ve * W_engine / 1000;
    % [m^3] Void volume as fraction of total volume
    V_void = vehicle.k_vv * V_total;
    % [m^3] Fixed volume for essential equipment
    V_fixed = vehicle.v_fix_m3;

    % [m^3] Total volume used by all components
    V_used = V_fuel + V_payload + V_crew + V_systems + V_engine + V_void + V_fixed;
    % [m^3] Available volume remaining
    V_available = V_total - V_used;

    % FIXED: Always calculate weight from volume budget (previous semester methodology)
    % This keeps the volume constraint ACTIVE even when V_available > 0
    % The if/else logic was causing F(2) = 0 when volume wasn't tight,
    % collapsing the dual-budget system to single-budget (weight only)
    % This created the artificial linear trend for medium-tau configs

    % Calculate maximum fuel volume (what's left after all other volumes)
    % [m^3] Maximum fuel volume available
    V_max_fuel = V_total - V_payload - V_crew - V_systems - V_engine - V_void - V_fixed;
    % [kg] Maximum fuel weight from available volume
    W_max_fuel = V_max_fuel * vehicle.rho_fuel_kg_m3;
    % [kg] Weight from volume budget - what TOGW would be if volume is exactly filled
    Wowe = (OEW + mission.payload_kg + W_crew + W_max_fuel);

    % Budget equations (from previous semester methodology)
    % F(1): Weight budget balance
    % [-] Weight budget equation with penalty
    F(1) = (Woew - OEW) + penalty_factor;

    % F(2): Volume budget balance - improved scaling
    % [kg] Weight difference for volume budget
    weight_diff = Wowe - TOGW_current;
    % [-] Volume budget equation with adaptive scaling and penalty
    F(2) = weight_diff / max(1000, abs(TOGW_current/10)) + penalty_factor;

    % Debug output if needed
    % [-] Debug output flag (disabled by default)
    if false
        % [-] Print budget equation debug information
        fprintf('  Budget Debug: OEW=%.0f, S_pln=%.1f, Woew=%.0f, Wowe=%.0f\n', ...
            OEW, S_pln, Woew, Wowe);
        % [-] Print equation values and available volume
        fprintf('    F(1)=%.3f, F(2)=%.3f, V_avail=%.1f\n', F(1), F(2), V_available);
    % [-] End debug output check
    end
% [-] End budget equations function
end

%% Supporting Functions
function T_W = calculate_thrust_constraints(wing_loading, M_cruise)
    % Calculate required T/W based on mission constraints
    % This is a simplified constraint analysis - can be enhanced

    % Takeoff constraint
    % [-] Minimum thrust-to-weight for reasonable takeoff performance
    T_W_takeoff = 0.3;

    % Supersonic cruise constraint
    % [-] Thrust-to-weight requirement increases with Mach number
    T_W_cruise = 0.4 + 0.1 * (M_cruise - 1.0);

    % Maneuverability constraint
    % [g] Maximum load factor for fighter aircraft
    n_max = 6.0;
    % [-] Thrust-to-weight for maneuverability (rough approximation)
    T_W_maneuver = (n_max * wing_loading / 5000) + 0.2;

    % Take most restrictive (highest) requirement
    % [-] Use maximum constraint value
    T_W = max([T_W_takeoff, T_W_cruise, T_W_maneuver]);
% [-] End thrust constraints function
end

function ff = calculate_fuel_fraction(LD, TSFC, t_cruise_min, M_cruise, TOGW_guess, opts)
    % Complete mission fuel fraction calculation including all flight phases
    % Now uses time-based TSFC method for cruise segment

    % [-] Check if options structure is provided
    if nargin < 5
        % [struct] Create default options structure
        opts = struct('units', 'SI');
    % [-] End input argument check
    end

    % Note: Old range-based calculation removed in favor of time-based method
    % Cruise conditions and atmospheric properties handled in time-based section below

    % Mission phase fuel fractions (from original compute_WR.m)
    % [-] Fuel fraction for engine start and warm-up
    W1_WTO = 0.990;

    % [-] Fuel fraction for taxi
    W2_W1 = 0.990;

    % [-] Fuel fraction for takeoff
    W3_W2 = 0.990;

    % [-] Fuel fraction for climb/acceleration to cruise altitude and speed
    W4_W3 = -0.036236*M_cruise + 0.998538;

    % [-] Calculate weight at start of cruise (W4)
    % W4 = TOGW - fuel burned during startup, taxi, takeoff, and climb
    W4_kg = TOGW_guess * W1_WTO * W2_W1 * W3_W2 * W4_W3;  % [kg]

    % [-] Fuel fraction for cruise segment using time-based TSFC method
    % This is more direct for fixed-time missions than Breguet range equation
    % F404 TSFC at cruise (M=1.4, 40k ft): ~1.84 lb/hr/lb
    TSFC_hr = 1.84;                    % [lb/hr/lb] Cruise TSFC
    TSFC_s = TSFC_hr / 3600;           % [lb/s/lb] Convert to per second
    cruise_time_s = t_cruise_min * 60;  % [s] Cruise time in seconds

    % Thrust required for steady cruise: T = Drag = Weight/L/D
    % Use weight at start of cruise (W4), not TOGW
    W_cruise_lbs = W4_kg * 2.20462;        % [lbs] Weight at cruise start
    T_req_lbs = W_cruise_lbs / LD;         % [lbf] Thrust = Drag at cruise

    % Cruise fuel consumed: m_fuel = TSFC × Thrust × time
    cruise_fuel_lbs = TSFC_s * T_req_lbs * cruise_time_s;  % [lbs]
    cruise_fuel_kg = cruise_fuel_lbs * 0.453592;           % [kg]

    % Cruise fuel fraction: W5/W4 = (W4 - fuel_burned) / W4
    W5_W4 = (W4_kg - cruise_fuel_kg) / W4_kg;

    % REMOVED: Artificial W5_W4 bounds - let physics determine actual values
    % W5_W4 = max(0.5, min(0.999, W5_W4));

    % [-] Fuel fraction for descent
    W6_W5 = 0.990;

    % [-] Fuel fraction for landing
    W7_W6 = 0.995;

    % [-] Calculate the total mission fuel fraction
    M_ff = W1_WTO * W2_W1 * W3_W2 * W4_W3 * W5_W4 * W6_W5 * W7_W6;

    % [-] Calculate mission fuel fraction
    mission_fuel = 1 - M_ff;

    % [-] Add 6% fuel reserve (fixed reserve as fraction of TOGW, per industry standard)
    % This is separate from mission fuel and represents reserves for:
    % - Holding patterns, diversions, weather delays
    % - Emergency situations
    % Note: 6% is typical for fighter aircraft per Raymer Ch. 3
    ff = mission_fuel + 0.06;

    % REMOVED: Artificial fuel fraction bounds - let physics determine actual values
    % ff = max(0.01, min(0.5, ff));

    % Debug output for fuel fraction calculation
    % [-] Debug output flag (disabled by default)
    if false
        % [-] Print fuel fraction debug information
        fprintf('  DEBUG FF: Time %.1f min, L/D=%.2f, TSFC=%.4f lb/hr/lb\n', t_cruise_min, LD, TSFC_hr);
        % [-] Print mission phase breakdown
        fprintf('    Mission phases: Start=%.3f, Taxi=%.3f, TO=%.3f, Climb=%.3f, Cruise=%.3f, Desc=%.3f, Land=%.3f\n', ...
                W1_WTO, W2_W1, W3_W2, W4_W3, W5_W4, W6_W5, W7_W6);
        % [-] Print total mission fuel fraction
        fprintf('    Total mission FF: M_ff=%.3f, Fuel fraction=%.3f\n', M_ff, ff);
    % [-] End debug output check
    end
% [-] End fuel fraction calculation function
end


function [I_str, K_w] = analyze_structure_enhanced(tau, S_pln, sweep_deg, TOGW_current, opts)
    % Returns structural parameters for X-29 class supersonic fighters
    % I_str: Fixed structural index (sweep effects handled in aerodynamics)
    % K_w: Wetted area ratio from geometry calculation

    % Get K_w from actual geometry
    % [struct] Get complete vehicle geometry using ELLIPTICAL FUSELAGE model
    geom = test_elliptical_fuselage(S_pln, tau, sweep_deg);
    % [-] Use real wetted area ratio from geometry
    K_w = geom.Kw;

    % Use fixed structural index
    % [kg/m^2] Structural index for baseline trade study
    I_str = 18.7;
end

function atm = get_atmosphere(h_altitude, units)
    % Enhanced atmosphere model with direct call to DD_atm_model

    % Check if units specified
    % [-] Check if units parameter provided
    if nargin < 2
        % [string] Default to SI units
        units = 'SI';
    % [-] End units check
    end

    % Use enhanced standard atmosphere directly (no external dependencies)
    % [struct] Calculate atmospheric properties using internal function
    atm = calculate_standard_atmosphere(h_altitude, units);
% [-] End atmosphere interface function
end

function atm = calculate_standard_atmosphere(h, units)
    % Calculate standard atmosphere properties

    % [-] Check if using SI units
    if strcmp(units, 'SI')
        % [m] Use altitude directly in meters
        h_m = h;
    % [-] Convert from imperial units
    else
        % [m] Convert altitude from feet to meters
        h_m = h * 0.3048;
    % [-] End units conversion
    end

    % Standard atmosphere calculation (SI units)
    % [-] Check if altitude is in troposphere
    if h_m <= 11000
        % [K] Temperature in troposphere with linear lapse rate
        T = 288.15 - 0.0065 * h_m;
        % [Pa] Pressure in troposphere
        P = 101325 * (T / 288.15)^5.256;
    % [-] Lower stratosphere (isothermal)
    else
        % [K] Constant temperature in lower stratosphere
        T = 216.65;
        % [Pa] Pressure in lower stratosphere with exponential decay
        P = 22632 * exp(-0.0001577 * (h_m - 11000));
    % [-] End atmospheric layer check
    end

    % Calculate other properties
    % [kg/m^3] Density from ideal gas law
    rho = P / (287.05 * T);
    % [m/s] Speed of sound from temperature
    a = sqrt(1.4 * 287.05 * T);
    % [Pa·s] Dynamic viscosity using Sutherland's law
    mu = 1.458e-6 * T^1.5 / (T + 110.4);
    % [m^2/s] Kinematic viscosity
    nu = mu / rho;

    % Ratios to sea level
    % [-] Temperature ratio to sea level
    theta = T / 288.15;
    % [-] Pressure ratio to sea level
    delta = P / 101325;
    % [-] Density ratio to sea level
    sigma = rho / 1.225;
    % [m] Geopotential height
    z = h_m;

    % Convert to requested units
    % [-] Check if imperial units requested
    if strcmp(units, 'EN') || strcmp(units, 'Imperial')
        % [R] Convert temperature to Rankine
        T = T * 1.8;
        % [lbf/ft^2] Convert pressure to imperial
        P = P * 0.0208854;
        % [slug/ft^3] Convert density to imperial
        rho = rho * 0.00194032;
        % [ft/s] Convert speed of sound to imperial
        a = a * 3.28084;
        % [lbf·s/ft^2] Convert dynamic viscosity to imperial
        mu = mu * 0.0208854;
        % [ft^2/s] Convert kinematic viscosity to imperial
        nu = nu * 10.7639;
        % [ft or m] Keep original height units
        z = h;
    % [-] End unit conversion check
    end

    % Package results
    % [struct] Initialize atmospheric properties structure
    atm = struct();
    % [m or ft] Altitude
    atm.h = h;
    % [K or R] Temperature
    atm.T = T;
    % [Pa or lbf/ft^2] Pressure
    atm.P = P;
    % [kg/m^3 or slug/ft^3] Density
    atm.rho = rho;
    % [m/s or ft/s] Speed of sound
    atm.a = a;
    % [Pa·s or lbf·s/ft^2] Dynamic viscosity
    atm.mu = mu;
    % [m^2/s or ft^2/s] Kinematic viscosity
    atm.nu = nu;
    % [-] Temperature ratio
    atm.theta = theta;
    % [-] Pressure ratio
    atm.delta = delta;
    % [-] Density ratio
    atm.sigma = sigma;
    % [m] Geopotential height
    atm.z = z;
% [-] End standard atmosphere calculation function
end
