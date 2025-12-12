%{
====================================================================================================
FUNCTION NAME: standard_atmosphere_1976.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 9/16/2025
LAST REVISION: 9/28/2025
====================================================================================================
FUNCTION DESCRIPTION:
This function implements the U.S. Standard Atmosphere 1976 model for atmospheric property
calculations in the MAE 4350 X-29 aircraft design project. It provides accurate atmospheric
conditions (temperature, pressure, density, viscosity, speed of sound) for altitudes from
sea level to 20 km, covering the troposphere and lower stratosphere regions relevant to
supersonic aircraft operations.

Key Features:
- Complete implementation of U.S. Standard Atmosphere 1976 model
- Support for both SI and Imperial unit systems
- Vectorized calculations for multiple altitudes
- Includes viscosity calculations using Sutherland's law
- Provides dimensionless ratios for aerodynamic calculations
- Robust error handling for out-of-range altitudes
- Optimized for aerospace vehicle design applications

Technical Methodology:
- Troposphere (0-11 km): Linear temperature lapse rate with hydrostatic equation
- Lower Stratosphere (11-20 km): Isothermal layer with exponential pressure decay
- Viscosity: Sutherland's law for temperature-dependent dynamic viscosity
- Speed of Sound: Calculated from temperature using ideal gas relations
- Unit Conversion: Automatic conversion between SI and Imperial systems
====================================================================================================
INPUTS:
h               - [m or ft] Geometric altitude (scalar or array)
                  Valid range: 0 ≤ h ≤ 20,000 m (0 ≤ h ≤ 65,617 ft)
                  Units depend on opts.units setting

opts            - [struct] Options structure (optional)
  .units        - [string] Unit system specification
                  'SI' (default): h in [m], outputs in SI units
                  'Imperial': h in [ft], outputs in Imperial units

OUTPUTS:
atm             - [struct] Atmospheric properties structure
  .h            - [m or ft] Geometric altitude (echo of input)
  .T            - [K or R] Temperature
  .P            - [Pa or lbf/ft²] Pressure
  .rho          - [kg/m³ or slug/ft³] Density
  .a            - [m/s or ft/s] Speed of sound
  .mu           - [Pa·s or lbf·s/ft²] Dynamic viscosity
  .nu           - [m²/s or ft²/s] Kinematic viscosity
  .T_ratio      - [-] Temperature ratio to sea level
  .P_ratio      - [-] Pressure ratio to sea level
  .rho_ratio    - [-] Density ratio to sea level
  .theta        - [-] Temperature ratio (common notation)
  .delta        - [-] Pressure ratio (common notation)
  .sigma        - [-] Density ratio (common notation)
  .T_sl         - [K or R] Sea level temperature (reference)
  .P_sl         - [Pa or lbf/ft²] Sea level pressure (reference)
  .rho_sl       - [kg/m³ or slug/ft³] Sea level density (reference)
  .a_sl         - [m/s or ft/s] Sea level speed of sound (reference)
  .Mach1        - [m/s or ft/s] Speed for Mach 1 (equals speed of sound)

DEPENDENCIES:
External Functions: None (self-contained implementation)

MATLAB Toolboxes Required: None (uses only base MATLAB functions)

REFERENCES:
- U.S. Standard Atmosphere, 1976, NASA-TM-X-74335, NOAA-S/T 76-1562
- Anderson, J. "Introduction to Flight" 8th Edition, McGraw-Hill, 2016
- Roskam, J. "Airplane Design Part I: Preliminary Sizing" DARcorporation, 1985
====================================================================================================
USAGE EXAMPLES:
% Single altitude in SI units
atm = standard_atmosphere_1976(10000);  % 10 km altitude

% Multiple altitudes in Imperial units
h_ft = [0, 10000, 20000, 30000];
atm = standard_atmosphere_1976(h_ft, struct('units','Imperial'));

% Cruise altitude for supersonic aircraft
atm_cruise = standard_atmosphere_1976(12200);  % X-29 cruise altitude
====================================================================================================
REVISION HISTORY:
9/16/2025 - D. Larin - Initial implementation with complete atmospheric model
9/20/2025 - D. Larin - Added Imperial units support and vectorization
9/24/2025 - D. Larin - Enhanced error handling and reference values
9/28/2025 - D. Larin - Comprehensive documentation and comment standardization
====================================================================================================
%}

function atm = standard_atmosphere_1976(h, opts)

    % [-] Check if options structure is provided
    if nargin < 2
        % [struct] Create default empty options structure
        opts = struct();
    % [-] End input argument check
    end
    % [-] Check if units field is specified in options
    if ~isfield(opts, 'units')
        % [string] Set default unit system to SI
        opts.units = 'SI';
    % [-] End units field check
    end

    % Physical constants (SI)
    % [m/s^2] Standard gravitational acceleration
    g0 = 9.80665;
    % [J/(kg·K)] Specific gas constant for dry air
    R = 287.05287;
    % [-] Ratio of specific heats for air
    gamma = 1.4;
    % [K] Sutherland's constant for viscosity calculations
    S = 110.4;
    % [kg/(m·s·K^0.5)] Viscosity coefficient for Sutherland's law
    beta_visc = 1.458e-6;

    % Sea level standard conditions (SI)
    % [K] Sea level standard temperature
    T0 = 288.15;
    % [Pa] Sea level standard pressure
    P0 = 101325;
    % [kg/m^3] Sea level standard density
    rho0 = 1.225;
    % [m/s] Sea level standard speed of sound
    a0 = 340.294;

    % Temperature gradient layers (SI)
    % Layer 1: Troposphere (0 to 11 km)
    % [m] Altitude of troposphere boundary
    h1 = 11000;
    % [K/m] Temperature lapse rate in troposphere
    L1 = -0.0065;

    % Layer 2: Lower Stratosphere (11 to 20 km)
    % [m] Altitude of model upper boundary
    h2 = 20000;
    % [K/m] Temperature lapse rate in stratosphere (isothermal)
    L2 = 0;

    % Convert input altitude to SI if needed
    % [-] Check if Imperial units are requested
    if strcmpi(opts.units, 'Imperial')
        % [m] Convert altitude from feet to meters
        h_m = h * 0.3048;
    % [-] Using SI units (meters)
    else
        % [m] Use altitude directly in meters
        h_m = h;
    % [-] End unit conversion check
    end

    % Ensure altitude is within valid range
    % [-] Check for negative altitudes
    if any(h_m < 0)
        % [-] Print warning for below sea level altitudes
        warning('Altitude below sea level. Using sea level values.');
        % [m] Clamp negative altitudes to sea level
        h_m(h_m < 0) = 0;
    % [-] End negative altitude check
    end
    % [-] Check for altitudes above model validity
    if any(h_m > 20000)
        % [-] Print warning for excessive altitudes
        warning('Altitude above 20 km. Model validity limited to 20 km.');
        % [m] Clamp excessive altitudes to model limit
        h_m(h_m > 20000) = 20000;
    % [-] End excessive altitude check
    end

    % Initialize output arrays
    % [K] Initialize temperature array
    T = zeros(size(h_m));
    % [Pa] Initialize pressure array
    P = zeros(size(h_m));
    % [kg/m^3] Initialize density array
    rho = zeros(size(h_m));

    % Troposphere (0 ≤ h < 11 km)
    % [-] Create logical index for troposphere altitudes
    idx1 = h_m < h1;
    % [-] Check if any altitudes are in troposphere
    if any(idx1)
        % [K] Calculate temperature with linear lapse rate
        T(idx1) = T0 + L1 * h_m(idx1);
        % [Pa] Calculate pressure using hydrostatic equation
        P(idx1) = P0 * (T(idx1) / T0).^(-g0 / (L1 * R));
        % [kg/m^3] Calculate density from ideal gas law
        rho(idx1) = P(idx1) ./ (R * T(idx1));
    % [-] End troposphere calculation
    end

    % Calculate conditions at 11 km for continuity
    % [K] Temperature at troposphere boundary
    T11 = T0 + L1 * h1;
    % [Pa] Pressure at troposphere boundary
    P11 = P0 * (T11 / T0)^(-g0 / (L1 * R));
    % [kg/m^3] Density at troposphere boundary
    rho11 = P11 / (R * T11);

    % Lower Stratosphere (11 km ≤ h ≤ 20 km) - Isothermal
    % [-] Create logical index for stratosphere altitudes
    idx2 = (h_m >= h1) & (h_m <= h2);
    % [-] Check if any altitudes are in stratosphere
    if any(idx2)
        % [K] Constant temperature in isothermal layer
        T(idx2) = T11;
        % [Pa] Calculate pressure with exponential decay
        P(idx2) = P11 * exp(-g0 * (h_m(idx2) - h1) / (R * T11));
        % [kg/m^3] Calculate density from ideal gas law
        rho(idx2) = P(idx2) / (R * T11);
    % [-] End stratosphere calculation
    end

    % Calculate derived properties
    % [m/s] Speed of sound from temperature
    a = sqrt(gamma * R * T);
    % [Pa·s] Dynamic viscosity using Sutherland's law
    mu = beta_visc * T.^1.5 ./ (T + S);
    % [m^2/s] Kinematic viscosity from dynamic viscosity and density
    nu = mu ./ rho;

    % Dimensionless ratios
    % [-] Temperature ratio to sea level standard
    T_ratio = T / T0;
    % [-] Pressure ratio to sea level standard
    P_ratio = P / P0;
    % [-] Density ratio to sea level standard
    rho_ratio = rho / rho0;

    % Build output structure in SI
    % [struct] Initialize atmospheric properties structure
    atm = struct();
    % [m] Store geometric altitude
    atm.h = h_m;
    % [K] Store temperature
    atm.T = T;
    % [Pa] Store pressure
    atm.P = P;
    % [kg/m^3] Store density
    atm.rho = rho;
    % [m/s] Store speed of sound
    atm.a = a;
    % [Pa·s] Store dynamic viscosity
    atm.mu = mu;
    % [m^2/s] Store kinematic viscosity
    atm.nu = nu;
    % [-] Store temperature ratio
    atm.T_ratio = T_ratio;
    % [-] Store pressure ratio
    atm.P_ratio = P_ratio;
    % [-] Store density ratio
    atm.rho_ratio = rho_ratio;
    % [-] Temperature ratio in common notation
    atm.theta = T_ratio;
    % [-] Pressure ratio in common notation
    atm.delta = P_ratio;
    % [-] Density ratio in common notation
    atm.sigma = rho_ratio;

    % Convert to Imperial units if requested
    % [-] Check if Imperial units are requested
    if strcmpi(opts.units, 'Imperial')
        % [ft] Keep original input altitude in feet
        atm.h = h;
        % [R] Convert temperature from Kelvin to Rankine
        atm.T = T * 1.8;
        % [lbf/ft^2] Convert pressure from Pa to lbf/ft²
        atm.P = P * 0.020885;
        % [slug/ft^3] Convert density from kg/m³ to slug/ft³
        atm.rho = rho * 0.00194032;
        % [ft/s] Convert speed of sound from m/s to ft/s
        atm.a = a * 3.28084;
        % [lbf·s/ft^2] Convert dynamic viscosity to Imperial units
        atm.mu = mu * 0.020885;
        % [ft^2/s] Convert kinematic viscosity from m²/s to ft²/s
        atm.nu = nu * 10.7639;
    % [-] End Imperial unit conversion
    end

    % Add sea level reference values for convenience
    % [-] Check if Imperial units for reference values
    if strcmpi(opts.units, 'Imperial')
        % [R] Sea level temperature in Rankine
        atm.T_sl = T0 * 1.8;
        % [lbf/ft^2] Sea level pressure in Imperial units
        atm.P_sl = P0 * 0.020885;
        % [slug/ft^3] Sea level density in Imperial units
        atm.rho_sl = rho0 * 0.00194032;
        % [ft/s] Sea level speed of sound in feet per second
        atm.a_sl = a0 * 3.28084;
    % [-] Using SI units for reference values
    else
        % [K] Sea level temperature in Kelvin
        atm.T_sl = T0;
        % [Pa] Sea level pressure in Pascals
        atm.P_sl = P0;
        % [kg/m^3] Sea level density in SI units
        atm.rho_sl = rho0;
        % [m/s] Sea level speed of sound in meters per second
        atm.a_sl = a0;
    % [-] End reference values unit check
    end

    % Add Mach 1 reference
    % [m/s or ft/s] Speed for Mach 1 (equals local speed of sound)
    atm.Mach1 = atm.a;
% [-] End standard_atmosphere_1976 function
end