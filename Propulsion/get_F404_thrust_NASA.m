%{
====================================================================================================
FUNCTION NAME: get_F404_thrust_NASA.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 11/26/2025
LAST REVISION: 11/26/2025
====================================================================================================
FUNCTION DESCRIPTION:
Empirical F404-GE-400 thrust model based on NASA TM-4140 flight test data from the X-29A.
"Measurement Effects on the Calculation of In-Flight Thrust for an F404 Turbofan Engine"
Timothy R. Conners, NASA Ames Research Center, Dryden Flight Research Facility, 1989.

This model uses multi-dimensional interpolation of actual flight test data to provide
accurate thrust estimates across the F404 operating envelope.

DATA SOURCE: NASA TM-4140, Table 2 - Net thrust values (WT method)
Flight Conditions Tested:
  - 10,000 ft, M = 0.4 and 0.8
  - 30,000 ft, M = 0.9 and 1.2
  - 40,000 ft, M = 0.8 and 1.6

PLA Settings:
  - 70°  = Part power (~60% thrust)
  - 87°  = Intermediate/Military power (max dry thrust)
  - 92.5° = Minimum afterburner
  - 109° = Mid afterburner
  - 130° = Maximum afterburner

F404-GE-400 Specifications (for reference):
  - Max AB thrust (SL): 16,000 lbf (71.2 kN)
  - Max MIL thrust (SL): 10,600 lbf (47.1 kN)
  - SFC Max AB: 1.84 lb/hr/lb
  - SFC Max MIL: 0.78 lb/hr/lb
====================================================================================================
INPUTS:
M               - [-] Free-stream Mach number (scalar)
                  Valid range: 0.3 ≤ M ≤ 2.0

h_ft            - [ft] Geometric altitude (scalar)
                  Valid range: 0 ≤ h ≤ 50,000 ft

PLA_deg         - [deg] Power Lever Angle (scalar, optional)
                  Valid range: 70° ≤ PLA ≤ 130°
                  Default: 87° (military power)
                  87° = Military power (max non-afterburner)
                  >87° = Afterburner engagement

OUTPUTS:
prop            - [struct] F404 performance results structure
  .thrust_lbf       - [lbf] Net thrust (Imperial)
  .thrust_N         - [N] Net thrust (SI)
  .thrust_kN        - [kN] Net thrust (SI)
  .TSFC_SI          - [kg/s/N] Thrust specific fuel consumption (SI)
  .TSFC_imperial    - [lb/hr/lb] Thrust specific fuel consumption (Imperial)
  .fuel_flow_kg_s   - [kg/s] Fuel mass flow rate
  .fuel_flow_lb_hr  - [lb/hr] Fuel mass flow rate
  .conditions       - [struct] Input conditions
====================================================================================================
%}

function prop = get_F404_thrust_NASA(M, h_ft, PLA_deg)

% Default to military power if PLA not specified
if nargin < 3 || isempty(PLA_deg)
    PLA_deg = 87;
end

% Clamp inputs to valid ranges
M = max(0.3, min(2.0, M));
h_ft = max(0, min(50000, h_ft));
PLA_deg = max(70, min(130, PLA_deg));

% =========================================================================
% NASA TM-4140 FLIGHT TEST DATA (WT Method - More Accurate)
% =========================================================================
% Net thrust in lbf from Table 2

% Flight condition grid points
alt_data = [10000, 10000, 30000, 30000, 40000, 40000];  % ft
mach_data = [0.4, 0.8, 0.9, 1.2, 0.8, 1.6];

% PLA settings
PLA_data = [70, 87, 92.5, 109, 130];  % degrees

% Thrust data matrix [6 conditions x 5 PLA settings] in lbf
% Rows: (10k,0.4), (10k,0.8), (30k,0.9), (30k,1.2), (40k,0.8), (40k,1.6)
% Cols: PLA 70, 87, 92.5, 109, 130
thrust_data = [
    4904, 7318, 7956, 9392, 11664;   % 10,000 ft, M=0.4
    4724, 7928, 8750, 10564, 13892;  % 10,000 ft, M=0.8
    3395, 4547, 5061, 6261, 8022;    % 30,000 ft, M=0.9
    5475, 5478, 6174, 7766, 10477;   % 30,000 ft, M=1.2
    2348, 2675, 2954, 3639, 4654;    % 40,000 ft, M=0.8
    4251, 4252, 4872, 6273, 9012;    % 40,000 ft, M=1.6
];

% =========================================================================
% INTERPOLATION STRATEGY
% =========================================================================
% Since we have sparse data points, use weighted distance interpolation
% For each query point (M, h_ft, PLA), find nearby data points and interpolate

% First, interpolate in PLA for each flight condition
thrust_at_PLA = zeros(6, 1);
for i = 1:6
    thrust_at_PLA(i) = interp1(PLA_data, thrust_data(i,:), PLA_deg, 'pchip', 'extrap');
end

% Now interpolate in altitude-Mach space
% Use inverse distance weighting with the 6 data points

% Calculate distances in normalized (alt, Mach) space
% Normalize altitude by 50000 ft, Mach by 2.0
alt_norm = h_ft / 50000;
mach_norm = M / 2.0;

alt_data_norm = alt_data / 50000;
mach_data_norm = mach_data / 2.0;

% Calculate distances
distances = sqrt((alt_data_norm - alt_norm).^2 + (mach_data_norm - mach_norm).^2);

% Check if we're very close to a data point
min_dist = min(distances);
if min_dist < 0.001
    % Use the closest point directly
    [~, idx] = min(distances);
    thrust_lbf = thrust_at_PLA(idx);
else
    % Inverse distance weighting
    % Use power of 2 for smoother interpolation
    weights = 1 ./ (distances.^2 + 0.001);  % Small epsilon to avoid division by zero
    weights = weights / sum(weights);
    thrust_lbf = sum(weights .* thrust_at_PLA');
end

% =========================================================================
% ALTITUDE CORRECTION FOR OUT-OF-RANGE CONDITIONS
% =========================================================================
% The data covers 10,000-40,000 ft. For sea level and higher altitudes,
% apply physics-based correction using density ratio

[~, ~, ~, rho_query] = atmosisa(h_ft * 0.3048);  % Convert to meters for atmosisa

% For altitudes below 10,000 ft, apply sea level correction
% This corrects for the fact that NASA data doesn't include sea level
% F404 specs: MIL = 10,600 lbf, Max AB = 16,000 lbf at SLS
if h_ft < 10000
    [~, ~, ~, rho_10k] = atmosisa(10000 * 0.3048);

    % First apply density correction
    density_correction = (rho_query / rho_10k)^0.7;  % Empirical exponent
    thrust_lbf = thrust_lbf * density_correction;

    % Then apply sea level calibration correction
    % This tapers from full correction at SL to no correction at 10,000 ft
    alt_factor = 1 - h_ft / 10000;  % 1 at SL, 0 at 10kft

    % Correction factors to match F404 specs at sea level
    % MIL: model gives ~8528 lbf, spec is 10,600 lbf -> factor = 1.243
    % AB:  model gives ~14206 lbf, spec is 16,000 lbf -> factor = 1.126
    if PLA_deg <= 87
        spec_correction = 1.243;
    else
        % Interpolate correction based on PLA
        AB_fraction = (PLA_deg - 87) / (130 - 87);
        spec_correction = 1.243 - AB_fraction * (1.243 - 1.126);
    end

    % Apply tapered correction
    correction = 1 + alt_factor * (spec_correction - 1);
    thrust_lbf = thrust_lbf * correction;
end

% For altitudes above 40,000 ft, scale down thrust
if h_ft > 40000
    [~, ~, ~, rho_40k] = atmosisa(40000 * 0.3048);
    density_correction = (rho_query / rho_40k)^0.7;
    thrust_lbf = thrust_lbf * density_correction;
end

% =========================================================================
% MACH NUMBER CORRECTION FOR OUT-OF-RANGE
% =========================================================================
% Data covers M=0.4 to M=1.6. Apply ram effect correction outside this range

if M < 0.4
    % Below M=0.4, thrust decreases slightly (less ram effect)
    mach_correction = 0.95 + 0.05 * (M / 0.4);
    thrust_lbf = thrust_lbf * mach_correction;
elseif M > 1.6
    % Above M=1.6, ram effect continues to help but diminishing returns
    % Use a mild correction based on ram pressure ratio
    ram_ratio_16 = 1 + 0.2 * 1.6^2;
    ram_ratio_M = 1 + 0.2 * M^2;
    mach_correction = (ram_ratio_M / ram_ratio_16)^0.3;  % Mild increase
    thrust_lbf = thrust_lbf * mach_correction;
end

% =========================================================================
% TSFC CALCULATION
% =========================================================================
% F404 documented TSFC values:
% - Military power (PLA ≤ 87°): 0.78 lb/hr/lb
% - Max Afterburner (PLA = 130°): 1.84 lb/hr/lb
% Interpolate between these based on PLA

if PLA_deg <= 87
    % Dry operation - TSFC varies slightly with throttle
    TSFC_imperial = 0.78 * (0.9 + 0.1 * (PLA_deg / 87));
else
    % Afterburner operation - TSFC increases with AB level
    % Linear interpolation from 87° to 130°
    AB_fraction = (PLA_deg - 87) / (130 - 87);
    TSFC_imperial = 0.78 + AB_fraction * (1.84 - 0.78);
end

% Convert TSFC to SI: lb/hr/lb -> kg/s/N
% 1 lb/hr/lb = 2.832545e-5 kg/s/N
TSFC_SI = TSFC_imperial * 2.832545e-5;

% =========================================================================
% UNIT CONVERSIONS AND OUTPUT
% =========================================================================

% Convert thrust to SI
thrust_N = thrust_lbf * 4.448222;
thrust_kN = thrust_N / 1000;

% Calculate fuel flow
fuel_flow_lb_hr = TSFC_imperial * thrust_lbf;
fuel_flow_kg_s = TSFC_SI * thrust_N;

% Package output
prop = struct();
prop.thrust_lbf = thrust_lbf;
prop.thrust_N = thrust_N;
prop.thrust_kN = thrust_kN;
prop.TSFC_SI = TSFC_SI;
prop.TSFC_imperial = TSFC_imperial;
prop.fuel_flow_kg_s = fuel_flow_kg_s;
prop.fuel_flow_lb_hr = fuel_flow_lb_hr;

% Store conditions
prop.conditions.Mach = M;
prop.conditions.altitude_ft = h_ft;
prop.conditions.altitude_m = h_ft * 0.3048;
prop.conditions.PLA_deg = PLA_deg;
prop.conditions.source = 'NASA TM-4140 (1989)';

end
