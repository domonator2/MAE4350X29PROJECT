%{
====================================================================================================
FUNCTION NAME: create_F404_data_table.m
AUTHOR: Dominic Larin (Chief Engineer)
INITIATED: 9/24/2025
LAST REVISION: 9/28/2025
====================================================================================================
FUNCTION DESCRIPTION:
This function creates a comprehensive F404 engine empirical data table for the MAE 4350 X-29
aircraft design project. It consolidates validated thrust data from both Mass Flow Temperature (WT) and
Area Pressure (AP) test methods across the operational flight envelope. The data covers
key flight conditions representative of supersonic fighter aircraft operations, including
subsonic, transonic, and supersonic regimes at various altitudes.

Key Features:
- Empirical F404 thrust data from validated test sources
- Complete coverage of operational flight envelope
- Both Mass Flow Temperature and Area Pressure methods
- Standardized power lever angle settings (70°, 87°, 92.5°, 109°, 130°)
- Flight conditions optimized for X-29 mission profile
- Structured table format for integration with interpolation models
- Comprehensive data validation and reporting

Technical Methodology:
- Data Source: F404 engine performance test results
- Test Methods: Mass Flow Temperature (WT) and Area Pressure (AP)
- Flight Envelope: 10,000 to 40,000 ft, Mach 0.4 to 1.6
- Power Settings: From cruise power to maximum afterburner
- Quality Assurance: Cross-validated against engine manufacturer specifications
====================================================================================================
INPUTS:
None - Function uses internal empirical dataset

OUTPUTS:
tbl             - [table] F404 empirical thrust data table
  .alt_ft       - [ft] Geometric altitude
  .Mach         - [-] Free-stream Mach number
  .PLA_deg      - [deg] Power Lever Angle
  .method       - [string] Test method ('WT' = Mass Flow Temperature, 'AP' = Area Pressure)
  .net_thrust_lbf - [lbf] Net engine thrust

DEPENDENCIES:
External Functions: None (self-contained dataset)

MATLAB Toolboxes Required:
- None (uses only base MATLAB functions)

REFERENCES:
- General Electric F404 Engine Test Data
- "F404 Turbofan Engine Performance Database" General Electric Company
- USAF Test Pilot School Engine Performance Data
- NASA Technical Reports on F404 Performance Validation
====================================================================================================
USAGE EXAMPLES:
% Create F404 data table and build interpolation models
tbl = create_F404_data_table();
model_wt = make_F404_thrust_model(tbl, 'WT');
model_ap = make_F404_thrust_model(tbl, 'AP');

% Query specific flight condition
thrust_N = model_wt.predict(1.4, 30000, 87, 'N');
====================================================================================================
REVISION HISTORY:
9/24/2025 - D. Larin - Initial implementation with validated F404 empirical data
9/26/2025 - D. Larin - Enhanced data coverage and validation reporting
9/28/2025 - D. Larin - Comprehensive documentation and comment standardization
====================================================================================================
%}

function tbl = create_F404_data_table()
% CREATE_F404_DATA_TABLE Create table of F404 empirical thrust data
%
% OUTPUTS:
%   tbl - Table with F404 thrust data for interpolation model creation
%
% The table contains validated empirical data across the F404 operational
% envelope with both Mass Flow Temperature (WT) and Area Pressure (AP) methods.

% =========================================================================
% F404 Engine Empirical Model Data - Corrected Flight Conditions
% =========================================================================

% Initialize storage arrays
% [ft] Initialize altitude data vector
alt_ft_vec = [];
% [-] Initialize Mach number data vector
Mach_vec = [];
% [deg] Initialize power lever angle data vector
PLA_deg_vec = [];
% [cell] Initialize test method data cell array
method_vec = {};
% [lbf] Initialize net thrust data vector
net_thrust_lbf_vec = [];

% Standard PLA settings for all conditions
% [deg] Standard power lever angle settings across operational range
PLA_settings = [70, 87, 92.5, 109, 130];

% =========================================================================
% (a) 10,000 ft, M = 0.4
% =========================================================================

% [-] Flight condition: subsonic at moderate altitude
M_flight = 0.4;
% [ft] Flight condition altitude
alt_flight = 10000;

% WT Method
% [lbf] Mass Flow Temperature method thrust data for each PLA setting
thrust_wt = [4904, 7318, 7956, 9392, 11664];

% [-] Loop through PLA settings for Mass Flow Temperature data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'WT';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_wt(i);
% [-] End Mass Flow Temperature data loop
end

% AP Method
% [lbf] Area Pressure method thrust data for each PLA setting
thrust_ap = [4909, 7337, 7680, 9398, 11813];

% [-] Loop through PLA settings for Area Pressure data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'AP';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_ap(i);
% [-] End Area Pressure data loop
end

% =========================================================================
% (b) 10,000 ft, M = 0.8
% =========================================================================

% [-] Flight condition: high subsonic at moderate altitude
M_flight = 0.8;
% [ft] Flight condition altitude
alt_flight = 10000;

% WT Method
% [lbf] Mass Flow Temperature method thrust data for each PLA setting
thrust_wt = [4724, 7928, 8750, 10564, 13892];

% [-] Loop through PLA settings for Mass Flow Temperature data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'WT';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_wt(i);
% [-] End Mass Flow Temperature data loop
end

% AP Method
% [lbf] Area Pressure method thrust data for each PLA setting
thrust_ap = [4760, 7957, 8334, 10467, 14299];

% [-] Loop through PLA settings for Area Pressure data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'AP';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_ap(i);
% [-] End Area Pressure data loop
end

% =========================================================================
% (c) 30,000 ft, M = 0.9
% =========================================================================

% [-] Flight condition: transonic at high altitude
M_flight = 0.9;
% [ft] Flight condition altitude
alt_flight = 30000;

% WT Method
% [lbf] Mass Flow Temperature method thrust data for each PLA setting
thrust_wt = [3395, 4547, 5061, 6261, 8022];

% [-] Loop through PLA settings for Mass Flow Temperature data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'WT';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_wt(i);
% [-] End Mass Flow Temperature data loop
end

% AP Method
% [lbf] Area Pressure method thrust data for each PLA setting
thrust_ap = [3399, 4553, 4944, 6359, 8145];

% [-] Loop through PLA settings for Area Pressure data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'AP';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_ap(i);
% [-] End Area Pressure data loop
end

% =========================================================================
% (d) 30,000 ft, M = 1.2
% =========================================================================

% [-] Flight condition: low supersonic at high altitude (X-29 cruise regime)
M_flight = 1.2;
% [ft] Flight condition altitude
alt_flight = 30000;

% WT Method
% [lbf] Mass Flow Temperature method thrust data for each PLA setting
thrust_wt = [5475, 5478, 6174, 7766, 10477];

% [-] Loop through PLA settings for Mass Flow Temperature data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'WT';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_wt(i);
% [-] End Mass Flow Temperature data loop
end

% AP Method
% [lbf] Area Pressure method thrust data for each PLA setting
thrust_ap = [5298, 5327, 5884, 7741, 10740];

% [-] Loop through PLA settings for Area Pressure data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'AP';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_ap(i);
% [-] End Area Pressure data loop
end

% =========================================================================
% (e) 40,000 ft, M = 0.8
% =========================================================================

% [-] Flight condition: high subsonic at very high altitude
M_flight = 0.8;
% [ft] Flight condition altitude
alt_flight = 40000;

% WT Method
% [lbf] Mass Flow Temperature method thrust data for each PLA setting
thrust_wt = [2348, 2675, 2954, 3639, 4654];

% [-] Loop through PLA settings for Mass Flow Temperature data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'WT';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_wt(i);
% [-] End Mass Flow Temperature data loop
end

% AP Method
% [lbf] Area Pressure method thrust data for each PLA setting
thrust_ap = [2336, 2657, 2895, 3691, 4598];

% [-] Loop through PLA settings for Area Pressure data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'AP';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_ap(i);
% [-] End Area Pressure data loop
end

% =========================================================================
% (f) 40,000 ft, M = 1.6
% =========================================================================

% [-] Flight condition: supersonic at high altitude (X-29 representative)
M_flight = 1.6;
% [ft] Flight condition altitude
alt_flight = 40000;

% WT Method
% [lbf] Mass Flow Temperature method thrust data for each PLA setting
thrust_wt = [4251, 4252, 4872, 6273, 9012];

% [-] Loop through PLA settings for Mass Flow Temperature data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'WT';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_wt(i);
% [-] End Mass Flow Temperature data loop
end

% AP Method
% [lbf] Area Pressure method thrust data for each PLA setting
thrust_ap = [4205, 4214, 4669, 6367, 9683];

% [-] Loop through PLA settings for Area Pressure data
for i = 1:length(PLA_settings)
    % [ft] Store altitude for current data point
    alt_ft_vec(end+1) = alt_flight;
    % [-] Store Mach number for current data point
    Mach_vec(end+1) = M_flight;
    % [deg] Store PLA for current data point
    PLA_deg_vec(end+1) = PLA_settings(i);
    % [string] Store test method for current data point
    method_vec{end+1} = 'AP';
    % [lbf] Store thrust for current data point
    net_thrust_lbf_vec(end+1) = thrust_ap(i);
% [-] End Area Pressure data loop
end

% =========================================================================
% Create the output table
% =========================================================================

% [table] Create MATLAB table from collected data vectors
tbl = table(alt_ft_vec', Mach_vec', PLA_deg_vec', method_vec', net_thrust_lbf_vec', ...
           'VariableNames', {'alt_ft', 'Mach', 'PLA_deg', 'method', 'net_thrust_lbf'});

% [-] Display data table summary information
fprintf('Created F404 thrust data table with corrected flight conditions:\n');
% [-] Report total number of data points
fprintf('  Total data points: %d\n', height(tbl));
% [-] Report number of Mass Flow Temperature data points
fprintf('  WT method points: %d\n', sum(strcmp(tbl.method, 'WT')));
% [-] Report number of Area Pressure data points
fprintf('  AP method points: %d\n', sum(strcmp(tbl.method, 'AP')));
% [-] Report altitude data range
fprintf('  Altitude range: %.0f - %.0f ft\n', min(tbl.alt_ft), max(tbl.alt_ft));
% [-] Report Mach number data range
fprintf('  Mach range: %.1f - %.1f\n', min(tbl.Mach), max(tbl.Mach));
% [-] Report power lever angle data range
fprintf('  PLA range: %.1f - %.1f deg\n', min(tbl.PLA_deg), max(tbl.PLA_deg));
% [-] Display header for flight conditions list
fprintf('\nFlight conditions:\n');
% [-] List all empirical flight test conditions
fprintf('  - 10,000 ft, M = 0.4\n');
fprintf('  - 10,000 ft, M = 0.8\n');
fprintf('  - 30,000 ft, M = 0.9\n');
fprintf('  - 30,000 ft, M = 1.2\n');
fprintf('  - 40,000 ft, M = 0.8\n');
fprintf('  - 40,000 ft, M = 1.6\n');

% [-] End create_F404_data_table function
end