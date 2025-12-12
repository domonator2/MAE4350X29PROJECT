function [T, SFC, fuel_flow] = get_F404_thrust_and_SFC(h, M, use_AB)
%GET_F404_THRUST_AND_SFC Get F404 engine thrust and SFC at given conditions
%
% Wrapper for NASA TM-4140 based propulsion model.
% Uses empirical data calibrated to F404-GE-400 specifications.
%
% Inputs:
%   h      - Altitude [m]
%   M      - Mach number
%   use_AB - Use afterburner (true/false)
%
% Outputs:
%   T         - Available thrust [N]
%   SFC       - Specific fuel consumption [kg/N/s]
%   fuel_flow - Fuel flow rate [kg/s]
%
% Author: Dominic Larin
% Date: December 2025

% Convert altitude to feet for NASA model
h_ft = h / 0.3048;

% Set PLA based on afterburner usage
if use_AB
    PLA = 130;  % Max afterburner
else
    PLA = 87;   % Military power (max dry)
end

% Call NASA TM-4140 based model
try
    prop = get_F404_thrust_NASA(M, h_ft, PLA);
    T = prop.thrust_N;
    SFC = prop.TSFC_SI;           % kg/N/s
    fuel_flow = prop.fuel_flow_kg_s;
catch
    % Fallback to simplified model if NASA model not available
    [T, SFC, fuel_flow] = simplified_F404(h, M, use_AB);
end

end

%% Simplified fallback model (only if NASA model unavailable)
function [T, SFC, fuel_flow] = simplified_F404(h, M, use_AB)
% Based on F404-GE-400 specifications
% Max AB thrust (SL): 16,000 lbf (71.2 kN)
% Max MIL thrust (SL): 10,600 lbf (47.1 kN)

T_sl_mil = 47100;   % Sea level MIL thrust [N]
T_sl_ab = 71200;    % Sea level AB thrust [N]

% TSFC from F404 specs (converted to SI: kg/s/N)
% MIL: 0.78 lb/hr/lb = 0.78 * 2.832545e-5 = 2.21e-5 kg/s/N
% AB:  1.84 lb/hr/lb = 1.84 * 2.832545e-5 = 5.21e-5 kg/s/N
SFC_mil = 0.78 * 2.832545e-5;
SFC_ab = 1.84 * 2.832545e-5;

% Altitude/Mach correction
[~, ~, ~, rho] = atmosisa(h);
rho_sl = 1.225;
sigma = rho / rho_sl;

% Thrust lapse with altitude
if M < 0.8
    T_ratio = sigma^0.7;
else
    % Ram effect helps at higher Mach
    ram = 1 + 0.2 * M^2;
    T_ratio = sigma^0.6 * (ram / (1 + 0.2 * 0.8^2))^0.3;
end

% Select thrust and SFC
if use_AB
    T = T_sl_ab * T_ratio;
    SFC = SFC_ab;
else
    T = T_sl_mil * T_ratio;
    SFC = SFC_mil;
end

% Fuel flow
fuel_flow = SFC * T;

end
