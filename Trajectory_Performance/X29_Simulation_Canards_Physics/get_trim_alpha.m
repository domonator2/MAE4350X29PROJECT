function [alpha_deg, CL_achieved, CD] = get_trim_alpha(W, V, h, S_ref, gamma_deg, n_load)
%GET_TRIM_ALPHA Solve for angle of attack to achieve trim (L = W*n*cos(gamma))
%
% Calculates the required angle of attack for trimmed flight given the
% flight conditions and load factor.
%
% Inputs:
%   W         - Weight [kg]
%   V         - True airspeed [m/s]
%   h         - Altitude [m]
%   S_ref     - Wing reference area [m²]
%   gamma_deg - Flight path angle [deg] (positive = climbing)
%   n_load    - Load factor (1.0 for level flight, > 1 for maneuvering)
%
% Outputs:
%   alpha_deg   - Required angle of attack [deg]
%   CL_achieved - Lift coefficient
%   CD          - Drag coefficient
%
% Author: Dominic Larin
% Date: December 2025

if nargin < 5, gamma_deg = 0; end
if nargin < 6, n_load = 1.0; end

g = 9.81;

% Atmosphere
[~, a, ~, rho] = atmosisa(h);
M = V / a;
q = 0.5 * rho * V^2;

% Required lift for trim: L = W * n * cos(gamma)
% For level flight (gamma=0, n=1): L = W
L_required = W * g * n_load * cosd(gamma_deg);

% Required CL
CL_required = L_required / (q * S_ref);

% Limit CL to achievable values
CL_max = 1.4;  % Maximum achievable CL
CL_required = min(CL_required, CL_max);
CL_required = max(CL_required, 0);

% Aerodynamic model parameters (Mach-dependent)
if M < 0.8
    % Subsonic
    CL_alpha = 0.08;   % 1/deg (typical for swept wing)
    CL_0 = 0.05;       % Zero-alpha lift
    CD0 = 0.020;
    K = 0.10;
elseif M < 1.0
    % Transonic rise
    frac = (M - 0.8) / 0.2;
    CL_alpha = 0.08 - 0.02 * frac;  % Decreases slightly
    CL_0 = 0.05 * (1 - frac);
    CD0 = 0.020 + 0.025 * frac;
    K = 0.10 + 0.04 * frac;
elseif M < 1.2
    % Transonic peak
    frac = (M - 1.0) / 0.2;
    CL_alpha = 0.06 - 0.01 * frac;
    CL_0 = 0.0;
    CD0 = 0.045 + 0.015 * frac;
    K = 0.14 + 0.03 * frac;
else
    % Supersonic
    beta = sqrt(M^2 - 1);
    % Supersonic lift curve slope: 4 / (beta * (1 + eta)) where eta is wing correction
    CL_alpha = 4 / beta / 2.5 * (180/pi);  % Convert to per-degree, ~0.05/deg at M=1.4
    CL_alpha = min(CL_alpha, 0.06);  % Cap
    CL_0 = 0.0;
    CD0 = 0.055 + 0.005 * min((M - 1.2) / 0.4, 1);  % Slowly decreases
    K = 0.17;
end

% Solve for alpha: CL = CL_0 + CL_alpha * alpha
% alpha = (CL - CL_0) / CL_alpha
if CL_alpha > 0.001
    alpha_deg = (CL_required - CL_0) / CL_alpha;
else
    alpha_deg = 5;  % Default
end

% Physical limits
alpha_deg = max(-5, min(alpha_deg, 15));

% Recalculate achieved CL and CD
CL_achieved = CL_0 + CL_alpha * alpha_deg;
CD = CD0 + K * CL_achieved^2;

end
