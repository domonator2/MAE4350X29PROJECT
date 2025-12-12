function [Cm_alpha, Cm_delta_c, Cm_delta_s, Cm_delta_f] = get_Cm_derivatives(M)
%GET_CM_DERIVATIVES Get X-29 pitching moment derivatives at given Mach
%
% Interpolates from flight test / analysis data.
% Handles transonic gap (M=0.8-1.1) via linear interpolation.
%
% Inputs:
%   M - Mach number
%
% Outputs:
%   Cm_alpha   - [1/rad] Pitching moment due to alpha (+ = unstable)
%   Cm_delta_c - [1/rad] Pitching moment due to canard deflection
%   Cm_delta_s - [1/rad] Pitching moment due to strake/flaperon
%   Cm_delta_f - [1/rad] Pitching moment due to flap
%
% Note: X-29 is statically unstable (Cm_alpha > 0) at subsonic speeds
%
% Author: Dominic Larin
% Date: December 2025

%% Data from stability analysis
% Subsonic (M = 0.1 - 0.7) and Supersonic (M = 1.2 - 1.6)
% Transonic gap filled by interpolation

M_data = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 1.2, 1.3, 1.4, 1.5, 1.6];

Cm_alpha_data = [1.28540519458483, 1.30048571063948, 1.31824220424665, ...
                 1.35973241715865, 1.4221085688453, 1.50630302031752, ...
                 1.62511347308201, 0.580808190518962, 0.397270749432178, ...
                 0.188673476607176, 0.0803785316997435, -0.0039486628960769];

Cm_delta_c_data = [0.787445647218162, 0.792326785354983, 0.800107422253404, ...
                   0.812789791504685, 0.830759742058472, 0.855128247997783, ...
                   0.888727877773245, 0.82881078236174, 0.762266502405749, ...
                   0.705410559187163, 0.662784643656045, 0.626510618662741];

Cm_delta_s_data = [-0.340632364003107, -0.340895961679933, -0.341337088253593, ...
                   -0.341958475606603, -0.342764022305544, -0.34375887114612, ...
                   -0.344949514304525, 0.259047687437136, -0.158303888702937, ...
                   -0.320998148810591, -0.394042728303987, -0.427722992054602];

Cm_delta_f_data = [-0.479935015882191, -0.484566141237985, -0.494951187304835, ...
                   -0.506722426939007, -0.522800599318584, -0.546784231566104, ...
                   -0.582656406765123, -0.753614313199735, -0.682443305433749, ...
                   -0.645592958257033, -0.606270622431548, -0.574173807862359];

%% Handle transonic gap (M = 0.8 - 1.1)
% Add interpolated points for smoother transition
M_transonic = [0.8, 0.9, 1.0, 1.1];

% Linear interpolation between M=0.7 and M=1.2 for transonic
for i = 1:length(M_transonic)
    Mt = M_transonic(i);
    frac = (Mt - 0.7) / (1.2 - 0.7);  % 0 at M=0.7, 1 at M=1.2

    Cm_alpha_transonic(i) = Cm_alpha_data(7) + frac * (Cm_alpha_data(8) - Cm_alpha_data(7));
    Cm_delta_c_transonic(i) = Cm_delta_c_data(7) + frac * (Cm_delta_c_data(8) - Cm_delta_c_data(7));
    Cm_delta_s_transonic(i) = Cm_delta_s_data(7) + frac * (Cm_delta_s_data(8) - Cm_delta_s_data(7));
    Cm_delta_f_transonic(i) = Cm_delta_f_data(7) + frac * (Cm_delta_f_data(8) - Cm_delta_f_data(7));
end

% Combine data (insert transonic between subsonic and supersonic)
M_full = [M_data(1:7), M_transonic, M_data(8:end)];
Cm_alpha_full = [Cm_alpha_data(1:7), Cm_alpha_transonic, Cm_alpha_data(8:end)];
Cm_delta_c_full = [Cm_delta_c_data(1:7), Cm_delta_c_transonic, Cm_delta_c_data(8:end)];
Cm_delta_s_full = [Cm_delta_s_data(1:7), Cm_delta_s_transonic, Cm_delta_s_data(8:end)];
Cm_delta_f_full = [Cm_delta_f_data(1:7), Cm_delta_f_transonic, Cm_delta_f_data(8:end)];

%% Interpolate to requested Mach
% Clamp to data range
M = max(0.1, min(1.6, M));

Cm_alpha = interp1(M_full, Cm_alpha_full, M, 'linear');
Cm_delta_c = interp1(M_full, Cm_delta_c_full, M, 'linear');
Cm_delta_s = interp1(M_full, Cm_delta_s_full, M, 'linear');
Cm_delta_f = interp1(M_full, Cm_delta_f_full, M, 'linear');

end
