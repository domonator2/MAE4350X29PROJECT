%{
====================================================================================================
SCRIPT NAME: get_geometry.m
AUTHOR: Jayden Martinez
INITIATED: 9/9/2025
LAST REVISION: 9/9/2025
====================================================================================================
SCRIPT DESCRIPTION:
This function defines the aircraft's geometry for use in multidisciplinary analysis, based on key
top-level geometric parameters. It is intended for the Parametric Sizing (PS) phase.
====================================================================================================
USER-DEFINED FUNCTIONS:
(diam_from_VL_spherocyl.m)|Solve for diameter D of a spherocylinder with volume V and cylindrical length L.
----------------------------------------------------------------------------------------------------
(trap_geom.m)|Calculates basic trapezoidal wing geometry parameters.
====================================================================================================
KEY ASSUMPTIONS & LOGIC:
- The vehicle is an X-29-like configuration (wing-canard-fuselage-vtail).
- The input 'Spln' is the total vehicle planform area. The wing planform
  area is assumed to be a fixed fraction of this total area.
- Canard and vertical tail areas are scaled relative to the wing area.
- Fuselage volume is a fraction of the total volume derived from 'tau'.
- Wetted areas are calculated for each component and summed. The 'Kw'
  (wetted area ratio) is an output of this calculation, not an input.

%}

% Inputs (SI):
%   Spln           - Total vehicle planform area [m^2].
%   tau            - Küchemann slenderness parameter, V_tot / Spln^1.5 [-].
%   LambdaLE_deg   - Wing leading-edge sweep angle [deg].
% Outputs:
%   geom           - A structure containing detailed vehicle geometry.
function geom = get_geometry(Spln, tau, LambdaLE_deg)

%% 1. Input Validation
% [] Checks for valid number of inputs and positive numeric values.
if nargin < 3 || ~isnumeric(Spln) || ~isnumeric(tau) || ~isnumeric(LambdaLE_deg) || Spln <= 0 || tau <= 0
    
    % [] Throws an error if inputs are invalid.
    error('get_geometry:InvalidInput', 'Inputs Spln, tau, and LambdaLE_deg must be positive numeric values.');
    
end

%% 2. Fixed Configuration Constants & Ratios
% [] Defines the units for the constants.
const.units = 'SI';

% Wing (X-29-like)
% [-] Defines the aspect ratio of the wing.
const.wing.AR      = 4.0;

% [-] Defines the taper ratio of the wing.
const.wing.lambda  = 0.40;

% [deg] Defines the dihedral angle of the wing.
const.wing.dihedral_deg = 0.0;

% [-] Defines the average thickness-to-chord ratio of the wing.
const.wing.tc_avg  = 0.05; 

% Canard
% [-] Defines the aspect ratio of the canard.
const.canard.AR     = 1.47;

% [-] Defines the taper ratio of the canard.
const.canard.lambda = 0.32;

% Vertical tail (single)  -- for S&C sizing handoffs
% [-] Defines the aspect ratio of the vertical tail.
const.vt.AR     = 2.64;

% [-] Defines the taper ratio of the vertical tail.
const.vt.lambda = 0.32;

% Fuselage proxy
% [m] Defines the vehicle length.
const.fuse.L   = 14.6558;

% [-] Defines the fraction of total volume assigned to fuselage.
const.Vfuse_frac = 0.55;

% Aerodynamic parameters
% [deg] Wing incidence angle (relative to fuselage reference line)
const.wing.iw = -2.25;

% [deg] Zero-lift angle of attack for the wing
const.wing.alphazero = -2;

% [-] Wave drag efficiency factor (1.0 for ideal Sears-Haack)
const.Ewd = 2;

% [-] Exposed wing area ratio (portion not covered by fuselage)
const.Sexposed_over_Sref = 0.85;

% Planform Area Ratios (based on original X-29 data: S_wing=185, S_canard=37, S_vtail=33.75)
% [-] Defines the assumed wing area as a fraction of total vehicle planform area.
const.S_wing_ratio  = 0.85;

% [-] Defines the canard area as a fraction of wing area.
const.S_canard_ratio = 37/185;

% [-] Defines the V-tail area as a fraction of wing area.
const.S_vtail_ratio = 33.75/185;

%% 3. Component Sizing based on Spln
% [m^2] Calculates the wing planform area.
S_w = Spln * const.S_wing_ratio;

% [m^2] Calculates the canard planform area.
S_c = S_w * const.S_canard_ratio;

% [m^2] Calculates the vertical tail planform area.
S_v = S_w * const.S_vtail_ratio;

% Reference area for aerodynamic coefficients (typically wing area)
Sref = S_w;

%% 4. Wing Geometry (Trapezoidal)
% [-] Extracts the wing aspect ratio from the constants structure.
AR = const.wing.AR;

% [-] Extracts the wing taper ratio from the constants structure.
lam = const.wing.lambda;

% [m] Calculates the wing span.
b_w   = sqrt(AR * S_w);

% [m] Calculates the wing root chord.
c_r_w = 2*S_w/(b_w*(1+lam));

% [m] Calculates the wing tip chord.
c_t_w = lam * c_r_w;

% [m] Calculates the wing mean geometric chord.
c_bar_w= (c_r_w + c_t_w) / 2;

% [m] Calculates the wing mean aerodynamic chord.
c_mac_w = (2/3)*c_r_w*((1 + lam + lam^2)/(1 + lam));

% [m] Calculates the spanwise location of the mean aerodynamic chord.
y_mac_w = (b_w/6)*((1 + 2*lam)/(1 + lam));

% Sweep conversions from LE for reporting
% [-] Calculates the sweep shift factor.
Kshift = (4/AR) * (1-lam)/(1+lam);

% [-] Calculates the tangent of the leading-edge sweep angle.
tanLE  = tand(LambdaLE_deg);

% [deg] Calculates the quarter-chord sweep angle.
Lam_c25_deg = atand(tanLE - 0.25*Kshift);

% [deg] Calculates the trailing-edge sweep angle.
Lam_TE_deg  = atand(tanLE - 1.00*Kshift);

% [deg] Calculates half-chord sweep angle (for transonic analysis)
Lam_c50_deg = atand(tanLE - 0.50*Kshift);

%% 5. System Volume & Fuselage Proxy
% [m^3] Calculates the total vehicle volume (Küchemann definition).
V_tot   = tau * Spln^(1.5);

% [m^3] Calculates the allocated fuselage volume.
V_fuse  = const.Vfuse_frac * V_tot;

% [m] Extracts the fuselage length from the constants structure.
L_fuse  = const.fuse.L;

% [m] Calculates the equivalent fuselage diameter.
% NOTE: Spherocylinder assumption produces L/D ≈ 14 (supersonic-optimized)
% but X-29 is a mixed subsonic/supersonic fighter requiring L/D ≈ 8-10
% per Nicolai & Carichner Fig 8.11. Current implementation uses simple
% spherocylinder; structural weight calibration factors (K_fus) compensate
% for the difference between spherocylinder and actual blended wing-body.
% See Geometry/FUSELAGE_FINENESS_RATIO.md for detailed analysis.
D_fuse  = diam_from_VL_spherocyl(V_fuse, L_fuse);

% [m^2] Maximum cross-sectional area for wave drag calculation
% Per Raymer "Aircraft Design: A Conceptual Approach" (6th ed., Section 12.5.2):
% "The maximum cross-sectional area should include wing thickness, nacelles, and
%  other components, not just the fuselage." For supersonic aircraft, this effective
%  area is typically 2.5-4.0 times the fuselage-only cross-section.
%
% Factor breakdown for fighter aircraft (Nicolai & Carichner, "Fundamentals of Aircraft
% and Airship Design", Ch. 9):
% - Wing root contribution: ~1.8x fuselage area
% - Engine inlet/nacelle: ~0.8x fuselage area
% - Non-circular cross-section factor: ~1.2x
% - Total effective multiplier: ~3.5x (conservative for transonic/supersonic analysis)
A_max_fuse = pi * (D_fuse/2)^2;
A_max = A_max_fuse * 3.5;  % Effective cross-sectional area multiplier

% [m] Reference length for wave drag (typically fuselage length)
L_ref = L_fuse;

%% 6. Canard & Vertical Tail Geometry (for handoffs)
% [m, m, m, m, m, m] Calculates the canard geometry parameters.
[b_c, cr_c, ct_c, c_bar_c, c_mac_c, y_mac_c] = trap_geom(S_c, const.canard.AR, const.canard.lambda);

% [m, m, m, m, m, m] Calculates the vertical tail geometry parameters.
[b_v, cr_v, ct_v, c_bar_v, c_mac_v, y_mac_v] = trap_geom(S_v, const.vt.AR, const.vt.lambda);

%% 7. Wetted Area Calculations (Bottom-Up)
% [m^2] Calculates the wetted area of the fuselage (spherocylinder).
S_wet_fuse = pi*D_fuse*(L_fuse - D_fuse) + pi*D_fuse^2;

% [m^2] Ensures a non-negative wetted area for the fuselage.
S_wet_fuse = max(S_wet_fuse, 0);

% [m^2] Calculates the wetted area of the wing.
S_wet_wing  = 2 * S_w * (1 + 0.5 * const.wing.tc_avg);

% [m^2] Calculates the wetted area of the canard (assuming thin surface).
S_wet_canard = 2 * S_c;

% [m^2] Calculates the wetted area of the vertical tail (assuming thin surface).
S_wet_vtail = 2 * S_v;

% [m^2] Calculates the total wetted area.
S_wet_total = S_wet_fuse + S_wet_wing + S_wet_canard + S_wet_vtail;

% [-] Calculates the wetted area ratio.
Kw = S_wet_total / Spln;

%% 7a. Create Components Array for Supersonic Drag Buildup
% Components structure for drag coefficient calculations
components = [];

% Wing component
components(1).name = 'wing';
components(1).S_wet = S_wet_wing;  % For subsonic drag calculation
components(1).Swet = S_wet_wing;   % For supersonic drag calculation
components(1).Lref = c_mac_w;  % Use MAC as reference length
components(1).FF = 1.0;  % Form factor (will be calculated in drag routine)
components(1).Q = 1.0;   % Interference factor
components(1).type = 'wing';
components(1).Cf_method.Re = 1e7;  % Placeholder Reynolds number
components(1).Cf_method.smooth = true;

% Fuselage component
components(2).name = 'fuselage';
components(2).S_wet = S_wet_fuse;  % For subsonic drag calculation
components(2).Swet = S_wet_fuse;   % For supersonic drag calculation
components(2).Lref = L_fuse;
components(2).FF = 1.0;
components(2).Q = 1.0;
components(2).type = 'fuselage';
components(2).Cf_method.Re = 1e8;  % Placeholder Reynolds number
components(2).Cf_method.smooth = true;

% Canard component
components(3).name = 'canard';
components(3).S_wet = S_wet_canard;  % For subsonic drag calculation
components(3).Swet = S_wet_canard;   % For supersonic drag calculation
components(3).Lref = c_mac_c;
components(3).FF = 1.0;
components(3).Q = 1.0;
components(3).type = 'wing';
components(3).Cf_method.Re = 1e6;  % Placeholder Reynolds number
components(3).Cf_method.smooth = true;

% Vertical tail component
components(4).name = 'vtail';
components(4).S_wet = S_wet_vtail;  % For subsonic drag calculation
components(4).Swet = S_wet_vtail;   % For supersonic drag calculation
components(4).Lref = c_mac_v;
components(4).FF = 1.0;
components(4).Q = 1.0;
components(4).type = 'tail';
components(4).Cf_method.Re = 1e6;  % Placeholder Reynolds number
components(4).Cf_method.smooth = true;

%% 8. Pack Outputs into 'geom' Structure
% [] Initializes the output structure.
geom = struct();

% inputs & global
% [] Stores the input values in the output structure.
geom.inputs = struct('Spln',Spln,'tau',tau,'LambdaLE_deg',LambdaLE_deg);

% [] Stores the constants in the output structure.
geom.const  = const;

% Top-level calculated values
% [m^3] Stores the total volume in the output structure.
geom.V_tot        = V_tot;

% [-] Stores the wetted area ratio in the output structure.
geom.Kw           = Kw;

% [m^2] Stores the total wetted area in the output structure.
geom.S_wet_total  = S_wet_total;

% Reference values for aerodynamics
geom.Sref = Sref;  % Reference area (wing area)
geom.Lref = L_ref;  % Reference length (fuselage length)
geom.Amax = A_max;  % Maximum cross-sectional area
geom.sweepLE = LambdaLE_deg;  % Leading edge sweep
geom.AR = AR;  % Wing aspect ratio

% Aerodynamic parameters
geom.iw = const.wing.iw;  % Wing incidence angle
geom.alphazero = const.wing.alphazero;  % Zero-lift angle
geom.Ewd = const.Ewd;  % Wave drag efficiency factor
geom.Sexposed_over_Sref = const.Sexposed_over_Sref;  % Exposed wing ratio

% Components array for drag buildup
geom.components = components;

% Total wetted area for functions that expect geom.Swet
geom.Swet = S_wet_total;

% wing block
% [] Stores the wing geometry in the output structure.
geom.wing = struct( ...
   'S', S_w, 'AR', AR, 'lambda', lam, 'dihedral_deg', const.wing.dihedral_deg, ...
   'span', b_w, 'b', b_w, 'c_r', c_r_w, 'c_t', c_t_w, 'c_bar', c_bar_w, ...
   'c_mac',c_mac_w,'y_mac',y_mac_w, ...
   'Lambda_LE_deg',LambdaLE_deg,'Lambda_c25_deg',Lam_c25_deg,'Lambda_c50_deg',Lam_c50_deg,'Lambda_TE_deg',Lam_TE_deg, ...
   'S_wet', S_wet_wing);

% fuselage block
% [] Stores the fuselage geometry in the output structure.
geom.fuse = struct('L',L_fuse,'D',D_fuse,'V',V_fuse,'S_wet',S_wet_fuse);

% canard block
% [] Stores the canard geometry in the output structure.
geom.canard = struct('S',S_c,'AR',const.canard.AR,'lambda',const.canard.lambda, ...
                     'span',b_c,'b',b_c,'c_r',cr_c,'c_t',ct_c,'c_bar',c_bar_c,'c_mac',c_mac_c,'y_mac',y_mac_c, ...
                     'S_wet', S_wet_canard);

% vertical tail block
% [] Stores the vertical tail geometry in the output structure.
geom.vtail = struct('S',S_v,'AR',const.vt.AR,'lambda',const.vt.lambda, ...
                    'span',b_v,'b',b_v,'c_r',cr_v,'c_t',ct_v,'c_bar',c_bar_v,'c_mac',c_mac_v,'y_mac',y_mac_v, ...
                    'S_wet', S_wet_vtail);

% convenience top-level fields (for quick access)
% [m^2] Stores the total planform area in the output structure.
geom.Spln = Spln;

% [-] Stores the Küchemann slenderness parameter in the output structure.
geom.tau  = tau;

% [m] Alias for wing span.
geom.span = b_w;

% [m] Alias for wing root chord.
geom.c_r  = c_r_w;

% [m] Alias for wing tip chord.
geom.c_t  = c_t_w;

% [m] Alias for wing mean aerodynamic chord.
geom.c_mac= c_mac_w;

% Additional reference values
geom.length = L_fuse;  % Vehicle length
geom.fuselage_d_over_b = D_fuse / b_w;  % Fuselage diameter to wing span ratio

end

%DIAM_FROM_VL_SPHEROCYL  Solve for diameter D of a spherocylinder
% with volume V and TOTAL length L.
%
% Model: V = (pi/4)*D^2*(L-D) + (pi/6)*D^3 = (pi/4)*L*D^2 - (pi/12)*D^3
% This is rearranged into a cubic polynomial to solve for D:
% (pi/12)*D^3 - (pi*L/4)*D^2 + V = 0
% Returns D > 0 (same units as L). If inputs are non-positive, returns 0.
function D = diam_from_VL_spherocyl(V, L)

    % [] Checks for non-positive inputs.
    if V <= 0 || L <= 0
        
        % [m] Sets diameter to 0 for invalid inputs.
        D = 0;
        
        % [] Exits the function.
        return
        
    end

    % Root function is a cubic polynomial in D.
    % (pi/12)*D^3 - (pi*L/4)*D^2 + 0*D + V = 0
    % [-] Defines the coefficients of the cubic polynomial.
    p = [pi/12, -pi*L/4, 0, V];
    
    % [-] Finds the roots of the polynomial.
    r = roots(p);

    % Try root find; if anything goes wrong, fall back to D0 (pure cylinder)
    try
        % Find the real, positive roots. There can be two, we need the smaller one
        % that is physically possible (D < L).
        D_real_pos = sort(r(imag(r) == 0 & real(r) > 0 & real(r) < L));
        
        % [] Checks if a valid root was found.
        if ~isempty(D_real_pos)
            
            % [m] Assigns the valid root to D.
            D = D_real_pos(1);
            
        else
            
            % Fallback to pure cylinder if no valid root found
            % [m] Calculates D assuming a pure cylinder.
            D = sqrt(4*V/(pi*L));
            
        end
        
    catch
        
        % Fallback for any other error
        % [m] Calculates D assuming a pure cylinder.
        D = sqrt(4*V/(pi*L));
        
    end
    
end

% Calculates basic trapezoidal wing geometry parameters.
function [b, cr, ct, cbar, cmac, ymac] = trap_geom(S, AR, lambda)

% [m] Calculates the span of the trapezoidal surface.
b     = sqrt(AR*S);

% [m] Calculates the root chord of the trapezoidal surface.
cr    = 2*S/(b*(1+lambda));

% [m] Calculates the tip chord of the trapezoidal surface.
ct    = lambda*cr;

% [m] Calculates the mean geometric chord of the trapezoidal surface.
cbar  = 2*S/b;

% [m] Calculates the mean aerodynamic chord of the trapezoidal surface.
cmac  = (2/3)*cr*((1+lambda+lambda^2)/(1+lambda));

% [m] Calculates the spanwise location of the mean aerodynamic chord.
ymac  = (b/6)*((1+2*lambda)/(1+lambda));

end
