%{
====================================================================================================
SCRIPT NAME: Raymer_aero_model.m
AUTHOR: Andrea Martinez
INITIATED: 9/14/2025
LAST REVISION: 9/22/2025
====================================================================================================
SCRIPT DESCRIPTION:
This function encompasses the aircraft's aerodynamic coefficients and parameters across all regimes.
====================================================================================================
%}

% Inputs (SI):
%   M              - Free-stream Mach Number. [-]
%   alpha_deg      - Angle of Attack. [deg]
%   h              - Altitude. [m]
%   geom           - A structure containing detailed vehicle geometry.

% Outputs:
%   L              - Lift. [N]
%   D              - Drag. [N]
%   LD_max         - Lift-to-Drag Ratio. [-]
%   CDi            - Induced Drag Coefficient. [-]
%   CL             - Lift Coefficient. [-]
%   CLopt             - Optimum Lift Coefficient. [-]
%   CLalpha        - Lift-Curve Slope Coefficient. [-]
%   CD             - Drag Coefficient. [-]
%   CD0            - Parasitic Drag Coefficient. [-]
%   CDi            - Induced Drag Coefficient. [-]
%   K              - Induced Drag Factor. [-]


function aero = Raymer_aero_model(M, alpha_deg, varargin)
    % Support both calling conventions:
    % raymer_aero_model(M, alpha, h, geom) - trajectory code
    % raymer_aero_model(M, alpha, geom, atm, opts) - synthesis code

    if nargin == 4
        % Old interface: (M, alpha, h, geom)
        h = varargin{1};
        geom = varargin{2};
        [~,atm.a,~,atm.rho,~,atm.mu] = atmosisa(h);
        opts = struct();
        opts.units = 'SI';
        opts.Ewave = 2.0;
    elseif nargin == 5
        % New interface: (M, alpha, geom, atm, opts)
        geom = varargin{1};
        atm = varargin{2};
        opts = varargin{3};
    else
        error('Invalid number of arguments');
    end

    % ----- Add Required Fields for Subsonic Model -----
    geom.Sref = geom.wing.S;  % Reference area (wing area)
    geom.AR = geom.wing.AR;   % Aspect ratio
    geom.Swet = geom.S_wet_total;  % Total wetted area
    geom.Cfe = 0.003;  % Equivalent skin friction coefficient (typical value)
    geom.sweep_tmax_deg = geom.wing.Lambda_c25_deg;  % Use quarter-chord sweep for thickness location
    geom.sweep_LE_deg = geom.wing.Lambda_LE_deg;  % Leading edge sweep
    
    % ----- Add Additional Fields for Transonic/Supersonic Models -----
    geom.Amax = pi * (geom.fuse.D/2)^2;  % Max cross-sectional area from fuselage diameter
    geom.Lref = geom.fuse.L;  % Use fuselage length as effective length
    geom.Ewd = 2.0;  % Wave drag efficiency factor (typical value)
    
    if ~isfield(geom, 'iw')
        geom.iw = 0;  % Wing incidence angle (degrees) - default: no incidence
    end
    if ~isfield(geom, 'alphazero')
        geom.alphazero = 0;  % Zero-lift angle (degrees) - default: zero-lift at 0 degrees
    end
    
    % ----- Drag Divergence Mach Estimate -----
    Lambda_c25 = deg2rad(geom.wing.Lambda_c25_deg);   % Wing sweep at quarter-chord
    t_c4 = geom.const.wing.tc_avg; % Average thickness-to-chord ratio
    CLopt_exp = 0.6;  % In-Flight Paper optimum CL at L/D max
    Ka = 0.95; % Airfoil technology factor (Supercritical Airfoil)

    Mdd = Ka/cos(Lambda_c25) - t_c4/(cos(Lambda_c25))^2 - CLopt_exp/(10*(cos(Lambda_c25))^3);  % Drag-Divergence Mach Number
    Mcr = Mdd - (0.1/80)^(1/3);  % Critical Mach
    Mss = 1.346;                   % Supersonic reference for Raymer

    % ----- Select Aerodynamic Regime -----
    if M <= Mcr
        % Subsonic
        out = Raymer_subsonic_model(M, alpha_deg, geom, opts);
        CL = out.CL; CLalpha = out.CL_alpha; CD = out.CD; CD0 = out.CD0_used; 
        % Calculate K from Oswald efficiency
        K = 1/(pi * geom.AR * out.e_oswald);
        CDi = K * CL^2;



     elseif M > Mcr && M < Mss
        % Transonic
        geom.MDD = Mdd;
        [CL, CLalpha, CD, CD0, CDi, K] = Raymer_transonic_model(M, alpha_deg, Mcr, Mss, geom, atm, opts);



    elseif M >= Mss
        % Supersonic
        out = Raymer_supersonic_model(M, alpha_deg, geom, atm, opts);
        CL = out.CL; CLalpha = out.CL_alpha; CD = out.CD; CD0 = out.CD0; 
        CDi = out.CDi; K = out.K;
    end

    % ----- Dynamic Pressure -----
    V = M * atm.a;
    q = 0.5 * atm.rho * V^2;
    
    Sref = geom.wing.S;

    % ---- Lift & Drag -----
    L    = q * CL * Sref;
    D    = q * CD * Sref;

    % ----- Performance Estimates -----
    LD_max = 1 / (2*sqrt(K*CD0));
    CLopt  = sqrt(CD0 / K);

    aero = struct();
    aero.L = L;
    aero.D = D;
    aero.L_over_D = L / D;  % Actual L/D at current conditions
    aero.LD_max = LD_max;
    aero.CL = CL;
    aero.CLopt = CLopt;
    aero.CL_alpha = CLalpha;  % Add both naming conventions
    aero.CLalpha = CLalpha;
    aero.CD = CD;
    aero.CD0 = CD0;
    aero.CDi = CDi;
    aero.K = K;
end
