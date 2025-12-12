
function out = raymer_aero_model(Mach, alpha_deg, geom, atm, opts)
% RAYMER_AERO_MODEL Unified aerodynamic model based on Raymer's methods.
%
% This function selects the appropriate aerodynamic model (subsonic, transonic,
% or supersonic) based on the input Mach number and computes the relevant
% aerodynamic coefficients.
%
% REFERENCES:
% - Raymer, D. "Aircraft Design: A Conceptual Approach" 6th Edition, AIAA, 2018
% - Nicolai, L. & Carichner, G. "Fundamentals of Aircraft and Airship Design", AIAA, 2010
% - Roskam, J. "Airplane Design Part VI: Preliminary Calculation of Aerodynamic..."
% - USAF Stability and Control DATCOM, McDonnell Douglas, 1978
%
% MODIFICATIONS for X-29 Project:
% - Enhanced supersonic drag model with realistic drag increments
% - Corrected maximum cross-sectional area calculation for wave drag
% - Added component-based drag buildup for improved accuracy
%
% Inputs:
%   Mach      - Free-stream Mach number(s).
%   alpha_deg - Angle(s) of attack in degrees.
%   geom      - Struct describing the airframe geometry.
%   atm       - Atmosphere structure (e.g., with rho, mu, a).
%   opts      - (Optional) struct of name-value pairs for configuration.
%
% Outputs:
%   out       - Structure containing the aerodynamic coefficients (CL, CD, etc.).

    arguments
        Mach (1,1) double {mustBeNonnegative}
        alpha_deg (:,1) double
        geom struct
        atm struct
        opts struct = struct()
    end

    if Mach < 0.8
        % --- Subsonic Regime ---
        out = raymer_subsonic_model_internal(Mach, alpha_deg, geom, opts);
    elseif Mach >= 0.8 && Mach <= 1.2
        % --- Transonic Regime ---
        out = raymer_transonic_model_internal(Mach, alpha_deg, geom, atm, opts);
    elseif Mach > 1.2 && Mach < 1.3
        % --- Transonic-Supersonic Blending Zone (M=1.2 to 1.3) ---
        out_trans = raymer_transonic_model_internal(1.2, alpha_deg, geom, atm, opts);
        out_super = raymer_supersonic_model_internal(Mach, alpha_deg, geom, atm, opts);

        % Blend factor: 0 at M=1.2 (pure transonic), 1 at M=1.3 (pure supersonic)
        blend = (Mach - 1.2) / (1.3 - 1.2);

        % Blend all coefficients
        out.CL = (1-blend)*out_trans.CL + blend*out_super.CL;
        out.CD = (1-blend)*out_trans.CD + blend*out_super.CD;
        out.CD0 = (1-blend)*out_trans.CD0 + blend*out_super.CD0;
        out.CDi = (1-blend)*out_trans.CDi + blend*out_super.CDi;
        out.CL_alpha = (1-blend)*out_trans.CL_alpha + blend*out_super.CL_alpha;
        out.K = (1-blend)*out_trans.K + blend*out_super.K;
        out.L_over_D = out.CL ./ out.CD;
    else
        % --- Supersonic Regime (M >= 1.3) ---
        out = raymer_supersonic_model_internal(Mach, alpha_deg, geom, atm, opts);
    end
end

%% ========================================================================
%                       INTERNAL MODEL IMPLEMENTATIONS
% =========================================================================

%% Subsonic Model
function out = raymer_subsonic_model_internal(Mach, alpha_deg, geom, opts)
    % Unpack geometry & options
    Sref  = need(geom.wing,'S');
    AR    = need(geom.wing,'AR');
    lam_t = deg2rad(abs(need(geom.wing,'Lambda_c25_deg')));
    lamLE = deg2rad(abs(need(geom.wing,'Lambda_LE_deg')));
    Sexposed = getOr(geom,'Sexposed_over_Sref',1.0);
    if isfield(geom,'Sexposed_over_Sref')
        Sexposed = geom.Sexposed_over_Sref;
    end
    d_over_b = 0.0;
    if isfield(geom,'fuselage_d_over_b')
        d_over_b = geom.fuselage_d_over_b;
    elseif isfield(geom,'fuse') && isfield(geom.fuse,'D') && isfield(geom.wing,'span')
        d_over_b = geom.fuse.D / geom.wing.span;
    end

    alpha_L0_deg = getOr(opts,'alpha_L0_deg',0.0);
    eta_airfoil  = getOr(opts,'eta_airfoil',0.95);

    % CL_alpha
    beta = sqrt(max(1e-9, 1 - Mach^2));
    F_fus = 1.07*(1 + d_over_b)^2;
    exposure = min(0.98, Sexposed * F_fus);
    term = (AR*beta/eta_airfoil)^2 * (1 + (tan(lam_t)^2)/(beta^2));
    CL_alpha = (2*pi*AR) ./ (2 + sqrt(4 + term));
    CL_alpha = CL_alpha * exposure;

    % Lift
    alpha_rad = deg2rad(alpha_deg - alpha_L0_deg);
    CL = CL_alpha .* alpha_rad;

    % Oswald e
    e_oswald = e_oswald_raymer(AR, lamLE);

    % Parasite drag CD0
    if isfield(opts,'CD0') && ~isempty(opts.CD0)
        CD0 = opts.CD0;
        CD0_detail = struct('path','override','CD0',CD0);
    elseif isfield(opts,'drag_buildup') && ~isempty(opts.drag_buildup)
        for k = 1:numel(opts.drag_buildup)
            if ~isfield(opts.drag_buildup(k),'type')
                opts.drag_buildup(k).type = opts.drag_buildup(k).name;
            end
        end
        [CD0, dtab] = cd0_component_buildup(opts.drag_buildup, Sref, Mach, opts);
        CD0_detail = struct('path','component_buildup','table',dtab);
    else
        % Default Cfe for simple drag estimation
        % NOTE: Trade study uses component-based drag buildup, not this path
        Cfe  = getOr(opts,'Cfe',getOr(geom,'Cfe',0.003));
        Swet = getOr(geom,'Swet',NaN);
        if isnan(Swet)
            error('Provide opts.CD0, or opts.drag_buildup, or geom.Swet.');
        end
        CD0 = Cfe * (Swet/Sref);
        CD0_detail = struct('path','Cfe_Swet','Cfe',Cfe,'Swet',Swet);
    end

    % Base drag
    CD_base = 0;
    if isfield(opts,'base_drag') && isfield(opts.base_drag,'A_base') && ~isempty(opts.base_drag.A_base)
        CD_base = cd_base_raymer(Mach, opts.base_drag.A_base, Sref);
    end

    % Excrescence drag with smooth taper
    % Calibrated to match X-29 subsonic L/D of 9.8 at M=0.6
    % Accounts for: interference drag, roughness, excrescences, leakage
    % Based on mach_sweep validation: increased from 5 to 10 counts
    CD_excrescence_max = 0.0010;  % 10 drag counts for subsonic interference/excrescences
    M_excr_start = 0.5;
    M_excr_end = 0.9;

    if Mach <= M_excr_start
        CD_excrescence = CD_excrescence_max;
    elseif Mach >= M_excr_end
        CD_excrescence = 0.0;
    else
        % Smooth cosine taper
        theta = pi * (Mach - M_excr_start) / (M_excr_end - M_excr_start);
        CD_excrescence = CD_excrescence_max * 0.5 * (1 + cos(theta));
    end

    % Induced drag
    K = 1 / (pi * AR * e_oswald);
    CDi_clean = K .* (CL.^2);

    % High-lift device increments (drag)
    dCD_flap = 0;
    dCDi_flap = 0;
    if isfield(opts,'flaps') && ~isempty(opts.flaps)
        fl = opts.flaps;
        if all(isfield(fl, {'delta_flap_deg','Cf_over_c','Sflap_over_Sref','typeTE'}))
            dCD_flap = dCD_flap_raymer(fl.delta_flap_deg, fl.Cf_over_c, fl.Sflap_over_Sref, fl.typeTE);
        end
        if isfield(fl,'deltas') && isfield(fl.deltas,'deltaCL_total') && ~isempty(fl.deltas.deltaCL_total)
            kf = 0.14; if isfield(fl,'span_is_full') && ~fl.span_is_full, kf = 0.28; end
            dCDi_flap = kf * (fl.deltas.deltaCL_total).^2 * cos(lamLE);
        end
    end

    % Total drag and L/D
    CD0_total = CD0 + CD_base + dCD_flap + CD_excrescence;
    CDi_total = CDi_clean + dCDi_flap;
    CD = CD0_total + CDi_total;
    L_over_D = CL ./ CD;

    % Package outputs
    out.CL        = CL(:);
    out.CD        = CD(:);
    out.CL_alpha  = CL_alpha;
    out.L_over_D  = L_over_D(:);
    out.K         = K;
    out.e_oswald  = e_oswald;
    out.CD0       = CD0_total;
    out.CDi       = CDi_total;
    out.components = struct('CD0_core',CD0_detail, 'CD_base',CD_base, ...
                            'dCD_flap',dCD_flap, 'dCDi_flap',dCDi_flap, ...
                            'CDi_clean',CDi_clean);
end


%% Transonic Model
function out = raymer_transonic_model_internal(Minf, alpha_deg, geom, atm, opts)
    Mcr = getOr(opts, 'Mcr', 0.75); % Critical Mach, example
    Mss = getOr(opts, 'Mss', 1.2);  % Supersonic reference Mach

    % Transonic CL
    out_sub = raymer_subsonic_model_internal(Mcr, alpha_deg, geom, opts);
    CLalpha_sub = out_sub.CL_alpha;
    out_sup = raymer_supersonic_model_internal(Mss, alpha_deg, geom, atm, opts);
    CLalpha_sup = out_sup.CL_alpha;
    
    if isnan(CLalpha_sup)
        CLalpha = CLalpha_sub;
    else
        CLalpha = pchip([Mcr Mss], [CLalpha_sub CLalpha_sup], Minf);
    end
    
    alpha_rad = deg2rad(alpha_deg - getOr(opts,'alpha_L0_deg',0.0));
    CL = CLalpha * alpha_rad;

    % Transonic CD
    % Get baseline CD0 at transonic entry (M=0.8) for smooth continuity
    M_trans_start = 0.8;
    out_sub_entry = raymer_subsonic_model_internal(M_trans_start, 0.0, geom, opts);
    CD0_at_entry = out_sub_entry.CD0;

    % Get CD0 at supersonic transition (M=1.2)
    out_sup_12 = raymer_supersonic_model_internal(1.2, 0.0, geom, atm, opts);
    CD0_at_12 = out_sup_12.CD0;
    if isnan(CD0_at_12)
        CD0_at_12 = CD0_at_entry;
    end

    % Wave drag rise using A-B-C-D-E construction
    CDw12 = raymer_CDw12(geom);
    MDD = Mcr + 0.08;  % Drag divergence Mach

    % Control points for smooth wave drag rise to peak at M=1.2
    % Remove plateau - drag should continuously rise through transonic
    M_pts = [M_trans_start, MDD, 0.95, 1.05, 1.15, 1.2];
    dCD_pts = [0, 0.002, 0.3*CDw12, 0.6*CDw12, 0.9*CDw12, CDw12];
    dCD0_wave = pchip(M_pts, dCD_pts, Minf);

    % Baseline CD0 interpolates from transonic entry to supersonic transition
    CD0_baseline = interp1([M_trans_start 1.2], [CD0_at_entry CD0_at_12], Minf, 'linear', 'extrap');
    CD0 = CD0_baseline + dCD0_wave;

    Ksub = out_sub.K;
    Ksup = out_sup.K;
    if isnan(Ksup)
        K = Ksub;
    else
        K = pchip([Mcr Mss], [Ksub Ksup], Minf);
    end
    
    CDi = K .* CL.^2;
    CD  = CD0 + CDi;
    L_over_D = CL ./ CD;

    out.CL        = CL(:);
    out.CD        = CD(:);
    out.CL_alpha  = CLalpha;
    out.L_over_D  = L_over_D(:);
    out.K         = K;
    out.CD0       = CD0;
    out.CDi       = CDi;
end

%% Supersonic Model
function out = raymer_supersonic_model_internal(Minf, alpha_deg, geom, atm, opts)
    % CL Supersonic
    cl_out = cl_supersonic_internal(Minf, alpha_deg, geom, opts);
    CL = cl_out.CL;
    CLalpha = cl_out.CLa;

    % CD Supersonic
    cd_out = cd_supersonic_internal(Minf, alpha_deg, geom, atm, opts, cl_out);
    CD = cd_out.CD_total;
    CD0 = cd_out.CD0_total;
    CDi = cd_out.CD_lift;
    K = cd_out.K;
    L_over_D = CL ./ CD;

    out.CL        = CL(:);
    out.CD        = CD(:);
    out.CL_alpha  = CLalpha;
    out.L_over_D  = L_over_D(:);
    out.K         = K;
    out.CD0       = CD0;
    out.CDi       = CDi;
end


%% ========================================================================
%                       HELPER FUNCTIONS
% =========================================================================

function val = need(s, field)
    if ~isfield(s,field) || isempty(s.(field)), error('Missing required field: %s', field); end
    val = s.(field);
end

function val = getOr(s, field, default)
    if isstruct(s) && isfield(s,field) && ~isempty(s.(field)), val = s.(field); else, val = default; end
end

function e = e_oswald_raymer(AR, lamLE)
    sweepLEdeg = rad2deg(lamLE);
    e_straight = 1.78*(1 - 0.045*AR^0.68) - 0.64;
    e_swept    = 4.61*(1 - 0.045*AR^0.68) * (cos(lamLE))^0.15 - 3.10;
    if sweepLEdeg <= 0, e = e_straight;
    elseif sweepLEdeg >= 30, e = e_swept;
    else
        t = sweepLEdeg/30; e = (1-t)*e_straight + t*e_swept;
    end
    e = max(0.1, e);
end

function CDw12 = raymer_CDw12(geom)
    Amax = geom.Amax;
    Lref = geom.Lref;
    Sref = geom.wing.S;
    Ewd  = getOr(geom, 'Ewd', 2.0);
    Dq_SH_12 = (9*pi/2) * (Amax / Lref)^2;
    CDw12 = (Ewd * Dq_SH_12) / Sref;
end

function cl_out = cl_supersonic_internal(M, alpha_deg, geom, opts)
    A = geom.wing.AR;
    LambdaLE_deg = geom.wing.Lambda_LE_deg;
    
    beta  = sqrt(M.^2 - 1.0);
    tanLE = tan(abs(LambdaLE_deg) * pi/180);
    Atan  = A .* tanLE;
    x     = tanLE ./ beta;
    Mmin  = sqrt(1 + tanLE.^2);

    xgrid = [1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.20 0.10];
    yA2   = [3.15 3.22 3.30 3.40 3.50 3.60 3.72 3.84 3.94 4.00];
    yA3   = [3.45 3.52 3.60 3.70 3.78 3.86 3.92 3.97 4.00 4.00];

    xinc = fliplr(xgrid);
    y20x = interp1(xinc, fliplr(yA2), min(max(x,xgrid(end)),xgrid(1)), 'pchip','extrap');
    y30x = interp1(xinc, fliplr(yA3), min(max(x,xgrid(end)),xgrid(1)), 'pchip','extrap');

    t = (Atan - 2.0) ./ (3.0 - 2.0);
    y_beta = (1 - t) .* y20x + t .* y30x;
    CLa_raw = y_beta ./ beta;

    Sratio    = getOr(opts, 'SexposedOverSref', 1.0);
    Ffuse_raw = 1.0; % Simplified
    scale_used = min(Sratio * Ffuse_raw, 0.98);
    CLa_supersonic_LE = CLa_raw .* scale_used;

    % Handle subsonic leading edge case (x > 1)
    % When LE is swept behind Mach cone, use subsonic-LE theory
    % Reference: USAF DATCOM Section 4.1.3.2, NACA RM-A53G08
    needsSubsonicLE = x > 1;
    invalidM  = M <= 1;

    if any(needsSubsonicLE)
        % Subsonic leading edge formula from DATCOM
        % For highly swept wings at supersonic speeds
        eta = getOr(opts, 'eta_airfoil', 0.95);  % Airfoil efficiency

        % DATCOM subsonic-LE formula for supersonic flight
        % CLa = (2*pi*AR) / sqrt(4 + (AR*beta/eta)^2)
        CLa_subsonic_LE = (2*pi*A) ./ sqrt(4 + (A*beta/eta).^2);
        CLa_subsonic_LE = CLa_subsonic_LE * scale_used;

        % Use subsonic-LE formula where needed
        CLa_supersonic_LE(needsSubsonicLE) = CLa_subsonic_LE;
    end

    % Set invalid Mach numbers to NaN
    if any(invalidM)
        CLa_supersonic_LE(invalidM) = NaN;
    end

    CLa = CLa_supersonic_LE;

    alpha_rad = deg2rad(alpha_deg);
    CL = CLa(:) .* alpha_rad(:);

    cl_out.CL = CL;
    cl_out.CLa = CLa;
end

function cd_out = cd_supersonic_internal(M, alpha_deg, geom, atm, opts, cl_in)
    Sref = geom.wing.S;
    AR = geom.wing.AR;
    sweepLE = geom.wing.Lambda_LE_deg;

    rho = atm.rho;
    mu = atm.mu;
    a = atm.a;

    CL = cl_in.CL;

    % Extract slenderness parameter for configuration-dependent drag
    tau = getOr(geom, 'tau', 0.05);  % Default to mid-range if not provided
    
    % Skin Friction Drag
    CD0_fric = 0;
    if isfield(geom, 'components')
        for k = 1:numel(geom.components)
            c = geom.components(k);
            Swet = c.Swet;
            Lref = c.Lref;
            Re = (rho * (a * M) * Lref) / mu;
            Cf = 0.455 / (log10(Re)^2.58 * (1 + 0.144*M^2)^0.65);
            CD0_fric = CD0_fric + Cf * (Swet / Sref);
        end
    else
        Cfe  = getOr(opts,'Cfe',getOr(geom,'Cfe',0.003));
        Swet = getOr(geom,'Swet', 2*Sref); % Approximation
        CD0_fric = Cfe * (Swet/Sref);
    end

    % Wave Drag
    if isfield(geom, 'Amax') && isfield(geom, 'Lref')
        Ewave = getOr(opts, 'Ewave', 2.0);
        Dq_SH = (9*pi/2) * (geom.Amax/geom.Lref)^2;
        CD_wave = Ewave * Dq_SH / Sref;
    else
        CD_wave = 0;
    end

    % Additional supersonic drag increments for fighter aircraft
    % Per Raymer "Aircraft Design: A Conceptual Approach" (6th ed., Table 12.3):
    % "Supersonic aircraft require additional drag increments beyond skin friction and wave drag"
    %
    % NOTE: These empirical drags were calibrated for spherocylinder fuselage model
    % Can be disabled via opts.disable_empirical_drags = true for elliptical fuselage
    %
    if getOr(opts, 'disable_empirical_drags', false)
        % Elliptical fuselage configuration - use only physics-based drags
        CD_trim = 0;
        CD_misc = 0;
    else
        % Standard configuration - use empirical corrections calibrated for spherocylinder
        % CD_trim: Control surface trim drag (Roskam "Airplane Design Part VI", Section 10.2.3)
        % X-29 destabilizing canard configuration requires significant trim authority
        % Scale with slenderness: slender configs (low tau) need more trim for stability
        % Calibrated to match X-29 (tau=0.07) L/D=4.2 at M=1.4
        CD_trim_base = 0.0028;  % Baseline for tau=0.07 (calibrated to X-29 L/D=4.2@M1.4)
        tau_effect_trim = -0.012 * (tau - 0.07);  % Minimal scaling, anchored at tau=0.07
        CD_trim = CD_trim_base + tau_effect_trim;
        CD_trim = max(0.0018, min(0.0045, CD_trim));  % Bound to reasonable range

        % CD_misc: Miscellaneous supersonic drag increments (USAF DATCOM, Section 4.1.5.2)
        % Includes: interference drag, antennas, probes, surface roughness, gaps, leakage, base drag
        % Scale with slenderness: slender configs have more interference drag
        % X-29 research aircraft has extensive instrumentation
        CD_misc_base = 0.0029;  % Baseline for tau=0.07 (calibrated to X-29 L/D=4.2@M1.4)
        tau_effect_misc = -0.018 * (tau - 0.07);  % Minimal scaling, anchored at tau=0.07
        CD_misc = CD_misc_base + tau_effect_misc;
        CD_misc = max(0.0022, min(0.0050, CD_misc));  % Bound to reasonable range
    end

    CD0_total = CD0_fric + CD_wave + CD_trim + CD_misc;

    % Drag due to lift
    beta = sqrt(M^2 - 1);
    LambdaLE_rad = deg2rad(abs(sweepLE));

    % Check if subsonic or supersonic leading edge
    tanLE = tan(LambdaLE_rad);
    x = tanLE / beta;  % LE parameter: x > 1 means subsonic LE

    if x > 1
        % Subsonic leading edge: Use modified induced drag formula
        % Reference: DATCOM Section 4.1.5.1(f), similar to subsonic theory
        e = e_oswald_raymer(AR, LambdaLE_rad);  % Oswald efficiency
        K = 1 / (pi * e * AR);
    else
        % Supersonic leading edge: Use standard supersonic formula
        den = 4*AR*beta - 2;
        K   = (AR*beta.^2 .* cos(LambdaLE_rad)) ./ den;
        if den <= 0, K = NaN; end
    end

    CD_lift  = K(:) .* (CL.^2);
    CD_total = CD0_total + CD_lift;

    cd_out.CD_total = CD_total;
    cd_out.CD0_total = CD0_total;
    cd_out.CD_lift = CD_lift;
    cd_out.K = K;
end

function CD_base = cd_base_raymer(M, A_base, Sref)
    D_over_q = (0.139 + 0.419*(M - 0.161)^2) * A_base;
    CD_base  = D_over_q / Sref;
end

function dCD = dCD_flap_raymer(delta_flap_deg, Cf_over_c, Sflap_over_Sref, typeTE)
    switch lower(strtrim(typeTE))
        case 'plain',   Fflap = 0.0144;
        case 'slotted', Fflap = 0.0074;
        otherwise, warning('Unknown flap type "%s"; using slotted constant.', typeTE); Fflap = 0.0074;
    end
    dCD = Fflap * (Cf_over_c) * (Sflap_over_Sref) * (delta_flap_deg - 10);
    dCD = max(0, dCD);
end

function [CD0, dtab] = cd0_component_buildup(components, Sref, Mach, opts)
    CD0 = 0;
    dtab = struct('name',{},'Cf',{},'FF',{},'Q',{},'CD0i',{});
    for i = 1:numel(components)
        c = components(i);
        Q = fieldOr(c,'Q',1.0);
        FF = 1.0;

        if isfield(c,'Cf_method') && isfield(c.Cf_method,'Cfe')
            Cf = c.Cf_method.Cfe;
        elseif isfield(c,'Cf_method') && isfield(c.Cf_method,'Re')
            Re = c.Cf_method.Re;
            smooth = fieldOr(c.Cf_method,'smooth',true);
            k_over_l = fieldOr(c.Cf_method,'k_over_l',0);
            Cf = cf_flatplate_turbulent(Re);
            if ~smooth && k_over_l>0
                Re_cut = 38.21*(1./k_over_l)^1.053;
                Cf = cf_flatplate_turbulent(min(Re,Re_cut));
            end
        else
            error('Component %d missing Cf_method.', i);
        end

        if isfield(c,'FF_method')
            if isfield(c.FF_method,'x_c_max')
                xcm  = c.FF_method.x_c_max;
                tc   = c.FF_method.t_c;
                LamM = deg2rad(c.FF_method.Lambda_m_deg);
                FF = (1 + (0.6/xcm)*tc + 100*tc^4) * (1.34*Mach^0.18 * cos(LamM)^0.28);
            elseif isfield(c.FF_method,'f') && isfield(c,'type') && strcmpi(c.type,'fuselage')
                f = c.FF_method.f;
                FF = 0.9 + 5/(f^1.5) + f/400;
            elseif isfield(c.FF_method,'f')
                f = c.FF_method.f;
                FF = 1 + 0.35/f;
            else
                FF = 1.0;
            end
        end

        CD0i = Cf * FF * Q * (c.S_wet / Sref);
        CD0  = CD0 + CD0i;

        dtab(end+1) = struct('name',fieldOr(c,'name',sprintf('comp%d',i)), ...
                             'Cf',Cf,'FF',FF,'Q',Q,'CD0i',CD0i);
    end
end

function Cf = cf_flatplate_turbulent(Re)
    if Re<=0, Cf = 0.003; return; end
    Cf = 0.455 / (log10(Re)^2.58);
end

function v = fieldOr(s,f,def)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = def; end
end
