function out = aero_subsonic_raymer(Mach, alpha_deg, geom, opts)
% AERO_SUBSONIC_RAYMER  Subsonic aerodynamics per Raymer Ch. 12 (wing-focused)
% Core outputs: CL, CD, CL_alpha (per rad), L_over_D
% ------------------------------------------------------------------------------
% Inputs
%   Mach      : scalar M (valid for subsonic M<1)
%   alpha_deg : scalar or vector α in degrees
%   geom      : struct with required fields (SI units)
%       Sref                 [m^2]  wing reference area
%       AR                   [-]    b^2/Sref
%       sweep_tmax_deg       [deg]  sweep at max thickness, Λ_tmax (for CLα)
%       sweep_LE_deg         [deg]  leading-edge sweep, Λ_LE (for e)
%       Sexposed_over_Sref   [-]    exposed wing area fraction (≤~1)
%       fuselage_d_over_b    [-]    fuselage dia/span (for F factor), set 0 if unknown
%       Swet                 [m^2]  wetted area (for Cfe·Swet/Sref path), optional if opts.CD0 supplied
%   opts      : struct with optional fields
%       alpha_L0_deg         [deg]  α_L=0 (default 0)
%       eta_airfoil          [-]    airfoil efficiency η (default 0.95)
%       CD0                  [-]    override parasite drag; else use Cfe*Swet/Sref
%       Cfe                  [-]    equivalent skin friction (if estimating CD0)
%       drag_buildup         struct array of components (optional, see helper below)
%       base_drag            struct('A_base', m^2)  -> adds CD_base
%       flaps                struct with high-lift device data (optional)
%         .deltas: struct for ΔCL due to devices (if known), else leave empty and we’ll only add ΔCD terms
%         .delta_flap_deg    [deg] trailing-edge flap deflection δ_flap (for ΔCD_flap)
%         .Cf_over_c         [-]   flap chord fraction (C_f/C) used in ΔCD_flap
%         .Sflap_over_Sref   [-]   flapped wing area fraction
%         .typeTE            char  'plain' or 'slotted' (ΔCD_flap uses Raymer 12.61 constants)
%         .span_is_full      logical true for full-span (k_f=0.14), false for half-span (k_f=0.28) in ΔCDi
%   NOTE: If you want component Cf via Reynolds number, add opts.fluid with rho, a, mu, V or qbar, and
%         per-component lengths/roughness; otherwise prefer Cfe path for quick conceptual runs.
% ------------------------------------------------------------------------------
% Returns (struct)
%   out.CL(:), out.CD(:), out.CL_alpha, out.L_over_D(:)
%   out.CD0_used, out.e_oswald, out.components (breakdown)
%   out.CLmax_clean (if provided/estimated), out.CLmax_highlift (if devices supplied)

    arguments
        Mach (1,1) double {mustBeNonnegative, mustBeLessThan(Mach,1)}
        alpha_deg (:,1) double
        geom struct
        opts struct = struct()
    end

    % ====== unpack geometry & options ======
    Sref  = need(geom,'Sref');
    AR    = need(geom,'AR');
    lam_t = deg2rad(need(geom,'sweep_tmax_deg'));
    lamLE = deg2rad(need(geom,'sweep_LE_deg'));
    Sexposed = getOr(geom,'Sexposed_over_Sref',1.0);
    d_over_b = getOr(geom,'fuselage_d_over_b',0.0);

    alpha_L0_deg = getOr(opts,'alpha_L0_deg',0.0);
    eta_airfoil  = getOr(opts,'eta_airfoil',0.95);

    % ====== CL_alpha (Raymer 12.6 with β, Λ_tmax, η) with exposure*F (12.9) ======
    beta = sqrt(max(1e-9, 1 - Mach^2));   % β = √(1-M^2)
    F_fus = 1.07*(1 + d_over_b)^2;        % Eq. 12.9
    exposure = min(0.98, Sexposed * F_fus);
    term = (AR*beta/eta_airfoil)^2 * (1 + (tan(lam_t)^2)/(beta^2));
    CL_alpha = (2*pi*AR) ./ (2 + sqrt(4 + term));
    CL_alpha = CL_alpha * exposure;        % suppress above ~1 per Raymer’s note

    % ====== Lift from α (linear) ======
    alpha_rad = deg2rad(alpha_deg - alpha_L0_deg);
    CL = CL_alpha .* alpha_rad;

    % ====== Oswald e (Raymer straight vs swept blend) ======
    e_oswald = e_oswald_raymer(AR, lamLE);

    % ====== Parasite drag CD0 ======
    if isfield(opts,'CD0') && ~isempty(opts.CD0)
        CD0 = opts.CD0;
        CD0_detail = struct('path','override','CD0',CD0);
    elseif isfield(opts,'drag_buildup') && ~isempty(opts.drag_buildup)
        [CD0, dtab] = cd0_component_buildup(opts.drag_buildup, Sref, Mach, opts);
        CD0_detail = struct('path','component_buildup','table',dtab);
    else
        Cfe  = getOr(opts,'Cfe',getOr(geom,'Cfe',NaN));
        Swet = getOr(geom,'Swet',NaN);
        if isnan(Cfe) || isnan(Swet)
            error('Provide opts.CD0, or opts.drag_buildup, or both geom.Swet and (opts.Cfe or geom.Cfe).');
        end
        CD0 = Cfe * (Swet/Sref);          % CD0 = Cfe·Swet/Sref
        CD0_detail = struct('path','Cfe_Swet','Cfe',Cfe,'Swet',Swet);
    end

    % ---- Base drag (Raymer 12.37) optional ----
    CD_base = 0;
    if isfield(opts,'base_drag') && isfield(opts.base_drag,'A_base') && ~isempty(opts.base_drag.A_base)
        CD_base = cd_base_raymer(Mach, opts.base_drag.A_base, Sref); % adds to CD0 path
    end

    % ====== Induced drag ======
    CDi_clean = (CL.^2) ./ (pi * AR * e_oswald);

    % ====== High-lift device increments (drag) ======
    dCD_flap = 0;
    dCDi_flap = 0;
    if isfield(opts,'flaps') && ~isempty(opts.flaps)
        fl = opts.flaps;
        % Parasitic drag of flaps (Raymer 12.61) – wing-referenced
        if all(isfield(fl, {'delta_flap_deg','Cf_over_c','Sflap_over_Sref','typeTE'}))
            dCD_flap = dCD_flap_raymer(fl.delta_flap_deg, fl.Cf_over_c, fl.Sflap_over_Sref, fl.typeTE);
        end
        % Induced drag increment from ΔCL due to flaps (Raymer 12.62)
        if isfield(fl,'deltas') && isfield(fl.deltas,'deltaCL_total') && ~isempty(fl.deltas.deltaCL_total)
            kf = 0.14; if isfield(fl,'span_is_full') && ~fl.span_is_full, kf = 0.28; end
            dCDi_flap = kf * (fl.deltas.deltaCL_total).^2 * cos(lamLE);
        end
    end

    % ====== Total drag and L/D ======
    CD = (CD0 + CD_base + dCD_flap) + (CDi_clean + dCDi_flap);
    L_over_D = CL ./ CD;

    % ====== Package outputs ======
    out.CL        = CL(:);
    out.CD        = CD(:);
    out.CL_alpha  = CL_alpha;
    out.L_over_D  = L_over_D(:);
    out.e_oswald  = e_oswald;
    out.CD0_used  = CD0 + CD_base + dCD_flap;
    out.components = struct('CD0_core',CD0_detail, 'CD_base',CD_base, ...
                            'dCD_flap',dCD_flap, 'dCDi_flap',dCDi_flap, ...
                            'CDi_clean',CDi_clean);

    % ====== Optional CLmax (clean & with devices) hooks ======
    out.CLmax_clean    = getOr(opts,'CLmax_clean',[]);
    out.CLmax_highlift = [];
    if isfield(opts,'flaps') && isfield(opts.flaps,'deltas') ...
       && isfield(opts.flaps.deltas,'deltaCL_max_total') && ~isempty(getOr(opts.flaps.deltas,'deltaCL_max_total',[]))
        if ~isempty(out.CLmax_clean)
            out.CLmax_highlift = out.CLmax_clean + opts.flaps.deltas.deltaCL_max_total;
        end
    end
end

% ===================== Helper functions =====================

function val = need(s, field)
    if ~isfield(s,field) || isempty(s.(field)), error('Missing required field: %s', field); end
    val = s.(field);
end

function val = getOr(s, field, default)
    if isstruct(s) && isfield(s,field) && ~isempty(s.(field)), val = s.(field); else, val = default; end
end

function e = e_oswald_raymer(AR, lamLE)
% Oswald e: blend straight vs swept (Raymer 12.48–12.49 pattern)
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

function CD_base = cd_base_raymer(M, A_base, Sref)
% Base drag (Raymer 12.37): (D/q)_base = [0.139 + 0.419(M-0.161)^2] A_base
    D_over_q = (0.139 + 0.419*(M - 0.161)^2) * A_base;
    CD_base  = D_over_q / Sref;
end

function dCD = dCD_flap_raymer(delta_flap_deg, Cf_over_c, Sflap_over_Sref, typeTE)
% Flap parasitic drag increment (Raymer 12.61):
% ΔCD_flap = F_flap * (C_f/C) * (S_flapped/Sref) * (δ_flap - 10)   (δ in degrees)
% F_flap: 0.0144 for plain flaps; 0.0074 for slotted flaps (per Raymer text figure)
    switch lower(strtrim(typeTE))
        case 'plain',   Fflap = 0.0144;
        case 'slotted', Fflap = 0.0074;
        otherwise, warning('Unknown flap type "%s"; using slotted constant.', typeTE); Fflap = 0.0074;
    end
    dCD = Fflap * (Cf_over_c) * (Sflap_over_Sref) * (delta_flap_deg - 10);
    dCD = max(0, dCD); % guard for small deflections
end

% ---------- Optional: Detailed component buildup (Raymer 12.30–12.33) ----------
function [CD0, dtab] = cd0_component_buildup(components, Sref, Mach, opts)
% components(i) fields expected (SI units unless noted):
%   .name        char
%   .S_wet       [m^2]
%   .Cf_method   struct  -> either:
%       A) 'equiv': struct('Cfe',value)                           -> Cf = Cfe
%       B) 'flatplate': struct('Re',value,'k_over_l',opt,'smooth',bool)
%   .FF_method   struct  -> one of:
%       'wing': struct('x_c_max',[],'t_c',[],'Lambda_m_deg',[])
%       'fuselage': struct('f',[])    where f = l/d
%       'nacelle' : struct('f',[])
%   .Q           [-] interference factor (default 1.0)
%
% Cf roughness handling: Subsonic roughness cutoff Re_cutoff = 38.21*(l/k)^1.053.
% If Re > Re_cutoff, Cf is frozen at Cf(Re_cutoff).
    CD0 = 0;
    dtab = struct('name',{},'Cf',{},'FF',{},'Q',{},'CD0i',{});
    for i = 1:numel(components)
        c = components(i);
        Q = fieldOr(c,'Q',1.0);
        FF = 1.0;

        % --- Cf
        if isfield(c,'Cf_method') && isfield(c.Cf_method,'Cfe')
            Cf = c.Cf_method.Cfe;   % equivalent skin friction path
        elseif isfield(c,'Cf_method') && isfield(c.Cf_method,'Re')
            Re = c.Cf_method.Re;
            smooth = fieldOr(c.Cf_method,'smooth',true);
            k_over_l = fieldOr(c.Cf_method,'k_over_l',0);  % surface roughness ratio
            Cf = cf_flatplate_turbulent(Re);
            if ~smooth && k_over_l>0
                Re_cut = 38.21*(1./k_over_l)^1.053; % Raymer roughness cutoff
                Cf = cf_flatplate_turbulent(min(Re,Re_cut)); % freeze Cf at Re_cutoff
            end
        else
            error('Component %d missing Cf_method.', i);
        end

        % --- FF (form factor): Raymer 12.30–12.33
        if isfield(c,'FF_method')
            if isfield(c.FF_method,'x_c_max') % wing/tail/strut/pylon
                xcm  = c.FF_method.x_c_max;   % (x/c)_m location of t/c(max)
                tc   = c.FF_method.t_c;       % thickness ratio
                LamM = deg2rad(c.FF_method.Lambda_m_deg); % sweep of max thickness line
                FF = (1 + (0.6/xcm)*tc + 100*tc^4) * (1.34*Mach^0.18 * cos(LamM)^0.28);
            elseif isfield(c.FF_method,'f') && isfield(c,'type') && strcmpi(c.type,'fuselage')
                f = c.FF_method.f;            % fineness ratio l/d
                FF = 0.9 + 5/(f^1.5) + f/400; % Raymer 12.31
            elseif isfield(c.FF_method,'f') % nacelle / smooth external store 12.32
                f = c.FF_method.f;
                FF = 1 + 0.35/f;
            else
                FF = 1.0;
            end
        end

        % --- contribution
        CD0i = Cf * FF * Q * (c.S_wet / Sref);
        CD0  = CD0 + CD0i;

        dtab(end+1) = struct('name',fieldOr(c,'name',sprintf('comp%d',i)), ...
                             'Cf',Cf,'FF',FF,'Q',Q,'CD0i',CD0i); %#ok<AGROW>
    end
end

function v = fieldOr(s,f,def)
    if isfield(s,f) && ~isempty(s.(f)), v = s.(f); else, v = def; end
end

function Cf = cf_flatplate_turbulent(Re)
% Fully turbulent flat-plate (e.g., Raymer): Cf = 0.455 / (log10 Re)^2.58
    if Re<=0, Cf = 0.003; return; end
    Cf = 0.455 / (log10(Re)^2.58);
end
