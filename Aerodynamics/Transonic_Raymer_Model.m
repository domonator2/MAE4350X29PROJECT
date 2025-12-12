function out = Transonic_Raymer_Model(Minf, alpha_deg, Mcr, Mss, geom, subFcn, supFcn, opts)
%{
====================================================================================================
SCRIPT NAME: Transonic_Raymer_Model.m
AUTHOR: Ethan Anderson
INITIATED: 09/13/2025
LAST REVISION: 09/20/2025
====================================================================================================
SCRIPT DESCRIPTION:
This function consolidates Raymer-based transonic aerodynamics for preliminary X-29 analysis. It computes
transonic lift and lift-curve slope by fairing between subsonic and supersonic endpoints, transonic drag
using Raymer's A–B–C–D–E construction for the zero-lift drag-rise increment with a baseline that connects
CD0 at the drag-rise Mach to CD0 at Mach 1.20, estimates the drag-rise Mach using the Douglas 0.10-slope
definition applied to the Raymer curve when requested, and evaluates the A-point wave-drag coefficient at
Mach 1.20 from Raymer's formulation.

USER-DEFINED FUNCTIONS:
(subFcn) -> subsonic endpoint function, signature: out = subFcn(M, alpha_deg, geom, opts)
(supFcn) -> supersonic endpoint function, signature: out = supFcn(M, alpha_deg, geom, opts)
----------------------------------------------------------------------------------------------------
(Raymer_compute_subsonic_CD.m) and (Raymer_K_subsonic.m) roles are provided by subFcn.
(Raymer_compute_supersonic_CD.m) and (Raymer_K_supersonic.m) roles are provided by supFcn.
====================================================================================================
ABBREVIATIONS:
A, B, C, D, and E denote Raymer's transonic drag-rise at Mach 1.20, Mach 1.05, Mach 1.00, M_D, and Mcr.
====================================================================================================
ADDITIONAL COMMENTS:
The Sears–Haack effective length must exclude constant cross-section segments. Area and length units must
be used consistently.
====================================================================================================
PERMISSION:
Use of this code, in part or in full, is authorized for Ethan Anderson only. No other use, redistribution,
or modification is permitted without Ethan Anderson's explicit consent.
====================================================================================================
%}

    % -------- defaults for optional inputs --------
    if nargin < 6 || isempty(subFcn), subFcn = @aero_subsonic_raymer; end
    % [-] If a subsonic function handle is not provided, default to aero_subsonic_raymer.

    if nargin < 7 || isempty(supFcn), supFcn = @CD_supersonic_raymer_rev_wrapper; end
    % [-] Default to wrapped CD_supersonic_raymer_rev so signature matches (M, alpha_deg, geom, opts).

    if nargin < 8 || isempty(opts), opts = struct(); end
    % [-] Create an empty options struct if not supplied.

    if ~isfield(opts,'useDouglasMDD'), opts.useDouglasMDD = true; end
    % [logical] Default to using the Douglas 0.10-slope criterion for M_D if not specified.

    % ===================== TRANSONIC LIFT =====================
    [CL_trans, CLalpha_trans] = Raymer_compute_transonic_CL(Minf, alpha_deg, Mcr, Mss, geom, subFcn, supFcn, opts);
    % [] Compute transonic CL and CLα by fairing endpoint slopes and applying (α+iw−α0).

    % ===================== TRANSONIC DRAG =====================
    [CD_trans, CD0_trans, CDi_trans, K_trans] = Raymer_compute_transonic_CD(CL_trans, Minf, Mcr, Mss, geom, subFcn, supFcn, opts);
    % [] Compute zero-lift CD0 via A–B–C–D–E plus baseline, fair K(M), then CDi and total CD.

    % ===================== PACKAGE OUTPUTS =====================
    out = struct();
    % [-] Initialize output struct container.

    out.CL        = CL_trans;
    % [-] Total lift coefficient at Minf based on blended slope and geometry angles.

    out.CL_alpha  = CLalpha_trans;
    % [1/rad] Effective transonic lift-curve slope at Minf.

    out.CD        = CD_trans;
    % [-] Total drag coefficient at Minf (CD0 + CDi).

    out.CD0       = CD0_trans;
    % [-] Zero-lift drag coefficient at Minf (baseline + wave increment).

    out.CDi       = CDi_trans;
    % [-] Induced drag coefficient at Minf (K·CL^2).

    out.K         = K_trans;
    % [-] Induced-drag factor at Minf blended between subsonic and supersonic endpoints.

    out.L_over_D  = out.CL ./ out.CD;
    % [-] Lift-to-drag ratio at Minf.
end

% =================================================================================================
% LOCAL FUNCTIONS
% =================================================================================================

function [CL, CLalpha] = Raymer_compute_transonic_CL(Minf, alpha_deg, Mcr, Mss, geom, subFcn, supFcn, opts)
% [] Transonic lift via endpoint fairing

    % Clamp Minf into the fairing bracket to avoid extrapolation.
    Minf = max(min(Minf, Mss), Mcr);
    % [-] Keep interpolation inside [Mcr, Mss].

    sub = subFcn(Mcr, alpha_deg, geom, opts);
    % [-] Subsonic endpoint evaluation at Mcr using the provided subsonic function.

    sup = supFcn(Mss, alpha_deg, geom, opts);
    % [-] Supersonic endpoint evaluation at Mss using the provided supersonic function.

    CLalpha_sub = fetch_CLalpha(sub, alpha_deg);
    % [1/rad] Extract subsonic slope; prefer explicit CL_alpha/CLa, fallback to CL/α.

    CLalpha_sup = fetch_CLalpha(sup, alpha_deg);
    % [1/rad] Extract supersonic slope; prefer explicit CL_alpha/CLa, fallback to CL/α.

    CLalpha = pchip([Mcr Mss], [CLalpha_sub CLalpha_sup], Minf);
    % [1/rad] Shape-preserving fairing of slope vs Mach; evaluated at Minf.

    alpha     = deg2rad(alpha_deg);
    % [rad] Convert geometric angle of attack to radians.

    iw        = deg2rad(geom.iw);
    % [rad] Wing incidence angle from geom.

    alphazero = deg2rad(geom.alphazero);
    % [rad] Zero-lift angle from geom.

    CL = CLalpha * (alpha + iw - alphazero);
    % [-] Linear lift relation with geometry angles applied.
end

function [CD, CD0, CDi, K] = Raymer_compute_transonic_CD(CL, Minf, Mcr, Mss, geom, subFcn, supFcn, opts)
% [] Transonic drag:
%    1) Choose M_D (given, Douglas 0.10-slope, or fallback Mcr+0.08).
%    2) Build Raymer ΔCD0(M) A–B–C–D–E and sample at Minf.
%    3) Build zero-lift baseline CD0(M) by connecting sub @M_D to sup @1.20.
%    4) Fair K(M) between endpoints; compute CDi = K·CL^2 and total CD.

    % Clamp Minf into the fairing bracket.
    Minf = max(min(Minf, Mss), Mcr);
    % [-] Keep interpolation inside [Mcr, Mss].

    % ---- choose M_D ----
    if isfield(geom,'MDD') && ~isempty(geom.MDD)
        M_D = geom.MDD;
        % [-] Use user-specified drag-rise Mach if provided.
    elseif opts.useDouglasMDD
        M_D = estimate_MDD_Douglas_from_wavecurve(Mcr, geom);
        % [-] Solve M_D via Douglas 0.10-slope on the ΔCD0(M) curve.
    else
        M_D = Mcr + 0.08;
        % [-] Boeing/Raymer fallback for drag-rise Mach when not solving.
    end
    M_D = max(min(M_D, 1.20), Mcr);
    % [-] Constrain M_D to [Mcr, 1.20].

    % ---- wave-rise A–B–C–D–E ----
    CDw12 = Raymer_CDw12(geom);
    % [-] A-point wave-drag coefficient at M=1.20 from Sears–Haack + efficiency.

    M_A = 1.20; dA = CDw12;
    % [-] Point A (1.20, A) sets the upper anchor for ΔCD0.

    M_B = 1.05; dB = dA;
    % [-] Point B (1.05, A) per Raymer.

    M_C = 1.00; dC = 0.5*dA;
    % [-] Point C (1.00, 0.5·A) per Raymer.

    M_E = Mcr;  dE = 0.0;
    % [-] Point E (Mcr, 0) starts the rise.

    dD  = 0.002;
    % [-] Fixed small increment at point D per Raymer guidance.

    if Minf >= M_C && Minf <= M_B
        dCD0_wave = interp1([M_C M_B], [dC dB], Minf, 'linear');
        % [-] Use exact linear segment between C and B in the middle of the rise.
    else
        dCD0_wave = pchip([M_E M_D M_C M_B M_A], [dE dD dC dB dA], Minf);
        % [-] Elsewhere, use a shape-preserving spline across all five points.
    end

    % ---- zero-lift baseline from endpoints ----
    sub_at_D  = subFcn(M_D,  0.0, geom, opts);
    % [-] Subsonic evaluation at M_D with α=0 for baseline CD0.

    sup_at_A  = supFcn(1.20, 0.0, geom, opts);
    % [-] Supersonic evaluation at M=1.20 with α=0 for baseline CD0.

    CD0_at_MDD = fetch_CD0_zero(sub_at_D);
    % [-] Extract zero-lift drag from subsonic result at M_D.

    CD0_at_120 = fetch_CD0_zero(sup_at_A);
    % [-] Extract zero-lift drag from supersonic result at M=1.20.

    CD0_baseline = pchip([M_D 1.20], [CD0_at_MDD CD0_at_120], Minf);
    % [-] Linear connection of the two baseline points (2-point PCHIP = line).

    CD0 = CD0_baseline + dCD0_wave;
    % [-] Total zero-lift CD in transonic is baseline plus wave increment.

    % ---- induced factor K fairing ----
    K_sub = fetch_K(subFcn(Mcr, alpha_deg_for_K(opts), geom, opts), geom);
    % [-] Induced-drag factor at Mcr from subsonic endpoint (direct, e-based, or CDi/CL²).

    K_sup = fetch_K(supFcn(Mss, alpha_deg_for_K(opts), geom, opts), geom);
    % [-] Induced-drag factor at Mss from supersonic endpoint.

    K = pchip([Mcr Mss], [K_sub K_sup], Minf);
    % [-] Shape-preserving fairing of K vs Mach; evaluated at Minf.

    % ---- induced and total ----
    CDi = K * CL^2;
    % [-] Induced drag using blended K and current total lift CL.

    CD  = CD0 + CDi;
    % [-] Total drag = zero-lift + induced.
end

function MDD_Doug = estimate_MDD_Douglas_from_wavecurve(Mcr, geom)
% [] Douglas drag-rise Mach solver:
%    Find the first Mach where d(ΔCD0)/dM = 0.10 on the A–B–C–E curve (no D in construction).
%    This uses a PCHIP spline and a simple bisection search between Mcr and 1.20.

    CDw12 = Raymer_CDw12(geom);
    % [-] A-point value sets curve scale.

    M_A = 1.20; dA = CDw12;
    % [-] A (1.20, A).

    M_B = 1.05; dB = dA;
    % [-] B (1.05, A).

    M_C = 1.00; dC = 0.5*dA;
    % [-] C (1.00, 0.5·A).

    M_lo = max(0.70, Mcr - 0.20);
    % [-] Add a low-M anchor to stabilize the spline.

    pp   = pchip([M_lo Mcr M_C M_B M_A], [0 0 dC dB dA]);
    % [-] Build ΔCD0(M) spline through low anchor, E, C, B, A.

    dpp  = fnder(pp, 1);
    % [-] First derivative of ΔCD0(M) with respect to Mach.

    a = Mcr; b = 1.20;
    % [-] Search bracket for bisection.

    if ppval(dpp, a) >= 0.10, MDD_Doug = a; return; end
    % [-] If slope already ≥0.10 at Mcr, set M_D = Mcr.

    if ppval(dpp, b) <  0.10, MDD_Doug = b; return; end
    % [-] If slope never reaches 0.10 by 1.20, set M_D = 1.20.

    for k = 1:50
        % [-] Fixed-iteration bisection for robustness.
        m = 0.5*(a+b);
        % [-] Midpoint Mach for test.
        if ppval(dpp, m) >= 0.10, b = m; else, a = m; end
        % [-] Tighten bracket towards the first crossing of 0.10.
    end
    MDD_Doug = 0.5*(a+b);
    % [-] Return midpoint as the converged M_D estimate.
end

function CDw12 = Raymer_CDw12(geom)
% [] A-point wave-drag at M=1.20:
%    CDw12 = Ewd * (9π/2) * (Amax/Lref)^2 / Sref
%    Uses Sears–Haack minimum-wave (D/q) term scaled by efficiency and referenced to Sref.

    Amax = geom.Amax;
    % [m^2] Maximum total cross-sectional area of configuration.

    Lref = geom.Lref;
    % [m] Sears–Haack effective length (exclude constant-area segments).

    Sref = geom.Sref;
    % [m^2] Reference area for nondimensionalization.

    Ewd  = geom.Ewd;
    % [-] Wave-drag efficiency factor (empirical).

    Dq_SH_12 = (9*pi/2) * (Amax / Lref)^2;
    % [m^2] Sears–Haack (D/q) at M=1.20, proportional to (A/L)^2.

    CDw12 = (Ewd * Dq_SH_12) / Sref;
    % [-] Convert (D/q) to coefficient by dividing by Sref and scaling by Ewd.
end

% ===================== tiny helpers (explained in full) =====================

function CLalpha = fetch_CLalpha(res, alpha_deg)
% [] Extract lift-curve slope from an endpoint result.
%    Priority: (1) res.CL_alpha or res.CLa
%              (2) CL/α (guarded)

    if (isfield(res,'CL_alpha') && ~isempty(res.CL_alpha))
        CLalpha = res.CL_alpha;                          % [1/rad]
    elseif (isfield(res,'CLa') && ~isempty(res.CLa))
        CLalpha = res.CLa;                               % [1/rad] supersonic field name
    else
        arad = deg2rad(max(1e-6, alpha_deg));            % [rad]
        CLalpha = res.CL / arad;                         % [1/rad]
    end
end

function CD0z = fetch_CD0_zero(res)
% [] Zero-lift CD0 at α≈0.
%    Priority: (1) CD0_used  (subsonic)
%              (2) CD0_total (supersonic)
%              (3) CD0
%              (4) CD (as last resort when α=0 was used)

    if isfield(res,'CD0_used') && ~isempty(res.CD0_used)
        CD0z = res.CD0_used;                             % [-]
    elseif isfield(res,'CD0_total') && ~isempty(res.CD0_total)
        CD0z = res.CD0_total;                            % [-]
    elseif isfield(res,'CD0') && ~isempty(res.CD0)
        CD0z = res.CD0;                                  % [-]
    else
        CD0z = res.CD;                                   % [-]
    end
end

function K = fetch_K(res, geom)
% [] Extract induced-drag factor K from an endpoint result.
%    Priority: (1) res.K                      -> direct value if provided
%              (2) 1/(π A R e)                -> compute from Oswald e if available
%              (3) mean(CDi_clean)/mean(CL^2) -> back-solve using induced and lift
%              (4) 1/(π A R * 0.8)            -> conservative fallback if nothing else exists

    if isfield(res,'K') && ~isempty(res.K)
        K = res.K;
        % [-] Best case: the endpoint already computed K explicitly.
        return
    end

    if isfield(res,'e_oswald') && ~isempty(res.e_oswald) && isfield(geom,'AR')
        K = 1/(pi * geom.AR * res.e_oswald);
        % [-] Compute K from aspect ratio and Oswald efficiency factor.
        return
    end

    if isfield(res,'components') && isfield(res.components,'CDi_clean') && ~isempty(res.components.CDi_clean)
        CL2 = max(1e-9, mean(res.CL.^2));
        % [-] Average CL^2 with guard to avoid divide-by-zero.
        K   = mean(max(0, res.components.CDi_clean)) / CL2;
        % [-] Back-solve K from average induced drag over average CL^2.
        return
    end

    K = 1/(pi * max(geom.AR,1) * 0.8);
    % [-] Fallback using a conservative Oswald e≈0.8 if nothing else is known.
end

function aK = alpha_deg_for_K(opts)
% [] Provide a small nonzero α to evaluate endpoints when extracting K.
%    This avoids 0/0 if an endpoint only exposes CDi via CL (since CDi ∝ CL^2).
%    Overridable via opts.alpha_for_K.

    if isfield(opts,'alpha_for_K') && ~isempty(opts.alpha_for_K)
        aK = opts.alpha_for_K;
        % [deg] Use user-provided small angle if specified.
    else
        aK = 0.5;
        % [deg] Default tiny angle; small enough to remain linear, nonzero to stabilize ratios.
    end
end

% ===================== supersonic adapter (signature wrapper) =====================

function out = CD_supersonic_raymer_rev_wrapper(M, alpha_deg, geom, opts)
% [-] Signature adapter so transonic code can call the supersonic routine as (M, α, geom, opts).
%     Internally forwards to CD_supersonic_raymer_rev(M, α, geom, atm, opts).

    if ~isfield(opts,'atm') || ~all(isfield(opts.atm, {'rho','mu','a'}))
        % [-] Ensure atmosphere is provided in options.
        % [units] opts.atm.rho [kg/m^3], opts.atm.mu [Pa·s], opts.atm.a [m/s] when opts.units='SI'
        %         If opts.units='Imperial', your supersonic code handles conversion;
        %         still require fields to exist here for robustness.
        error('opts.atm with fields rho, mu, a is required for CD_supersonic_raymer_rev.');
        % [-] Clear, early failure if atmosphere data is missing or incomplete.
    end

    out = CD_supersonic_raymer_rev(M, alpha_deg, geom, opts.atm, opts);
    % [-] Forward call to the original supersonic solver with atmosphere split out of opts.
    % [units] M [-], alpha_deg [deg], geom fields in the same units family indicated by opts.units.
end
