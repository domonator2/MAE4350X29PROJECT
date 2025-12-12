function out = Raymer_supersonic_model(M, alpha_deg, geom, atm, opts)

    % Create geometry structure compatible with CD_supersonic_raymer_rev
    geom_raw = geom; % Use the provided geometry structure
    geom = struct();
    geom.Sref = geom_raw.wing.S;  % Reference area (wing area)
    geom.AR = geom_raw.wing.AR;   % Aspect ratio
    geom.sweepLE = geom_raw.wing.Lambda_LE_deg;  % Leading edge sweep in degrees
    geom.length = geom_raw.fuse.L;  % Overall aircraft length for wave drag
    geom.Amax = pi * (geom_raw.fuse.D/2)^2;  % Maximum cross-sectional area (fuselage)

    % Create components array for skin friction calculations
    geom.components = [];
    % Wing component
    geom.components(1).name = 'wing';
    geom.components(1).Swet = geom_raw.wing.S_wet;
    geom.components(1).Lref = geom_raw.wing.c_mac;  % Use MAC as reference length

    % Fuselage component  
    geom.components(2).name = 'fuselage';
    geom.components(2).Swet = geom_raw.fuse.S_wet;
    geom.components(2).Lref = geom_raw.fuse.L;

    % Canard component
    geom.components(3).name = 'canard';
    geom.components(3).Swet = geom_raw.canard.S_wet;
    geom.components(3).Lref = geom_raw.canard.c_mac;

    % Vertical tail component
    geom.components(4).name = 'vtail';
    geom.components(4).Swet = geom_raw.vtail.S_wet;
    geom.components(4).Lref = geom_raw.vtail.c_mac;
    
    out_CD = CD_supersonic_raymer_rev(M, alpha_deg, geom, atm, opts);
    out_CL = CL_supersonic_raymer_rev(M, alpha_deg, geom, opts);
    
    % Combine into unified output structure
    out = struct();
    out.CL = out_CL.CL;
    out.CL_alpha = out_CL.CLa;
    out.CD = out_CD.CD_total;
    out.CD0 = out_CD.CD0_total;
    out.K = out_CD.K;
    out.CDi = out.K * out.CL^2;

end

function out = CL_supersonic_raymer_rev(M, alpha_deg, geom, opts)
% CL_supersonic_raymer_rev  Unified supersonic lift coefficient and slope
%
%   This function returns both the lift‑curve slope ``CLa`` (per rad)
%   and the lift coefficient ``CL`` for a thin, trapezoidal wing in
%   fully supersonic flow.  It digitises Raymer's Fig. 12.7d for
%   taper ratio λ≈1/3 (right‑hand half only) and applies optional
%   scaling factors for exposed area and fuselage lift.  The same
%   inputs and units conventions used here are mirrored in the drag
%   function ``CD_supersonic_raymer_rev`` for a consistent workflow.
%
%   The implementation is valid only when ``tan(|Λ_LE|)/β ≤ 1``, i.e.
%   for Mach numbers ``M ≥ sqrt(1+tan^2(|Λ_LE|))``.  Inputs below
%   this threshold will result in NaN outputs and may trigger a
%   warning if requested.
%
% Inputs
%   M        – Mach number(s).  Numeric scalar or row vector.  Values
%              should be >1 for supersonic flow.  Non‑supersonic
%              inputs yield NaN outputs.
%   alpha_deg – Angle(s) of attack in degrees.  Scalar or row vector.
%              If a vector, a separate lift coefficient is computed
%              for each angle and Mach combination.  Small‑angle
%              assumptions apply.
%   geom     – Structure describing the wing geometry.  Required
%              fields:
%                * ``AR``      – Aspect ratio (dimensionless).
%                * ``sweepLE`` – Leading‑edge sweep angle in degrees
%                                (sign OK; magnitude is used).
%              Any additional fields are ignored by this routine but
%              may be present for compatibility with the drag
%              function.
%   opts     – (Optional) structure of name–value pairs controlling
%              scaling and behaviour.  Recognised fields:
%                * ``SexposedOverSref`` – Exposed area ratio (default 1.0).
%                * ``Ffuselage``        – Fuselage lift factor ``F``.  If
%                    omitted, it may be computed from diameter and span.
%                * ``FuselageDiameter`` – Fuselage diameter ``d`` (same
%                    length units as ``WingSpan``).
%                * ``WingSpan``         – Wing span ``b`` (same units as
%                    ``FuselageDiameter``).
%                * ``ScaleProductCap``  – Cap on ``SexposedOverSref*F``
%                    (default 0.98).  Set to ``Inf`` to disable.
%                * ``WarnIfOutOfRange`` – If true (default), warn when
%                    ``M < sqrt(1+tan^2(|Λ_LE|))``.
%                * ``units``            – Unit system ('SI' [default]
%                    or 'Imperial').  Controls interpretation of
%                    ``FuselageDiameter`` and ``WingSpan``.
%
% Outputs
%   out – Structure containing:
%         ``CL``          : Lift coefficient matrix (nMach × nAlpha).
%         ``CL_raw``      : Unscaled lift coefficient before applying
%                            area and fuselage factors.
%         ``CLa``         : Scaled lift‑curve slope (1 × nMach).
%         ``CLa_raw``     : Unscaled lift‑curve slope (1 × nMach).
%         ``beta``        : sqrt(M^2−1) compressibility parameter.
%         ``tanLE``       : tan(|Λ_LE|).
%         ``Atan``        : A * tan(|Λ_LE|).
%         ``x``           : tan(|Λ_LE|)/β (chart abscissa).
%         ``side``        : Cell array marking 'right' or 'left-needed'.
%         ``M_min_right_half`` : Minimum Mach for right‑half validity.
%         ``SexposedOverSref`` : Area ratio actually applied.
%         ``Ffuselage_used``   : Fuselage lift factor after capping.
%         ``Ffuselage_raw``    : Fuselage lift factor before capping.
%         ``d_over_b``         : Ratio used for fuselage factor (NaN if
%                                not computed).
%         ``note``        : Descriptive string summarising the method.
%
% This routine is designed to be readable and maintainable.  It
% normalises inputs, applies scaling consistently and returns
% diagnostic fields for plotting and debugging.  See Raymer
% Chapter 12 for theoretical background.

% Default options
if nargin < 4, opts = struct(); end
Sratio    = getfield_or_default(opts, 'SexposedOverSref', 1.0);
Ffuse_in  = getfield_or_default(opts, 'Ffuselage',       1.0);
d_in      = getfield_or_default(opts, 'FuselageDiameter', []);
b_in      = getfield_or_default(opts, 'WingSpan',        []);
cap       = getfield_or_default(opts, 'ScaleProductCap', 0.98);
warnFlag  = getfield_or_default(opts, 'WarnIfOutOfRange', true);
unitsFlag = getfield_or_default(opts, 'units', 'SI');

% Validate geometry
if ~isstruct(geom) || ~isfield(geom,'AR') || ~isfield(geom,'sweepLE')
    error('geom must contain fields AR and sweepLE.');
end
A           = geom.AR;
LambdaLE_deg= geom.sweepLE;

% Force row vectors for broadcasting
M         = M(:).';
alpha_deg = alpha_deg(:).';

% Compute compressibility and geometry terms
beta  = sqrt(M.^2 - 1.0);
tanLE = tan(abs(LambdaLE_deg) * pi/180);
Atan  = A .* tanLE;
x     = tanLE ./ beta;
Mmin  = sqrt(1 + tanLE.^2);

% Digitised Raymer Fig. 12.7d (λ≈1/3) right half for Atan rails [2,3]
xgrid = [1.00 0.90 0.80 0.70 0.60 0.50 0.40 0.30 0.20 0.10];
yA2   = [3.15 3.22 3.30 3.40 3.50 3.60 3.72 3.84 3.94 4.00];
yA3   = [3.45 3.52 3.60 3.70 3.78 3.86 3.92 3.97 4.00 4.00];

% Interpolate along x for each rail (reverse for increasing‑x interpolation)
xinc = fliplr(xgrid);
y20x = interp1(xinc, fliplr(yA2), min(max(x,xgrid(end)),xgrid(1)), 'pchip','extrap');
y30x = interp1(xinc, fliplr(yA3), min(max(x,xgrid(end)),xgrid(1)), 'pchip','extrap');

% Interpolate between Atan rails 2→3
t = (Atan - 2.0) ./ (3.0 - 2.0);
y_beta = (1 - t) .* y20x + t .* y30x;

% Unscaled slope (per rad)
CLa_raw = y_beta ./ beta;

% Compute fuselage lift factor if not explicitly given
d_over_b = NaN;
Ffuse_raw = Ffuse_in;
if isempty(Ffuse_in) || Ffuse_in == 1.0
    if ~isempty(d_in) && ~isempty(b_in) && all(isfinite([d_in(:); b_in(:)])) && all([d_in(:); b_in(:)]>0)
        if strcmpi(unitsFlag,'Imperial')
            ft2m = 0.3048;
            diam = d_in .* ft2m;
            span = b_in .* ft2m;
        else
            diam = d_in;
            span = b_in;
        end
        d_over_b = diam ./ span;
        Ffuse_raw = 1.07 .* (1 + d_over_b).^2;
    end
end

% Apply exposed area and fuselage scaling with cap
scale_raw  = Sratio .* Ffuse_raw;
if isfinite(cap)
    scale_used = min(scale_raw, cap);
else
    scale_used = scale_raw;
end
CLa = CLa_raw .* scale_used;

% Identify left‑half and invalid Mach points
needsLeft = x > 1;
invalidM  = M <= 1;
if any(needsLeft | invalidM)
    CLa_raw(needsLeft | invalidM) = NaN;
    CLa    (needsLeft | invalidM) = NaN;
    if warnFlag && any(needsLeft)
        % Use plain underscores in the identifier to avoid invalid escape
        warning(['CL_supersonic_raymer_rev: M below right‑half validity. ', ...
            'Some values have M < sqrt(1 + tan^2(|Λ_LE|)) ≈ %.3f. ', ...
            'Left half of Fig. 12.7 not digitised.'], Mmin);
    end
end

% Broadcast alpha to match M and compute lift coefficients
alpha_rad = deg2rad(alpha_deg);
nM = numel(M);
nA = numel(alpha_rad);
% Expand alpha for each Mach
if isscalar(alpha_rad)
    alpha_mat = alpha_rad .* ones(nM,1);
else
    if nA ~= nM
        % alpha vector of length nM or scalar only allowed
        error('alpha_deg must be scalar or have the same length as M.');
    end
    alpha_mat = alpha_rad(:); % column vector (nM×1)
end

% Compute raw and scaled CL (nM×1)
CL_raw_vec = CLa_raw(:) .* alpha_mat;
CL_vec     = CLa(:)     .* alpha_mat;

% Reshape outputs: always return nM×1 column vectors
CL_raw = CL_raw_vec(:);
CL     = CL_vec(:);


% Build output structure
out.CL          = CL;
out.CL_raw      = CL_raw;
out.CLa         = CLa;
out.CLa_raw     = CLa_raw;
out.beta        = beta;
out.tanLE       = tanLE;
out.Atan        = Atan;
out.x           = x;
out.side        = repmat({'right'}, size(x));
out.side(needsLeft | invalidM) = {'left-needed'};
out.M_min_right_half = Mmin;
out.SexposedOverSref = Sratio;
out.Ffuselage_used   = scale_used ./ max(Sratio, eps);
out.Ffuselage_raw    = scale_raw  ./ max(Sratio, eps);
out.d_over_b         = d_over_b;
out.note = sprintf(['Raymer Fig. 12.7d (λ≈1/3) right‑half used. ', ...
    'Atan=%.3f. Valid for M≥%.3f (|Λ_LE|=%.2f°).'], Atan, Mmin, abs(LambdaLE_deg));
end

%% Helper function: return field value or default
function v = getfield_or_default(s, fld, def)
    if isstruct(s) && isfield(s, fld) && ~isempty(s.(fld))
        v = s.(fld);
    else
        v = def;
    end
end

function out = CD_supersonic_raymer_rev(M, alpha_deg, geom, atm, opts)
% CD_supersonic_raymer_rev  Unified supersonic drag coefficient
%
%   This routine estimates the total drag coefficient ``CD`` for a
%   thin, trapezoidal wing in fully supersonic flow by summing
%   zero‑lift parasite drag, wave drag and drag due to lift.  It
%   implements Daniel Raymer’s quick methods with consistent units
%   handling and calls the lift function ``CL_supersonic_raymer_rev``
%   internally to obtain the lift coefficient ``CL(M,α)`` and slope
%   ``CLa``.  If desired, an alternative lift function can be passed
%   via ``opts.CLFcn`` along with parameters in ``opts.CLParams``.
%
%   Parasite drag is built up from turbulent flat‑plate skin friction
%   (Raymer Eq. 12.27) with an optional Mach‑dependent Reynolds
%   cutoff (Eq. 12.29).  Wave drag is estimated using the
%   Sears–Haack minimum‑wave drag formula scaled by ``opts.Ewave``.
%   Drag due to lift uses Raymer’s quick approximation (Eq. 12.51).
%   All inputs are converted to SI internally if ``opts.units`` is
%   'Imperial'.
%
% Inputs
%   M         – Mach number(s).  Numeric row vector.  All values
%               should be >1 for supersonic validity.  Any value ≤1
%               produces NaN outputs in the result.
%   alpha_deg – Angle(s) of attack in degrees.  Scalar or row vector
%               with the same length as ``M``.  Each Mach uses the
%               corresponding angle.  Small‑angle assumptions apply.
%   geom      – Struct describing the airframe geometry.  Required
%               fields:
%                 * ``Sref``    – Reference wing area [m² or ft²].
%                 * ``AR``      – Aspect ratio (dimensionless).
%                 * ``sweepLE`` – Leading‑edge sweep in degrees (sign OK).
%                 * ``components`` – Array of structs for skin friction.
%                   Each component must contain ``name``, ``Swet`` and
%                   ``Lref``.  Optional ``Cf_override`` and ``k_rough``
%                   override the default friction model and roughness.
%               Optional fields:
%                 * ``length`` – Overall aircraft length ℓ for wave drag.
%                 * ``Amax``   – Maximum cross‑section area for wave drag.
%               Additional fields are ignored.
%   atm       – Atmosphere structure with fields ``rho`` (density),
%               ``mu`` (dynamic viscosity) and ``a`` (speed of sound).
%               Units must be consistent with ``opts.units``.
%   opts      – (Optional) struct of name–value pairs controlling
%               calculations.  Recognised fields:
%                 * ``units``      – 'SI' (default) or 'Imperial'.
%                 * ``Ewave``      – Wave‑drag efficiency factor (default 2.0).
%                 * ``CfModel``    – Skin friction model ('RaymerTurbulent').
%                 * ``Re_cutoff_enable`` – Apply Reynolds cutoff (default true).
%                 * ``Re_cutoff_k_default`` – Default roughness height
%                       for cutoff when component lacks ``k_rough``.
%                       Given in metres (or feet) according to ``units``.
%                 * ``CLFcn``      – Function handle for computing lift.
%                       Default: @CL_supersonic_raymer_rev.
%                 * ``CLParams``   – Parameters structure passed to
%                       ``CLFcn``.  Fields depend on the chosen lift
%                       function.  For the default lift function these
%                       include ``SexposedOverSref``, ``Ffuselage``,
%                       ``FuselageDiameter``, ``WingSpan``,
%                       ``ScaleProductCap``, ``WarnIfOutOfRange`` and
%                       ``units``.
%
% Outputs
%   out – Structure containing:
%     ``CD_total``       : Total drag (nM × 1).
%     ``CD0_total``      : Zero‑lift drag (friction + wave) (nM × 1).
%     ``CD0_fric``       : Friction component of zero‑lift drag (nM × 1).
%     ``CD_wave``        : Wave drag portion (nM × 1).
%     ``K``              : Drag‑due‑to‑lift factor (1 × nM).
%     ``CL``             : Lift coefficient (nM × 1).
%     ``CLa``            : Lift‑curve slope (1 × nM).
%     ``Re``, ``Re_eff`` : Reynolds numbers and cut‑off values per
%                           component (nM × nComp).
%     ``Cf``             : Skin‑friction coefficients per component
%                           (nM × nComp).
%     ``CD0_components`` : Component contributions to friction drag
%                           (nM × nComp).
%     ``component_names``: Cell array of component names.
%     ``units_input``    : Echo of the ``units`` string.
%     ``notes``          : Descriptive string summarising methods used.
%

% Handle optional opts
if nargin < 5, opts = struct(); end
units      = getfield_or_default(opts, 'units', 'SI');
Ewave      = getfield_or_default(opts, 'Ewave', 2.0);
CfModel    = getfield_or_default(opts, 'CfModel', 'RaymerTurbulent');
ReCutEnable= getfield_or_default(opts, 'Re_cutoff_enable', true);
k_default  = getfield_or_default(opts, 'Re_cutoff_k_default', 1.5e-6);
CLFcn      = getfield_or_default(opts, 'CLFcn', @CL_supersonic_raymer_rev);
CLParams   = getfield_or_default(opts, 'CLParams', struct());

% Validate basic inputs
if ~isfield(geom,'Sref') || ~isfield(geom,'AR') || ~isfield(geom,'sweepLE')
    error('geom must contain Sref, AR and sweepLE.');
end
if ~isfield(geom,'components') || isempty(geom.components)
    error('geom.components must be a non‑empty array.');
end
if ~isfield(atm,'rho') || ~isfield(atm,'mu') || ~isfield(atm,'a')
    error('atm must contain rho, mu and a.');
end

% Enforce row vectors
M        = M(:).';
alpha_deg= alpha_deg(:).';
nM       = numel(M);

% Convert geometry and flow to SI if needed
[rho_SI, mu_SI, a_SI, geom_SI, k_SI] = convertToSI(units, atm, geom, k_default);
Sref = geom_SI.Sref;
A    = geom_SI.AR;
sweepLE = geom_SI.sweepLE;
len  = getfield_or_default(geom_SI,'length', NaN);
Amax = getfield_or_default(geom_SI,'Amax',   NaN);

% Mach‑dependent parameters
beta = sqrt(M.^2 - 1);


% --- Lift coefficient and slope ---
% Prepare geometry for lift function (AR and sweep only)
geomCL = struct('AR', A, 'sweepLE', sweepLE);

% Build options for lift function: ensure units and user CLParams
cl_opts = struct('units', units);
% Merge CLParams fields into cl_opts
fnCP = fieldnames(CLParams);
for ii = 1:numel(fnCP)
    cl_opts.(fnCP{ii}) = CLParams.(fnCP{ii});
end

% Call the lift function (assumed signature M, alpha_deg, geom, opts)
% The lift function returns a struct with CL and CLa fields
resCL = CLFcn(M, alpha_deg, geomCL, cl_opts);
if ~isstruct(resCL) || ~isfield(resCL,'CL')
    error('Lift function must return a struct with a field ''CL''.');
end
CL = resCL.CL;
if isfield(resCL,'CLa') && ~isempty(resCL.CLa)
    CLa = resCL.CLa;
else
    % Estimate slope from nonzero alpha
    alpha_rad = deg2rad(alpha_deg);
    [~, jBest] = min(abs(alpha_rad) + 1e-12);
    if abs(alpha_rad(jBest)) > 0
        CLa = CL ./ alpha_rad(jBest);
    else
        CLa = NaN(1,nM);
    end
end

% Ensure CL is nM×1 (alpha vector same length as M or scalar)
CL = reshape(CL, nM, []);
if size(CL,2) ~= 1
    error('Lift function should return CL of size nMach×1 when alpha has length nMach or is scalar.');
end

% --- Skin‑friction drag build‑up ---
nc      = numel(geom_SI.components);
Re      = nan(nM,nc);
Re_eff  = nan(nM,nc);
Cf_mat  = nan(nM,nc);
CD0_comp= nan(nM,nc);
names   = strings(1,nc);

for k = 1:nc
    c   = geom_SI.components(k);
    names(k) = string(getfield_or_default(c,'name',sprintf('comp%d',k)));
    Swet = requireField(c,'Swet', sprintf('components(%d).Swet', k));
    Lref = requireField(c,'Lref', sprintf('components(%d).Lref', k));
    Cf_override = getfield_or_default(c,'Cf_override',[]);
    k_rough     = getfield_or_default(c,'k_rough', k_SI);
    % Reynolds number
    Re(:,k) = (rho_SI .* (a_SI .* M) .* Lref) ./ mu_SI;
    % Cutoff if enabled
    if ReCutEnable && isfinite(k_rough) && k_rough > 0
        Rcut = 44.62 * ((Lref / k_rough) ^ 1.053) .* (M .^ 1.16);
        Re_eff(:,k) = min(Re(:,k), Rcut.');
    else
        Re_eff(:,k) = Re(:,k);
    end
    % Skin friction
    if ~isempty(Cf_override)
        Cf_mat(:,k) = Cf_override + zeros(nM,1);
    else
        switch lower(CfModel)
            case 'raymerturbulent'
                Re_eff_pos = max(Re_eff(:,k), 1e3);
                denomM = (1 + 0.144*(M.^2)).^0.65;
                Cf_mat(:,k) = 0.455 ./ ((log10(Re_eff_pos)).^2.58 .* denomM.');
                Cf_mat(:,k) = max(Cf_mat(:,k), 0);
            otherwise
                error('Unsupported CfModel: %s', CfModel);
        end
    end
    % Contribution to friction drag
    CD0_comp(:,k) = Cf_mat(:,k) * (Swet / Sref);
end

% Sum friction drag across components
CD0_fric = sum(CD0_comp, 2);

% --- Wave drag using Sears–Haack ---
if isfinite(len) && isfinite(Amax)
    Dq_SH   = (9*pi/2) * (Amax/len)^2;
    % Fixed: Wave drag now increases with Mach (not decreases)
    % Gentle linear increase: ~15% rise per 0.1 Mach above M=1.2
    Dq_wave = Ewave * (1 + 0.15*(M - 1.2)) * Dq_SH;
    CD_wave = Dq_wave / Sref;
else
    CD_wave = NaN(nM,1);
end
CD0_total = CD0_fric + CD_wave;

% --- Drag due to lift: Raymer quick K ---
LambdaLE_rad = deg2rad(abs(sweepLE));
den = 4*A*beta - 2;
K   = (A*beta.^2 .* cos(LambdaLE_rad)) ./ den;
invalid = (den <= 0) | (M <= 1);
if any(invalid)
    K(invalid)         = NaN;
    CD0_total(invalid) = NaN;
end

% --- Total drag ---
CD_lift  = K(:) .* (CL.^2);
CD_total = CD0_total + CD_lift;

% Populate output structure
out.CD_total       = CD_total;
out.CD0_total      = CD0_total;
out.CD0_fric       = CD0_fric;
out.CD_wave        = CD_wave;
out.K              = K;
out.CL             = CL;
out.CLa            = CLa;
out.Re             = Re;
out.Re_eff         = Re_eff;
out.Cf             = Cf_mat;
out.CD0_components = CD0_comp;
out.component_names= cellstr(names);
out.units_input    = units;
out.notes = sprintf(['Supersonic drag via Raymer quick method using lift from %s. ', ...
    'units=%s; CfModel=%s; Re_cutoff=%s; Ewave=%.3g; |Λ_LE|=%.2f°'], ...
    func2str(CLFcn), units, CfModel, ternary(ReCutEnable,'on','off'), Ewave, abs(sweepLE));
end

%% -----------------------------------------------------------------------
% Helper functions (getfield_or_default already defined above)

function x = requireField(s, fld, context)
    if ~isfield(s, fld) || isempty(s.(fld))
        error('Missing required field %s.', context);
    end
    x = s.(fld);
end

function [rho_SI, mu_SI, a_SI, geom_SI, k_SI] = convertToSI(unitsFlag, atm, geom, k_default)
    switch lower(unitsFlag)
        case 'si'
            rho_SI = atm.rho;
            mu_SI  = atm.mu;
            a_SI   = atm.a;
            geom_SI = geom;
            k_SI   = k_default;
        case 'imperial'
            ft2m = 0.3048;
            ft2m2 = ft2m^2;
            rho_SI = atm.rho * (14.59390294 / 0.028316846592);
            mu_SI  = atm.mu  * (14.59390294 / 0.3048);
            a_SI   = atm.a   * ft2m;
            geom_SI = geom;
            geom_SI.Sref = geom.Sref * ft2m2;
            if isfield(geom,'length') && ~isempty(geom.length)
                geom_SI.length = geom.length * ft2m;
            end
            if isfield(geom,'Amax') && ~isempty(geom.Amax)
                geom_SI.Amax   = geom.Amax   * ft2m2;
            end
            % Convert each component
            for kk = 1:numel(geom.components)
                geom_SI.components(kk) = geom.components(kk);
                if isfield(geom.components(kk),'Swet') && ~isempty(geom.components(kk).Swet)
                    geom_SI.components(kk).Swet = geom.components(kk).Swet * ft2m2;
                end
                if isfield(geom.components(kk),'Lref') && ~isempty(geom.components(kk).Lref)
                    geom_SI.components(kk).Lref = geom.components(kk).Lref * ft2m;
                end
                if isfield(geom.components(kk),'k_rough') && ~isempty(geom.components(kk).k_rough)
                    geom_SI.components(kk).k_rough = geom.components(kk).k_rough * ft2m;
                end
            end
            k_SI = k_default * ft2m;
        otherwise
            error('Unknown units flag: %s. Use ''SI'' or ''Imperial''.', unitsFlag);
    end
end

function y = ternary(cond, a, b)
    if cond, y = a; else, y = b; end
end
