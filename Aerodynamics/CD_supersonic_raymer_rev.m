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
    CD_wave = Ewave * Dq_SH / Sref + zeros(nM,1);
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
out.CD_lift        = CD_lift;
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
% Helper functions
function v = getfield_or_default(s, fld, def)
    if isstruct(s) && isfield(s, fld) && ~isempty(s.(fld))
        v = s.(fld);
    else
        v = def;
    end
end

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