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