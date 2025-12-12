%{
====================================================================================================
SCRIPT NAME: make_F404_tsfc_model.m
AUTHOR: Carlos Gonzalez
INITIATED: 09/07/2025
LAST REVISION: 09/26/2025
====================================================================================================
SCRIPT DESCRIPTION:
Builds a thrust specific fuel consumption (TSFC) interpolator for the GE F404-GE-400 engine:
    c = c(Mach, Altitude_ft, PLA_deg)  ->  TSFC in kg/(N*s)

The function accepts any of the following input column combinations in table `tbl`:
  • Direct TSFC:    tsfc_kg_per_N_s | tsfc_lbm_lbf_hr
  • Fuel + Thrust:  Wf_lbm_per_hr + net_thrust_lbf  |  Wf_kg_per_s + net_thrust_N
  • NASA columns:   enginecoreflowLbhr [+ afterburnerFuelFlowLbhr] + (net_thrust_lbf | net_thrust_N)

Robust features:
  • Column aliasing for Mach / altitude / PLA / method (see “ALIASES” below).
  • Automatic source selection (chooses the candidate with the most valid rows).
  • Optional “auto-clean” gating (PLA, Thrust, Fuel flow, TSFC range).
  • Duplicate-knot averaging, scatteredInterpolant with clamped queries.
  • SI outputs: TSFC [kg/(N*s)], fuel flow [kg/s], convenience unit [kg/(kN*hr)].

ALIASES HANDLED:
  - alt_ft | altitude_ft | alt_m | altitude_m   (meters auto-converted to feet)
  - PLA_deg | PLA | Throttle_deg
  - method  (filtered to 'WT' or 'AP' if present; default = 'WT')

----------------------------------------------------------------------------------------------------
SECTION MAP (what each section does):
  • CLEANING DEFAULTS      : Defines/tunes gates & thresholds (opts) with sensible defaults.
  • CONSTANTS              : Declares unit conversion constants (SI-first).
  • NORMALIZE INDEPENDENTS : Finds aliased columns, converts altitude to ft, filters by method.
  • BUILD CANDIDATES       : Produces all possible TSFC-in-SI vectors from available columns.
  • OPTIONAL AUTO-CLEAN    : Removes rows failing PLA/Thrust/Fuel/TSFC gates.
  • ASSEMBLE VECTORS       : Flattens data (after clean) into column vectors for interpolation.
  • CLEANING               : Drops non-finite/≤0 points; verifies non-empty dataset.
  • DEDUPLICATE KNOTS      : Averages duplicates at identical (M,h,PLA) nodes.
  • INTERPOLANT & BOUNDS   : Builds scatteredInterpolant and records query limits.
  • PREDICTORS & HELPERS   : Defines clamped predictor & unit helpers (fuel flow, kg/(kN*hr)).
  • REPORT                 : Summarizes source, counts, bounds, TSFC range, sanity notes.
  • PACKAGE                : Returns the model struct, units, methods, and report.
%}

function model = make_F404_tsfc_model(tbl, method, opts)
% make_F404_tsfc_model  Build GE F404-GE-400 TSFC interpolator c(M,h,PLA)
%
% OUTPUT STRUCT (SI):
%   model.F                       : scatteredInterpolant over [Mach, alt_ft, PLA_deg] -> TSFC [kg/(N*s)]
%   model.predict(M,h_ft,PLA_deg) : TSFC [kg/(N*s)]   (queries are safely clamped to bounds)
%   model.fuel_flow_kgps(M,h,PLA,T_N) : fuel flow [kg/s] for provided thrust [N]
%   model.tsfc_kg_per_kN_hr(M,h,PLA)  : TSFC in [kg/(kN*hr)] for reporting tables
%   model.bounds, model.units, model.method, model.report, model.info

    if nargin < 2 || isempty(method), method = 'WT'; end
    if nargin < 3 || isempty(opts),   opts = struct(); end

    %% =========================[ CLEANING DEFAULTS ]========================
    % Define default gates/thresholds; merge with user-supplied opts.
    def = struct('auto_clean', false, ...
                 'min_PLA_deg', 70, ...
                 'min_T_N', 2.0e4, ...     % [N]  (~4500 lbf)
                 'min_Wf_kgps', 0.15, ...  % [kg/s] (~1200 lbm/hr)
                 'tsfc_low', 7.0e-6, ...
                 'tsfc_high', 1.5e-4);
    fn = fieldnames(def);
    for k = 1:numel(fn)
        f = fn{k};
        if ~isfield(opts,f) || isempty(opts.(f)), opts.(f) = def.(f); end
    end

    %% =============================[ CONSTANTS ]=============================
    % Unit conversions (SI-first); used for candidate construction and reporting.
    lbm_to_kg = 0.45359237;
    lbf_to_N  = 4.4482216152605;
    hr_to_s   = 3600;
    m_to_ft   = 3.28083989501312;
    lbm_lbf_hr_to_kg_N_s = lbm_to_kg / (lbf_to_N * hr_to_s); % ~2.832e-5

    %% ===================[ NORMALIZE INDEPENDENT VARS ]=====================
    % Resolve column aliases, convert meters->feet, and (optionally) filter by method.
    V = string(tbl.Properties.VariableNames);
    find1 = @(names) find(ismember(lower(V), lower(string(names))), 1, 'first');

    iMach = find1("Mach");                         assert(~isempty(iMach), 'Table must contain "Mach".');
    iAlt  = find1(["alt_ft","altitude_ft","alt_m","altitude_m"]);  assert(~isempty(iAlt),  'Need alt_ft/alt_m.');
    iPLA  = find1(["PLA_deg","PLA","Throttle_deg"]);                assert(~isempty(iPLA), 'Need PLA_deg/PLA/Throttle_deg.');
    iMeth = find1("method");

    Mach = double(tbl.(V(iMach)));
    alt  = double(tbl.(V(iAlt)));
    if any(strcmpi(V(iAlt), ["alt_m","altitude_m"])), alt_ft = alt.*m_to_ft; else, alt_ft = alt; end
    PLA_deg = double(tbl.(V(iPLA)));

    if ~isempty(iMeth)
        keep = strcmpi(string(tbl.(V(iMeth))), string(method));
        if any(keep)
            tbl = tbl(keep,:);  Mach = Mach(keep,:);  alt_ft = alt_ft(keep,:);  PLA_deg = PLA_deg(keep,:);
        end
    end
    V = string(tbl.Properties.VariableNames);

    %% =====================[ BUILD TSFC CANDIDATE SET ]=====================
    % Build all viable TSFC-in-SI vectors from the available columns (and name them).
    cands = struct('name', {}, 'val', {});

    % (A) NASA core/AB + thrust
    has_core = ismember('enginecoreflowLbhr', V);
    has_AB   = ismember('afterburnerFuelFlowLbhr', V);
    has_T_lbf = ismember('net_thrust_lbf', V);
    has_T_N   = ismember('net_thrust_N',   V);

    if has_core && (has_T_lbf || has_T_N)
        core = double(tbl.enginecoreflowLbhr);
        AB   = zeros(height(tbl),1);
        if has_AB
            AB = double(tbl.afterburnerFuelFlowLbhr); AB(~isfinite(AB)) = 0;
        end
        core(~isfinite(core)) = NaN;
        Wf_lbm_per_hr = core + AB;
        if has_T_N
            cands(end+1).name = 'enginecoreflowLbhr[+AB] + net_thrust_N';
            cands(end  ).val  = (Wf_lbm_per_hr * lbm_to_kg / hr_to_s) ./ double(tbl.net_thrust_N);
        end
        if has_T_lbf
            cands(end+1).name = 'enginecoreflowLbhr[+AB] + net_thrust_lbf';
            cands(end  ).val  = (Wf_lbm_per_hr ./ double(tbl.net_thrust_lbf)) * lbm_lbf_hr_to_kg_N_s;
        end
    end

    % (B) Fuel + thrust (prefer SI pair first)
    if all(ismember({'Wf_kg_per_s','net_thrust_N'}, V))
        cands(end+1).name = 'Wf_kg_per_s + net_thrust_N';
        cands(end  ).val  = (double(tbl.Wf_kg_per_s) ./ double(tbl.net_thrust_N));
    end
    if all(ismember({'Wf_lbm_per_hr','net_thrust_lbf'}, V))
        cands(end+1).name = 'Wf_lbm_per_hr + net_thrust_lbf';
        cands(end  ).val  = (double(tbl.Wf_lbm_per_hr) ./ double(tbl.net_thrust_lbf)) * lbm_lbf_hr_to_kg_N_s;
    end

    % (C) Direct TSFC
    if ismember('tsfc_kg_per_N_s', V)
        cands(end+1).name = 'tsfc_kg_per_N_s';
        cands(end  ).val  = double(tbl.tsfc_kg_per_N_s);
    end
    if ismember('tsfc_lbm_lbf_hr', V)
        cands(end+1).name = 'tsfc_lbm_lbf_hr -> SI';
        cands(end  ).val  = double(tbl.tsfc_lbm_lbf_hr) * lbm_lbf_hr_to_kg_N_s;
    end

    % pick best by count of finite positive values
    best_idx = 0; best_count = -1;
    for i = 1:numel(cands)
        v = cands(i).val; ok = isfinite(v) & v > 0; cnt = sum(ok);
        if cnt > best_count, best_count = cnt; best_idx = i; end
    end
    if best_idx==0 || best_count==0
        error(['No usable TSFC source found. Ensure at least one exists: ', ...
               'tsfc_kg_per_N_s | tsfc_lbm_lbf_hr | (Wf_lbm_per_hr & net_thrust_lbf) | ', ...
               '(Wf_kg_per_s & net_thrust_N) | enginecoreflowLbhr[+afterburnerFuelFlowLbhr]+thrust.']);
    end
    tsfc_SI = cands(best_idx).val;
    src     = cands(best_idx).name;

    %% =========================[ OPTIONAL AUTO-CLEAN ]=======================
    % Apply PLA/Thrust/Fuel/TSFC gates to remove suspect rows when enabled.
    if opts.auto_clean
        % thrust (N)
        if ismember('net_thrust_N', V)
            TN = double(tbl.net_thrust_N);
        elseif ismember('net_thrust_lbf', V)
            TN = double(tbl.net_thrust_lbf)*lbf_to_N;
        else
            TN = NaN(height(tbl),1);
        end
        % fuel (kg/s)
        if ismember('Wf_kg_per_s', V)
            Wf_kgps = double(tbl.Wf_kg_per_s);
        elseif ismember('enginecoreflowLbhr', V)
            ABc = zeros(height(tbl),1);
            if ismember('afterburnerFuelFlowLbhr', V)
                ABc = double(tbl.afterburnerFuelFlowLbhr); ABc(~isfinite(ABc))=0;
            end
            Wf_kgps = (double(tbl.enginecoreflowLbhr)+ABc)*lbm_to_kg/hr_to_s;
        elseif ismember('Wf_lbm_per_hr', V)
            Wf_kgps = double(tbl.Wf_lbm_per_hr)*lbm_to_kg/hr_to_s;
        else
            Wf_kgps = (any(isfinite(TN))) .* tsfc_SI .* TN + (~any(isfinite(TN)))*NaN;
        end

        % gates
        gPLA = PLA_deg >= opts.min_PLA_deg;
        gT   = true(height(tbl),1);  if any(isfinite(TN)),      gT  = TN      > opts.min_T_N;     end
        gWF  = true(height(tbl),1);  if any(isfinite(Wf_kgps)), gWF = Wf_kgps >= opts.min_Wf_kgps; end
        gC   = isfinite(tsfc_SI) & tsfc_SI >= opts.tsfc_low & tsfc_SI <= opts.tsfc_high;

        idx_keep_tbl = gPLA & gT & gWF & gC;
        dropped_auto = sum(~idx_keep_tbl);
    else
        idx_keep_tbl = true(height(tbl),1);
        dropped_auto = 0;
    end

    %% ===========================[ ASSEMBLE VECTORS ]=======================
    % Flatten surviving rows into column vectors for interpolation.
    M   = Mach(:);      M   = M(idx_keep_tbl);
    hft = alt_ft(:);    hft = hft(idx_keep_tbl);
    PLA = PLA_deg(:);   PLA = PLA(idx_keep_tbl);
    cSI = tsfc_SI(:);   cSI = cSI(idx_keep_tbl);

    %% =============================[ CLEANING ]=============================
    % Remove non-finite/≤0 TSFC or independent values; fail fast if empty.
    good = isfinite(M) & isfinite(hft) & isfinite(PLA) & isfinite(cSI) & (cSI > 0);
    n_dropped_general = sum(~good);
    M = M(good); hft = hft(good); PLA = PLA(good); cSI = cSI(good);
    assert(~isempty(cSI), 'All TSFC rows became invalid after cleaning for source: %s', src);

    %% =====================[ DEDUPLICATE EXACT KNOTS ]=====================
    % Average multiple entries at identical (Mach, Alt_ft, PLA_deg) nodes.
    K = [M, hft, PLA];
    [Ku, ~, idxu] = unique(K, 'rows', 'stable');
    cSIu = accumarray(idxu, cSI, [], @mean);
    M = Ku(:,1); hft = Ku(:,2); PLA = Ku(:,3); cSI = cSIu;

    %% ==========================[ INTERPOLANT & BOUNDS ]===================
    % Build the scatteredInterpolant and record the valid query bounds.
    F = scatteredInterpolant(M, hft, PLA, cSI, 'linear', 'nearest');
    bounds.Mach    = [min(M) max(M)];
    bounds.alt_ft  = [min(hft) max(hft)];
    bounds.PLA_deg = [min(PLA) max(PLA)];

    %% ==========================[ PREDICTORS & HELPERS ]===================
    % Clamped predictor and convenience unit helpers for downstream users.
    predict = @(Mq,hq,PLAQ) F( clamp(Mq,bounds.Mach), ...
                               clamp(hq,bounds.alt_ft), ...
                               clamp(PLAQ,bounds.PLA_deg) );
    fuel_flow_kgps    = @(Mq,hq,PLAQ,thrust_N) predict(Mq,hq,PLAQ) .* thrust_N;      % [kg/s]
    tsfc_kg_per_kN_hr = @(Mq,hq,PLAQ) predict(Mq,hq,PLAQ) * (1000*hr_to_s);          % [kg/(kN*hr)]

    %% ===============================[ REPORT ]============================
    % Compose a human-readable summary (source, counts, bounds, TSFC range).
    cmin = min(cSI); cmax = max(cSI);
    txt  = {};
    txt{end+1} = sprintf('Source: %s', src);
    txt{end+1} = sprintf('Rows kept: %d  (dropped: %d)', numel(cSI), dropped_auto + n_dropped_general);
    txt{end+1} = sprintf('Unique knots: %d', size(Ku,1));
    txt{end+1} = sprintf('Bounds  Mach:[%.2f %.2f]  Alt(ft):[%.0f %.0f]  PLA:[%.1f %.1f]', ...
        bounds.Mach, bounds.alt_ft, bounds.PLA_deg);
    txt{end+1} = sprintf('TSFC range (kg/(N*s)): [%.3g, %.3g]', cmin, cmax);

    % Graded plausibility check using configured thresholds.
    LOW_TOL  = opts.tsfc_low;  HIGH_TOL = opts.tsfc_high;
    margin_lo = cmin - LOW_TOL;  margin_hi = HIGH_TOL - cmax;
    if (margin_lo >= 0) && (margin_hi >= 0)
        rng_msg = 'OK';
    elseif (margin_lo > -2e-6) && (margin_hi >= 0)
        rng_msg = 'OK (dry-lean edge)';
    else
        rng_msg = 'Check source units';
    end
    txt{end+1} = sprintf('Range sanity: %s', rng_msg);
    if opts.auto_clean
        txt{end+1} = sprintf('Auto-clean: PLA≥%.0f°, T>%.0f N, Wf≥%.2f kg/s, TSFC∈[%.2g, %.2g]', ...
            opts.min_PLA_deg, opts.min_T_N, opts.min_Wf_kgps, opts.tsfc_low, opts.tsfc_high);
    end

    %% =============================[ PACKAGE ]=============================
    % Return fully-formed model struct with predictors, bounds, and report.
    model = struct();
    model.F = F;
    model.bounds = bounds;
    model.units = struct('tsfc','kg/(N*s)','also','kg/(kN*hr)','fuel_flow','kg/s');
    model.method = string(method);
    model.predict = predict;
    model.predict_fast = predict;
    model.fuel_flow_kgps = fuel_flow_kgps;
    model.tsfc_kg_per_kN_hr = tsfc_kg_per_kN_hr;
    model.report = strjoin(txt, newline);
    model.info = struct('source',src,'dropped_rows',dropped_auto + n_dropped_general,'unique_knots',size(Ku,1));

    fprintf('\nTSFC MODEL [%s]\n%s\n', upper(model.method), model.report);
end

%% --------------------------------- helpers --------------------------------
function y = clamp(x, lims)
% clamp  Clamp numeric array x to [lims(1), lims(2)]
    y = min(max(x, lims(1)), lims(2));
end

