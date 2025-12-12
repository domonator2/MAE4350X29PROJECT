%{
====================================================================================================
FUNCTION NAME: make_F404_thrust_model.m
AUTHORS: Carlos Gonzalez (Propulsion Team), Dominic Larin (Chief Engineer)
INITIATED: 9/7/2025
LAST REVISION: 9/28/2025
====================================================================================================
SCRIPT DESCRIPTION:
Builds a thrust interpolator for the GE F404-GE-400 engine:
    T = T(Mach, Altitude_ft, PLA_deg)  ->  thrust in lbf (or N when requested)

The function ingests a table `tbl` containing the F404 data, filters by test method
('WT' or 'AP'), and constructs a scatteredInterpolant surface over (Mach, alt_ft, PLA_deg).
Queries are clamped to the valid bounds. An optional “2.5-D” surrogate (bilinear in M,h
per PLA plane + PCHIP across PLA) is provided for quick plots — use the exact path
for simulations.

EXPECTED INPUT COLUMNS (after any upstream aliasing):
  • Mach, alt_ft, PLA_deg, method, net_thrust_lbf
  • method may be categorical/string/char: 'WT' or 'AP'

OUTPUT STRUCT:
  • .F                : scatteredInterpolant over [Mach, alt_ft, PLA_deg] -> thrust_lbf
  • .bounds           : query clamping limits for Mach / alt_ft / PLA_deg
  • .units            : struct('thrust','lbf','altitude','ft','PLA','deg')
  • .method           : 'WT' or 'AP'
  • .predict          : exact evaluator (scatteredInterpolant), units 'lbf' (default) or 'N'
  • .predict_fast     : alias to .predict (exact; robust and fast for this dataset)
  • .predict_approx   : optional 2.5-D surrogate (bilinear in M,h per PLA + PCHIP across PLA)

NOTES:
  - Exact path uses 'linear' scatteredInterpolant (identity at data points).
  - The 2.5-D surrogate is for visualization/quick looks; do not rely on it for <1% accuracy
    over sparse grids. Prefer .predict / .predict_fast for analysis and sims.

----------------------------------------------------------------------------------------------------
SECTION MAP (what each section does):
  • METHOD FILTER        : Normalizes method labels and subsets table to requested 'WT'/'AP'.
  • BASIC CHECKS         : Removes missing rows and validates non-empty subset.
  • EXACT INTERPOLANT    : Builds scatteredInterpolant and records query bounds.
  • PACKAGE MODEL        : Assembles the returned struct and function handles.
  • HELPERS (exact)      : Clamping, broadcasting, and unit conversion for the exact path.
  • HELPERS (approx 2.5D): Lightweight surrogate for plots (bilinear per-PLA + PCHIP across PLA).
%}

function model = make_F404_thrust_model(tbl, method)
% MAKE_F404_THRUST_MODEL Build GE F404-GE-400 thrust interpolator T(M,h,PLA)
%
% INPUTS:
%   tbl    - Table with F404 empirical data
%   method - Test method 'WT' or 'AP' (default: 'WT')
%
% OUTPUTS:
%   model  - Structure with interpolation functions and metadata

    % [-] Check if test method is specified, default to Mass Flow Temperature
    if nargin < 2 || isempty(method), method = 'WT'; end

    %% ===========================[ METHOD FILTER ]==========================
    % Normalize/validate METHOD column (handles categorical/string/char) and subset.
    % [string] Normalize method column to uppercase strings
    meth_col = upper(strtrim(string(tbl.method)));
    % [string] Normalize requested method to uppercase string
    want     = upper(strtrim(string(method)));
    % [-] Check if requested method exists in data
    if ~any(meth_col == want)
        % [-] Throw error for unavailable test method
        error('Requested method "%s" not found. Available: %s', ...
              method, strjoin(unique(meth_col).', ', '));
    % [-] End method availability check
    end

    %% ============================[ BASIC CHECKS ]==========================
    % Subset required columns and drop missing; ensure non-empty dataset.
    % [table] Extract required columns for specified test method
    sub = tbl(meth_col == want, {'Mach','alt_ft','PLA_deg','net_thrust_lbf'});
    % [table] Remove rows with missing data
    sub = rmmissing(sub);
    % [-] Ensure dataset is non-empty after filtering
    assert(~isempty(sub), 'No rows after filtering; check your table contents.');

    %% =========================[ EXACT INTERPOLANT ]========================
    % Build scatteredInterpolant in (Mach, alt_ft, PLA_deg) -> thrust_lbf and set bounds.
    % [scatteredInterpolant] Create 3D thrust interpolation function
    F = scatteredInterpolant( ...
            sub.Mach, sub.alt_ft, sub.PLA_deg, sub.net_thrust_lbf, ...
            'linear','nearest');

    % [-] Define Mach number interpolation bounds
    bounds.Mach   = [min(sub.Mach)    max(sub.Mach)];
    % [ft] Define altitude interpolation bounds
    bounds.alt_ft = [min(sub.alt_ft)  max(sub.alt_ft)];
    % [deg] Define power lever angle interpolation bounds
    bounds.PLA    = [min(sub.PLA_deg) max(sub.PLA_deg)];

    %% =============================[ PACKAGE MODEL ]========================
    % Return the model struct with evaluators and meta info.
    % [scatteredInterpolant] Store main interpolation function
    model.F       = F;
    % [struct] Store interpolation bounds for input validation
    model.bounds  = bounds;
    % [struct] Store unit specifications for outputs
    model.units   = struct('thrust','lbf','altitude','ft','PLA','deg');
    % [string] Store selected test method
    model.method  = want;

    % Exact evaluator (vectorized). Units: 'lbf' (default) or 'N'.
    % [function_handle] Create exact thrust prediction function
    model.predict = @(M,h_ft,PLA_deg,units) ...
                      local_predict_scattered(F, bounds, M, h_ft, PLA_deg, units);

    % Fast evaluator = exact (robust & fast for this problem size).
    % [function_handle] Create fast prediction function (alias to exact)
    model.predict_fast = model.predict;

    % Optional: 2.5-D surrogate (bilinear in M,h @ each PLA plane + PCHIP across PLA).
    % [function_handle] Create approximate prediction function for visualization
    model.predict_approx = @(M,h_ft,PLA_deg,units) ...
                      local_predict_approx_25D(sub, bounds, M, h_ft, PLA_deg, units);
% [-] End make_F404_thrust_model function
end

%% ====================== helpers (exact path) ==============================
function T = local_predict_scattered(F, bounds, M, h_ft, PLA_deg, units)
% LOCAL_PREDICT_SCATTERED Exact scatteredInterpolant evaluation with clamping & units
%
% INPUTS:
%   F       - scatteredInterpolant function
%   bounds  - Structure with interpolation bounds
%   M       - Mach number(s)
%   h_ft    - Altitude(s) in feet
%   PLA_deg - Power lever angle(s) in degrees
%   units   - Output units 'lbf' or 'N' (optional)
%
% OUTPUTS:
%   T       - Thrust in requested units

    % [-] Check if units are specified, default to pounds-force
    if nargin < 6 || isempty(units), units = 'lbf'; end
    % [varies] Clamp inputs to bounds and broadcast to common size
    [Mc,hc,PLAc] = local_clamp_and_broadcast(bounds, M, h_ft, PLA_deg);
    % [lbf] Evaluate thrust using scatteredInterpolant
    Tlbf = F(Mc, hc, PLAc);
    % [varies] Convert thrust to requested units
    T    = local_units(Tlbf, units);
% [-] End local_predict_scattered function
end

function [Mc, hc, PLAc] = local_clamp_and_broadcast(bounds, M, h_ft, PLA_deg)
% LOCAL_CLAMP_AND_BROADCAST Clamp inputs to bounds and broadcast to common size
%
% INPUTS:
%   bounds  - Structure with min/max bounds for each variable
%   M       - Mach number(s)
%   h_ft    - Altitude(s) in feet
%   PLA_deg - Power lever angle(s) in degrees
%
% OUTPUTS:
%   Mc      - Clamped and broadcasted Mach numbers
%   hc      - Clamped and broadcasted altitudes
%   PLAc    - Clamped and broadcasted power lever angles

    % [-] Clamp Mach number to valid interpolation range
    Mc   = max(min(M,      bounds.Mach(2)),   bounds.Mach(1));
    % [ft] Clamp altitude to valid interpolation range
    hc   = max(min(h_ft,   bounds.alt_ft(2)), bounds.alt_ft(1));
    % [deg] Clamp power lever angle to valid interpolation range
    PLAc = max(min(PLA_deg,bounds.PLA(2)),    bounds.PLA(1));

    % [cell] Store sizes of clamped inputs
    sizes  = {size(Mc), size(hc), size(PLAc)};
    % [-] Store number of elements in each input
    numelv = [numel(Mc), numel(hc), numel(PLAc)];
    % [-] Find indices of non-scalar inputs
    idxNS  = find(numelv>1);

    % [-] Determine target size for broadcasting
    if isempty(idxNS)
        % [-] All inputs are scalar
        target = [1 1];
    % [-] At least one input is non-scalar
    else
        % [-] Use size of first non-scalar input as target
        target = sizes{idxNS(1)};
        % [-] Verify all non-scalar inputs have same size
        for k = idxNS
            % [-] Check size compatibility
            if ~isequal(sizes{k}, target)
                % [-] Throw error for incompatible array sizes
                error('Inputs must be scalars or same-sized arrays.');
            % [-] End size compatibility check
            end
        % [-] End non-scalar input loop
        end
    % [-] End target size determination
    end

    % [-] Broadcast scalar inputs to target size
    if isscalar(Mc),   Mc   = repmat(Mc,   target); end
    % [ft] Broadcast scalar altitude to target size
    if isscalar(hc),   hc   = repmat(hc,   target); end
    % [deg] Broadcast scalar PLA to target size
    if isscalar(PLAc), PLAc = repmat(PLAc, target); end
% [-] End local_clamp_and_broadcast function
end

function T = local_units(Tlbf, units)
% LOCAL_UNITS Convert thrust from lbf (native) to requested units
%
% INPUTS:
%   Tlbf  - Thrust in pounds-force
%   units - Requested output units 'lbf' or 'N'
%
% OUTPUTS:
%   T     - Thrust in requested units

    % [string] Check if units are specified, default to pounds-force
    if nargin<2||isempty(units), units='lbf'; end
    % [string] Convert thrust based on requested units
    switch lower(units)
      % [lbf] Return thrust in pounds-force (no conversion)
      case 'lbf', T = Tlbf;
      % [N] Convert thrust from pounds-force to Newtons
      case 'n',   T = Tlbf * 4.4482216152605;
      % [-] Error for unsupported units
      otherwise,  error('units must be ''lbf'' or ''N''.');
    % [-] End unit conversion
    end
% [-] End local_units function
end

%% ================== helpers (optional approx 2.5-D) =======================
function T = local_predict_approx_25D(sub, bounds, M, h_ft, PLA_deg, units)
% LOCAL_PREDICT_APPROX_25D Bilinear (M,h) per PLA plane + PCHIP across PLA
%
% This function provides a 2.5D approximation method for visualization and
% quick analysis. It creates bilinear interpolation surfaces for each PLA
% setting, then uses PCHIP interpolation across PLA values.
%
% WARNING: For plots/quick checks only - not recommended for high-accuracy analysis
%
% INPUTS:
%   sub     - Filtered data table for specific test method
%   bounds  - Structure with interpolation bounds
%   M       - Mach number(s)
%   h_ft    - Altitude(s) in feet
%   PLA_deg - Power lever angle(s) in degrees
%   units   - Output units 'lbf' or 'N' (optional)
%
% OUTPUTS:
%   T       - Thrust in requested units

    % [string] Check if units are specified, default to pounds-force
    if nargin<6||isempty(units), units='lbf'; end

    % ---- Clamp queries to bounds
    % [-] Clamp Mach number to interpolation bounds
    M  = min(max(M,  bounds.Mach(1)),   bounds.Mach(2));
    % [ft] Clamp altitude to interpolation bounds
    h  = min(max(h_ft,bounds.alt_ft(1)),bounds.alt_ft(2));
    % [deg] Clamp power lever angle to interpolation bounds
    Pq = min(max(PLA_deg,bounds.PLA(1)),bounds.PLA(2));

    % ---- Remember original shape; vectorize queries
    % [-] Store original shape of input arrays
    shp = size(M);
    % [-] Check that all inputs have compatible sizes
    if ~isequal(size(h), shp) || ~isequal(size(Pq), shp)
        % [-] Throw error for incompatible input sizes
        error('Inputs M, h_ft, PLA_deg must be the same size or scalars.');
    % [-] End size compatibility check
    end
    % [-] Count total number of query points
    N = numel(M);
    % [-] Vectorize Mach number array
    Mv  = M(:);
    % [ft] Vectorize altitude array
    hv  = h(:);
    % [deg] Vectorize power lever angle array
    Pv  = Pq(:);

    % ---- Build per-PLA bilinear interpolants (over regular H×M grids)
    % [-] Number of Mach grid points for bilinear interpolation
    nM = 101;
    % [-] Number of altitude grid points for bilinear interpolation
    nH = 101;
    % [-] Create regular Mach number grid
    Mvec = linspace(bounds.Mach(1),   bounds.Mach(2),   nM);
    % [ft] Create regular altitude grid
    Hvec = linspace(bounds.alt_ft(1), bounds.alt_ft(2), nH);
    % [deg] Get unique PLA values from data (sorted)
    PLA_set = unique(sub.PLA_deg, 'sorted')';

    % [-] Count number of unique PLA values
    nPLA = numel(PLA_set);
    % [cell] Initialize cell array for bilinear interpolants per PLA
    G_planes = cell(1, nPLA);

% [cell] Create bilinear interpolants for each PLA plane
G_planes = cell(1, numel(PLA_set));
% [-] Loop through each unique PLA value
for ip = 1:numel(PLA_set)
    % [-] Create mask for current PLA value (with tolerance)
    mask = abs(sub.PLA_deg - PLA_set(ip)) < 1e-9;
    % [-] Extract Mach data for current PLA
    Mk = sub.Mach(mask);
    % [ft] Extract altitude data for current PLA
    Hk = sub.alt_ft(mask);
    % [lbf] Extract thrust data for current PLA
    Tk = sub.net_thrust_lbf(mask);
    % per-plane scattered fit (2-D) then sample to regular (H,M) grid
    % [scatteredInterpolant] Create 2D interpolant for current PLA plane
    Fk = scatteredInterpolant(Mk, Hk, Tk, 'linear', 'nearest');
    % [-] Create meshgrid for regular sampling
    [Mg, Hg] = meshgrid(Mvec, Hvec);
    % [lbf] Sample thrust on regular grid
    Tg = Fk(Mg, Hg);
    % [griddedInterpolant] Create bilinear interpolant for current PLA plane
    G_planes{ip} = griddedInterpolant({Hvec, Mvec}, Tg, 'linear', 'nearest');
% [-] End PLA plane loop
end

% Clamp + broadcast
% [-] Re-clamp Mach number to bounds (redundant but safe)
M  = min(max(M,  bounds.Mach(1)),   bounds.Mach(2));
% [ft] Re-clamp altitude to bounds (redundant but safe)
h  = min(max(h_ft,bounds.alt_ft(1)),bounds.alt_ft(2));
% [deg] Re-clamp PLA to bounds (redundant but safe)
Pq = min(max(PLA_deg,bounds.PLA(1)),bounds.PLA(2));
% [-] Get sizes of input arrays
[sM,sH,sP] = deal(size(M),size(h),size(Pq));
% [-] Determine target size for broadcasting
target = [max([sM(1),sH(1),sP(1)]), max([sM(2),sH(2),sP(2)])];
% [-] Broadcast scalar Mach to target size
if isscalar(M),  M  = repmat(M,  target); end
% [ft] Broadcast scalar altitude to target size
if isscalar(h),  h  = repmat(h,  target); end
% [deg] Broadcast scalar PLA to target size
if isscalar(Pq), Pq = repmat(Pq, target); end

% Evaluate bilinear on each PLA plane at (h,M)
% [-] Count number of PLA planes
np = numel(PLA_set);
% [lbf] Initialize thrust array for all PLA planes
Tplanes = zeros([size(M), np]);
% [-] Loop through each PLA plane
for ip = 1:np
    % [lbf] Evaluate bilinear interpolant at query points
    Tplanes(:,:,ip) = G_planes{ip}(h, M);
% [-] End PLA plane evaluation loop
end

% PCHIP across PLA at each (i,j)
% [deg] Convert PLA set to row vector for interpolation
PLA_vec = PLA_set(:).';
% [-] Get dimensions of query arrays
[nr,nc] = size(M);
% [lbf] Initialize final thrust array
Tq = zeros(nr,nc);
% [-] Loop through each row of query points
for i = 1:nr
    % [-] Loop through each column of query points
    for j = 1:nc
        % [lbf] Extract thrust values across all PLA planes at (i,j)
        Tij = squeeze(Tplanes(i,j,:));
        % [lbf] Interpolate across PLA using PCHIP
        Tq(i,j) = interp1(PLA_vec, Tij, Pq(i,j), 'pchip');
    % [-] End column loop
    end
% [-] End row loop
end

% Convert units
% [string] Convert thrust based on requested units
switch lower(units)
  % [lbf] Return thrust in pounds-force (no conversion)
  case 'lbf', T = Tq;
  % [N] Convert thrust from pounds-force to Newtons
  case 'n',   T = Tq * 4.4482216152605;
  % [-] Error for unsupported units
  otherwise,  error('units must be ''lbf'' or ''N''.');
% [-] End unit conversion
end
% [-] End local_predict_approx_25D function
end
