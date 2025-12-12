%% plot_CL_CD_supersonic.m
% Demonstration script for unified supersonic lift and drag functions
%
%   This script shows how to compute and plot the coefficient of lift
%   ``C_L`` and coefficient of drag ``C_D`` as functions of Mach
%   number and angle of attack using the revised functions
%   ``CL_supersonic_raymer_rev`` and ``CD_supersonic_raymer_rev``.  It
%   sets up a representative flight condition, defines the airframe
%   geometry in SI units and computes lift and drag for a range of
%   angles of attack.  The results are plotted with clearly labelled
%   axes and legends and a quick parasite drag breakdown is printed.
%
%   The geometry and options sections are easy to edit so that users
%   can substitute their own aircraft data.  See the documentation
%   within the functions for parameter descriptions.


clear; clc; close all

%% ----- Flight condition (SI units) -----
% These values correspond to an altitude of roughly 11 km.  Modify as
% needed for your specific case.
atm.rho = 0.364;       % kg/m^3
atm.mu  = 1.46e-5;     % Pa·s
atm.a   = 295.0;       % m/s

%% ----- Analysis ranges -----
Mvec       = linspace(1.15, 2.00, 150);  % Mach numbers (row vector)
alphas_deg = 0:5:20;                     % Angles of attack [deg] (row vector)

%% ----- Airframe geometry (SI units) -----
% Core planform parameters for the X‑29 (converted from ft to m).
geom.Sref    = 185.0 * 0.09290304;  % m^2 (reference area)
geom.AR      = 4.0;                % aspect ratio
geom.sweepLE = -29.27;             % deg (forward sweep allowed)

% Skin‑friction component list.  Replace with your best estimates.
% Each entry must include Swet [m^2] and Lref [m].  Names are for
% labelling only.  Roughness height k_rough and Cf_override are
% optional.
geom.components = struct( ...
    'name',   {'fuselage','wing','canard','vtail'}, ...
    'Swet',   {29.7,       24.2,   6.5,     7.4}, ...    % m^2 (≈ 320,260,70,80 ft^2)
    'Lref',   {13.7,        2.20,  0.91,    1.52}, ...   % m   (≈ 45,7.22,3.0,5.0 ft)
    'k_rough',{[],          [],    [],      []} ...
);

% Wave drag inputs.  Provide aircraft length and max cross‑section area.
geom.length = 14.7;  % m (overall length; adjust as needed)
geom.Amax   = 0.70;  % m^2 (max cross‑section area; adjust as needed)

%% ----- Options for lift and drag -----
opts.units             = 'SI';        % use SI units throughout
opts.Ewave             = 2.0;         % wave‑drag efficiency factor
opts.CfModel           = 'RaymerTurbulent';
opts.Re_cutoff_enable  = true;        % apply Mach‑dependent Re cutoff
opts.Re_cutoff_k_default = 1.5e-6;     % default roughness height [m]

% Parameters for the lift function (passed via opts.CLParams)
opts.CLFcn    = @CL_supersonic_raymer_rev;
opts.CLParams = struct();
opts.CLParams.SexposedOverSref = 188.84/185.0;  % ≈1.0207 (dimensionless)
opts.CLParams.Ffuselage        = 1.0;           % if unknown, leave 1.0
opts.CLParams.ScaleProductCap  = 0.98;          % default cap
opts.CLParams.WarnIfOutOfRange = true;
opts.CLParams.units            = opts.units;    % match units

%% ----- Preallocate result arrays -----
nM = numel(Mvec);
nA = numel(alphas_deg);
CLmat = nan(nM, nA);
CDmat = nan(nM, nA);
CD0_table = []; % to capture parasite drag at alpha=0

%% ----- Compute lift and drag for each angle of attack -----
for j = 1:nA
    alpha_j = alphas_deg(j);
    % Compute lift coefficient using the new lift function.  The lift
    % function requires alpha vector same length as M or scalar.
    resCL = CL_supersonic_raymer_rev(Mvec, alpha_j, geom, opts.CLParams);
    CLmat(:,j) = resCL.CL;
    % Compute drag coefficient using the drag function.  This will call
    % the lift function internally for consistency.
    outCD = CD_supersonic_raymer_rev(Mvec, alpha_j, geom, atm, opts);
    CDmat(:,j) = outCD.CD_total;
    if j == 1
        % Record parasite breakdown for alpha=0
        CD0_table = [Mvec(:), outCD.CD0_fric(:), outCD.CD_wave(:), outCD.CD0_total(:)];
    end
end

%% ----- Plot C_L vs Mach for each angle -----
figure('Color','w'); hold on; grid on; box on;
colours = lines(nA);
for j = 1:nA
    plot(Mvec, CLmat(:,j), 'LineWidth', 1.8, 'Color', colours(j,:));
end
% Draw right‑half validity guideline
tanLE = tan(abs(geom.sweepLE)*pi/180);
M_min_right = sqrt(1 + tanLE^2);
xline(M_min_right, '--k', 'M_{min}', 'LabelOrientation','horizontal', ...
    'LabelVerticalAlignment','bottom');
xlabel('Mach number, M');
ylabel('C_L');
title('Supersonic C_L vs Mach (unified lift function)');
% Use TeX interpreter for Greek letter alpha and degree symbol.
legend(arrayfun(@(a) sprintf('\\alpha = %d^{\\circ}', a), alphas_deg, ...
    'UniformOutput', false), 'Location','northwest', 'Interpreter','tex');

%% ----- Plot C_D vs Mach for each angle -----
figure('Color','w'); hold on; grid on; box on;
for j = 1:nA
    plot(Mvec, CDmat(:,j), 'LineWidth', 1.8, 'Color', colours(j,:));
end
% Also plot parasite drag only (alpha=0) as dashed black
plot(CD0_table(:,1), CD0_table(:,4), 'k--', 'LineWidth', 1.5);
xline(M_min_right, '--k', 'M_{min}', 'LabelOrientation','horizontal', ...
    'LabelVerticalAlignment','bottom');
xlabel('Mach number, M');
ylabel('C_D');
title('Supersonic C_D vs Mach (unified lift & drag functions)');
% For drag legend, include parasite-only curve and use TeX interpreter.
legend([arrayfun(@(a) sprintf('\\alpha = %d^{\\circ}', a), alphas_deg, ...
    'UniformOutput', false), {'C_{D0} (\alpha=0)'}], ...
    'Location','northeast', 'Interpreter','tex');

%% ----- Console summary of parasite drag at alpha=0 -----
fprintf('\nParasite drag breakdown at \n alpha = 0 deg:\n');
T = array2table(CD0_table, 'VariableNames', {'M','CD0_fric','CD_wave','CD0_total'});
disp(T(1:5,:));
Swet_total = sum([geom.components.Swet]);
fprintf('Swet/Sref = %.3f\n', Swet_total/geom.Sref);
fprintf('Wave inputs: length = %.3f m, Amax = %.3f m^2, Ewave = %.2f\n', ...
    geom.length, geom.Amax, opts.Ewave);
fprintf('Right‑half validity M_min ≈ %.3f (|sweepLE| = %.2f deg)\n\n', ...
    M_min_right, abs(geom.sweepLE));