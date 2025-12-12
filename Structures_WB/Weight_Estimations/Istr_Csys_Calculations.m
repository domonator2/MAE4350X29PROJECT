% Structural Index and Constant System Weight Ranges (optimization)
clear;clc;

% Weight Component Constants
C = Constants();

% Weight Component Calculator for Structures and Crew Provisions
[W_wing, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Weight_Calculator(C);


% Structural Index Range
S_wet = linspace(800,1300,100); %[ft^2] range of the wet surface area

Istr = zeros(size(S_wet));

for i = 1:length(S_wet)
    Istr(i) = W_str / S_wet(i);
end

S_wet_metric = S_wet*0.092903; %[m^2]
Istr_metric = Istr*4.88; % [kg/m^2]

% Constant System Weight Convergence
tolerance = 1e-10; % Convergence tolerance
max_iter = 1000; % Max iterations to prevent infinite loops

Csys_guess = 2000; %[lbs] Initial system weight guess

% Iterative convergence
iter = 0;
TOGW = 0;

while abs(TOGW - C.W_TO) > tolerance && iter < max_iter
    % Compute system weight
    W_sys = Csys_guess + C.fsys*C.W_dry;
    
    % Compute OEW
    OEW = W_str + W_sys + C.W_eng + W_cp;
    
    % Compute TOGW
    TOGW = OEW + C.W_ppl + C.W_crew; %[lbs]
    
    % Csys convergence to TOGW
    Csys = Csys_guess - (TOGW - C.W_TO);
    
    iter = iter + 1;
end

% Constant System Weight Range
Csys_range = linspace(1500,2500,100); %[lbs] range of constant system weight

Wsys_range = zeros(size(Csys_range));
W_OEW = zeros(size(Csys_range));
TOGW_range = zeros(size(Csys_range));

for i = 1:length(Csys_range)
    Wsys_range(i) = Csys_range(i) + C.fsys*C.W_dry;
    W_OEW(i) = W_str + Wsys_range(i) + C.W_eng + W_cp;
    TOGW_range(i) = W_OEW(i) + C.W_ppl + C.W_crew;
end

fprintf('Structure Weight Components (lbs)\n');
fprintf('---------------------------------\n');
fprintf('W_wing = %.4f\n', W_wing);
fprintf('W_emp = %.4f\n', W_emp);
fprintf('W_fus = %.4f\n', W_fus);
fprintf('W_n = %.4f\n', W_n);
fprintf('W_g = %.4f\n\n', W_g);

fprintf('Structure Weight (lbs)\n');
fprintf('----------------------\n');
fprintf('W_str = %.4f\n\n', W_str);

fprintf('Engine Weight (lbs)\n');
fprintf('-------------------\n');
fprintf('W_eng = %.2f\n\n', C.W_eng);

fprintf('Crew and Crew Provisions Weight (lbs)\n');
fprintf('-------------------------------------\n');
fprintf('W_crew = %.2f\n', C.W_crew);
fprintf('W_cp = %.4f\n\n', W_cp);

fprintf('System Weight (lbs)\n');
fprintf('-------------------\n');
fprintf('W_sys = %.2f\n\n', W_sys);

fprintf('Constant System Weight (lbs)\n');
fprintf('----------------------------\n');
fprintf('C_sys = %.2f\n\n', Csys);

Plot_Csys_vs_TOGW(TOGW_range, Csys_range);

fprintf('Operating Empty Weight (lbs)\n');
fprintf('----------------------------\n');
fprintf('W_OEW = %.2f\n\n', OEW);

fprintf('Take Off Gross Weight (lbs)\n');
fprintf('----------------------------\n');
fprintf('W_TOGW = %.2f\n\n', TOGW);

fprintf('Structural Index Range (kg/m^2)\n');
fprintf('-------------------------------\n');
fprintf('%.2f <= Istr <= %.2f\n', min(Istr)*4.88, max(Istr)*4.88);
fprintf('Mean Istr = %.2f\n\n', mean(Istr)*4.88)

Plot_Istr_vs_Swet(S_wet_metric, Istr_metric);