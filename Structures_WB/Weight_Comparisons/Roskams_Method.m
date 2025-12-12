function Ro = Roskams_Method()

% Raymers Method Istr and Csys Calculations
CRo = Constants_Roskam();

% Structure Weight Calculations
[W_wing, W_v, W_c, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Weight_Calculator_Roskam(CRo);

% Structural Index
Istr = W_str/CRo.S_wet;

% Constant System Weight Convergence
tolerance = 1e-10; % Convergence tolerance
max_iter = 100000; % Max iterations to prevent infinite loops

Csys_guess = 2640; %[lbs] Initial system weight guess

% Iterative convergence
iter = 0;
TOGW = 0;

while abs(TOGW - CRo.W_TO) > tolerance && iter < max_iter
    % Compute system weight
    W_sys = Csys_guess + CRo.fsys*CRo.W_dry;
    
    % Compute OEW
    OEW = W_str + W_sys + CRo.W_eng + W_cp;
    
    % Compute TOGW
    TOGW = OEW + CRo.W_ppl + CRo.W_crew; %[lbs]
    
    % Csys convergence to TOGW
    Csys = Csys_guess - (TOGW - CRo.W_TO);
    
    iter = iter + 1;
end

OEW = W_str + W_sys + W_cp + CRo.W_eng;
TOGW = OEW + CRo.W_ppl + CRo.W_crew;
 
%% Structure Weight

Ro.W_wing = W_wing;
Ro.W_v = W_v;
Ro.W_c = W_c;
Ro.W_emp = W_emp;
Ro.W_fus = W_fus;
Ro.W_n = W_n;
Ro.W_g = W_g;
Ro.W_str = W_str;
Ro.W_cp = W_cp;
Ro.W_sys = W_sys;
Ro.Istr = Istr;
Ro.Csys = Csys;
Ro.OEW = OEW;
Ro.TOGW = TOGW;

end