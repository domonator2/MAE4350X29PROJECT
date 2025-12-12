function Ra = Raymers_Method()

% Raymers Method Istr and Csys Calculations
CRa = Constants_Raymer();

% Structure Weight Calculations
[W_wing, W_vt, W_c, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Weight_Calculator_Raymer(CRa);

% Structural Index
Istr = W_str/CRa.S_wet;

% Constant System Weight Convergence
tolerance = 1e-10; % Convergence tolerance
max_iter = 100000; % Max iterations to prevent infinite loops

Csys_guess = 2765; %[lbs] Initial system weight guess

% Iterative convergence
iter = 0;
TOGW = 0;

while abs(TOGW - CRa.W_TO) > tolerance && iter < max_iter
    % Compute system weight
    W_sys = Csys_guess + CRa.fsys*CRa.W_dry;
    
    % Compute OEW
    OEW = W_str + W_sys + CRa.W_en + W_cp;
    
    % Compute TOGW
    TOGW = OEW + CRa.W_ppl + CRa.W_crew; %[lbs]
    
    % Csys convergence to TOGW
    Csys = Csys_guess - (TOGW - CRa.W_TO);
    
    iter = iter + 1;
end

OEW = W_str + W_sys + W_cp + CRa.W_en;
TOGW = OEW + CRa.W_ppl + CRa.W_crew;
 
%% Structure Weight

Ra.W_wing = W_wing;
Ra.W_v = W_vt;
Ra.W_c = W_c;
Ra.W_emp = W_emp;
Ra.W_fus = W_fus;
Ra.W_n = W_n;
Ra.W_g = W_g;
Ra.W_str = W_str;
Ra.W_cp = W_cp;
Ra.W_sys = W_sys;
Ra.Istr = Istr;
Ra.Csys = Csys;
Ra.OEW = OEW;
Ra.TOGW = TOGW;

end