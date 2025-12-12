function [W_wing, W_v, W_c, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Weight_Calculator_Roskam(CRo)

% Wing Weight Calculation
W_wing = (3.08*(((CRo.Kw * CRo.n_ult * CRo.W_TO)/CRo.tc_max) * ( (tan(CRo.sweep_LE) - 2*(1-CRo.lambda_wing)/(CRo.AR_wing*(1+CRo.lambda_wing)) )^2 +1)*(10^-6) )^0.593 * (CRo.AR_wing*(1+CRo.lambda_wing))^0.89 *(CRo.Sw)^0.741);

% Empennage Weight Calculation
W_v = CRo.Kv*CRo.Sv*(3.81 * ( (CRo.Sv^0.2)*CRo.VD/(1000*cos(CRo.sweep_c2_v)^0.5) ) - 0.287);
W_c = CRo.Kc*CRo.Sc*(3.81 * ( (CRo.Sc^0.2)*CRo.VD/(1000*cos(CRo.sweep_c2_c)^0.5) ) - 0.287);
W_emp = W_v + W_c;

% Fuselage Weight Calculation
W_fus = 10.43*(CRo.Kinl^1.42) * (CRo.qD/100)^0.283 * (CRo.W_TO/1000)^0.95 * (CRo.lf/CRo.hf)^0.71;

% Nacelles Weight Calculation
W_n = 0.055*CRo.T_TO;

% Landing Gear Weight Calculation
W_g = (62.61*(CRo.W_TO/1000)^0.84);

% Structure Weight Calculation
W_str = W_wing + W_emp + W_fus + W_n + W_g;

% Crew Provision Weight Calculation
W_ox = 16.9*(CRo.N_crew)^1.494;
W_fur = 22.9*(CRo.N_crew*CRo.qD/100)^0.743 + 107*(1*CRo.W_TO/100000)^0.585;

W_cp = W_ox + W_fur;

end 