function [W_wing, W_v, W_c, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Weight_Calculator(C)

% Wing Weight Calculation
W_wing = (3.08*(((C.Kw * C.n_ult * C.W_OEW)/C.tc_max) * ( (tan(C.sweep_LE) - 2*(1-C.lambda_wing)/(C.AR_wing*(1+C.lambda_wing)) )^2 +1)*(10^-6) )^0.593 * (C.AR_wing*(1+C.lambda_wing))^0.89 *(C.Sw)^0.741);

% Empennage Weight Calculation
W_v = C.Kv*C.Sv*(3.81 * ( (C.Sv^0.2)*C.VD/(1000*cos(C.sweep_c2_v)^0.5) ) - 0.287);
W_c = C.Kc*C.Sc*(3.81 * ( (C.Sc^0.2)*C.VD/(1000*cos(C.sweep_c2_c)^0.5) ) - 0.287);
W_emp = W_v + W_c;

% Fuselage Weight Calculation
W_fus = 10.43*(C.Kinl^1.42) * (C.qD/100)^0.283 * (C.W_TO/1000)^0.95 * (C.lf/C.hf)^0.71;

% Nacelles Weight Calculation
W_n = 0.055*C.T_TO;

% Landing Gear Weight Calculation
W_g = (62.61*(C.W_TO/1000)^0.84);

% Structure Weight Calculation
W_str = W_wing + W_emp + W_fus + W_n + W_g;

% Crew Provision Weight Calculation
W_ox = 16.9*(C.N_crew)^1.494;
W_fur = 22.9*(C.N_crew*C.qD/100)^0.743 + 107*(1*C.W_TO/100000)^0.585;

W_cp = W_ox + W_fur;

end 