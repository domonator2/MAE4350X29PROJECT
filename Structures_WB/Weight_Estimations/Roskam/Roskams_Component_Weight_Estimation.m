% Roskam's Class II Method for Component Weight Estimation
% Aircraft - X29 
% Type - Fighter Jet

clear;clc;

%% Wing Weight Estimation

Kw = 1.00; %[] for fixed wing aircraft
n_ult = 6.5; %[G] ultimate load factor
W_GW = 13906; %[lbs] Gross Weight
sweep_LE = deg2rad(-29.27); %[rad] leading edge sweep angle
AR_wing = 4; %[] aspect ratio
lambda_wing = 0.4; % [] taper ratio
Sw = 185; %[ft^2] wing area
tc_max = 0.07; %[] max thickness-to-chord ratio

W_wing = (3.08 *( ((Kw*n_ult*W_GW)/tc_max) * ( (tan(sweep_LE) - 2*(1-lambda_wing)/(AR_wing*(1+lambda_wing)) )^2 +1)*1E-6 )^0.593 * (AR_wing*(1+lambda_wing))^0.89 *(Sw)^0.741);
% [] 0.80 is the fudge factor to account for composite materials used

fprintf('Wing Weight: %.2f lbs\n', W_wing);

%% Empennage Weight Estimation

% Vertical Tail
Kv = 1.0;
Sv = 33.75; %[ft^2] vertical tail area 
AR_v = 2.64; %[] vertical tail aspect ratio
[~, V, qD, ~] = flight_condition_units(40000, 'M', 0.9, 'EN'); % Vertical Dive Speed
VD = V/1.6878;
lambda_v = 0.306;
sweep_LE_v = deg2rad(47);
sweep_c2_v = tan(sweep_LE_v) - (1/AR_v)*(0.5)*(1-lambda_v)/(1+lambda_v);

W_v = Kv*Sv*(3.81 * ( (Sv^0.2)*VD/(1000*cos(sweep_c2_v)^0.5) ) - 0.287);
fprintf('Vertical Tail Weight: %.2f lbs\n', W_v);

% Canard
Sc = 37.0; %[ft^2]
Kc = 1.1; 
AR_c = 1.47;
lambda_c = 0.318;
sweep_LE_c = deg2rad(42);
sweep_c2_c = tan(sweep_LE_c) - (1/AR_c)*(0.5)*(1-lambda_c)/(1+lambda_c);

W_c = Kc*Sc*(3.81 * ( (Sc^0.2)*VD/(1000*cos(sweep_c2_c)^0.5) ) - 0.287);
fprintf('Canard Weight: %.2f lbs\n', W_c);

% Empennage
W_emp = W_v + W_c;
fprintf('Empennage Weight: %.2f lbs\n', W_emp);

%% Fuselage Weight Estimation
Kinl = 1.25; %[] for buried engine designs
lf = 48.0833; %[ft] fuselage max length
hf = 6; %[ft] fuselage max height
W_TO = 17800;

W_fus = 10.43*(Kinl^1.42) * (qD/100)^0.283 * (W_TO/1000)^0.95 * (lf/hf)^0.71;
fprintf('Fuselage Weight: %.2f lbs\n', W_fus);

%% Nacelles Weight Estimation
T_TO = 16000; %[lbs] Take off thrust

W_n = 0.055*T_TO;
fprintf('Nacelles Weight: %.2f lbs\n', W_n);

%% Landing Gear Weight Estimation
W_g = 62.61*(W_TO/1000)^0.84;
fprintf('Landing Gear Weight: %.2f lbs\n\n', W_g);

%% Structure Weight Estimation
W_str = W_wing + W_emp + W_fus + W_n + W_g;
fprintf('Structure Weight: %.2f lbs\n', W_str);

%% Engine Weight (known)
W_eng = 2282; %[lbs] engine weight
fprintf('Engine Weight: %.2f lbs\n', W_eng);

%% Systems Weight Estimation
Csys = 2000; %[lbs] Constant System Weight
fsys = 0.2; %[lbs/lbs] 
W_dry = 13906; % [lbs] dry weight

W_sys = Csys + fsys*W_dry;
fprintf('Systems Weight: %.2f lbs\n', W_sys);

%% Crew Provisions
W_ox = 16.9*(1)^1.494;
W_fur = 22.9*(1*qD/100)^0.743 + 107*(1*W_TO/100000)^0.585;

W_cp = W_ox + W_fur;
fprintf('Crew Provisions Weight: %.2f lbs\n\n', W_cp);

%% Operating Empty Weight Estimation
mu_a = 0;
OEW = (1 + mu_a)*(W_str + W_eng + W_sys + W_cp);
fprintf('Operating Empty Weight: %.2f lbs\n', OEW);

%% Take Off Gross Weight Estimation
W_crew = 180; %[lbf] Crew Weight
W_ppl = 3980; %[lbf] Fuel Weight
TOGW = OEW + W_crew + W_ppl;
fprintf('Take Off Gross Weight: %.2f lbs\n\n', TOGW);

%% Structural Index Calculation
S_wet = 10.7639*(120.1); % [ft^2] Wetted Surface
Istr = W_str/S_wet; % [lbs/ft^2] Structural Index
fprintf('Istr =  %.4f kg/m^2 \n', Istr*4.88); % [kg/m^2]

