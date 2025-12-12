function CRo = Constants_Roskam ()
%% Wing Constants
CRo.Kw = 1.00; %[] for fixed wing aircraft
CRo.n_ult = 1.5*8; %[G] ultimate load factor
CRo.W_OEW = 13906; %[lbs] take off gross weight
CRo.W_TO = 17800; %[lbs] take off gross weight
CRo.sweep_LE = deg2rad(-29.27); %[rad] leading edge sweep angle
CRo.AR_wing = 4; %[] aspect ratio
CRo.lambda_wing = 0.4; % [] taper ratio
CRo.Sw = 185; %[ft^2] wing area
CRo.tc_max = 0.07; %[] max thickness-to-chord ratio

%% Empennage Constants
% Vertical Tail
CRo.Kv = 1.0;
CRo.Sv = 33.75; %[ft^2] vertical tail area 
CRo.AR_v = 2.64; %[] vertical tail aspect ratio
CRo.MD = 1.6; % []Max Mach Value
[~, V, CRo.qD, ~] = flight_condition_units(40000, 'M', CRo.MD, 'EN'); % Vertical Dive Speed
CRo.VD = V/1.6878;
CRo.lambda_v = 0.306;
CRo.sweep_LE_v = deg2rad(47);
CRo.sweep_c2_v = tan(CRo.sweep_LE_v) - (2/CRo.AR_v)*(1-CRo.lambda_v)/(1+CRo.lambda_v);

% Canard
CRo.Sc = 37.0; %[ft^2]
CRo.Kc = 1.1; 
CRo.AR_c = 1.47;
CRo.lambda_c = 0.318;
CRo.sweep_LE_c = deg2rad(42);
CRo.sweep_c2_c = tan(CRo.sweep_LE_c) - (2/CRo.AR_c)*(1-CRo.lambda_c)/(1+CRo.lambda_c);

%% Fuselage Constants
CRo.Kinl = 1.25; %[] for buried engine designs
CRo.lf = 48.0833; %[ft] fuselage max length
CRo.hf = 5.02; %[ft] fuselage max height

%% Nacelles Constants
CRo.T_TO = 16000; %[lbs] Take off thrust

%% Engine Constants
CRo.W_eng = 2282;
CRo.W_ppl = 3980;

%% Systems Constants
CRo.fsys = 0.2; %[lbs/lbs] 
CRo.W_dry = 13906; % [lbs] dry weight (equivalent to empty weight)

%% Crew Provisions
CRo.N_crew = 1; %[persons] just the pilot

%% Take off gross Constants
CRo.W_crew = 180; %[lbs] pilot and gear
CRo.W_ppl = 3980; %[lbs] fuel weight

%% Structural Index
CRo.S_wet = 1200; %[ft^2] total wetted area

end