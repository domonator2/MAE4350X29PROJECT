function C = Constants ()
%% Wing Constants
C.Kw = 1.00; %[] for fixed wing aircraft
C.n_ult = 1.5*8; %[G] ultimate load factor
C.sweep_LE = deg2rad(-29.27); %[rad] leading edge sweep angle
C.AR_wing = 4; %[] aspect ratio
C.lambda_wing = 0.4; % [] taper ratio
C.Sw = 185; %[ft^2] wing area
C.tc_max = 0.07; %[] max thickness-to-chord ratio

%% Empennage Constants
% Vertical Tail
C.Kv = 1.0;
C.Sv = 33.75; %[ft^2] vertical tail area 
C.AR_v = 2.64; %[] vertical tail aspect ratio
C.MD = 1.6; % []Max Mach Value
[~, V, C.qD, ~] = flight_condition_units(40000, 'M', C.MD, 'EN'); % Vertical Dive Speed
C.VD = V/1.6878;
C.lambda_v = 0.306;
C.sweep_LE_v = deg2rad(47);
C.sweep_c2_v = tan(C.sweep_LE_v) - (1/C.AR_v)*(0.5)*(1-C.lambda_v)/(1+C.lambda_v);

% Canard
C.Sc = 37.0; %[ft^2]
C.Kc = 1.1; 
C.AR_c = 1.47;
C.lambda_c = 0.318;
C.sweep_LE_c = deg2rad(42);
C.sweep_c2_c = tan(C.sweep_LE_c) - (1/C.AR_c)*(0.5)*(1-C.lambda_c)/(1+C.lambda_c);

%% Fuselage Constants
C.Kinl = 1.25; %[] for buried engine designs
C.lf = 48.0833; %[ft] fuselage max length
C.hf = 5.02; %[ft] fuselage max height
C.W_OEW = 13906;
C.W_TO = 17800; %[lbs] take off gross weight

%% Nacelles Constants
C.T_TO = 16000; %[lbs] Take off thrust

%% Engine Constants
C.W_eng = 2282;
C.W_ppl = 3980;

%% Systems Constants
C.Csys = 2000; %[lbs] constant systems weight
C.fsys = 0.2; %[lbs/lbs] 
C.W_dry = 13906; % [lbs] dry weight (equivalent to empty weight)

%% Crew Provisions
C.N_crew = 1; %[persons] just the pilot

%% Take off gross Constants
C.W_crew = 180; %[lbs] pilot and gear
C.W_ppl = 3980; %[lbs] fuel weight

end