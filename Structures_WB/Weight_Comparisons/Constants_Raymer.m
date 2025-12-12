function CRa = Constants_Raymer()
    %% Wing Constants
    CRa.Kdw = 1.0; %[] for non-delta wings
    CRa.Kvs = 1.0; %[] for non-variable sweep wings
    CRa.W_TO = 17800; %[lbs] take-off gross weight
    CRa.Nz = 8*1.5; %[g] ultimate load factor
    CRa.Sw = 185; %[ft^2] wing area
    CRa.A_wing = 4; %[] aspect ratio
    CRa.tc_root = 0.07; %[] root thickness over chord ratio
    CRa.lambda_wing = 0.4; % [] taper ratio
    CRa.sweep_wing_LE = deg2rad(-29.27); %[rad] leading edge sweep angle
    CRa.sweep_wing = tan(CRa.sweep_wing_LE) - (1/(CRa.A_wing)) * (((1-CRa.lambda_wing)/(1+CRa.lambda_wing)) );
    CRa.Scsw = 63.15; %[ft^2] control system area on the wing (used webdigitizer)

    %% Emennage Constants
    CRa.Krht = 1.0; % [] for non-rolling horizontal tails
    CRa.Ht_Hv = 0; % [] for conventional vertical tails
    CRa.Svt = 33.75; %[ft^2] vertical tail area 
    CRa.M = 1.6; % [] Mach number
    CRa.Lt = 11.66; %[ft] wing 25% MAC to tail 25% MAC (used webdigitizer)
    CRa.Sr = 15.02; %[ft^2] rudder area (used webdigitizer)
    CRa.A_vt = 2.64; %[] vertical tail aspect ratio
    CRa.lambda_vt = 0.306; %[] vertical tail taper ratio
    CRa.sweep_vt_LE = deg2rad(47); %[rad] vertical tail LE sweep angle
    CRa.sweep_vt = tan(CRa.sweep_vt_LE) - (1/(CRa.A_vt))*((1-CRa.lambda_vt)/(1+CRa.lambda_vt));
    
    CRa.Fw = 7.06; %[ft] fuselage length at canard locations (used webdigitizer)
    CRa.Bc = 3.84; %[ft] canard span (used webdigitier)
    CRa.Sc = 37.0; %[ft^2] canard area

    %% Fuselage Constants
    CRa.Kdwf = 1.0; %[] for non-delta wings
    CRa.L = 48.083; %[ft] fuselage structural length
    CRa.D = 5.02; %[ft] fuselage structural depth (used webdigitizer)
    CRa.W = 4.60; %[ft] fuselage structural width (used webdigitizer - excluded the inlets)

    %% Landing Gear Constants
    CRa.Kcb = 1.0; %[] for non-cross-beam gear
    CRa.Ktpg = 1.0; %[] for non-tripod landing gear
    CRa.Wl = 17800; %[lbs] landing design gross weight
    CRa.Nl = 1.2*1.5; %[] landing gear ultimate load factor
    CRa.Lm = 3.41; %[ft] length of main landing gear (used webdigitizer)
    CRa.Ln = 3.22; %[ft] length of nose landing gear (used webdigitizer)
    CRa.Nnw = 1; %[] number of nose wheels

    %% Nacelles Constants
    CRa.Nen = 1; %[] number of engines
    CRa.T = 16000; %[lbs] total engine thrust
    CRa.W_en = 2282; %[lbs] engine weight
    CRa.Sfw = 250.20; %[ft^2] firewall surface area (used webdigitizer and assumption that it was a cylinder)
    CRa.Kvg = 1; %[] for non-variable geometry 
    CRa.Ld = 17.3; %[ft] duct length (used webdigitizer)
    CRa.Kd = 2.75; %[] for rectangular air ducts
    CRa.Ls = 5.5; %[ft] single duct length (used webdigitizer)
    CRa.De = 35/12; %[ft] engine diameter
    CRa.Ltp = 1.95; %[ft] lenth of the the tailpipe (used webdigitizer)

    %% Crew Provisions Constants
    CRa.Np = 1; % 
    CRa.Wc = 500; % [lbs] max cargo weight (assumed number)
    CRa.Sf = 1200 - 185 - 33.75 - 37.0;
    CRa.W_crew = 180; %[lbs] crew weight

    %% Systems Constants
    CRa.fsys = 0.2; %[lbs/lbs] system weight percentage of empty weight
    CRa.W_dry = 13906; % [lbs] dry weight = OEW

    %% Engine Constants
    CRa.W_en = 2282; %[lbs] engine weight
    CRa.W_ppl = 3980; %[lbs] fuel weight
    
    %% Structural Index
    CRa.S_wet = 1200; %[ft^2] total wetted area

end