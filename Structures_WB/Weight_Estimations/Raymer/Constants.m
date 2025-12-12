function C = Constants()

    % Wing Variables
    C.Kdw = 1.0; %[] for non-delta wings
    C.Kvs = 1.0; %[] for non-variable sweep wings
    C.W_TO = 13906; %[lbs] take-off gross weight
    C.Nz = 8*1.5; %[g] ultimate load factor
    C.Sw = 185; %[ft^2] wing area
    C.A_wing = 4; %[] aspect ratio
    C.tc_root = 0.07; %[] root thickness over chord ratio
    C.lambda_wing = 0.4; % [] taper ratio
    C.sweep_wing_LE = deg2rad(-29.27); %[rad] leading edge sweep angle
    C.sweep_wing = tan(C.sweep_wing_LE) - (4/(C.A_wing)) * ((1-0.25)*((1-C.lambda_wing)/(1+C.lambda_wing)) );
    C.Scsw = 31.57; %[ft^2] control system area on the wing

    % Emennage Variables
    C.Krht = 1.0; % [] for non-rolling horizontal tails
    C.Ht_Hv = 0; % [] for conventional vertical tails
    C.Svt = 33.75; %[ft^2] vertical tail area 
    C.M = 1.6; % [] Mach number
    C.Lt = 6.78; %[ft] tail length
    C.Sr = 66.74; %[ft^2] rudder area
    C.A_vt = 2.64; %[] vertical tail aspect ratio
    C.lambda_vt = 0.306; %[] vertical tail taper ratio
    C.sweep_vt_LE = deg2rad(47); %[rad] vertical tail LE sweep angle
    C.sweep_vt = tan(C.sweep_vt_LE) - (4/(C.A_vt))*((1-0.25)*((1-C.lambda_vt)/(1+C.lambda_vt)) );
    C.Fw = 6.34; %[ft] fuselage length at canard locations
    C.Bc = 3.93; %[ft] canard span
    C.Sc = 37.0; %[ft^2] canard area

    % Fuselage Variables
    C.Kdwf = 1.0; %[] for non-delta wings
    C.L = 48.083; %[ft] fuselage structural length
    C.D = 5.02; %[ft] fuselage structural depth
    C.W = 3.77; %[ft] fuselage structural width 

    % Landing Gear Variables
    C.Kcb = 1.0; %[] for non-cross-beam gear
    C.Ktpg = 1.0; %[] for non-tripod landing gear
    C.Wl = 13906; %[lbs] landing design gross weight ~ OEW
    C.Nl = 1.2*1.5; %[] landing gear ultimate load factor
    C.Lm = 5.03; %[ft] length of main landing gear
    C.Ln = 3.67; %[ft] length of nose landing gear
    C.Nnw = 1; %[] number of nose wheels

    % Nacelles Variables
    C.Nen = 1; %[] number of engines
    C.T = 16000; %[lbs] otal engine thrust
    C.W_en = 2282; %[lbs] engine weight
    C.Sfw = 276; %[ft^2] firewall surface area
    C.Kvg = 1; %[] for non-variable geometry
    C.Ld = 17.3; %[ft] duct length
    C.Kd = 2.75; %[] for rectangular air ducts
    C.Ls = 5.5; %[ft] single duct length
    C.De = 35/12; %[ft] engine diameter
    C.Ltp = 1.95; %[ft] lenth of the the tailpipe

end