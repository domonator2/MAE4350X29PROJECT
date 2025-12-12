function [W_wing, W_vt, W_c, W_emp, W_fus, W_n, W_g, W_str] = Weight_Calculator(C)

    % Wing Calculation
    W_wing = 0.0103*C.Kdw*C.Kvs* (C.W_TO*C.Nz)^0.5 * C.Sw^0.622 * C.A_wing^0.785 * (C.tc_root)^-0.4 * (1+C.lambda_wing)^0.05 * (cos(C.sweep_wing))^-1 * C.Scsw^0.04;

    % Empennage Calculation
    W_vt = 0.452*C.Krht*(1 + C.Ht_Hv)^0.5 *(C.W_TO*C.Nz)^0.488 * C.Svt^0.718 * C.M^0.341 * C.Lt^-1 * (1 + C.Sr/C.Svt)^0.348 * C.A_vt^0.223 * (1 + C.lambda_vt)^0.25 *(cos(C.sweep_vt))^-0.323;
    W_c = 3.316*(1 + C.Fw/C.Bc)^-2.0 * (C.W_TO*C.Nz/1000)^0.260 * C.Sc^0.806;

    W_emp = W_vt + W_c;

    % Fuselage Calculation
    W_fus = 0.499 * C.Kdwf * C.W_TO^0.35 * C.Nz^0.25 * C.L^0.5 * C.D^0.849 * C.W^0.685;

    % Landing Gear Calculation
    
    W_main = C.Kcb*C.Ktpg*(C.Wl*C.Nl)^0.25 * (C.Lm*12)^0.973;
    W_nose = (C.Wl*C.Nl)^0.290 * (C.Ln*12)^0.5 * C.Nnw^0.525;

    W_g = W_main + W_nose;

    % Nacelles Calculation
    W_engine_mounts = 0.013 * C.Nen^0.795 * C.T^0.579 * C.Nz;
    W_engine_section = 0.01 * C.W_en^0.717 * C.Nen * C.Nz;
    W_firewall = 1.13*C.Sfw;
    W_airduct = 13.29*C.Kvg * C.Ld^0.643 * C.Kd^0.182 * C.Nen^1.498 * (C.Ls/C.Ld)^-0.373 * C.De;
    W_tailpipe = 3.5 * C.De * C.Ltp * C.Nen;

    W_n = W_engine_mounts + W_engine_section + W_firewall + W_airduct + W_tailpipe;

    % Structures Calculation
    W_str = W_wing + W_emp + W_fus + W_g + W_n;

    % Crew Provisions
    

end