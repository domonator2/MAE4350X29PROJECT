function [W_wing, W_vt, W_c, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Weight_Calculator_Raymer(CRa)

    % Wing Calculation
    W_wing = 0.0103*CRa.Kdw*CRa.Kvs* (CRa.W_TO*CRa.Nz)^0.5 * CRa.Sw^0.622 * CRa.A_wing^0.785 * (CRa.tc_root)^-0.4 * (1+CRa.lambda_wing)^0.05 * (cos(CRa.sweep_wing))^-1 * CRa.Scsw^0.04;

    % Empennage Calculation
    W_vt = 0.452*CRa.Krht*(1 + CRa.Ht_Hv)^0.5 *(CRa.W_TO*CRa.Nz)^0.488 * CRa.Svt^0.718 * CRa.M^0.341 * CRa.Lt^-1 * (1 + CRa.Sr/CRa.Svt)^0.348 * CRa.A_vt^0.223 * (1 + CRa.lambda_vt)^0.25 *(cos(CRa.sweep_vt))^-0.323;
    W_c = 3.316*(1 + CRa.Fw/CRa.Bc)^-2.0 * (CRa.W_TO*CRa.Nz/1000)^0.260 * CRa.Sc^0.806;

    W_emp = W_vt + W_c;

    % Fuselage Calculation
    W_fus = 0.499 * CRa.Kdwf * CRa.W_TO^0.35 * CRa.Nz^0.25 * CRa.L^0.5 * CRa.D^0.849 * CRa.W^0.685;

    % Landing Gear Calculation
    
    W_main = CRa.Kcb*CRa.Ktpg*(CRa.Wl*CRa.Nl)^0.25 * (CRa.Lm*12)^0.973;
    W_nose = (CRa.Wl*CRa.Nl)^0.290 * (CRa.Ln*12)^0.5 * CRa.Nnw^0.525;

    W_g = W_main + W_nose;

    % Nacelles Calculation
    W_engine_mounts = 0.013 * CRa.Nen^0.795 * CRa.T^0.579 * CRa.Nz;
    W_engine_section = 0.01 * CRa.W_en^0.717 * CRa.Nen * CRa.Nz;
    W_firewall = 1.13*CRa.Sfw;
    W_airduct = 13.29*CRa.Kvg * CRa.Ld^0.643 * CRa.Kd^0.182 * CRa.Nen^1.498 * (CRa.Ls/CRa.Ld)^-0.373 * CRa.De;
    W_tailpipe = 3.5 * CRa.De * CRa.Ltp * CRa.Nen;

    W_n = W_engine_mounts + W_engine_section + W_firewall + W_airduct + W_tailpipe;

    % Structures Calculation
    W_str = W_wing + W_emp + W_fus + W_g + W_n;

    % Crew Provisions
    W_furnishing = 0.0577 * CRa.Np^0.1 * CRa.Wc^0.393 * CRa.Sf^0.75;
    W_handling_gear = 3.0E-4 * CRa.W_TO;

    W_cp = W_furnishing + W_handling_gear;

end