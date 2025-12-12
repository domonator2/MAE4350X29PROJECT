function [Minf, V, q, ReL] = flight_condition_units(h, option, value, units)

% Compute Mach number and velocity    
isMach = strcmpi(option, 'M');
isVelocity = strcmpi(option, 'V');
    
    if isMach
        Minf = value; % Mach number remains unchanged
        [~, ~, rho, a, mu] = atmospheric_model_units(h, units);
        V = Minf * a; % Compute velocity (using speed of sound in m/s)
        q = (1/2) * rho * (V^2); % dynamic pressure
        ReL = (rho * V) / mu; % Reynolds number * L
    
    elseif isVelocity
        V = value; % Velocity is already converted to m/s
        [~, ~, rho, a, mu] = atmospheric_model_units(h, units);
        q = (1/2) * rho * (V^2);
        ReL = (rho * V) / mu;
        Minf = V / a; % Compute Mach number (using speed of sound in m/s)
    else
        error('Invalid option. Enter M for Mach number or V for velocity.');
    end
end