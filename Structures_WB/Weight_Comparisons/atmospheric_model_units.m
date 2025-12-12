function [T, P, rho, a, mu] = atmospheric_model_units(h, units) 
    % Constants
    g0 = 9.80665; % gravity acceleration (m/s^2)
    R = 287.0528; % Gas Constant (J/kg·K)

    % Convert altitude to meters if input is in English units
    if strcmpi(units, 'EN')
        h = h * 0.3048; % Convert feet to meters
    end

    % Compute atmospheric properties in SI units
    [T, P, rho, a, mu] = calculate_atmosphere(h, g0, R);

    % Convert outputs to English units if requested
    if strcmpi(units, 'EN')
        T = T * 1.8; % Convert temperature to Rankine
        P = P * 0.020885; % Convert pressure to lbf/ft^2
        rho = rho * 0.00194032; % Convert density to slug/ft^3
        a = a * (3.28084); % Convert speed of sound to ft/s
        mu = mu * (0.020885); % Convert dynamic viscosity to lbf·s/ft^2
    end

end

function [T, P, rho, a, mu] = calculate_atmosphere(h, g0, R)
    % Compute atmospheric properties based on altitude

    % Convert to geopotential altitude
    z = (h * 6371e3) / (h + 6371e3); % Geopotential altitude approximation

    % Atmospheric layer conditions
    if (0 <= z) && (z <= 11000)
        T0 = 288.15; % Base temperature (K)
        p0 = 101325; % Base pressure (Pa)
        Ti = -0.0065; % Temperature gradient (K/m)
        T = T0 + Ti * (z - 0);
        P = p0 * (T / T0)^(-g0 / (R * Ti));

    elseif (11000 < z) && (z <= 20000)
        T = 216.65; % Constant temperature (K)
        p0 = 22632.049; % Base pressure (Pa)
        P = p0 * exp(-g0 * (z - 11000) / (R * T));

    elseif (20000 < z) && (z <= 32000)
        T0 = 216.65;
        p0 = 5474.882;
        Ti = 0.001;
        T = T0 + Ti * (z - 20000);
        P = p0 * (T / T0)^(-g0 / (R * Ti));

    elseif (32000 < z) && (z <= 47000)
        T0 = 228.65;
        p0 = 868.017;
        Ti = 0.0028;
        T = T0 + Ti * (z - 32000);
        P = p0 * (T / T0)^(-g0 / (R * Ti));

    elseif (47000 < z) && (z <= 52000)
        T = 270.65;
        p0 = 110.906;
        P = p0 * exp(-g0 * (z - 47000) / (R * T));

    elseif (52000 < z) && (z <= 61000)
        T0 = 270.65;
        p0 = 59.001;
        Ti = -0.002;
        T = T0 + Ti * (z - 52000);
        P = p0 * (T / T0)^(-g0 / (R * Ti));

    elseif (61000 < z) && (z <= 79000)
        T0 = 252.65;
        p0 = 18.210;
        Ti = -0.004;
        T = T0 + Ti * (z - 61000);
        P = p0 * (T / T0)^(-g0 / (R * Ti));

    elseif (79000 < z) && (z <= 90000)
        T = 180.65;
        p0 = 1.038;
        P = p0 * exp(-g0 * (z - 79000) / (R * T));

    else
        error('Altitude is out of range for this atmospheric model.');
    end

    % Compute density, speed of sound, and dynamic viscosity
    rho = P / (R * T);
    a = sqrt(1.4 * R * T);
    mu = (1.458e-6 * T^(3/2)) / (T + 110.4);
end
