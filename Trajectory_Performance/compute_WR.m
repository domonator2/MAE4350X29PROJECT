%{
====================================================================================================
SCRIPT NAME: compute_WR.m
AUTHOR: Damian Do
INITIATED: 9/7/2025
LAST REVISION: 9/13/2025
====================================================================================================
SCRIPT DESCRIPTION:
This function will calculate the weight ratio for a given mission segment.
====================================================================================================
USER-DEFINED FUNCTIONS:
None.
====================================================================================================
%}
function WR = compute_WR(L_D_max, SFC, Mission, M_cruise)

% [m] Extracts the cruise altitude from the Mission structure.
h_cruise = Mission.h_cruise;

% [s] Extracts the cruise time from the Mission structure.
t_cruise = Mission.t_cruise;

% [m/s, kg/m^3, etc.] Computes atmospheric properties at cruise altitude.
[a_cruise, ~, ~, ~, ~, ~, ~, ~, ~, ~] = DD_atm_model(h_cruise, "SI");

% [m/s] Calculates the cruise velocity.
V_cruise = M_cruise*a_cruise;

% [m] Calculates the supersonic cruise range.
R_cruise = V_cruise*t_cruise;

% [1/s] Converts the specific fuel consumption from 1/hr to 1/s.
SFC_si = SFC / 3600;

% [-] Fuel fraction for engine start and warm-up.
W1_WTO = 0.990;

% [-] Fuel fraction for taxi.
W2_W1 = 0.990;

% [-] Fuel fraction for takeoff.
W3_W2 = 0.990;

% [-] Fuel fraction for climb/acceleration to cruise altitude and speed.
W4_W3 = -0.036236*M_cruise + 0.998538;

% [-] Fuel fraction for the cruise segment (from Breguet Range equation).
W5_W4 = exp(-R_cruise*SFC_si*9.81/(V_cruise*L_D_max));

% [-] Fuel fraction for descent.
W6_W5 = 0.990;

% [-] Fuel fraction for landing.
W7_W6 = 0.995;

% [-] Calculates the total mission fuel fraction.
M_ff = W1_WTO*W2_W1*W3_W2*W4_W3*W5_W4*W6_W5*W7_W6;

% [-] Calculates the ratio of fuel weight to TOGW.
ff = (1 - M_ff) * 1.06;

% [-] Calculates the weight ratio.
WR = 1/(1 - ff);

end

%------------------------------ Helper atm model function -----------------
function [a, rho, mu, nu, theta, delta, sigma, z, P, T] = DD_atm_model(h,units)

% [J/(kg*K)] Defines the air gas constant in SI units.
R = 287.0528;

% [-] Defines the ratio of specific heats for air.
gamma = 1.4;

% [m/s^2] Defines the gravitational constant in SI units.
g0 = 9.805545;

% [K] Defines the sea level temperature reference.
T0 = 288.15;

% [Pa] Defines the sea level pressure reference.
P0 = 101325;

% [kg/m^3] Calculates the sea level density reference.
rho0 = P0/(R*T0);

% [m] Defines the radius of the Earth.
RE = 6378137;

% [] Checks if the units are given in English units.
if strcmp(units, "EN")

   % [m] Converts altitude from ft to m.
   h=h*0.3048;

end

% [m] Calculates the geopotential altitude.
z = (RE*h)/(RE+h);

% [] Computes temperature and pressure depending on input altitude.
if 0<=z && z<=11000
    
    % [K] Defines the temperature at the bottom of the layer.
    Ti = 288.15;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = -6.5E-3;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 101325;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 0;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*((Ti + delTi*(z-zi))/Ti)^(-g0/(R*delTi));

elseif 11000 < z && z <= 20000
    
    % [K] Defines the temperature at the bottom of the layer.
  Ti = 216.65;
  
  % [K/m] Defines the temperature lapse rate.
  delTi = 0; 
  
  % [Pa] Defines the pressure at the bottom of the layer.
  Pi = 22632.049;
  
  % [m] Defines the altitude at the bottom of the layer.
  zi = 11000;

  % [K] Calculates the temperature at the given altitude.
  T = Ti + delTi*(z - zi);  
  
  % [Pa] Calculates the pressure at the given altitude.
  P = Pi * exp((-g0 * (z - zi)) / (R * Ti)); 

elseif 20000<z && z<=32000
    
    % [K] Defines the temperature at the bottom of the layer.
    Ti = 216.65;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = 1E-3;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 5474.882;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 20000;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*((Ti + delTi*(z-zi))/Ti)^(-g0/(R*delTi));

elseif 32000<z && z<=47000

    % [K] Defines the temperature at the bottom of the layer.
    Ti = 228.65;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = 2.8E-3;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 868.017;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 32000;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*((Ti + delTi*(z-zi))/Ti)^(-g0/(R*delTi));

elseif 47000<z && z<=52000

    % [K] Defines the temperature at the bottom of the layer.
    Ti = 270.65;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = 0;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 110.906;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 47000;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*exp((-g0*(z-zi))/(R*Ti));

elseif 52000<z && z<=61000

    % [K] Defines the temperature at the bottom of the layer.
    Ti = 270.65;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = -2E-3;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 59.001;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 52000;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*((Ti + delTi*(z-zi))/Ti)^(-g0/(R*delTi));

elseif 61000<z && z<=79000

    % [K] Defines the temperature at the bottom of the layer.
    Ti = 252.65;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = -4E-3;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 18.210;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 61000;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*((Ti + delTi*(z-zi))/Ti)^(-g0/(R*delTi));

elseif 79000<z && z<=90000

    % [K] Defines the temperature at the bottom of the layer.
    Ti = 180.65;
    
    % [K/m] Defines the temperature lapse rate.
    delTi = 0;
    
    % [Pa] Defines the pressure at the bottom of the layer.
    Pi = 1.038;
    
    % [m] Defines the altitude at the bottom of the layer.
    zi = 79000;

    % [K] Calculates the temperature at the given altitude.
    T = Ti + delTi*(z-zi);
    
    % [Pa] Calculates the pressure at the given altitude.
    P = Pi*exp((-g0*(z-zi))/(R*Ti));   
    
end

% [m/s] Calculates the speed of sound.
a = sqrt(gamma*R*T);

% [kg/m^3] Calculates the air density.
rho = P/(R*T);

% [Ns/m^2] Calculates the dynamic viscosity.
mu = (1.458E-6)*(T^1.5)*(1/(T+110.4));

% [m^2/s] Calculates the kinematic viscosity.
nu = mu/rho;

% [-] Calculates the temperature ratio.
theta = T/T0;

% [-] Calculates the pressure ratio.
delta = P/P0;

% [-] Calculates the density ratio.
sigma = rho/rho0;

% [] Checks if the units are given in English units.
if strcmp(units, "EN")

   % [ft/s] Converts speed of sound from m/s to ft/s.
   a=a*3.28084;
   
   % [slug/ft^3] Converts air density from kg/m^3 to slug/ft^3.
   rho=rho*0.00194032;
   
   % [lbfs/ft^2] Converts dynamic viscosity from Ns/m^2 to lbfs/ft^2.
   mu=mu*0.02088543;
   
   % [ft^2/s] Converts kinematic viscosity from m^2/s to ft^2/s.
   nu=nu*10.7639;
   
   % [ft] Converts geopotential altitude from m to ft.
   z=z*3.28084;
   
   % [lbf/ft^2] Converts pressure from Pa to lbf/ft^2.
   P=P*0.020885;
   
   % [R] Converts temperature from K to degrees Rankine.
   T=T*1.8;

end

end