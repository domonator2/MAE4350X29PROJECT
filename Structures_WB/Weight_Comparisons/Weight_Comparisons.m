% This code will be used to compare the methods used as well as determine the best output for Istr and Csys

clear;clc;

% Call Constants for Raymer and Roskam;
CRo = Constants_Roskam();
CRa = Constants_Raymer();

% Use Methods to Estimate Weight Components
Ra = Raymers_Method();
Ro = Roskams_Method();

% --- 2. Create Structural Component Weight Table ---
fprintf('====================================================\n');
fprintf('  TABLE 1: STRUCTURAL COMPONENT WEIGHTS \n');
fprintf('====================================================\n');

% Define row names for structural components
StrRowNames = {
    'Wing Weight (W_wing)';
    'Vertical Tail (W_v)';
    'Canard Weight (W_c)';
    'Empennage (W_emp)';
    'Fuselage Weight (W_fus)';
    'Landing Gear Weight (W_g)';
    'Nacelle Weight (W_n)';
};

% Create Roskam data vector for structure
RoskamStrValues = [
    Ro.W_wing*0.453592;
    Ro.W_v*0.453592;
    Ro.W_c*0.453592;
    Ro.W_emp*0.453592;
    Ro.W_fus*0.453592;
    Ro.W_g*0.453592;
    Ro.W_n*0.453592;
];

% Create Roskam data vector for structure
RoskamSrtPercent = [
    100*(Ro.W_wing/Ro.W_str);
    100*(Ro.W_v/Ro.W_str);
    100*(Ro.W_c/Ro.W_str);
    100*(Ro.W_emp/Ro.W_str);
    100*(Ro.W_fus/Ro.W_str);
    100*(Ro.W_g/Ro.W_str);
    100*(Ro.W_n/Ro.W_str);
];

% Create Raymer data vector for structure
RaymerStrValues = [
    Ra.W_wing*0.453592;
    Ra.W_v*0.453592;
    Ra.W_c*0.453592;
    Ra.W_emp*0.453592;
    Ra.W_fus*0.453592;
    Ra.W_g*0.453592;
    Ra.W_n*0.453592;
];

% Create Raymer data vector for structure
RaymerStrPercent = [
    100*(Ra.W_wing/Ra.W_str);
    100*(Ra.W_v/Ra.W_str);
    100*(Ra.W_c/Ra.W_str);
    100*(Ra.W_emp/Ra.W_str);
    100*(Ra.W_fus/Ra.W_str);
    100*(Ra.W_g/Ra.W_str);
    100*(Ra.W_n/Ra.W_str);
];

% Create and display the Structural Weight Table
StrTable = table(RoskamStrValues, RaymerStrValues, RoskamSrtPercent, RaymerStrPercent, ...
    'VariableNames', {'Roskam [kg]', 'Raymer [kg]', 'Roskam W/W_str [%]', 'Raymer W/W_str [%]'}, ...
    'RowNames', StrRowNames);
disp(StrTable);

% --- 3. Create Operating Empty Weight (OEW) Table ---
fprintf('\n====================================================\n');
fprintf('  TABLE 2: OPERATING EMPTY WEIGHT COMPONENTS \n');
fprintf('====================================================\n');

% Define row names for OEW components
OewRowNames = {
    'Total Structural Weight (W_str)';
    'Total System Weight (W_sys)';
    'Fixed Equipment/CP (W_cp)';
    'Engine Weight (W_eng)';
    'OPERATING EMPTY WEIGHT (OEW)';
};

% Create Roskam data vector for OEW
RoskamOewValues = [
    Ro.W_str*0.453592;
    Ro.W_sys*0.453592;
    Ro.W_cp*0.453592;
    CRo.W_eng*0.453592;
    Ro.OEW*0.453592
];

% Create Roskam data vector for OEW
RoskamOewPercent = [
    Ro.W_str/Ro.OEW * 100;
    Ro.W_sys/Ro.OEW * 100;
    Ro.W_cp/Ro.OEW * 100;
    CRo.W_eng/Ro.OEW * 100;
    Ro.OEW/Ro.OEW * 100
];

% Create Raymer data vector for OEW
RaymerOewValues = [
    Ra.W_str*0.453592;
    Ra.W_sys*0.453592;
    Ra.W_cp*0.453592;
    CRa.W_en*0.453592;
    Ra.OEW*0.453592
];

% Create Raymer data vector for OEW
RaymerOewPercent = [
    Ra.W_str/Ra.OEW * 100;
    Ra.W_sys/Ra.OEW * 100;
    Ra.W_cp/Ra.OEW * 100;
    CRa.W_en/Ra.OEW * 100;
    Ra.OEW/Ra.OEW * 100
];


% Create and display the OEW Table
OewTable = table(RoskamOewValues, RaymerOewValues, RoskamOewPercent, RaymerOewPercent, ...
    'VariableNames', {'Roskam [kg]', 'Raymer [kg]', 'Roskam W/OEW [%]', 'Raymer W/OWE [%]'}, ...
    'RowNames', OewRowNames);
disp(OewTable);


% --- 4. Create Key Metrics and Final Weights Table ---
fprintf('\n====================================================\n');
fprintf('  TABLE 3: KEY METRICS & FINAL WEIGHTS\n');
fprintf('====================================================\n');

% Define row names for Key Metrics
MetricsRowNames = {
    'Structural Index (Istr) [kg/m^2]';
    'Constant System Weight (Csys) [kg]';
    'Take Off Gross Weight [kg]';
};

% Create Roskam data vector for metrics
RoskamMetricsValues = [
    Ro.Istr*4.88;
    Ro.Csys*0.453592;
    Ro.TOGW*0.453592;

];

% Create Raymer data vector for metrics
RaymerMetricsValues = [
    Ra.Istr*4.88;
    Ra.Csys*0.453592;
    Ra.TOGW*0.453592;
];

% Create and display the Metrics Table
MetricsTable = table(RoskamMetricsValues, RaymerMetricsValues, ...
    'VariableNames', {'Roskam', 'Raymer'}, ...
    'RowNames', MetricsRowNames);
disp(MetricsTable);

% --- 5. Comparison of Key Differences ---
fprintf('\n====================================================\n');
fprintf('  KEY COMPARISON METRICS (Differences)\n');
fprintf('====================================================\n');

% Calculate differences
Csys_diff = abs(Ro.Csys*0.453592 - Ra.Csys*0.453592);
Istr_diff = abs(4.88*Ro.Istr - 4.88*Ra.Istr);

fprintf('  1. Structural Index (Istr) Difference: %.4f kg/m^2\n', Istr_diff);
fprintf('  2. Constant System Weight (Csys) Difference: %.2f kg\n', Csys_diff);
fprintf('  3. Take Off Gross Weight (TOGW) Difference: %.2f kg\n', abs(Ro.TOGW*0.453592 - Ra.TOGW*0.453592));
fprintf('====================================================\n');