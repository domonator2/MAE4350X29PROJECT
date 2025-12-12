% Raymers Method Istr and Csys Calculations

clear;clc;

% Raymer Constants
C = Constants();

% Structure Weight Calculations
[W_wing, W_vt, W_c, W_emp, W_fus, W_n, W_g, W_str, W_cp] = Structures_Weight_Calculator(C);

fprintf('W_wing: %.2f\n', W_wing);
fprintf('W_emp: %.2f\n', W_emp);
fprintf('W_fus: %.2f\n', W_fus);
fprintf('W_n: %.2f\n', W_n);
fprintf('W_g: %.2f\n\n', W_g);

fprintf('W_str: %.2f', W_str);
