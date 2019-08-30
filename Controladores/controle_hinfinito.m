clear all; close all; clc
q11 = 1;
q22 = 0.5;
q12 = 0;
ru = 10.8;

gamma = 100000000;
cont = 0;
while (1)
cont = cont+0.1;
t11_a =(1+ (100 - cont))*sqrt(q11*ru);

gamma_a = sqrt(t11_a^4/(-q11 + t11_a^2*ru^(-1)));

t12_a = sqrt(-q22/(gamma_a^(-2)*t11_a^2 - ru^(-1)));

k_a = -q12 + t11_a*t12_a*ru^-1 - gamma_a^(-2)*t11_a^3*t12_a;

if k_a > 0 && gamma_a < gamma
    gamma = gamma_a;
    t11 = t11_a;
    t12 = t12_a;
    k   =   k_a;
end

if gamma_a > gamma
    break;
end
end

[{'gamma'} gamma;
 {'t11'}   t11;
 {'t12'}   t12;
 {'k'}     k;
 {'q11'}   q11;
 {'q22'}   q22;
 {'q12'}   q12;
 {'ru'}    ru]

fprintf('data.hinfty_controller.gamma = %f; \n',gamma)
fprintf('data.hinfty_controller.t11 = %f;\n', t11)
fprintf('data.hinfty_controller.t12 = %f;\n', t12)
fprintf('data.hinfty_controller.k = %f;\n', k)
fprintf('data.hinfty_controller.q11 = %f;\n', q11)
fprintf('data.hinfty_controller.q22 = %f;\n', q22)
fprintf('data.hinfty_controller.q12 = %f;\n', q12)
fprintf('data.hinfty_controller.ru = %f;\n', ru)

disp(' ')
disp(' ')
disp(' ')

fprintf('ru = %f; \n',ru)
fprintf('gama = %f; \n',gamma)
fprintf('t_11 = %f; \n',t11)
fprintf('t_12 = %f; \n',t12)


