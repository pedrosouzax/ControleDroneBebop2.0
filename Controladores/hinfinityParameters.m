function [t11,t12,ru] = hinfinityParameters(Hinfparam)

gamma = 100000000;
cont = 0;
while (1)
cont = cont+0.1;
t11_a =(1+ (100 - cont))*sqrt(Hinfparam(1)*Hinfparam(4));

gamma_a = sqrt(t11_a^4/(-Hinfparam(1) + t11_a^2*Hinfparam(4)^(-1)));

t12_a = sqrt(-Hinfparam(2)/(gamma_a^(-2)*t11_a^2 - Hinfparam(4)^(-1)));

k_a = -Hinfparam(3) + t11_a*t12_a*Hinfparam(4)^-1 - gamma_a^(-2)*t11_a^3*t12_a;

if k_a > 0 && gamma_a < gamma
    gamma = gamma_a
    t11 = t11_a
    t12 = t12_a
    k   =   k_a
end

if gamma_a > gamma
    break;
end
end



