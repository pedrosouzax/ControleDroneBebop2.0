function [ru, gama, t11, t12] = parameters_Hinf(ru, q11, q12, q22)
% gama = input('Escolha um valor para gama > 0 de valor elevado: '); %Escolha um valor para gama
% gama_anterior = gama+1;
% ru = input('Escolha um valor para ru: ');
% q11 = input('Escolha um valor para q11: ');
% q12 = input('Escolha um valor para q12: ');
% q22 = input('Escolha um valor para q22: ');
% t11 = input('Escolha um valor para t11 (t11>>sqrt(q11*ru)): ');%Escolha t11 de modo que t11>>sqrt(q11*ru)

gama = 10e8; %input('Escolha um valor para gama > 0 de valor elevado: '); %Escolha um valor para gama
% ru = 0.7; %input('Escolha um valor para ru: ');
% q11 = 0.3; %input('Escolha um valor para q11: ');
% q12 = 0.0; %input('Escolha um valor para q12: ');
% q22 = 1.2; %input('Escolha um valor para q22: ');
% %textt11 = sprintf('Escolha um valor para t11 (t11>>sqrt(q11*ru) ou (t11 >> %f)): ', sqrt(q11*ru));
t11_a = 100;%100*sqrt(q11*ru); %input(textt11);%Escolha t11 de modo que t11>>sqrt(q11*ru)

while(true)
    gama_a = (t11_a^2)/ (sqrt(-q11 + (t11_a^2)/ru));
    t12_a = sqrt((-q22)/((-1/ru)+((t11_a^2)/gama_a^2)));
    k_a = -q12 + ((t11_a*t12_a)/ru) - (((t11_a^3)*t12_a)/gama_a^2);
    
    if k_a>0 && gama_a<gama
        gama = gama_a;
        t11 = t11_a;
        t12 = t12_a;
        k = k_a;
    end
    
    if gama < gama_a
        break;
    
    else 
        t11_a = t11_a - 0.001;
    end
    
    if t11_a <= sqrt(q11*ru)
        break;
    end
end
    