function [v,u,F0, deltaF0, THETA_GP, THETA, us, Zt,yt] = controller_v2 (control, Kp,Kd,k, q,dqd,d2qd, R_inv,T_zero, T11, T11_inv,q_til,dq_til,dq,qd, THETA, qe, qed,X_ga, F0_gp, sigma_quad,l, THETA_GP, Zt,yt, etol, p_max, omega_quad,i,dt)

%% Modelo nominal        
F0 = droneNominalModel(k, d2qd, dq_til, dqd, q_til, q, X_ga(3), X_ga(1), X_ga(2));


if(control == 1)% FeedBack Linearization
    
    %Lei de controle PD
    u = -Kp*q_til - Kd*dq_til;
    %u(1:4)= tanh(u(1:4));

    %Feedback Linearization
    v = Feedback_linearization(k,q,d2qd,dqd,u);
    
    v(1:4) = tanh(v(1:4));
    
    % Retorna com zero pois não existe estas variáveis para esse
    % controlador
    us = zeros(4,1);
    deltaF0 = zeros(4,1);
    Zt=0;
    yt = 0;
    THETA = 0;
    THETA_GP = zeros(4,1);
    

elseif(control == 2)% H infinito
    
    %Lei de Controle H infinito
    u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til; q_til];
    
    %Vetor de controle H infinito
    v = F0 + T11_inv*u;
    
    v(1:4)= tanh(v(1:4));
    
    % Retorna com zero pois não existe estas variáveis para esse
    % controlador
    us = zeros(4,1);
    deltaF0 = zeros(4,1); 
    Zt=0;
    yt = 0;
    THETA = 0;
    THETA_GP = zeros(4,1);
    
elseif(control == 3)% H infinito com Rede Neural
    
    %Lei de Controle H infinito
    u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til; q_til];
    
    %Rede Neural
    [E, THETA, us] = Rede_Neural(THETA, dt, i, q_til, dq_til, T_zero, T11,  q, dq, qd, dqd, d2qd);
    
    %Vetor de Controle H infinito com Rede Neural 
    v = F0 +  E*THETA + T11_inv * u + us;
    
    v(1:4)= tanh(v(1:4));
    
    deltaF0 = E*THETA; %Delta F0 da Rede Neural
    
    % Retorna com zero pois não existe estas variáveis para esse
    % controlador
    Zt=0;
    yt = 0;
    THETA_GP = zeros(4,1);
    
elseif(control == 4)% H infinito com GP THETA Adaptativo
    
    %Lei de Controle H infinito
    u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til; q_til];
    u(1:4)= tanh(u(1:4));

    %GP com THETA adaptativo
    [mu_GPtheta, THETA_GP, us] = GP_theta_adaptativo(sigma_quad,l, THETA_GP, dt, i, q_til, dq_til, T_zero,T11, qe', qed);
    
    %Vetor de Controle com H infinito GP adaptativo
    v = F0 + mu_GPtheta + T11_inv*u + us;
    
    deltaF0 = mu_GPtheta; %Delta F0 do GP adaptativo THETA    
    
    
    % Retorna com zero pois não existe estas variáveis para esse
    % controlador
    Zt=0;
    yt = 0;
    THETA = 0;
    
elseif(control == 5)% H infinito com GP Offline
    
    %Lei de Controle H infinito
    u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til; q_til];
    u(1:4)= tanh(u(1:4));

    %Vetor de Controle com H infinito GP offline
    v = F0 + F0_gp(i,:)' + T11_inv*u; % Estimando a incerteza
    
    deltaF0 = F0_gp(i,:);
    
    % Retorna com zero pois não existe estas variáveis para esse
    % controlador
    us = zeros(4,1);
    Zt=0;
    yt = 0;
    THETA = 0;
    THETA_GP = zeros(4,1);
    
elseif(control == 6)% H infinito com GP Adaptativo Artigo
    
    %Lei de Controle H infinito
    u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til; q_til];
    u(1:4)= tanh(u(1:4));
    
    zt_1 = qed;
    
    %GP adaptativo - artigo 
    [mu,cov,Zt,zt_1] = GP_adaptativo_art(sigma_quad,l,zt_1,Zt,yt,i,omega_quad,etol, p_max);
    
    %Vetor de Controle com H infinito GP adaptativo do artigo
    v = F0 + mu(:,i) + T11_inv*u;
    
    %Zt(:,length(Zt)+1) = zt_1; 
    yt(:,i+1) = mu(:,i);%Salva o valor da média calculado
    
    Zt(:,length(Zt(1,:))+1) = qe; %Vetor que salva os valores conhecidos a cada passo
    stdv = sqrt(diag(cov)); %Desvio padrão
    
    deltaF0(:,i) = mu(:,i); %Delta F0 do GP do artigo
    
    % Retorna com zero pois não existe estas variáveis para esse
    % controlador
    us = zeros(4,1);
    THETA = 0;
    THETA_GP = zeros(4,1);
    
else
   disp('Controle inexistente');
   return
end

% if (control ~= 3 && control ~= 5)
%     us(4,i) = 0;
%     deltaF0(4,i) = 0;
%     
% elseif(control ~= 4)
%     deltaF0(4,i) = 0;   
%     
% elseif (control ~= 6)
%     Zt=0;
%     yt = 0;
%     deltaF0(4,i) = 0;
% end