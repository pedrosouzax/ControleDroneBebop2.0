% function[norm_v] = Controlador_Main(THETA_GP)
clc
clear all
close all

%Control Settings

%Escolha o controle:

%Feedback Linearization ---------------------> control = 1
%H infinito ---------------------------------> control = 2
%H infinito com Rede Neural Adaptativa ------> control = 3
%H infinito com GP Offline ------------------> control = 4
%H infinito com GP Adaptativo ---------------> control = 5
%H infinito com GP Adaptativo do Artigo -----> control = 6

control = 1;

qtd_exp = 1;

for control=1:5
    for exp=1:qtd_exp

        
format long
path(path,'.\trajectory')
path(path,'.\Plots')


%% Variaveis de alteração

%Trajactory Settings

traj_des = 1; % Escolha: 1- Trajetória Elipse; 2 - Trajetória ponto final

%%%%%%%%%%%%% Trajetória Elipse %%%%%%%%%%%%%%%%%%
vel_d = 0.3; % Velocidade desejada em m/s

A = 0.5; % Amplitude da trajetória
omg = vel_d/A; % Velocidade desenada em rad/s 
tempo_voo = 12.6;%tempo de voo em segundos

%Sample Time
dt = 0.01;
N = tempo_voo/dt;


        
        [X_ga, q, dq, d2q, v, THETA_GP, qf,dqf, d2qf,traj_opt, to, tf, ax, ay, az, afi, Kp, Kd, k, kuncertain, T11, T12, T11_inv, T_zero, R, R_inv, THETA, sigma_quad, l, F0_gp, cov_gp,  K_xx,  K_sx,  K_ss, valores_conhecidos, yt, omega_quad, qe, qed] = initialization();
      
    for i=1:N
        %% Planejamento de trajetória
        t = i*dt; 
        t_v(i)= t;
        
       [qd,dqd,d2qd,d3qd,to,tf,ax,ay,az,afi] = trajectory_desired(traj_des,q,dq,d2q, qf, dqf, d2qf, t, dt, traj_opt, vel_d, A, omg, to,tf,ax,ay,az,afi);
         

        %% Calculo do erro da posição, velocidade e aceleração
        q_til(1:3) = qd(1:3) - q(1:3); %diferença entre a posição desejada e atual
        q_til(4) = normalize_angle_f(qd(4) - q(4), -pi);
        dq_til = dqd - dq;
        d2q_til = d2qd - d2q;


        %% Modelo nominal        
        F0 = droneNominalModel(kuncertain, d2qd, dq_til, dqd, q_til, q, X_ga(3), X_ga(1), X_ga(2));

        %% Controlador

        [v,u,deltaF0, THETA_GP, ~, us, valores_conhecidos,yt] = controller (control, Kp,Kd,k, q,dqd,d2qd, R_inv,T_zero, T11,T11_inv,q_til,dq_til,dq,qd, THETA, qe,qed, F0, F0_gp, sigma_quad,l, THETA_GP, valores_conhecidos,yt,omega_quad,i,dt);

        %% Drone Model
        [d2q,dq,q,d,B] = droneDynamic(k,t,dt,q,dq,v);
        
        qe = q' + dq'+ d2q';
        qed = qd' + dqd' + d2qd';
 
        
        %% Dados salvos para os plots
        
        %[u_v, v_v, us_v, F0_v, deltaF0_v, d_v, B_v, dB_v, q_v, dq_v, d2q_v, qe_v, qd_v, dqd_v, d2qd_v, qed_v, q_til_v,dq_til_v, d2q_til_v] = save_dataPlot(u, v, us, F0, deltaF0, d, B, q, dq, d2q, qe, qd, dqd, d2qd, qed, q_til,dq_til, d2q_til, i);

        %%%%%%%% Controlador %%%%%%%        
        u_v(:,i) = u; %matriz do vetor da lei de controle
        v_v(:,i) = v; %matriz do vetor de controle
        us_v(:,i) = us; %Controle variavel da Rede Neural
        
        %%%%%%% Modelo %%%%%%%%
        F0_v(:,i) = F0; %matriz do vetor F0
        deltaF0_v(:,i) = deltaF0; %matriz do vetor de incerteza de F0
        
        %%%%%% Dinamica %%%%%2%%
        d_v(:,i) = d;
        B_v(:,:,i) = B;
        dB_v(:,i) = B*d; %Aceleraï¿½ï¿½o do distï¿½rbio
        
        %%%%% Dados de posição, velocidade, aceleração e qe = q+dq+d2q  %%%%%%%%%%%%
        q_v(:,i) = q;
        dq_v(:,i) = dq;
        d2q_v(:,i) = d2q;        
        qe_v(:,i) =  qe;        

        
        %%%%%% Posição desejada, Velocidade desejada, Aceleração desejada e qed = qd+dqd+d2qd  %%%%%%%
        qd_v(:,i) = qd;
        dqd_v(:,i) = dqd;
        d2qd_v(:,i) = d2qd;
        qed_v(:,i) =  qed;
        
        %%%%% Erro de posição e Velocidade %%%%%%%%
        
        q_til_v(:,i) = q_til';
        dq_til_v(:,i) = dq_til';
        d2q_til_v(:,i) = d2q_til';

         %drawnow;
         %pause(0.1);            
    

    normEPM_til(exp,control,i) =  norm(q_til_v(1:3,i));
    normEPsiM_til(exp,control,i) = norm(q_til_v(4,i));
    normU(exp,control,i) = norm(v_v(i));    
   
    t
    
    end


   clear_data(exp, control, qtd_exp, u_v, v_v, us_v, F0_v, deltaF0_v, d_v, B_v, dB_v, q_v, dq_v, d2q_v, qe_v, qd_v, dqd_v, d2qd_v, qed_v, q_til_v,dq_til_v, d2q_til_v, normEPM_til, normEPsiM_til, normU);


    % Dados para calcular o indice de desempenho u médio
    u_barra(exp, control, :,:) = abs(v_v);
 

    end
    
    if(control == 2)
        %save_data(qe_v, qed_v, F0_v, v_v, q_v, qd_v);
    end
    
        
    %%%%%%%%%%Plota o grafico 3D para cada controle %%%%%%%%%%%%%%%%%

    plot3dim(qd_v,q_v,control)
    
    %% Calculo dos Indices de desmpenho para os gráficos
    somaEPM = sum(normEPM_til,1)/qtd_exp;
    somaEPsiM = sum(normEPsiM_til,1)/qtd_exp;
    somau_norm_media = sum(normU,1)/qtd_exp;
    
    EPM(control,:) = somaEPM(1,control,:);
    EPsiM(control,:) = somaEPsiM(1,control,:);
    u_norm_media(control,:) = somau_norm_media(1,control,:);
    
    
 %%%%%%%%%%Plota os graficos de desempenho %%%%%%%%%%%%%%%%%

    %plotDesempenhoL22(control,EPM,EPsiM,u_norm_media, t_v);
    
    %% Calculo dos Indice de desempenho 
    L22(control) = norm(q_til_v(1,:) + q_til_v(2,:) + q_til_v(3,:))
    L2psi(control) = norm(q_til_v(4,:)) 
    % Calculo do u médio
    soma1 = sum(u_barra,1)/qtd_exp;    
    ubarra_med1(control,:,:) = soma1(1,control,:,:);
    soma2(control,:) = sum(ubarra_med1(control,:,:),2);
    ubarra_med(control) = sum(soma2(control,:),2)
    
end % Fim do for de controle

    %% Calculando os Indices de Desempenho
    %L22 = norm(q_til_v(1,:) + q_til_v(2,:) + q_til_v(3,:))
    %L2psi = norm(q_til_v(4,:)) 
    %L2 = norm(q_til_v(1:3,:))
    
%     for i=1:tempo_voo*10
%         for qpos=1:length(q)
%             u_barra(qpos) = sum(u_barra_med(:,qpos,i));
%             u_barra(qpos) = u_barra(qpos)/qtd_exp;
%         end
%         u_barra_i(i) = sum(u_barra(:,:));
%     end
%     
%    u_barra_sum = sum(u_barra_i)
  

save_data(qe_v, qed_v, F0_v, v_v, q_v, qd_v);


%plot_graphics(t_v, q_v, qd_v,dq_v,dqd_v,d2q_v,d2qd_v, u_v, v_v, F0_v, deltaF0_v, d_v, dB_v, us_v, control)