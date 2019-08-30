% function[norm_v] = experimento(THETA_GP)
clc
clear all
close all

%Control Settings

%Escolha o controle:

%Feedback Linearization ---------------------> control = 1
%H infinito ---------------------------------> control = 2
%H infinito com Rede Neural Adaptativa ------> control = 3
%H infinito com GP Adaptativo ---------------> control = 4
%H infinito com GP Offline ------------------> control = 5
%H infinito com GP Adaptativo do Artigo -----> control = 6

control = 1;

qtd_exp = 1;

for exp=1:qtd_exp
    
            
    format long
    path(path,'.\Controladores')
    path(path,'.\Dados')
    path(path,'.\Drone')
    path(path,'.\Plots')
    path(path,'.\Resultados')
    path(path,'.\trajectory')



    %% Variaveis de alteração

    %Trajactory Settings

    traj_des = 1; % Escolha: 1- Trajetória Elipse; 2 - Trajetória ponto final

    %%%%%%%%%%%%% Trajetória Elipse %%%%%%%%%%%%%%%%%%
    vel_d = 0.3; % Velocidade desejada em m/s

    A = 0.5; % Amplitude da trajetória
    omg = vel_d/A; % Velocidade desenada em rad/s 
    tempo_voo = 20;%tempo de voo em segundos

    %Sample Time
    dt = 0.01;
    N = tempo_voo/dt;

    %% Valor de incerteza do modelo do drone

    up =20; %incerteza em porcentagem
    a = -up/100;
    b = up/100;
    uncertain = (-a +((b-a)*rand()));


        for control=1:5

            [X_ga, q, dq, d2q, v, THETA_GP, qf,dqf, d2qf,traj_opt, to, tf, ax, ay, az, afi, k, kuncertain, Kp, Kd, T11, T12, T11_inv, T_zero, R, R_inv, THETA, sigma_quad, l, F0_gp, cov_gp, Zt, yt, etol, p_max, omega_quad, qe, qed] = initialization(uncertain);

            for i=1:N
                %% Planejamento de trajetória
                t = i*dt; 
                t_v(i)= t;

               [qd,dqd,d2qd,d3qd,to,tf,ax,ay,az,afi] = trajectory_desired(traj_des,q,dq,d2q, qf, dqf, d2qf, t, dt, traj_opt, vel_d, A, omg, to,tf,ax,ay,az,afi);


            %        q_til = [1;1;1;1];
            %     dq = [1;1;1;1];
            %     
            %     dq_til = [1;1;1;1];
            %     dqd = [2;2;2;2];
            %     
            %     d2qd = [0;0;0;0];
            % 
            % q(4) = pi;

            %     %% Calculo do erro da posição, velocidade e aceleração
                q_til = q(1:3) - qd(1:3); %diferença entre a posição desejada e atual
                q_til(4) = normalize_angle_f( q(4) - qd(4), -pi);
                dq_til = dq - dqd;
                d2q_til = d2q - d2qd;

            %     q_til =  qd(1:3) - q(1:3); %diferença entre a posição desejada e atual
            %     q_til(4) = normalize_angle_f( qd(4)-q(4), -pi);
            %     dq_til = dqd - dq;
            %     d2q_til =  d2qd- d2q;




                %% Controlador

                %[v,u,deltaF0, THETA_GP, ~, us, Zt,yt] = controller (control, Kp,Kd,k, q,dqd,d2qd, R_inv,T_zero, T11,T11_inv,q_til,dq_til,dq,qd, THETA, qe,qed, F0, F0_gp, sigma_quad,l, K_xx, K_ss, K_sx, THETA_GP, Zt,yt,omega_quad,i,dt);

                [v,u,F0, deltaF0, THETA_GP, THETA, us, Zt,yt] = controller_v2 (control, Kp,Kd,k, q,dqd,d2qd, R_inv,T_zero, T11, T11_inv,q_til,dq_til,dq,qd, THETA, qe, qed, X_ga, F0_gp, sigma_quad,l, THETA_GP, Zt,yt,etol, p_max,omega_quad,i,dt);


                %% Drone Model
                [d2q,dq,q,d,B] = droneDynamic(kuncertain,t,dt,q,dq,v,tempo_voo);

                qe = q + dq+ d2q;
                qed = qd' + dqd' + d2qd';


                %% Dados salvos para os plots

                %[u_v, v_v, us_v, F0_v, deltaF0_v, d_v, B_v, dB_v, q_v, dq_v, d2q_v, qe_v, qd_v, dqd_v, d2qd_v, qed_v, q_til_v,dq_til_v, d2q_til_v] = save_dataPlot(u, v, us, F0, deltaF0, d, B, q, dq, d2q, qe, qd, dqd, d2qd, qed, q_til,dq_til, d2q_til, i);

                %%%%%%% Controlador %%%%%%%        
                u_v(:,i) = u; %matriz do vetor da lei de controle
                v_v(:,i) = v; %matriz do vetor de controle
                us_v(:,i) = us; %Controle variavel da Rede Neural

                %%%%%% Modelo %%%%%%%%
                F0_v(:,i) = F0; %matriz do vetor F0
                deltaF0_v(:,i) = deltaF0; %matriz do vetor de incerteza de F0

                %%%%% Dinamica %%%%%2%%
                d_v(:,i) = d;
                B_v(:,:,i) = B;
                dB_v(:,i) = B*d; %Aceleraï¿½ï¿½o do distï¿½rbio

                %%%% Dados de posição, velocidade, aceleração e qe = q+dq+d2q  %%%%%%%%%%%%
                q_v(:,i) = q;
                dq_v(:,i) = dq;
                d2q_v(:,i) = d2q;        
                qe_v(:,i) =  qe;        



                %%%%% Posição desejada, Velocidade desejada, Aceleração desejada e qed = qd+dqd+d2qd  %%%%%%%
                qd_v(:,i) = qd;
                dqd_v(:,i) = dqd;
                d2qd_v(:,i) = d2qd;
                qed_v(:,i) =  qed;

                %%%% Erro de posição e Velocidade %%%%%%%%

                q_til_v(:,i) = q_til';
                dq_til_v(:,i) = dq_til';
                d2q_til_v(:,i) = d2q_til';

                %if(t>4)
                    qmed_til_v(exp,control,:,i) = q_til';
                    dqmed_til_v(exp,control,:,i) = dq_til';
                    d2qmed_til_v(exp,control,:,i) = d2q_til';

                     %drawnow;
                     %pause(0.1);            


                    normEPM_til(exp,control,i) =  norm(q_til_v(1:3,i));
                    normEPsiM_til(exp,control,i) = norm(q_til_v(4,i));
                    normU(exp,control,i) = norm(v_v(1:4,i));    
                %end
                t

            end


            if(control == 3)
                save_data(qe_v, qed_v, deltaF0_v, v_v, q_v, qd_v);
            end

           % clear_data(exp, control, qtd_exp, u_v, v_v, us_v, F0_v, deltaF0_v, d_v, B_v, dB_v, q_v, dq_v, d2q_v, qe_v, qd_v, dqd_v, d2qd_v, qed_v, q_til_v,dq_til_v, d2q_til_v, normEPM_til, normEPsiM_til, normU);

            clearvars -except N control exp qtd_exp dt t traj_des vel_d A omg tempo_voo uncertain t_v u_v v_v us_v F0_v deltaF0_v d_v B_v dB_v q_v dq_v d2q_v qe_v qd_v dqd_v d2qd_v qed_v q_til_v dq_til_v d2q_til_v qmed_til_v normEPM_til normEPsiM_til normU u_barra

            % Dados para calcular o indice de desempenho u médio
            u_barra(exp, control, :,:) = abs(v_v);

            if(exp == qtd_exp)
                %%%%%%%%%%Plota o grafico 3D para cada controle %%%%%%%%%%%%%%%%%
                plot3dim_exp(qd_v,q_v,control);
                plot_graphics_exp_x(t_v,q_v, qd_v, control);
                plot_graphics_exp_y(t_v,q_v, qd_v, control);
                plot_graphics_exp_z(t_v,q_v, qd_v, control);
                plot_graphics_exp_psi(t_v,q_v, qd_v, control);
            end


        end% Fim do for de controle
    

   
end % Fim do for de experimentos


  
[L22, L2psi, ubarra_med] = Desempenho(qtd_exp, qmed_til_v, u_barra, normEPM_til, normEPsiM_til, normU, t_v);
    

    
%save_data(qe_v, qed_v, F0_v, v_v, q_v, qd_v);


%plot_graphics(t_v, q_v, qd_v,dq_v,dqd_v,d2q_v,d2qd_v, u_v, v_v, F0_v, deltaF0_v, d_v, dB_v, us_v, control)