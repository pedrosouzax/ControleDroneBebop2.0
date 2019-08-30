clc
clear all

path(path,'.\trajectory')

%Sample Time
dt = 0.1;

%Posições atuais
% x = str2double(get(gui.posAtual_x,'string'));
% y = str2double(get(gui.posAtual_y,'string'));
% z = str2double(get(gui.posAtual_z,'string'));
% fi = str2double(get(gui.posAtual_yaw,'string'));

x = 0;
y=0;
z=0;
fi= 0;

%Velodidades atuais
dx = 0;
dy = 0;
dz = 0;
dfi = 0; 

%Define posição desejada  
% xf = str2double(get(gui.posDes_x,'string'));
% yf = str2double(get(gui.posDes_y,'string'));
% zf = str2double(get(gui.posDes_z,'string'));
% fif = (pi/180)*(str2double(get(gui.posDes_yaw,'string')));
xf = 3;
yf = 3;
zf = 3;
fif = pi/2;
q = [x ; y; z; fi];
 
dq = [dx; dy; dz; dfi];
d2q = zeros(4,1);

qd = [xf; yf; zf; fif];

qf=[xf;yf;zf;fif];
dqf=[0;0;0;0];
d2qf=[0;0;0;0];

Kp = [50 0 0 0 
      0 50 0 0
      0  0 50 0 
      0 0 0 50];

Kd = [10 0 0 0  
      0 10 0 0
      0 0 10 0 
      0 0 0 10];

%MATRIZ DE CONTROLE VIRTUAL
v = [0; 0; 0; 0]; 
 
%CONSTANTES PARA O MODELO UTILIZADAS NAS MATRIZES M E N
k1 = 1.74199;
k2 = 0.94016;
k3 = 1.54413;
k4 = 0.89628;
k5 = 3.34885;
k6 = 3.29467;
k7 = 6.51209;
k8 = 3.92187;

k = [k1; k2; k3; k4; k5; k6; k7; k8];

% DEFINIÇÕES DE MATRIZES E CONSTANTES
% ru = 0.7;
% gama = 0.77;
% t_11 = -0.35;
% t_12 = 1.3;

q11 = 0.5;
q22 = 0;
q12 = 0.5;
ru = 1;

%[gamma,t_11,t_12,~] = hinfinityParameters(q11,q22,q12,ru)
t_11 = 1
t_12 = 1
T11 = t_11 * eye(4);

T12 = t_12 * eye(4);
T11_inv = inv(T11);

R = ru*eye(4);

T_zero = [ T11          T12
          zeros(4)       eye(4)];
             
%Valor de THETA conhecido para rede neural
THETA = zeros(28,1);

%Trajectory variables
% vel_d = str2double(get(gui.speed,'string')); %Drone Speed
vel_d = 0.1;
traj_opt = 1;

    for i=1:1000
    %% Planejamento de trajetória
    t = i*dt; 
    t_v(i)= t;

        if traj_opt == 1
            to = t;
            td = norm(qf-q)/vel_d;
            tf = to+td;
            [ax,ay,az,afi] = trajectory_par_pol_f(dt,to,tf,q,dq,d2q,qf,dqf,d2qf);
            traj_opt = 0;
        end

        if t<=tf
            [xd,dxd,d2xd,d3xd] = traj_pol_f(ax,t,to);
            [yd,dyd,d2yd,d3yd] = traj_pol_f(ay,t,to);
            [zd,dzd,d2zd,d3zd] = traj_pol_f(az,t,to);
            [fid,dfid,d2fid,d3fid] = traj_pol_f(afi,t,to);
        end

        qd = [xd; yd; zd; fid];
        dqd = [dxd; dyd; dzd; dfid];
        d2qd = [d2xd; d2yd; d2zd; d2fid];

        xd_v(i)=xd;
        yd_v(i)=yd;
        zd_v(i)= zd;
        fid_v(i) =fid;
        dxd_v(i)=dxd;
        dyd_v(i)=dyd;
        dzd_v(i)= dzd;
        dfid_v(i) =dfid;
        d2xd_v(i)=d2xd;
        d2yd_v(i)=d2yd;
        d2zd_v(i)=d2zd;
        d2fid_v(i) =d2fid;

    %% Calculo do erro da posição, velocidade e aceleração
    q_til(1:3,1) = qd(1:3) - q(1:3); %diferença entre a posiçao desejada e atual
    q_til(4,1) = normalize_angle_f(qd(4) - q(4), -pi);
    
    
    dq_til = dqd - dq;
   
    d2q_til = d2qd - d2q;

    
%     dq_til(1:3) = dqd(1:3) - dq(1:3);
%     dq_til(4) =  normalize_angle_f(dqd(4) - dqd(4),-pi);
%     
%     d2q_til(1:3) = d2qd(1:3) - d2q(1:3);
%     d2q_til(4) =  normalize_angle_f(d2qd(4) - d2q(4),-pi);
    
   
    
    R_inv = inv(R);

    q_til_v(:,i) = q_til;

%     %% FeedBack Linearization + PD Control    
%     u = Kp*q_til + Kd*dq_til;
%     v = Feedback_linearization(k,q,d2qd,dqd,u);

    %% Modelo nominal        
    F0 = droneNominalModel(k, d2qd, dq_til, dqd, q_til, q, ru, gamma, t_11, t_12);

    %% H infinito 
     u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til; q_til];
     v =  T11_inv*u + F0;
    % Saturação da lei de controle
%     u(1:4)= tanh(u(1:4));

    u_v(:,i) = u;
    F0_v(:,i) = F0;
    v_v(:,i) = v;

    %% Gaussian Processes
    % [v_gp] = GP(q_d, q, u, T11_inv);
    % 
    % v = v_gp + F0 ;

    %% Controle robusto adaptativo

    % [Matriz_E, THETA, us] = Rede_Neural(THETA, dt, i, q_til, dq_til, T_zero, T11,  q, dq, q_d, dq_d, d2q_d);

    %v = F0 + Matriz_E*THETA + T11_inv * u + us; 


    %% Drone Model
    [d2q,dq,q] = droneDynamic(k,dt,q,dq,v);

    x = q(1);
    y = q(2);
    z = q(3);
    fi = q(4);

    qv(:,i) = q;    
    
%     set(gcf,'CurrentAxes',gui.graph_xyz)
%     move_object_f(quad1_obj,[0;0;fi],[x;y;z]);
%     drawnow
%     
%     set(gcf,'CurrentAxes',gui.graph_xy)
%     move_object_f(quad2_obj,[0;0;fi],[x;y;z]);
%     drawnow
%     
%     set(gui.posAtual_x,'string',x);
%     set(gui.posAtual_y,'string',y);
%     set(gui.posAtual_z,'string',z);
%     set(gui.posAtual_yaw,'string',fi*180/pi);
%     
%     pause(0.01);


    end




%% Plot da posição no tempo
figure(1)
subplot(2,2,1), plot(t_v, qv(1,:),'g');
title('Gráfico posição X');
xlabel('Tempo(s)');
ylabel('X(m)');
%figure(2)
subplot(2,2,2), plot(t_v, qv(2,:),'r');
title('Gráfico posição Y');
xlabel('Tempo (s)');
ylabel('Y (m)');
%figure(3)
subplot(2,2,3), plot(t_v, qv(3,:),'b');
title('Gráfico posição Z');
xlabel('Tempo (s)');
ylabel('Z(m)');
%figure(4)
subplot(2,2,4), plot(t_v, qv(4,:),'m');
title('Gráfico posição PSI');
xlabel('Tempo (s)');
ylabel('PSI (rad)');

%% Plot da Lei de controle no tempo
figure(2)
subplot(2,2,1), plot(t_v, u_v(1,:), 'g');
title('Lei de Controle (u)')
xlabel('Tempo (s)');
ylabel('u (x)');

subplot(2,2,2), plot(t_v, u_v(2,:), 'r');
title('Lei de Controle (u)')
xlabel('Tempo (s)');
ylabel('u (y)');

subplot(2,2,3), plot(t_v, u_v(3,:), 'b');
title('Lei de Controle (u)')
xlabel('Tempo (s)');
ylabel('u (z)');

subplot(2,2,4), plot(t_v, u_v(4,:), 'm');
title('Lei de Controle (u)')
xlabel('Tempo (s)');
ylabel('u (psi)');

% %% Plot do Vetor de Controle no tempo
% figure(3)
% subplot(2,2,1), plot(t_v, v_v(1,:), 'g');
% title('Vetor de Controle (v)')
% xlabel('Tempo (s)');
% ylabel('v (x)');
% 
% subplot(2,2,2), plot(t_v, v_v(2,:), 'r');
% title('Vetor de Controle (v)')
% xlabel('Tempo (s)');
% ylabel('v (y)');
% 
% subplot(2,2,3), plot(t_v, v_v(3,:), 'b');
% title('Vetor de Controle(v)')
% xlabel('Tempo (s)');
% ylabel('v (z)');
% 
% subplot(2,2,4), plot(t_v, v_v(4,:), 'm');
% title('Vetor de Controle (v)')
% xlabel('Tempo (s)');
% ylabel('v (psi)');
% 
% %% Plot do Modelo Nominal no tempo
% figure(4)
% subplot(2,2,1), plot(t_v, F0_v(1,:), 'g');
% title('F0 no tempo')
% xlabel('Tempo (s)');
% ylabel('F0 (x)');
% 
% subplot(2,2,2), plot(t_v, F0_v(2,:), 'r');
% title('F0 no tempo');
% xlabel('Tempo (s)');
% ylabel('F0 (y)');
% 
% subplot(2,2,3), plot(t_v, F0_v(3,:), 'b');
% title('F0 no tempo');
% xlabel('Tempo (s)');
% ylabel('F0 (z)');
% 
% subplot(2,2,4), plot(t_v, F0_v(4,:), 'm');
% title('F0 no tempo');
% xlabel('Tempo (s)');
% ylabel('F0 (psi)');

%% Plot erro de posição
figure (5)
subplot(2,2,1), plot(t_v, q_til_v(1,:), 'g');
title('Erro de posição')
xlabel('Tempo(s)');
ylabel('Cordenada x (m)');

subplot(2,2,2), plot(t_v, q_til_v(2,:), 'r');
title('Erro de posição')
xlabel('Tempo(s)');
ylabel('Cordenada y (m)');

subplot(2,2,3), plot(t_v, q_til_v(3,:), 'b');
title('Erro de posição')
xlabel('Tempo(s)');
ylabel('Cordenada z (m)');

subplot(2,2,4), plot(t_v, q_til_v(4,:), 'm');
title('Erro de posição')
xlabel('Tempo(s)');
ylabel('Cordenada psi (rad)');

