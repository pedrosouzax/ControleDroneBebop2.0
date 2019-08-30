% function[norm_v] = Controlador_Main(THETA_GP)
clc
clear all
close all

format long
path(path,'.\trajectory')

%Sample Time
dt = 0.1;

%Posições atuais
% x = str2double(get(gui.posAtual_x,'string'));
% y = str2double(get(gui.posAtual_y,'string'));
% z = str2double(get(gui.posA   tual_z,'string'));
% fi = str2double(get(gui.posAtual_yaw,'string'));

X_ga = [0.8021,  -12.0236,    0.7315];%melhor Matriz de transforma��o sem transposta

%X_ga = [X_ga(1), X_ga(2), X_ga(3)];

x = 0;
y=0;
z=0;
fi=0;

%Velodidades atuais
dx = 0;
dy = 0;
dz = 0;
dfi = 0; 

%Define posi��o desejada  
% xf = str2double(get(gui.posDes_x,'string'));
% yf = str2double(get(gui.posDes_y,'string'));
% zf = str2double(get(gui.posDes_z,'string'));
% fif = (pi/180)*(str2double(get(gui.posDes_yaw,'string')));
% xf = 3;
% yf = 1.4;
% zf = 2;
% fif = pi/2;

xf = 2.7;
yf = 2;
zf = 3.1;
fif = pi/3;

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

up =20; %incerteza em porcentagem
a = -up/100;
b = up/100;
kuncertain = k + (-a +((b-a).*rand(8,1))).*k;

% DEFINI��ES DE MATRIZES E CONSTANTES
%[ru, gama, t_11, t_12] = parameters_Hinf(0.7, 1.2, 0, 0.3) %parameters_Hinf(ru, q11, q12, q22)
% ru = 0.7;
% gama = 0.77;
% t_11 = -0.26; 
% t_12 = 3.5;
% ru = 0.70;
% gama = 0.767;
% t_11 = -0.3499;
% t_12 = 3.633;


T11 = X_ga(1) * eye(4);
T12 = X_ga(2) * eye(4);
T11_inv = inv(T11);

R = X_ga(3)*eye(4);

T_zero = [ T11          T12
          zeros(4)       eye(4)];
             
%Valor de THETA conhecido para rede neural
THETA = zeros(28,1);

%Valor de THETA conhecido para rede neural
% THETA_GP = [THETA_GP(1); THETA_GP(2); THETA_GP(3); THETA_GP(4)];

THETA_GP = [0.159832873135406;0.237741112701634;0.004570340314342;0.336590141950529] %0.039760438957400

%Trajectory variables
% vel_d = str2double(get(gui.speed,'string')); %Drone Speed
vel_d = 0.1;
traj_opt = 1;


% input_train = load('qe_v.txt');
% output_train = load('F0_v.txt');

% for(i=1:length(input_train))
%     input_train1(:,i) = (input_train(:,i) - min(input_train(:,i)))/(max(input_train(:,i))-min(input_train(:,i)));
% end
%  
% for(i=1:length(output_train))
%     output_train(:,i) = (output_train(:,i) - min(output_train(:,i)))/(max(output_train(:,i))-min(output_train(:,i)));
% end

% network = feedforwardnet([7 7 7 7]);
% net = create_train_feedforward(network,input_train(:,1:1000),output_train(1,1:1000));

% W1 = net.IW{1,1};
% b1 = net.b{1};
% 
% W2 = net.LW{2,1};
% b2 = net.b{2};
% 
% W3 = net.LW{3,2};
% b3 = net.b{3};

% network.trainParam.showWindow = false;
% 
% vet_saida = sim(net,input_train);
% 
% qe = [q' dq' qd' dq' dq'];

target = load('F0_v.txt');
Xtrain = load('qe_v.txt');
Xtest  = load('qed_v.txt');

sigma_quad = 1e-5;
l=10;

for j=1:4

K_xx(j,:,:) = diag(kernel(sigma_quad, Xtest(j,:), Xtest(j,:), l));
K_ss(j,:,:) = diag(kernel(sigma_quad, Xtrain(j,:), Xtrain(j,:), l));
K_sx(j,:,:) = diag(kernel(sigma_quad, Xtrain(j,:), Xtest(j,:), l));

end

% [F0_gp,cov_gp] = GPoff(sigma_quad, l, target', Xtrain, Xtest);
% 
% stdv = sqrt(diag(cov_gp(:,:,1)));
% 
% F0_gp_v = F0_gp;

valores_conhecidos = Xtest(:,2);
yt = target(:,2);

qe = zeros(4,1);

qed = zeros(4,1);

    for i=1:1000
    %% Planejamento de trajet�ria
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
        
        qd_v(:,i) = qd;
        dqd_v(:,i) = dqd;
        d2qd_v(:,i) = d2qd;

    %% Calculo do erro da posi��o, velocidade e acelera��o
    q_til(1:3) = qd(1:3) - q(1:3); %diferen�a entre a posi�ao desejada e atual
    q_til(4) = normalize_angle_f(qd(4) - q(4), -pi);
    %q_til(4) = wrapTo2Pi(qd(4) - q(4));
    
    dq_til = dqd - dq;
   
    d2q_til = d2qd - d2q;

    
%     dq_til(1:3) = dqd(1:3) - dq(1:3);
%     dq_til(4) =  normalize_angle_f(dqd(4) - dqd(4),-pi);
%     
%     d2q_til(1:3) = d2qd(1:3) - d2q(1:3);
%     d2q_til(4) =  normalize_angle_f(d2qd(4) - d2q(4),-pi);
    
   
    
    R_inv = inv(R);

    q_til_v(:,i) = q_til';
    dq_til_v(:,i) = dq_til';
    d2q_til_v(:,i) = d2q_til';
    
    

    %% PD Control    
    %u = Kp*q_til.' + Kd*dq_til;

    %% FeedBack Linearization
    %v = Feedback_linearization(k,q,d2qd,dqd,u);

    %% Modelo nominal        
    F0 = droneNominalModel(kuncertain, d2qd, dq_til, dqd, q_til, q, X_ga(3), X_ga(1), X_ga(2));

    %% H infinito
    u = -R_inv * [eye(4) zeros(4)] * T_zero * [dq_til
                                              q_til'];
                                          
    u(1:4)= tanh(u(1:4));
    %u= tanh(u);
    %u(1:3) = u(1:3)/norm(u(1:3));
    
    
   v = F0 + T11_inv*u;
%      if(i <=  999)
%      v = F0 + T11_inv*u;
%      else
%         vet_saida = net(qe');
%         v = vet_saida + T11_inv*u;
%      end
     
    
    u_v(:,i) = u;
    F0_v(:,i) = F0;
    v_v(:,i) = v;
    %vet_saida(:,i) = vet_saida;



    %% Controle robusto com rede reural adaptativo

     [E, THETA, us] = Rede_Neural(THETA, dt, i, q_til, dq_til, T_zero, T11,  q, dq, qd, dqd, d2qd);
%      
%      v = F0 +  E*THETA + T11_inv * u + us; 
%      
%      v = F0 +  T11_inv * u;
%      %v = E*THETA + T11_inv * u + us; 


% 
     F0_est_v(:,i) = E*THETA;
% 
%       Theta_v(:,i) = THETA;
% 
       us_v(:,i) = us;
       
%% Gaussian Processes com controle robusto
    
      %v = F0_gp(i) + T11_inv*u;
%      v = F0 + F0_gp(:,i) + T11_inv*u + us; % Estimando a incerteza

       
       
%% Gaussian Processes com controle Robusto Adaptativo com THETA_GP

% [mu, THETA_GP, us] = GP_theta_adaptativo(sigma_quad, K_xx, K_ss, K_sx, THETA_GP, dt, i, q_til, dq_til, T_zero,T11, q, dq, qd, dqd, d2qd)
% 
%  v = F0 + mu*THETA_GP + T11_inv*u + us;
%  
%  F0_gp_v(:,i) = mu*THETA_GP; %Delta F0
 
 
%% Gaussian Processes com controle Robusto Adaptativo do artigo

% valor_desejado = qed; 
% 
% [mu,cov,valores_conhecidos,valor_desejado] = GP_adaptativo_art(sigma_quad,l,valor_desejado,valores_conhecidos,yt,i);
% 
% v = F0 + mu(:,i) + T11_inv*u;
% 
% % valores_conhecidos(:,i) = valor_desejado; 
% yt(:,i+1) = mu(:,i);
% 
% 
% valores_conhecidos(:,i+1) = qe; 
% 
% F0_gp(:,i) = mu(:,i);
% stdv = sqrt(diag(cov));
% 

    %% Drone Model
    [d2q,dq,q,d,B] = droneDynamic(k,t,dt,q,dq,v);
    t
    
    d_v(:,i) = d;
    B_v(:,:,i) = B;
    nedB_v(:,i) = B*d; %Acelera��o do dist�rbio

    x = q(1);
    y = q(2);
    z = q(3);
    fi = q(4);

    qv(:,i) = q;
    dq_v(:,i) = dq;
    d2q_v(:,i) = d2q;
    
%     qe = [q' dq' qd' dqd' d2qd'];
%     qe_v(:,i) = [q' dq' qd' dqd' d2qd'];
%   

    qe = q' + dq'+ d2q';
    qe_v(:,i) =  qe;
    
    qed = qd' + dqd' + d2qd';
    qed_v(:,i) =  qed;
    
    
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
% 

    end
    
    norm_v = norm(q_til_v)
    %pause(0.1);
    
%      end
    
save('qe_v.txt','qe_v','-ascii')
save('qed_v.txt','qed_v','-ascii')
save('F0_v.txt','F0_v','-ascii')
% save('xx_v.txt','xx_v','-ascii')
% save('v_v.txt','v_v','-ascii')
save('qv.txt', 'qv','-ascii')

save('qd_v.txt', 'qd_v','-ascii')

    
% % fileID = fopen('abc.txt','w'); 
% % formatspec = '%4.2f,%4.2f,%4.2f,%4.2f-'
% % fprintf(fileID,formatspec,qv);
% % fclose(fileID);
% 
% norm_v = norm(q_til_v)
% % erro_relativo = (1 - ((vet_saida)./F0_v(1,:)))
% % error_square = norm(erro_relativo)
% % norm_v
% 
% % figure

% % plot(t_v,F0_v(1,:));
% % hold on
% % plot(t_v,F0_gp(1,:));

% figure (1)
% 
% 
% plot(t_v, F0_gp(1,:)-2*stdv', '-.r');
% 
% hold on
% 
% plot(t_v, F0_gp(1,:)+2*stdv', '-.r');
% x2 = [t_v, fliplr(t_v)];
% inBetween = [F0_gp(1,:)+2*stdv', fliplr(F0_gp(1,:)-2*stdv')];
% fill(x2, inBetween, 'g');
% 
% plot(t_v, F0_gp(1,:), '-.b');
% title('F0 GP')
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_xtil = ylabel('$\tilde{x}(m)$','Interpreter','latex');
% set(y_xtil,'FontSize',17);
% 
% figure(2)
% plot(t_v, F0_gp(2,:), 'r');
% title('F0 GP')
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_ytil = ylabel('$\tilde{y}(m)$','Interpreter','latex');
% set(y_ytil,'FontSize',17);
% 
% figure(3)
% plot(t_v, F0_gp(3,:), 'b');
% title('F0 GP');
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_ztil = ylabel('$\tilde{z}(m)$','Interpreter','latex');
% set(y_ztil,'FontSize',17);
% 
% figure(4)
% plot(t_v, F0_gp(4,:), 'm');
% title('F0 GP')
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_psitil = ylabel('$\tilde{\psi}(rad)$','Interpreter','latex');
% set(y_psitil,'FontSize',17);

% figure(5)
% plot(t_v, F0_gp(1,:)-2*stdv', '-.r'); 

%plot_graphics(t_v, qv, qd_v,dq_v,dqd_v,d2q_v,d2qd_v, u_v, v_v, F0_v, F0_est_v, F0_gp_v, q_til_v, d_v, B_v, d_v, us_v);
