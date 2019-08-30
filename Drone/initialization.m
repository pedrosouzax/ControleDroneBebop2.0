function [X_ga, q, dq, d2q, v, THETA_GP, qf,dqf, d2qf,traj_opt, to, tf, ax, ay, az, afi, k, kuncertain, Kp, Kd, T11, T12, T11_inv, T_zero, R, R_inv, THETA, sigma_quad, l, F0_gp, cov_gp, Zt, yt, etol, p_max, omega_quad, qe, qed] = inicialization(uncertain)

%X_ga = [0.8021,  -12.0236,    0.7315];%melhor Matriz de transformaï¿½ï¿½o sem transposta
%X_ga = [0.187042762538844, 0.199239778549170, 0.289458041143591];

%X_ga = [ X_ga(1), X_ga(2), X_ga(3)];
%X_ga = [2.213594, 2.969848, 5];
%X_ga = [-3.752886082944024 -14.517295616416320 -0.218471681096723];
X_ga = [-1.088275985782010,-5.523761210177261,0.944310599276061];

x = 0;
y=0;
z=0;
fi=0;

%Velodidades atuais
dx = 0;
dy = 0;
dz = 0;
dfi = 0; 

% Posições Finais desejada- NÃO ALTERE AQUI - ALTERE MAIS ABAIXO
xf = 0;
yf = 0;
zf = 0;
fif = 0;

q = [x ; y; z; fi]; 
dq = [dx; dy; dz; dfi];
d2q = zeros(4,1);


%MATRIZ DE CONTROLE VIRTUAL
v = [0; 0; 0; 0]; 
 


% DEFINIï¿½ï¿½ES DE MATRIZES E CONSTANTES
%[ru, gama, t_11, t_12] = parameters_Hinf(0.7, 1.2, 0, 0.3) %parameters_Hinf(ru, q11, q12, q22)
% ru = 0.7;
% gama = 0.77;
% t_11 = -0.26; 
% t_12 = 3.5;
% ru = 0.70;
% gama = 0.767;
% t_11 = -0.3499;
% t_12 = 3.633;



            
%Valor de THETA conhecido para rede neural
% THETA_GP = [THETA_GP(1); THETA_GP(2); THETA_GP(3); THETA_GP(4)];

%THETA_GP = [0.159832873135406;0.237741112701634;0.004570340314342;0.336590141950529] %0.039760438957400
THETA_GP =   [0.612879500245117;
               0.182171338091438;
               0.260273000706138;
               0.457826272500228];

sum_qtil = 0;
sum_psitil = 0;
sum_u = 0;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ATENÇÃO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Definições pré-definidas das trajetórias (Altere apenas se souber o que está fazendo)

qf=[xf;yf;zf;fif];
dqf=[0;0;0;0];
d2qf=[0;0;0;0];

traj_opt = 1;
to = 0;
tf = 0;
ax=0;
ay=0;
az=0;
afi=0;

%Definições pré-definidas dos controladores (Altere apenas se souber o que está fazendo)


%------------------------------- Feedback Linearization -----------------%

%MATRIZES DE CONSTANSTE PROPORCIAL E DERIVATIVO PARA O PD
Kp = [50 0 0 0 
      0 50 0 0
      0  0 50 0 
      0 0 0 50];

Kd = [10 0 0 0  
      0 10 0 0
      0 0 10 0 
      0 0 0 10];

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

kuncertain = k + uncertain*k;

% kuncertain = 1.2*k;

%--------------------------------- H infinito ----------------------------%

T11 = X_ga(1) * eye(4);
T12 = X_ga(2) * eye(4);
T11_inv = inv(T11);

T_zero = [ T11          T12
          zeros(4)       eye(4)];


R = X_ga(3)*eye(4);
R_inv = inv(R);

%------------------------- Rede Neural Adaptativa ------------------------%

%Valor de THETA conhecido para rede neural
THETA =   [0.049424539731510;
   0.061707295451636;
   0.000461925433276;
   0.084625949075640;
  -0.016795187320992;
   0.033830358165730;
   0.194154634138103;
   0.182778602748309;
   0.012818117397264;
   0.092727957537253;
   0.149779259989359;
   0.033765763435889;
   0.140992308309683;
   0.137906576754882;
   0.124481892922731;
   0.011478649662710;
   0.075406168428600;
   0.158245482247798;
   0.110618521134542;
   0.055980525947054;
   0.101501855994308;
   0.032678259167206;
   0.091563683432064;
   0.165588211089541;
   0.146271708852501;
   0.055041114197478;
   0.131968142153905;
   0.074970661079984];

%--------------------------GP configurações Gerais------------------------%
sigma_quad = 1e-5;
l=10;

target = load('F0_v.txt');
Xtrain = load('qe_v.txt');
Xtest  = load('qed_v.txt');

%------------------------------ GP Offline -------------------------------%

[F0_gp,cov_gp] = GPoff(sigma_quad, l, target', Xtrain, Xtest);

% stdv = sqrt(diag(cov_gp(:,:,1)));

%------------------------- GP com THETA Adaptativo -----------------------%

% for j=1:4
% 
%     K_xx(j,:,:) = kernel(sigma_quad, Xtest(j,:), Xtest(j,:), l, length(Xtest));
%     K_ss(j,:,:) = kernel(sigma_quad, Xtrain(j,:), Xtrain(j,:), l, length(Xtrain));
%     K_sx(j,:,:) = kernel(sigma_quad, Xtrain(j,:), Xtest(j,:), l, length(Xtrain));
% 
% end

%------------------------- GP com THETA Adaptativo - Artigo -----------------------%
    Zt = Xtest(:,2);
    yt = target(:,2);
    
    qe = zeros(4,1);
    qed = zeros(4,1);
    
    etol = 0.1;
    p_max = 30;

    omega_quad = 1; %ruido branco        


end
