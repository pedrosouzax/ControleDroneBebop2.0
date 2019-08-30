close all; clear all; clc;

Y = load('F0_v.txt'); 
% Xtrain = load('qv.txt');
% Xtest  = load('qd_v.txt');
Xtrain = load('q_v.txt');
Xtest  = load('qd_v.txt');
tempo = 1:1:2000;



Ytrain = Y(1,:);

Ztrain  = Xtrain(1,:);
Ztest = Xtest(1,:);


sigma_quad = 0.1;
l = 100;

% [x,y] = meshgrid(Ztest,Ztest);
% K_zz = kernel_v2(sigma_quad, x, y, l);
% 
% [x,y] = meshgrid(Ztrain,Ztrain);
% K_ss = kernel_v2(sigma_quad, x, y, l);
% 
% [x,y] = meshgrid(Ztrain,Ztest);
% K_sz = kernel_v2(sigma_quad, x, y, l);


K_zz = kernel(sigma_quad, Ztest, Ztest, l, 2000);
 
K_ss = kernel(sigma_quad, Ztrain, Ztrain, l, 2000);

K_sz = kernel(sigma_quad, Ztrain, Ztest, l, 2000);



 L = chol(K_ss + (1e-9)*eye(2000));
% J  = L\Ytrain';
% alpha = L'\J;
% mu = K_sz'*alpha;
mu = (K_sz/L)/L'*Ytrain';
%mu = K_sz*L*Ytrain';
v = L\K_sz;
cov = K_ss - (v*v');
cov_n = abs(cov');
% f_post = mvnrnd(mu,cov_n);
f_post = mu + diag(cov_n);




% plot(Ztest, f_post,'r')
figure
hold on
plot(tempo, Ytrain,'b');
plot(tempo, mu, '-.r');
%plot(tempo, f_post,'y');
%plot(Ztest, diag(cov_n), 'g');


