Xtrain = [0:0.5:2*pi]';
Ytrain = sin(Xtrain);
Xtest = [0:0.5:2*pi]';

sigma_quad = 1;
l=1;

K_xx = diag(kernel(sigma_quad, Xtest, Xtest, l));
K_ss = diag(kernel(sigma_quad, Xtrain, Xtrain, l));
K_sx = diag(kernel(sigma_quad, Xtrain, Xtest, l));


L = chol(K_xx + sigma_quad*eye(1));
J  = L\Ytrain;
alpha = L'\J;
mu = K_sx'*alpha;
v = L\K_sx';
cov = K_ss - (v*v');
%cov = K_ss - K_sx*L*K_sx';

figure (1)
plot(Xtrain,Ytrain);

figure (2)
plot(Xtest,mu);