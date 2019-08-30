%% Gaussian Processes offline: Calcula o modelo antes da execu��o

%Par�metros: 
%Y � o target, 
%Xtest s�o os pontos que eu quero estimar a fun��o
%Xtrain s�o os pontos conhecidos

%sigma_quad � 
%l �

function [mu, cov] = GPoff(sigma_quad, l, Ytrain, Xtest, Xtrain)
    for j=1:4

        K_xx = kernel(sigma_quad, Xtrain(j,:), Xtrain(j,:), l, length(Xtrain));
        K_ss = kernel(sigma_quad, Xtest(j,:), Xtest(j,:), l, length(Xtrain));
        K_sx = kernel(sigma_quad, Xtrain(j,:), Xtest(j,:), l, length(Xtrain));


        L = chol(K_xx + sigma_quad*eye(length(Xtrain)));
%         J  = L\Ytrain(:,j);
%         alpha = L'\J;
%         mu(j,:) = K_sx'*alpha;
        mu = (K_sx/L)/L'*Ytrain;
        v = L\K_sx';
        cov(:,:,j) = K_ss - (v*v');
        
        %cov(:,:,j) = K_ss - K_sx*L*K_sx';
        
        
        %cov_n(:,:,j) = abs(cov');
        % f_post = mvnrnd(mu,cov_n);
        %f_post = mu + diag(cov);

%         sigma_quad = 0.1;
%         l = 100;
% 
%         K= kernel(sigma_quad,q_d(j),q_d(j),l)
%         L = chol(K+sigma_quad*eye(1))
%         J = (L\q_d(j))'
%         alfa = L'\J';
%         f_media(:,j) = kernel(sigma_quad,q(j),q_d(j),l)'*alfa;
%         v = L\kernel(sigma_quad,q(j),0,1);
%         var = kernel(sigma_quad,q(j),q(j),1) - v'*v;

    end
    
end



