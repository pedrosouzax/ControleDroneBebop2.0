function [ mu, THETA_GP, us] = GP_theta_adaptativo( sigma_quad, l,THETA_GP, T, i, q_til, dq_til, T_zero,T11, qe, qed)
 
    
    K_sx = zeros(4,1); 
    
    for j=1:4 
        
             %L = chol(K_xx(j,i,i) + sigma_quad*eye(1));
             %J  = L\Ytrain(:,j);
             %alpha = L'\J;
             %K_sx(:,j) = kernel(sigma_quad,qed(j),qe,l,length(qe));
             K_sx(j) = kernel(sigma_quad,qed(j),qe(j),l,1);
             %v = L\K_sx(j,i,i)';
             %cov(:,:,j) = K_ss(j,i,i) - (v*v');
    end
    
    
    
    %mu = diag(mu);
  
    
    
    Z = 30*eye(4);
    
    K_sx = diag(K_sx);
    
    N1 = 0.00005; % É o Ka do Roberto
    delta = 0.002;
%     
%    N1 = 5; % É o Ka do Roberto
%     delta = 2;    

    B2 = [eye(4);zeros(4)];
    aux_THETAd = -inv(Z)*K_sx*T11*B2'*T_zero*[dq_til; q_til];
        
    % Algoritmo 3.31da dissertação do Roberto                                                        
    if ((THETA_GP'*THETA_GP <= N1) | ((THETA_GP'*THETA_GP > N1) & (THETA_GP'*aux_THETAd <= 0)))
        THETA_d = aux_THETAd;
    else
        THETA_d = (aux_THETAd - (((((THETA_GP'*THETA_GP - N1)*THETA_GP'*aux_THETAd))./(delta*THETA_GP')*THETA_GP)*THETA_GP));
    end
        
    % Atualização dos parâmetros
    if i>2
    THETA_GP = THETA_GP + THETA_d*T;
    end
    
    q_til
    mult_q_til = (q_til')*(q_til);
    kxe=(1/0.2)*sqrt(abs(mult_q_til)); %Não sei porquê essa função é definida assim
    
    sgn = sign(T11*[eye(4) zeros(4)]*T_zero*[dq_til; q_til]); 
    sgn
    
    us = -kxe*sgn;
    
    us
    
    i
    
    mu = K_sx*THETA_GP;
    
    THETA_GP
    
    
%     if (i <2)
%         us = 0;
%     end
    
    %us = -inv(T11)*kxe*[B2'*T_zero*[dq_til; q_til.']]./[B2'*T_zero*[dq_til; q_til.']];
    
    %us(:,cont)=-inv(T11)*kxe*[B2'*To*xtil(:,cont)]./[B2'*To*xtil(:,cont)];
    
    
                                           
   %us=-inv(T11)*kxe*[[eye(4) zeros(4)]'*T_zero*[dq_til; q_til.']]./[[eye(4) zeros(4)]'*T_zero*[dq_til
   %                                                                                             q_til.']];
    

end

