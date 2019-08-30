function [ E, THETA, us] = Rede_Neural( THETA, T, i, q_til, dq_til, T_zero,T11, q, dq, qd, dqd, d2qd)
 
    pk=7;%número de neurônios
    n=4;%número de redes neurais
    
    %Pesos
    for k=1:n
        for i=1:pk
            W(k,i,:) = [-1, -1, -1, -1, -1,-1, -1,-1, 1,1,1,1,1,1,1,1,1,1,1,1];
        end
    end
    
    %Bias
    for k=1:n
       B(k,:) = [ -3, -2, -1, 0, 1, 2, 3]; 
    end
   
    qe = [q' dq' qd' dqd' d2qd']; 
    
    E =zeros(n,pk*n);
    
    for k=1:n
        for i=1:pk
            sumwq=0;
            for j=1:5*n
               sumwq = W(k,i,j)*qe(j) + sumwq;
            end
            H = tanh(sumwq + B(k,i));
            E(k,i+pk*(k-1)) = H;
        end
    end
    
    
    %Z = 30*eye(28);
    
    Z = 3000*eye(28);
    
    N1 = 0.00005; % É o Ka do Roberto
    delta = 0.002;

    B2 = [eye(4);zeros(4)];
    aux_THETAd = -inv(Z)*E'*T11*B2'*T_zero*[dq_til; q_til];
        
    % Algoritmo 3.31da dissertação do Roberto                                                        
    if ((THETA'*THETA <= N1) | ((THETA'*THETA > N1) & (THETA'*aux_THETAd <= 0)))
        THETA_d = aux_THETAd;
    else
        THETA_d = (aux_THETAd - (((((THETA'*THETA - N1)*THETA'*aux_THETAd))./(delta*THETA'*THETA))*THETA));
    end
        
    % Atualização dos parâmetros
    if i>2
    THETA = THETA + THETA_d*T;
    end
    
    q_til
    mult_q_til = (q_til')*(q_til);
    kxe=(1/0.2)*sqrt(abs(mult_q_til)); %Não sei porquê essa função é definida assim
    
    sgn = sign(T11*[eye(4) zeros(4)]*T_zero*[dq_til; q_til]); 
    sgn
    
    us = -kxe*sgn;
    
    us
    
    i
    
%     if (i <2)
%         us = 0;
%     end
%     
    %us = -inv(T11)*kxe*[B2'*T_zero*[dq_til; q_til.']]./[B2'*T_zero*[dq_til; q_til.']];
    
    %us(:,cont)=-inv(T11)*kxe*[B2'*To*xtil(:,cont)]./[B2'*To*xtil(:,cont)];
    
    
                                           
   %us=-inv(T11)*kxe*[[eye(4) zeros(4)]'*T_zero*[dq_til; q_til.']]./[[eye(4) zeros(4)]'*T_zero*[dq_til
   %                                                                                             q_til.']];
    

end

