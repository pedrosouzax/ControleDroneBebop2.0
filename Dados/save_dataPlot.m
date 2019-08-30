function [u_v, v_v, us_v, F0_v, deltaF0_v, d_v, B_v, dB_v, q_v, dq_v, d2q_v, qe_v, qd_v, dqd_v, d2qd_v, qed_v, q_til_v,dq_til_v, d2q_til_v] = save_dataPlot(u, v, us, F0, deltaF0, d, B, q, dq, d2q, qe, qd, dqd, d2qd, qed, q_til,dq_til, d2q_til, i)
        
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
end
        