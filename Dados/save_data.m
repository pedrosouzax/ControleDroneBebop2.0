function [] = save_data(qe_v, qed_v, F0_v, v_v, q_v, qd_v)
    save('.\Dados\qe_v.txt','qe_v','-ascii')
    save('.\Dados\qed_v.txt','qed_v','-ascii')
    save('.\Dados\F0_v.txt','F0_v','-ascii')
    %save('xx_v.txt','xx_v','-ascii')
    save('.\Dados\v_v.txt','v_v','-ascii')
    save('.\Dados\q_v.txt', 'q_v','-ascii')
    
    save('.\Dados\qd_v.txt', 'qd_v','-ascii')
end