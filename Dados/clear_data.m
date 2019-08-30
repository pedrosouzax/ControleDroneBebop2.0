function [] = clear_data(exp, control, qtd_exp, u_v, v_v, us_v, F0_v, deltaF0_v, d_v, B_v, dB_v, q_v, dq_v, d2q_v, qe_v, qd_v, dqd_v, d2qd_v, qed_v, q_til_v,dq_til_v, d2q_til_v, normEPM_til, normEPsiM_til, normU);
    clearvars -except control exp qtd_exp u_v v_v us_v F0_v deltaF0_v d_v B_v dB_v q_v dq_v d2q_v qe_v qd_v dqd_v d2qd_v qed_v q_til_v dq_til_v d2q_til_v normEPM_til normEPsiM_til normU
end

% %Limpando variaveis para o proximo experimento como era feito antes 
% q = zeros(4,1);
% dq = zeros(4,1);
% d2q = zeros(4,1);
% d3q = zeros(4,1);
% qe = zeros(4,1);
% dqe = zeros(4,1);
% qd = zeros(4,1);
% dqd = zeros(4,1);
% d2qd = zeros(4,1);
% d3qd = zeros(4,1);
% qe = zeros(4,1);
% qed = zeros(4,1);
% q_til = zeros(1,4);
% dq_til = zeros(1,4);
% d2q_til = zeros(1,4);
% clear us
% clear u; %matriz do vetor da lei de controle
% clear v; %matriz do vetor de controle
% clear F0; %matriz do vetor F0
% clear deltaF0; %matriz do vetor de incerteza de F0