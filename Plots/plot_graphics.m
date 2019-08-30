function [] = plot_graphics(t_v, q_v, qd_v,dq_v,dqd_v,d2q_v,d2qd_v, u_v, v_v, F0_v, deltaF0_v, d_v, dB_v, us_v, control)
%% Plot da posição no tempo
figure
plot(t_v, q_v(1,:),'g');
%title('$Grafico posicao X$','Interpreter','latex');
%title('Gráfico posição X');
xlabel('$Tempo(s)$','Interpreter','latex');
y_x = ylabel('$x(m)$','Interpreter','latex');
set(y_x, 'FontSize', 17);

figure
plot(t_v, q_v(2,:),'r');
%title('Gráfico posição Y');
xlabel('$Tempo(s)$','Interpreter','latex');
y_y = ylabel('$y(m)$','Interpreter','latex');
set(y_y, 'FontSize', 17);

figure
plot(t_v, q_v(3,:),'b');
%title('Gráfico posição Z');
xlabel('$Tempo(s)$','Interpreter','latex');
y_z = ylabel('$z(m)$','Interpreter','latex');
set(y_z, 'FontSize', 17);

figure
plot(t_v, q_v(4,:),'m');
%title('Gráfico posição PSI');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psi = ylabel('$\psi (rad)$','Interpreter','latex');
set(y_psi, 'FontSize', 17);

%% Plot da Lei de controle no tempo
figure
plot(t_v, u_v(1,:), 'g');
%title('Lei de Controle (u)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_u_x = ylabel('$u_x$','Interpreter','latex');
set(y_u_x, 'FontSize', 20);

figure
plot(t_v, u_v(2,:), 'r');
%title('Lei de Controle (u)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_u_y = ylabel('$u_y$','Interpreter','latex');
set(y_u_y, 'FontSize', 20);

figure
plot(t_v, u_v(3,:), 'b');
%title('Lei de Controle (u)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_u_z = ylabel('$u_z$','Interpreter','latex');
set(y_u_z,'FontSize',20);

figure
plot(t_v, u_v(4,:), 'm');
%title('Lei de Controle (u)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_u_psi = ylabel('$u_{\psi}$','Interpreter','latex');
set(y_u_psi,'FontSize',20);

%% Plot do Vetor de Controle no tempo
figure
plot(t_v, v_v(1,:), 'g');
%title('Vetor de Controle (v)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_v_x = ylabel('$v_x (m/s)$','Interpreter','latex');
set(y_v_x,'FontSize',20);

figure
plot(t_v, v_v(2,:), 'r');
%title('Vetor de Controle (v)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_v_y = ylabel('$v_y (m/s)$','Interpreter','latex');
set(y_v_y,'FontSize',20);

figure
plot(t_v, v_v(3,:), 'b');
%title('Vetor de Controle(v)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_v_z = ylabel('$v_z (m/s)$','Interpreter','latex');
set(y_v_z,'FontSize',20);

figure
plot(t_v, v_v(4,:), 'm');
%title('Vetor de Controle (v)')
xlabel('$Tempo(s)$','Interpreter','latex');
y_v_psi = ylabel('$v_{psi} (rad/s)$','Interpreter','latex');
set(y_v_psi,'FontSize',20);

%% Plot do Vetor de Controle Variável (us)

if(control == 3 || control == 5)
    figure
    plot(t_v, us_v(1,:), 'g');
%    title('Controle Variável (us)');
    xlabel('$Tempo(s)$','Interpreter','latex');
    y_v_x = ylabel('$us_x$','Interpreter','latex');
    set(y_v_x,'FontSize',20);

    figure
    plot(t_v, us_v(2,:), 'r');
%   title('Controle Variável (us)');
    xlabel('$Tempo(s)$','Interpreter','latex');
    y_v_y = ylabel('$us_y$','Interpreter','latex');
    set(y_v_y,'FontSize',20);

    figure
    plot(t_v, us_v(3,:), 'b');
%    title('Controle Variável (us)');
    xlabel('$Tempo(s)$','Interpreter','latex');
    y_v_z = ylabel('$us_z$','Interpreter','latex');
    set(y_v_z,'FontSize',20);

    figure
    plot(t_v, us_v(4,:), 'm');
%   title('Controle Variável (us)');
    xlabel('$Tempo(s)$','Interpreter','latex');
    y_v_psi = ylabel('$us_{psi}$','Interpreter','latex');
    set(y_v_psi,'FontSize',20);

end

%% Plot do Modelo Nominal no tempo
figure
plot(t_v, F0_v(1,:), 'g');
%title('F0 no tempo')
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_x = ylabel('$F0_x$','Interpreter','latex');
set(y_F0_x,'FontSize',14);

figure
plot(t_v, F0_v(2,:), 'r');
%title('F0 no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_y = ylabel('$F0_y$','Interpreter','latex');
set(y_F0_y,'FontSize',14);

figure
plot(t_v, F0_v(3,:), 'b');
%title('F0 no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_z = ylabel('$F0_z$','Interpreter','latex');
set(y_F0_z,'FontSize',14);

figure
plot(t_v, F0_v(4,:), 'm');
%title('F0 no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_psi = ylabel('$F0_{psi}$','Interpreter','latex');
set(y_F0_psi,'FontSize',14);

%% Plot do Delta F0 no tempo
figure
plot(t_v, deltaF0_v(1,:), 'g');
%title('$\Delta F0_x$ Estimado no tempo','Interpreter','latex')
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_x = ylabel('$\Delta F0_x$','Interpreter','latex');
set(y_F0_x,'FontSize',14);

figure
plot(t_v, deltaF0_v(2,:), 'r');
%title('$\Delta F0_y$ Estimado no tempo','Interpreter','latex');
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_y = ylabel('$\Delta F0_y$','Interpreter','latex');
set(y_F0_y,'FontSize',14);

figure
plot(t_v, deltaF0_v(3,:), 'b');
%title('$\Delta F0_z$ Estimado no tempo','Interpreter','latex');
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_z = ylabel('$\Delta F0_z$','Interpreter','latex');
set(y_F0_z,'FontSize',14);

figure
plot(t_v, deltaF0_v(4,:), 'm');
%title('$\Delta F0_{psi}$ Estimado no tempo','Interpreter','latex');
xlabel('$Tempo(s)$','Interpreter','latex');
y_F0_psi = ylabel('$\Delta F0_{psi}$','Interpreter','latex');
set(y_F0_psi,'FontSize',14);

%% Posição real vs posição desejada
figure
plot(t_v, q_v(1,:), 'g', t_v, qd_v(1,:), 'r-');
%title('Pos Real vs Pos Desejada');
legend('real', 'desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_xtil = ylabel('$x(m)$','Interpreter','latex');
set(y_xtil,'FontSize',17);

figure
plot(t_v, q_v(2,:), 'g', t_v, qd_v(2,:), 'r-');
%title('Pos Real vs Pos Desejada');
legend('real', 'desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_ytil = ylabel('$y(m)$','Interpreter','latex');
set(y_ytil,'FontSize',17);

figure
plot(t_v, q_v(3,:), 'g', t_v, qd_v(3,:), 'r-');
%title('Pos Real vs Pos Desejada');
legend('real', 'desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_ztil = ylabel('$z(m)$','Interpreter','latex');
set(y_ztil,'FontSize',17);

figure
plot(t_v, q_v(4,:), 'g', t_v, qd_v(4,:), 'r-');
%title('Pos Real vs Pos Desejada');
legend('real', 'desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$\psi(rad)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

%% Velocidade real vs Velocidade desejada
figure
plot(t_v, dq_v(1,:), 'g', t_v, dqd_v(1,:), 'r-');
%title('Vel Real vs Vel Desejada');
legend('vel real', ' vel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_xtil = ylabel('$x(m/s)$','Interpreter','latex');
set(y_xtil,'FontSize',17);

figure
plot(t_v, dq_v(2,:), 'g', t_v, dqd_v(2,:), 'r-');
%title('Vel Real vs Vel Desejada');
legend('vel real', 'vel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_ytil = ylabel('$y(m/s)$','Interpreter','latex');
set(y_ytil,'FontSize',17);

figure
plot(t_v, dq_v(3,:), 'g', t_v, dqd_v(3,:), 'r-');
%title('Vel Real vs Vel Desejada');
legend('vel real', 'vel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_ztil = ylabel('$z(m/s)$','Interpreter','latex');
set(y_ztil,'FontSize',17);

figure
plot(t_v, dq_v(4,:), 'g', t_v, dqd_v(4,:), 'r-');
%title('Vel Real vs Vel Desejada');
legend('vel real', 'vel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$\psi(rad/s)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

%% Aceleração real vs Aceleração desejada
figure
plot(t_v, d2q_v(1,:), 'g', t_v, d2qd_v(1,:), 'r-');
%title('Acel Real vs Acel Desejada');
legend('acel real', ' acel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_xtil = ylabel('$x(m/s^2)$','Interpreter','latex');
set(y_xtil,'FontSize',17);

figure
plot(t_v, d2q_v(2,:), 'g', t_v, d2qd_v(2,:), 'r-');
%title('Acel Real vs Acel Desejada');
legend('acel real', 'acel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_ytil = ylabel('$y(m/s^2)$','Interpreter','latex');
set(y_ytil,'FontSize',17);

figure
plot(t_v, d2q_v(3,:), 'g', t_v, d2qd_v(3,:), 'r-');
%title('Acel Real vs Acel Desejada');
legend('acel real', 'acel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_ztil = ylabel('$z(m/s^2)$','Interpreter','latex');
set(y_ztil,'FontSize',17);

figure
plot(t_v, d2q_v(4,:), 'g', t_v, d2qd_v(4,:), 'r-');
%title('Acel Real vs Acel Desejada');
legend('acel real', 'acel desejada');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$\psi(rad/s^2)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

%% Disturbio
figure
plot(t_v, d_v(1,:), 'b');
%title('Distúrbio no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d(N)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

figure
plot(t_v, d_v(2,:), 'b');
%title('Distúrbio no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d(N)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

figure
plot(t_v, d_v(3,:), 'b');
%title('Distúrbio no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d(N)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

figure
plot(t_v, d_v(4,:), 'b');
%title('Distúrbio no tempo');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d(N)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

%% Aceleração do Distúrbio
figure
plot(t_v, dB_v(1,:), 'r');
%title('Aceleração do Distúrbio no tempo x');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d*B(m/s^2)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

figure
plot(t_v, dB_v(2,:), 'r');
%title('Aceleração do Distúrbio no tempo y');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d*B(m/s^2)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

figure
plot(t_v, dB_v(3,:), 'r');
%title('Aceleração do Distúrbio no tempo z');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d*B(m/s^2)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

figure
plot(t_v, dB_v(4,:), 'r');
%title('Aceleração do Distúrbio no tempo \Psi');
xlabel('$Tempo(s)$','Interpreter','latex');
y_psitil = ylabel('$d*B(m/s^2)$','Interpreter','latex');
set(y_psitil,'FontSize',17);

end
