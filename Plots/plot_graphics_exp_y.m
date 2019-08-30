function [] = plot_graphics_exp_y(t_v, q_v, qd_v, control)

if control == 1
    cordalinha = 'g';    
    figure (52)
    plot(t_v, qd_v(2,:),'b--', 'LineWidth', 2);
    x_x = xlabel('$Tempo(s)$','Interpreter','latex');
    set(x_x, 'FontSize', 25);
    y_x = ylabel('$y(m)$','Interpreter','latex');
    set(y_x, 'FontSize', 25);
    grid on
    hold on

    plot(t_v, q_v(2,:),cordalinha, 'LineWidth', 2);
    set(gca, 'FontSize',25);
    hold on
elseif control == 2
    cordalinha = 'r';
    figure (52)
    plot(t_v, q_v(2,:),cordalinha, 'LineWidth', 2);
    set(gca, 'FontSize',25);
    hold on
elseif control == 3
    cordalinha = 'c';
    figure (52)
    plot(t_v, q_v(2,:),cordalinha, 'LineWidth', 2);
    set(gca, 'FontSize',25);
    hold on
elseif control == 4
    cordalinha = 'm';
    figure (52)
    plot(t_v, q_v(2,:),cordalinha, 'LineWidth', 2);
    set(gca, 'FontSize',25);
    hold on
else
    cordalinha = 'y';
    figure (52)
    plot(t_v, q_v(2,:),cordalinha, 'LineWidth', 2);
    set(gca, 'FontSize',25);
end
    
    
%% Plot da posição no tempo

% 
% figure (51)
% plot(t_v, q_v(1,:),cordalinha, 'LineWidth', 2);
% %title('$Grafico posicao X$','Interpreter','latex');
% %title('Gráfico posição X');
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_x = ylabel('$x(m)$','Interpreter','latex');
% set(y_x, 'FontSize', 25);
% grid on
% hold on
% 
% figure (52)
% plot(t_v, qd_v(2,:),'b--', 'LineWidth', 2);
% hold on
% 
% figure (52)
% plot(t_v, q_v(2,:),cordalinha, 'LineWidth', 2);
% %title('Gráfico posição Y');
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_y = ylabel('$y(m)$','Interpreter','latex');
% set(y_y, 'FontSize', 25);
% grid on
% hold on
% 
% 
% figure (53)
% plot(t_v, qd_v(3,:),'b--', 'LineWidth', 2);
% hold on
% 
% figure (53)
% plot(t_v, q_v(3,:),cordalinha, 'LineWidth', 2);
% %title('Gráfico posição Z');
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_z = ylabel('$z(m)$','Interpreter','latex');
% set(y_z, 'FontSize', 25);
% grid on
% hold on
% 
% figure (54)
% plot(t_v, qd_v(4,:),'b--', 'LineWidth', 2);
% hold on
% 
% figure (54)
% plot(t_v, q_v(4,:),cordalinha, 'LineWidth', 2);
% %title('Gráfico posição PSI');
% xlabel('$Tempo(s)$','Interpreter','latex');
% y_psi = ylabel('$\psi (rad)$','Interpreter','latex');
% set(y_psi, 'FontSize', 25);
% grid on
% hold on
% 
% 
% 
% end
% 
% figure (51)
% legend('Referência', 'FL', 'H\infty', 'H\infty RN', 'H\infty GP THETA', 'H\infty GP OFF','Location','northeast');


figure (52)
legend({'Referência', 'FL', 'H\infty', 'H\infty RN',  'H\infty GP THETA','H\infty GP OFF'}, 'FontSize',17, 'Location','northeast');

% figure (53)
% legend('Referência', 'FL', 'H\infty', 'H\infty RN', 'H\infty GP THETA', 'H\infty GP OFF','Location','northeast');
% 
% figure (54)
% legend('Referência', 'FL', 'H\infty', 'H\infty RN', 'H\infty GP THETA', 'H\infty GP OFF','Location','northeast');
% 
