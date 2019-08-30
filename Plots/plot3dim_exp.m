function [] = plot3dim_exp(qd_v,q_v,control)

%     figure
%     traj_desejada = plot3(qd_v(1,:),qd_v(2,:),qd_v(3,:),'b--');
%     set (traj_desejada, 'LineWidth', 2)
%     title('Gráfico trajeto real vs desejado');
%     xl = xlabel('$x(m)$','Interpreter','latex');
%     yl = ylabel('$y(m)$','Interpreter','latex');
%     zl = zlabel('$z(m)$','Interpreter','latex');
%     set(xl, 'FontSize', 25);
%     set(yl, 'FontSize', 25);
%     set(zl, 'FontSize', 25);
%     grid on
%     hold on
    
 if(control == 1)
     
    figure(1)
    traj_desejada = plot3(qd_v(1,:),qd_v(2,:),qd_v(3,:),'b--');
    set (traj_desejada, 'LineWidth', 2)
    set(gca, 'FontSize',25);
    %title('Gráfico trajeto real vs desejado');
    xl = xlabel('$x(m)$','Interpreter','latex');
    yl = ylabel('$y(m)$','Interpreter','latex');
    zl = zlabel('$z(m)$','Interpreter','latex');
    set(xl, 'FontSize', 25);
    set(yl, 'FontSize', 25);
    set(zl, 'FontSize', 25);
    grid on
    hold on
    
    cordalinha = 'g';
    traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
    set (traj_real, 'LineWidth', 2);
    hold on
    
 elseif(control == 2)
     cordalinha = 'r';
     figure(1)
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     hold on
 elseif(control == 3)
     cordalinha = 'c';
     figure(1)
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     hold on
 elseif(control == 4)
     cordalinha = 'm';
     figure(1)
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     hold on
 elseif(control == 5)
     cordalinha = 'y';
     figure(1)
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     hold on
%  else
%      cordalinha = 'k';
%      traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
%      set (traj_real, 'LineWidth', 2);
 end
legend({'Referência', 'FL', 'H\infty', 'H\infty RN',  'H\infty GP THETA','H\infty GP OFF'}, 'FontSize',16, 'Location','northeast');


