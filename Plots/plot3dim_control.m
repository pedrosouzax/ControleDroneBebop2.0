function [] = plot3dim_control(qd_v,q_v,control)

    figure
    traj_desejada = plot3(qd_v(1,:),qd_v(2,:),qd_v(3,:),'b--');
    set (traj_desejada, 'LineWidth', 2)
    title('Gráfico trajeto real vs desejado');
    xl = xlabel('$x(m)$','Interpreter','latex');
    yl = ylabel('$y(m)$','Interpreter','latex');
    zl = zlabel('$z(m)$','Interpreter','latex');
    set(xl, 'FontSize', 17);
    set(yl, 'FontSize', 17);
    set(zl, 'FontSize', 17);
    grid on
    hold on
    
 if(control == 1)
     
%     figure
%     traj_desejada = plot3(qd_v(1,:),qd_v(2,:),qd_v(3,:),'b--');
%     set (traj_desejada, 'LineWidth', 2)
%     %title('Gráfico trajeto real vs desejado');
%     xl = xlabel('$x(m)$','Interpreter','latex');
%     yl = ylabel('$y(m)$','Interpreter','latex');
%     zl = zlabel('$z(m)$','Interpreter','latex');
%     set(xl, 'FontSize', 17);
%     set(yl, 'FontSize', 17);
%     set(zl, 'FontSize', 17);
%     grid on
%     hold on
    
    cordalinha = 'g';
    traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
    set (traj_real, 'LineWidth', 2);
    %hold on
    legend('Referência', 'FL','Location','northeast');


    
 elseif(control == 2)
     cordalinha = 'r';
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     %hold on
     legend('Referência', 'H\infty','Location','northeast');

 elseif(control == 3)
     cordalinha = 'c';
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     %hold on
     legend('Referência','H\infty RN','Location','northeast');


 elseif(control == 4)
     cordalinha = 'm';
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     %hold on
     legend('Referência','H\infty GP THETA','Location','northeast');
     

 elseif(control == 5)
     cordalinha = 'y';
     traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
     set (traj_real, 'LineWidth', 2);
     %hold on
     legend('Referência','H\infty GP OFF','Location','northeast');

     


%  else
%      cordalinha = 'k';
%      traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
%      set (traj_real, 'LineWidth', 2);
 end
%legend('Referência', 'FL', 'H\infty', 'H\infty RN', 'H\infty GP OFF', 'H\infty GP THETA','Location','northeast');


