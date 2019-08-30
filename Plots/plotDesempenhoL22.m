function [] = plotDesempenhoL22(control,EPM, EPsiM,u_norm_media, t_v)
    
if(control == 1)
    figure(2)
    plotEPM = semilogy(t_v,(EPM(control,:)),'g');
    set (plotEPM, 'LineWidth', 2)
    %title('norma L2 do vetor de erros de posição');
    xl = xlabel('$Tempo(s)$','Interpreter','latex');
    set(xl, 'FontSize', 25);
    yl = ylabel('$\varepsilon(m)$','Interpreter','latex');
    set(yl, 'FontSize', 25);
    set(gca, 'FontSize',25);
    grid on
    hold on
    
    figure(3)
    plotEPsiM = semilogy(t_v,(EPsiM(control,:)),'g');
    set (plotEPsiM, 'LineWidth', 2)
    %title('norma L2 do vetor de erros de orientação');
    xl = xlabel('$Tempo(s)$','Interpreter','latex');
    set(xl, 'FontSize', 25);
    yl = ylabel('$\varepsilon_\psi(rad)$','Interpreter','latex');
    set(yl, 'FontSize', 25);
    set(gca, 'FontSize',25);
    grid on
    hold on
    
    figure(4)
    plotu_norm_media = semilogy(t_v,(u_norm_media(control,:)),'g');
    set (plotu_norm_media, 'LineWidth', 2)
    %title('norma L2 do vetor de controle');
    xl = xlabel('$Tempo(s)$','Interpreter','latex');
    set(xl, 'FontSize', 25);
    yl = ylabel('$u$','Interpreter','latex');
    set(yl, 'FontSize', 25);
    set(gca, 'FontSize',25);
    grid on
    hold on
    
 elseif(control == 2)
     figure(2)
     cordalinha = 'r';
     plotEPM = semilogy(t_v,(EPM(control,:)),cordalinha);
     set (plotEPM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
     figure(3)
     cordalinha = 'r';
     plotEPsiM = semilogy(t_v,(EPsiM(control,:)),cordalinha);
     set (plotEPsiM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
     figure(4)
     cordalinha = 'r';
     plotu_norm_media = semilogy(t_v,(u_norm_media(control,:)),cordalinha);
     set (plotu_norm_media, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
 elseif(control == 3)
     figure(2)
     cordalinha = 'c';
     plotEPM = semilogy(t_v,(EPM(control,:)),cordalinha);
     set (plotEPM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
     figure(3)
     cordalinha = 'c';
     plotEPsiM = semilogy(t_v,(EPsiM(control,:)),cordalinha);
     set (plotEPsiM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
      figure(4)
     cordalinha = 'c';
     plotu_norm_media = semilogy(t_v,(u_norm_media(control,:)),cordalinha);
     set (plotu_norm_media, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
 elseif(control == 4)
     figure(2)
     cordalinha = 'm';
     plotEPM = semilogy(t_v,(EPM(control,:)),cordalinha);
     set (plotEPM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
     figure(3)
     cordalinha = 'm';
     plotEPsiM = semilogy(t_v,(EPsiM(control,:)),cordalinha);
     set (plotEPsiM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
     figure(4)
     cordalinha = 'm';
     plotu_norm_media = semilogy(t_v,(u_norm_media(control,:)),cordalinha);
     set (plotu_norm_media, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
 elseif(control == 5)
     figure(2)
     cordalinha = 'y';
     plotEPM = semilogy(t_v,(EPM(control,:)),cordalinha);
     set (plotEPM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
     figure(3)
     cordalinha = 'y';
     plotEPsiM = semilogy(t_v,(EPsiM(control,:)),cordalinha);
     set (plotEPsiM, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
      figure(4)
     cordalinha = 'y';
     plotu_norm_media = semilogy(t_v,(u_norm_media(control,:)),cordalinha);
     set (plotu_norm_media, 'LineWidth', 2);
     set(gca, 'FontSize',25);
     %hold on
     
%  else
%      cordalinha = 'k';
%      traj_real = plot3(q_v(1,:),q_v(2,:),q_v(3,:),cordalinha);
%      set (traj_real, 'LineWidth', 2);
end
figure (2)
legend({'FL', 'H\infty', 'H\infty RN',  'H\infty GP THETA','H\infty GP OFF'}, 'FontSize',17, 'Location','northeast');


figure (3)
legend({'FL', 'H\infty', 'H\infty RN',  'H\infty GP THETA','H\infty GP OFF'}, 'FontSize',17, 'Location','northeast');

figure (4)
legend({'FL', 'H\infty', 'H\infty RN',  'H\infty GP THETA','H\infty GP OFF'}, 'FontSize',17, 'Location','northeast');

