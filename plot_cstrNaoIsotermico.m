%Estados / Vari�veis Controladas
% x1 = h - n�vel dentro do tanque; 
% x2 = Ca - Concentra��o de sa�da do produto A;
% x3 = T - Temperatura dentro do reator; 
% x4 = R - Velocidade de rea��o?
tamLetra = 10;
tamTitulo = 12;
espes = 3;
saidas=x;

h = figure();
h.WindowState = 'maximized';

subplot(3,2,1);
plot(saidas(1,1:iteracoes),'r','linewidth',espes);
hold on
plot(x_pred_vect(1,1:iteracoes),'k--','linewidth',espes);
if controle==1 
    plot(ref_h(1:iteracoes),'b--','linewidth',espes);
    legend({'Real','Predi��o','Refer�ncia'},'FontSize',tamLetra);
else
   legend({'Real','Predi��o'},'FontSize',tamLetra); 
end
xlim([0 iteracoes])
% ylim([0.92 1.2])
xlabel('Itera��es (-)','FontSize',tamLetra);
ylabel('Altura (m)','FontSize',tamLetra);
title('$h$ - Altura do tanque','interpreter','latex','FontSize',tamTitulo)
grid


subplot(3,2,3);
plot(saidas(2,1:iteracoes),'r','linewidth',espes);
hold on
plot(x_pred_vect(2,1:iteracoes),'k--','linewidth',espes)%,'marker','x');
if controle==1 
    plot(ref_ca(1:iteracoes),'b--','linewidth',espes);
    legend({'Real','Predi��o','Refer�ncia'},'FontSize',tamLetra);
else
   legend({'Real','Predi��o'},'FontSize',tamLetra); 
end
xlim([0 iteracoes])
% ylim([0.98 1.08])
xlabel('Itera��es (-)','FontSize',tamLetra);
ylabel('Concentra��o (kmol/m^3)','FontSize',tamLetra);
title('$C_a$ - Concentra\c{c}\~{a}o do produto A','interpreter','latex','FontSize',tamTitulo);
grid


subplot(3,2,5);
plot(saidas(3,1:iteracoes),'r','linewidth',espes);
hold on
plot(x_pred_vect(3,1:iteracoes),'k--','linewidth',espes);
if controle==1 
    plot(ref_T(1:iteracoes),'b--','linewidth',espes);
    legend({'Real','Predi��o','Refer�ncia'},'FontSize',tamLetra);
else
   legend({'Real','Predi��o'},'FontSize',tamLetra); 
end
xlim([0 iteracoes])
% ylim([390 410])
xlabel('Itera��es (-)','FontSize',tamLetra);
ylabel('Temperatura Interna (K)','FontSize',tamLetra);
title('$T$ - Temperatura dentro do tanque','interpreter','latex','FontSize',tamTitulo);
grid


% figure
%----- Entradas / Vari�veis Manipuladas -----
% u1 = qo - Vaz�o de sa�da;
% u2 = Caf - Concentra��o do produto A na alimenta��o do tanque; 
% u3 = Qh/pCp - taxa de remo��o de calor normalizada;
subplot(3,2,2);
% plot(entradas_atrasadas_vect(1,1:iteracoes),'k--','linewidth',2);
% hold on
ax=gca;
xlim([0 iteracoes])
ylim([4.8e-3 5.2e-3])
ax.YAxis.Exponent = -3;
% plot(entradas_comp_vect(1,1:iteracoes),'b','linewidth',2);
plot(entradas(1,1:iteracoes),'r','linewidth',espes);
title('Entradas','FontSize',tamTitulo)
legend({'$q_0$ - Vaz\~{a}o de sa\''ida'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlabel('Itera��es (-)','FontSize',tamLetra);
ylabel('Vaz�o (m^3/s^{-1})','FontSize',tamLetra)
grid

subplot(3,2,4);
% plot(entradas_atrasadas_vect(2,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(2,1:iteracoes),'b','linewidth',2);
plot(entradas(2,1:iteracoes),'r','linewidth',espes);
legend({'$C_{af}$ - Concentra\c{c}\~{a}o produto A na alimenta\c{c}\~{a}o do tanque'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlim([0 iteracoes])
ylim([4.9 5.2])
xlabel('Itera��es (-)','FontSize',tamLetra);
ylabel('Concentra��o (kmol/m^3)','FontSize',tamLetra)
grid

subplot(3,2,6);
% plot(entradas_atrasadas_vect(3,1:iteracoes),'k--','linewidth',2);
% hold on
% plot(entradas_comp_vect(3,1:iteracoes),'b','linewidth',2);
plot(entradas(3,1:iteracoes),'r','linewidth',espes);
legend({'${Qh}/{pc_p}$ - Taxa de remo\c{c}\~{a}o de calor normalizada'},'interpreter','latex','Location','best','FontSize',tamLetra);
xlim([0 iteracoes])
ylim([0.73 0.8])
xlabel('Itera��es (-)','FontSize',tamLetra);
ylabel({'Taxa de remo��o de'; 'calor (Km^3/s^{-1})'},'FontSize',tamLetra)
grid

sgtitle('FSP N�o-Linear')

% h2 = figure();
% h2.WindowState = 'maximized';
% %----- Perturba��es ------
% % q = qi - vazao de entrada; 
% % Ti = Ti - Temeratura externa;
% subplot(2,1,1)
% plot(perturbacoes(1,1:iteracoes),'r','linewidth',espes);
% hold on
% plot(x_pred_vect(4,1:iteracoes),'k--','linewidth',espes);
% title('$q_i$ - Vaz\~{a}o de entrada','interpreter','latex','FontSize',tamTitulo);
% legend({'Real','Predi��o'},'Location','best','FontSize',tamLetra);
% ax=gca;
% xlim([0 iteracoes])
% if max(perturbacoes(1,:))>5e-3
%    ylim([4.9e-3 5.2e-3])
% else
%    ylim([4.8e-3 5.1e-3])
% end
% ax.YAxis.Exponent = -3;
% ylabel('Vaz�o (m^3/s^{-1})','FontSize',tamLetra)
% xlabel('Itera��es (-)','FontSize',tamLetra);
% grid
% 
% 
% subplot(2,1,2);
% plot(perturbacoes(2,1:iteracoes),'r','linewidth',espes);
% hold on
% plot(x_pred_vect(5,1:iteracoes),'k--','linewidth',espes);
% title('$T_i$ - Temperatura externa','interpreter','latex','FontSize',tamTitulo)
% legend({'Real','Predi��o'},'Location','best','FontSize',tamLetra);
% ax=gca;
% xlim([0 iteracoes])
% if max(perturbacoes(2,:))>350
%     ylim([348 360])
% else
%     ylim([340 352])
% end
% ylabel('Temperatura Externa (K)','FontSize',tamLetra)
% xlabel('Itera��es (-)','FontSize',tamLetra);
% grid
