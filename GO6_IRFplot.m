%--------------------------------------------------------------------------%
% This codes plots some of the results from GO5_IRF                        %
%                                                                          %
%                                                                          %
%--------------------------------------------------------------------------%
clc
disp('Plot IRFs...')

if saveplt

%% 1. IRF TO FACTOR f AND g: ALL VOLs and GDPs + WEIGHTED
%--------------------------------------------------------------------------
FigSize(20,16)
subplot(2,2,1); plot(IR.IRf_VOL'); hold on; h2 = plot(IRw.IRf_wVOL','Color','k','LineWidth',3);  grid on; set(gca,'GridLineStyle',':');...
    title('(A) Volatility response to $\zeta_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top');
subplot(2,2,2); plot(IR.IRf_GDP'); hold on; h2 = plot(IRw.IRf_wGDP','Color','k','LineWidth',3);  grid on; set(gca,'GridLineStyle',':');...
    title('(B) GDP growth response to $\zeta_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top');
subplot(2,2,3); plot(IR.IRg_VOL'); hold on; h2 = plot(IRw.IRg_wVOL','Color','k','LineWidth',3);  grid on; set(gca,'GridLineStyle',':');...
    title('(C) Volatility response to $\xi_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top'); xlabel('Quarters','Fontsize',9)
subplot(2,2,4); plot(IR.IRg_GDP'); hold on; h2 = plot(IRw.IRg_wGDP','Color','k','LineWidth',3);  grid on; set(gca,'GridLineStyle',':');...
    title('(D) GDP growth response to $\xi_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top');
xlabel('Quarters')
SaveFigure('Output\IR_FACTOR',0)
clf('reset')

%% 2. IRF TO FACTOR f AND g: SWATHE
FigSize(20,8)
subplot(1,2,1); PlotSwathe(IRw.IRf_wGDP,[IRw.IRf_sup_wGDP; IRw.IRf_inf_wGDP]);  grid on; set(gca,'GridLineStyle',':'); ...
    title('(A) GDP growth response to $\zeta_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top');
subplot(1,2,2); PlotSwathe(IRw.IRf_wVOL,[IRw.IRf_sup_wVOL; IRw.IRf_inf_wVOL]);  grid on; set(gca,'GridLineStyle',':'); ...
    title('(B) Volatility response to $\zeta_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top');
SaveFigure('Output\IRmg_FACTOR_A',0)
clf('reset')

subplot(1,2,1); PlotSwathe(IRw.IRg_wVOL,[IRw.IRg_sup_wVOL; IRw.IRg_inf_wVOL]);  grid on; set(gca,'GridLineStyle',':'); ...
    title('(A) Volatility response to $\xi_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top'); xlabel('Quarters','Fontsize',9)
subplot(1,2,2); PlotSwathe(IRw.IRg_wGDP,[IRw.IRg_sup_wGDP; IRw.IRg_inf_wGDP]);  grid on; set(gca,'GridLineStyle',':'); ...
    title('(B) GDP growth response to $\xi_t$','Interpreter','Latex'); axis tight; set(gca,'Layer','top');
xlabel('Quarters','Fontsize',9)
SaveFigure('Output\IRmg_FACTOR_B',0)
clf('reset')


%% 3. COMPARE GLOBAL, COUNTRY-SPECIFIC, BIVARIATE (SWATHE)
FigSize(20,16)
clear h1 h2
subplot(2,2,1); 
    h1 = PlotSwathe(IRw.IRn_wGDP,[IRw.IRn_sup_wGDP; IRw.IRn_inf_wGDP]); hold on;
    h2 = plot(IR.IRn_GDP(32,:),'LineStyle','-','Marker','*','Color',cmap(2),'LineWidth',1); hold on
    title({'\textbf{(A) GDP growth response to $\eta_{it}$}';'\hspace{1cm}(Multi-country model)'},'Interpreter','Latex'); 
    grid on; set(gca,'GridLineStyle',':'); axis tight; set(gca,'Layer','top'); ylim([-0.3 0]);
    legend([h1.bar h2],{'Average economy';'United States'},'Interpreter','Latex','Location','SouthEast')
subplot(2,2,2); 
    h1 = PlotSwathe(IRw.IRbiv_wGDP,[IRw.IRbiv_sup_wGDP; IRw.IRbiv_inf_wGDP],rgb('very dark green')); hold on;
    h2 = plot(IR.IRbiv_GDP(32,:),'LineStyle','-','Marker','*','Color',cmap(2),'LineWidth',1); hold on
    title({'\textbf{(B) GDP growth response to $\eta_{it}$}';'\hspace{1cm}(Single-country model)'},'Interpreter','Latex'); 
    grid on; set(gca,'GridLineStyle',':'); axis tight; set(gca,'Layer','top'); ylim([-0.3 0]);
    legend([h1.bar h2],{'Average economy';'United States'},'Interpreter','Latex','Location','SouthEast')
subplot(2,2,3); 
    PlotSwathe(IRw.IRn_wVOL,[IRw.IRn_sup_wVOL; IRw.IRn_inf_wVOL]); hold on; 
    plot(IR.IRn_VOL(32,:),'LineStyle','-','Marker','*','Color',cmap(2),'LineWidth',1); hold on
    title({'\textbf{(C) Volatility response to $\eta_{it}$}';'\hspace{1cm}(Multi-country model)'},'Interpreter','Latex'); 
    grid on; set(gca,'GridLineStyle',':'); axis tight; set(gca,'Layer','top'); ylim([0 0.3]);
subplot(2,2,4); 
    PlotSwathe(IRw.IRbiv_wVOL,[IRw.IRbiv_sup_wVOL; IRw.IRbiv_inf_wVOL],rgb('very dark green')); hold on;
    plot(IR.IRbiv_VOL(32,:),'LineStyle','-','Marker','*','Color',cmap(2),'LineWidth',1); hold on
    title({'\textbf{(D) Volatility response to $\eta_{it}$}';'\hspace{1cm}(Single-country model)'},'Interpreter','Latex'); 
    grid on; set(gca,'GridLineStyle',':'); axis tight; set(gca,'Layer','top'); ylim([0 0.3]);
SaveFigure('Output\IRbiv1',0)
clf('reset')
            
end

close all
disp('Done!')

