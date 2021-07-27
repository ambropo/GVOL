%--------------------------------------------------------------------------%
% This codes plots some of the results from GO7_FEVD                       %
%                                                                          %
%                                                                          %
%--------------------------------------------------------------------------%
clc
disp('Plot FEVDs (block diagonal VCV with Cholesky)...')

if saveplt

%% 1. FEVD -- WEIGHTED AVERAGE OF VOLs and GDPs decompositions
%==========================================================================
FigSize(20,8)
% subplot(1,3,1:2)
hplot = bar(100.*[v2g_w' v2eta_i_w' v2eta_other_w' v2f_w' v2eps_i_w' v2eps_other_w'],'stacked');
title('Volatility ($v_{it}$), average','Interpreter','Latex'); axis tight; %ylim([0 100])
hcol = colormap(lines(6)); 
    for ii=1:6; hplot(ii).FaceColor = hcol(ii,:); end
    for ii=1:6; hplot(ii).EdgeColor = 'w'; end
    for ii=1:6; hplot(ii).LineWidth = 0.001; end
    hatchfill2(hplot(1),'single','HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(2),'cross' ,'HatchAngle',45 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(3),'single','HatchAngle',180,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(4),'single','HatchAngle',135,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(5),'cross' ,'HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    set(gca,'Layer','top');
% Legend hatch
    [~,hleg,~,~] = legendflex(hplot,{'Global financial $\hat{\xi}$','Volatility $\hat{\eta}_i$','Volatility (other) $\sum \hat{\eta}_j$',...
        'Global growth $\hat{\zeta}$','GDP growth $\hat{\varepsilon}_i$','GDP growth (other) $\sum \hat{\varepsilon}_j$'},...
        'anchor',[5 5],'buffer',[-10 10],'fontsize',9,'interpreter','Latex');
    hatchfill2(hleg(7) ,'single','HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(8) ,'cross' ,'HatchAngle',45 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(9) ,'single','HatchAngle',180,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(10),'single','HatchAngle',135,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(11),'cross' ,'HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
SaveFigure('Output\FEVD_BDChol_WEIGHT_A',0)
clf('reset')

hplot = bar(100.*[y2g_w' y2eta_i_w' y2eta_other_w' y2f_w' y2eps_i_w' y2eps_other_w'],'stacked');
title('GDP Growth ($\Delta y_{it}$), average','Interpreter','Latex'); axis tight; %ylim([0 100])
    for ii=1:6; hplot(ii).FaceColor = hcol(ii,:); end
    for ii=1:6; hplot(ii).EdgeColor = 'w'; end
    for ii=1:6; hplot(ii).LineWidth = 0.001; end
    hatchfill2(hplot(1),'single','HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(2),'cross' ,'HatchAngle',45 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(3),'single','HatchAngle',180,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(4),'single','HatchAngle',135,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hplot(5),'cross' ,'HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    set(gca,'Layer','top');
% Legend hatch
    [~,hleg,~,~] = legendflex(hplot,{'Global financial $\hat{\xi}$','Volatility $\hat{\eta}_i$','Volatility (other) $\sum \hat{\eta}_j$',...
        'Global growth $\hat{\zeta}$','GDP growth $\hat{\varepsilon}_i$','GDP growth (other) $\sum \hat{\varepsilon}_j$'},...
        'anchor',[3 3],'buffer',[-10 -10],'fontsize',9,'interpreter','Latex');
    hatchfill2(hleg(7) ,'single','HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(8) ,'cross' ,'HatchAngle',45 ,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(9) ,'single','HatchAngle',180,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(10),'single','HatchAngle',135,'HatchLineWidth',.5,'HatchLineStyle','-');
    hatchfill2(hleg(11),'cross' ,'HatchAngle',90 ,'HatchLineWidth',.5,'HatchLineStyle','-');
SaveFigure('Output\FEVD_BDChol_WEIGHT_B',0)
clf('reset')

%% 2. FEVD - ALL COUNTRIES
%==========================================================================
% GFEVD VOL
FigSize
for jj=1:ncountry
    subplot(4,8,jj)
    hplot = bar(100.*[v2g(jj,:); v2eta_i(jj,:); v2eta_other(jj,:); v2f(jj,:); v2eps_i(jj,:); v2eps_other(jj,:)]','stacked');
    hcol = colormap(lines(6)); 
    for ii=1:6; hplot(ii).FaceColor = hcol(ii,:); end
    for ii=1:6; hplot(ii).EdgeColor = 'w'; end
    for ii=1:6; hplot(ii).LineWidth = 0.001; end
    set(gca,'Layer','top');
    title(['VOL - ' cnames{jj}])
    axis tight
end
legopt=LegOption; legopt.interpreter='Latex';
LegSubplot({'$\hat{\xi}$','$\hat{\eta}_i$','$\sum \hat{\eta}_j$','$\hat{\zeta}$','$\hat{\varepsilon}_i$','$\sum \hat{\varepsilon}_j$'},legopt)
SaveFigure('Output\FEVD_BDChol_VOL',0)
clf('reset')

% GFEVD GDP
FigSize
for jj=1:ncountry
    subplot(4,8,jj)
    hplot = bar(100.*[y2g(jj,:); y2eta_i(jj,:); y2eta_other(jj,:); y2f(jj,:); y2eps_i(jj,:); y2eps_other(jj,:)]','stacked');
    hcol = colormap(lines(6)); 
    for ii=1:6; hplot(ii).FaceColor = hcol(ii,:); end
    for ii=1:6; hplot(ii).EdgeColor = 'w'; end
    for ii=1:6; hplot(ii).LineWidth = 0.001; end
    set(gca,'Layer','top');
    title(['GDP - ' cnames{jj}])
    axis tight
end
legopt=LegOption; legopt.interpreter='Latex';
LegSubplot({'$\hat{\xi}$','$\hat{\eta}_i$','$\sum \hat{\eta}_j$','$\hat{\zeta}$','$\hat{\varepsilon}_i$','$\sum \hat{\varepsilon}_j$'},legopt)
SaveFigure('Output\FEVD_BDChol_GDP',0)
clf('reset')

end

close all
disp('Done!')

