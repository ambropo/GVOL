%--------------------------------------------------------------------------%
% This codes plots some of the results from GO2_REGS                       %
%                                                                          %
%                                                                          %
%--------------------------------------------------------------------------%
clc
disp('Run FIGURES...')

if saveplt

%% 1. Unconditional
%==========================================================================
% Boxplot
FigSize(20,7)
errorbar(1:ncountry,UNC.corrvy(:,1),UNC.corrvy(:,2),'o','Color',cmap(1),'Linewidth',1.2,'MarkerFaceColor',cmap(2))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long(UNC.corrvy(:,3)),'yLim',[-0.8 0.5],'XTickLabelRotation',45,'Layer','top');
hold on; plot(zeros(ncountry,1),'k:','Linewidth',0.5)
ylabel('Correlation');
SaveFigure('Output\UNC_BOXvy',0)
clf('reset')
% Pairwise correlation: VOL only 
FigSize(20,7)
h = bar(UNC.PCv_num(:,1),'BarWidth',0.4); 
ylimm = [-0.2 .75];
set(h(1),'FaceColor',cmap(2),'Edgecolor',cmap(2)); hold on; 
bar1(1:ncountry) = mean(UNC.PCv_num(:,1)); plot(bar1,':','LineWidth',2,'Color',cmap(2))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long,'yLim',ylimm,'XTickLabelRotation',45,'Layer','top');
ylabel('Correlation');
SaveFigure('Output/PC_V',0)
clf('reset')
% Pairwise correlation: VOL and GDP
FigSize(20,7)
h = bar([UNC.PCv_num(:,1) UNC.PCy_num(:,2)],'BarWidth',1.2); hold on;
bar1(1:ncountry) = mean(UNC.PCv_num(:,1)); plot(bar1,':','LineWidth',2,'Color',cmap(2)); hold on
bar2(1:ncountry) = mean(UNC.PCy_num(:,2)); plot(bar2,':','LineWidth',2,'Color',cmap(1)); hold on
ylimm = [-0.2 .75];
set(h(1),'FaceColor',cmap(2),'Edgecolor',cmap(2))
set(h(2),'FaceColor',cmap(1),'Edgecolor',cmap(1))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long,'yLim',ylimm,'XTickLabelRotation',45,'Layer','top');
ylabel('Correlation');
SaveFigure('Output/PC',0)
clf('reset')

%% 2. Dynamic f and g
%==========================================================================
% Boxplot
FigSize(20,7)
errorbar(1:ncountry,DYN.corrne(:,1),DYN.corrne(:,2),'o','Color',cmap(1),'Linewidth',1.2,'MarkerFaceColor',cmap(2))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long(DYN.corrne(:,3)),'yLim',[-0.8 0.5],'XTickLabelRotation',45,'Layer','top');
hold on; plot(zeros(ncountry,1),'k:','Linewidth',0.5)
ylabel('Correlation');
SaveFigure('Output\DYN_BOXne',0)
clf('reset')
% Pairwise correlation
FigSize(20,7)
h = bar([DYN.PCn_num(:,1) DYN.PCe_num(:,1)],'BarWidth',1.2); hold on;
bar1(1:ncountry) = mean(DYN.PCn_num(:,1)); plot(bar1,':','LineWidth',2,'Color',cmap(2)); hold on
bar2(1:ncountry) = mean(DYN.PCe_num(:,1)); plot(bar2,':','LineWidth',2,'Color',cmap(1)); hold on
ylimm = [-0.2 .75];
set(h(1),'FaceColor',cmap(2),'Edgecolor',cmap(2))
set(h(2),'FaceColor',cmap(1),'Edgecolor',cmap(1))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long,'yLim',ylimm,'XTickLabelRotation',45,'Layer','top');
ylabel('Correlation');
SaveFigure('Output/PC_F_G',0)
clf('reset')

%% 3. Dynamic f only (ie, gnot)
%==========================================================================
%  Boxplot
FigSize(20,7)
errorbar(1:ncountry,DYNgnot.corrne(:,1),DYNgnot.corrne(:,2),'o','Color',cmap(1),'Linewidth',1.2,'MarkerFaceColor',cmap(2))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long(DYNgnot.corrne(:,3)),'yLim',[-0.8 0.5],'XTickLabelRotation',45,'Layer','top');
hold on; plot(zeros(ncountry,1),'k:','Linewidth',0.5)
ylabel('Correlation');
SaveFigure('Output\DYNgnot_BOXne',0)
clf('reset')
% Pairwise correlation
FigSize(20,7)
h = bar([DYNgnot.PCn_num(:,1) DYNgnot.PCe_num(:,1)],'BarWidth',1.2); hold on;
bar1(1:ncountry) = mean(DYNgnot.PCn_num(:,1)); plot(bar1,':','LineWidth',2,'Color',cmap(2)); hold on
bar2(1:ncountry) = mean(DYNgnot.PCe_num(:,1)); plot(bar2,':','LineWidth',2,'Color',cmap(1)); hold on
ylimm = [-0.2 .75];
set(h(1),'FaceColor',cmap(2),'Edgecolor',cmap(2))
set(h(2),'FaceColor',cmap(1),'Edgecolor',cmap(1))
set(gca,'xLim',[0 ncountry+1],'xTick',(1:ncountry),'xTickLabel',cnames_long,'yLim',ylimm,'XTickLabelRotation',45,'Layer','top');
ylabel('Correlation');
SaveFigure('Output/PC_F',0)
clf('reset')
end
close all

disp('Done!')
