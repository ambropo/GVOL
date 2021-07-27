%--------------------------------------------------------------------------%
% This code estimates the factors f (zeta) and g (xi)                      % 
%                                                                          %
%                                                                          %
%--------------------------------------------------------------------------%
clc
disp('Estimate factors...')

%% 1. FACTORS 
%==========================================================================
% 1.1 Estimate the factors
%--------------------------------------------------------------------------
% Prelims
approx = 6; nlag=1; const=1; % code does not allow for nlag endo to be different from nlag exog

% Compute weighted averages (equal weights)
[vwbar, ~] = CrossWeightAverage(v,weights');
[ywbar, ~] = CrossWeightAverage(y,weights');
zvwbar = zscore(vwbar); 
zywbar = [NaN; zscore(CommonSample(ywbar))]; % add one row of NaNs in y because of first diff
date_num = Date2Num(date);
nobs = size(date,1);

% Estimate a VAR in: [vwbar ywbar]
VARxnames = {'vwbar','ywbar'}; 
[DATA, fo, lo] = CommonSample([vwbar ywbar]);
[FACTORf, FACTORfopt] = VARmodel(DATA,nlag,const);
FACTORfopt.vnames = VARxnames; 
[table, beta] = VARprint(FACTORf,FACTORfopt,approx);
f = [nan(fo+nlag,1); FACTORf.resid(:,2)];
fn = [nan(fo+nlag,1); zscore(FACTORf.resid(:,2))];

% Estimate a VAR in: [vwbar ywbar] + f
VARxnames = {'vwbar','ywbar'}; VARxnames_ex = {'f'};
[DATA, fo, lo] = CommonSample([vwbar ywbar f]);
[FACTORg, FACTORgopt] = VARmodel(DATA(:,1:2),nlag,const,DATA(:,3));
FACTORgopt.vnames = VARxnames; FACTORgopt.vnames_ex = VARxnames_ex;
[table, beta] = VARprint(FACTORg,FACTORgopt,approx);
VARxnames = {'vwbar','ywbar'}; 
g = [nan(fo+nlag,1); FACTORg.resid(:,1)];
gn = [nan(fo+nlag,1); zscore(FACTORg.resid(:,1))];

% 1.1 Export time series of the factors 
%--------------------------------------------------------------------------
if savexls
    xlswrite(fout,TabPrint([zvwbar zywbar fn gn],{'vwbar','ywbar','fn','gn'},date,approx),'Factors'); 
end  

% Correlation between f, g, ywbar, and vwbar; serial correlation
[~, TABLE] = CorrTable(CommonSample([f g ywbar vwbar MACRO.EBP MACRO.VIX MACRO.LMN]),{'f';'g';'ywbar';'vwbar';'EBP';'VIX';'LMN';});
if savexls
    xlswrite(fout,TABLE,'FactorsCorr'); 
end    

% Stats of the factors
[~, Stats] = SummStats([fn gn],1,{'fn','gn'});
[ACF.acf_f, ACF.lag_f, ACF.ci_f] = autocorr(CommonSample(f));
[ACF.acf_g, ACF.lag_g, ACF.ci_g] = autocorr(CommonSample(g));


%% 3. FIGURES
%==========================================================================
if saveplt 

% 3.1 Plot factors
%----------------------------------------------------------------------
FigSize(20,14)
subplot(2,1,1)
plot(fn,'LineStyle','-','Color',cmap(1),'LineWidth',2); hold on
axis tight; grid on; set(gca,'GridLineStyle',':')
plot(ones(length(fn)),'LineStyle',':','Color',cmap(1),'LineWidth',1); hold on
plot(-ones(length(fn)),'LineStyle',':','Color',cmap(1),'LineWidth',1); hold on
DatesPlot(date_num(1),nobs,12)
title('Panel A: Common Growth Shock $\hat{\zeta}_t$','Interpreter','Latex')
subplot(2,1,2)
plot(gn,'LineStyle','-','Color',cmap(1),'LineWidth',2); hold on
plot(ones(length(fn)),'LineStyle',':','Color',cmap(1),'LineWidth',1); hold on
plot(-ones(length(fn)),'LineStyle',':','Color',cmap(1),'LineWidth',1); hold on
axis tight; grid on; set(gca,'GridLineStyle',':')
DatesPlot(date_num(1),nobs,12)
title('Panel B: Common Financial Shock $\hat{\xi}_t$','Interpreter','Latex')
SaveFigure('Output\F_plot',0)
clf('reset')

end

disp('Done!')



