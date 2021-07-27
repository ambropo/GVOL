%--------------------------------------------------------------------------%
% This code estimates the country specific models: 1) unconditional        %
% analysis; 2) conditional analysis on both f (zeta) and g (xi);           %
% 3) conditional on f (zeta) only                                          %
%--------------------------------------------------------------------------%
clc
disp('Run country specific models...')

%% 1. UNCONDITIONAL
%==========================================================================
% Results are stored in UNC
%--------------------------------------------------------------------------

% 1.1 Regression of v on y for correlations
%--------------------------------------------------------------------------
for ii=1:length(cnames)
    [aux, fo, lo] = CommonSample([v(:,ii) y(:,ii)]);
    aux_norm = zscore(aux);
    OLSout = OLSmodel(aux_norm(:,1),aux_norm(:,2),0);
    OUTvy(1,ii) = OLSout.beta;
    OUTvy(2,ii) = OLSout.bstd;
    OUTvy(3,ii) = OLSout.tstat;
    OUTvy(4,ii) = OLSout.tprob;
    OUTvy(5,ii) = OLSout.rbar;
    temp = corrcoef(aux(:,1),aux(:,2));
    OUTvy(6,ii) = temp(1,2);
    conf(1,ii) = OLSout.beta-OLSout.bint(1);
end
% Sort correlation coefficients in a vector (cell)
[aux1, aux2] = sort(OUTvy(6,:));
UNC.CORRvy = TabPrint(aux1,cnames(aux2)',{'Corr'})';
% Correlation and confidence bands (double)
UNC.corrvy = [aux1; conf(aux2); aux2]';
% Correlation with p-values, t-stat, etc (cell)
UNC.OLSvy = TabPrint(OUTvy,cnames',{'beta';'std err';'tstat';'p';'adj R2';'Corr'},approx);

% Compute pair corr of VOL/GDP
[UNC.PCv_num, UNC.PCv] = PairCorrUnbalanced(v,cnames_long);
[UNC.PCy_num, UNC.PCy] = PairCorrUnbalanced(Y,cnames_long);

% Compute unbalanced corr matrix of VOL/GDP
[UNC.CORv, OUTtemp] = CorrUnbalanced(v,cnames); UNC.CORv_N = OUTtemp.N;
[UNC.CORy, OUTtemp] = CorrUnbalanced(y,cnames); UNC.CORy_N = OUTtemp.N;

% Compute the CD statistics of VOL/GDP
UNC.CDv = CDstat(UNC.CORv,UNC.CORv_N);
UNC.CDy = CDstat(UNC.CORy,UNC.CORy_N);

% Compute summ stats of LEVEL of VOL
[~, UNC.STATSv] = SummStats(v,0,cnames');
[~, UNC.STATSvwbar] = SummStats(CrossWeightAverage(v, weights'),0,{'vwbar'});
% Compute summ stats of GROWTH RATE of GDP
[~, UNC.STATSy] = SummStats(y,0,cnames');
[~, UNC.STATSywbar] = SummStats(CrossWeightAverage(y, weights'),0,{'ywbar'});
% Compute summ stats of LEVEL of GDP
[~, UNC.STATSy_lev] = SummStats(Y,0,cnames');
[~, UNC.STATSy_levwbar] = SummStats(CrossWeightAverage(Y, weights'),0,{'Ywbar'});

% Write
if savexls 
    xlswrite(fout, UNC.CORRvy,'UNC_CORRvy'); 
    xlswrite(fout, UNC.OLSvy,'UNC_OLSvy'); 
    xlswrite(fout, UNC.PCy,'UNC_PCy'); 
	xlswrite(fout, UNC.PCv,'UNC_PCv'); 	
    xlswrite(fout, UNC.STATSv,'v_stats'); 
	xlswrite(fout, UNC.STATSvwbar,'vwbar_stats'); 
    xlswrite(fout, UNC.STATSy,'y_stats'); 
	xlswrite(fout, UNC.STATSywbar,'ywbar_stats'); 
    xlswrite(fout, UNC.STATSy_lev,'y_lev_stats'); 
	xlswrite(fout, UNC.STATSy_levwbar,'y_levwbar_stats'); 
end

% 1.2 Lead-lag corralation between y and v
%--------------------------------------------------------------------------
lead_lag = 4; steps = lead_lag*2+1;
StoreCorr = nan(ncountry,steps); %CORR(y,v_t+j)
for ii=1:ncountry 
    aux = CommonSample([y(:,ii),v(:,ii)]);
    for j=1:steps
        tempnum = corrcoef(aux(1+lead_lag:end-lead_lag,1),aux(1+j-1:end+j-steps,2));
        if isnan(tempnum)==1
            StoreCorr(ii,j) = NaN;
        else
            StoreCorr(ii,j) = tempnum(1,2);
        end
    end
end
UNC.XCORRvy.avg = nanmean(StoreCorr,1); 
UNC.XCORRvy.med = nanmedian(StoreCorr,1); 
UNC.XCORRvy.stdev = sqrt(nanvar(StoreCorr,1)/(ncountry));
UNC.XCORRvy.low = prctile(StoreCorr,10,1); UNC.XCORRvy.upp = prctile(StoreCorr,90,1);
UNC.XCORRvy.CORR = StoreCorr;


%% 2. DYNAMIC MODEL (f and g)
%==========================================================================
% Results are stored in DYN and VAR
%--------------------------------------------------------------------------

% Initialize
approx = 6;
VARnames = {'v','y'}; 
% Create data for zlag
[~, aux] = VARmakexy([vwbar ywbar],nlag,0);
zlag = [nan(nlag,2*nlag); aux];
% Create labels for zlag
index = 1; clear aux;
for jj=1:nlag
    for ii=1:length(VARxnames)
        aux{index} = [VARxnames{ii} '(-' num2str(jj) ')'];
        index = index + 1;
    end
end

% 2.1 Volatility equation
%--------------------------------------------------------------------------
% Create labels for full list of exogenous [zlag f g]
VARnames_ex = [aux,'fn','gn'];
% Initialize
DYN.RESv = nan(nobs,ncountry);
TABLEv = [];
% Estimate VOL equation
for ii=1:length(cnames)
    % Common sample
    [DATA, fo, lo] = CommonSample([v(:,ii) y(:,ii) zlag fn gn]);
    % Estimate VARout.(cnames{ii})
    [VAR.(cnames{ii}), VARopt.(cnames{ii})] = VARmodel(DATA(:,1:2),nlag,const,DATA(:,3:end));
    VARopt.(cnames{ii}).vnames = VARnames; 
    VARopt.(cnames{ii}).vnames_ex = VARnames_ex;
    [table, beta] = VARprint(VAR.(cnames{ii}),VARopt.(cnames{ii}),approx);
    % Select 2nd column only (=VOL) and add country name
    table(1,2) = cnames(ii);
    TABLEv = [TABLEv table(:,2)]; % 2nd column, because first column is labels!
    % Save residuals from equation 1 (VOL)
    DYN.RESv(fo+max(nlag)+1:end-lo,ii) = VAR.(cnames{ii}).eq1.resid; % "eq1" is because first equation is VOL
    DYN.fo_v(ii) = fo; DYN.lo_v(ii) = lo;
    % Save Durbin-Watson test from equation 1 (VOL)
    DYN.DWv(ii,1) = VAR.(cnames{ii}).eq1.dw; % "eq1" is because first equation is VOL
    % Save VARout.(cnames{ii})
    DYN.VARv.(cnames{ii}) = VAR.(cnames{ii});
end
% Save VARX table
DYN.OLSv = [table(:,1) TABLEv];
% Save impact of factors (f and g, these are 8 rows at the end of DYN.OLSv, before nobs and R2)
DYN.IMPv = [DYN.OLSv(1,:); DYN.OLSv(end-10:end-3,:)]';
% Save Pairwise Correlation of volatility innovations (n)
[DYN.PCn_num, DYN.PCn] = PairCorrUnbalanced(DYN.RESv,cnames_long);
% Compute unbalanced corr matrix of VOL innovations (n)
[DYN.CORn, OUTtemp] = CorrUnbalanced(DYN.RESv,cnames);
DYN.CORn_N = OUTtemp.N;
% Compute the CD statistics of VOL
DYN.CDn = CDstat(DYN.CORn,DYN.CORn_N);
if savexls
    % xlswrite(fout,DYN.OLSv,'DYN_OLSv');
    xlswrite(fout,DYN.IMPv,'DYN_IMPv');
    xlswrite(fout,DYN.PCn, 'DYN_PCn');
    xlswrite(fout,TabPrint(DYN.RESv,cnames',date,approx),'DYN_RESv');
end

% 2.2 Output equation
%--------------------------------------------------------------------------
% Create labels for full list of exogenous [zlag f g]
VARnames_ex = [aux,'fn'];
% Initialize
DYN.RESy = nan(nobs,ncountry);
TABLEy = [];
% Estimate VOL equation
for ii=1:length(cnames)
    % Common sample
    [DATA, fo, lo] = CommonSample([v(:,ii) y(:,ii) zlag fn]);
    % Estimate VARout.(cnames{ii})
    [VAR.(cnames{ii}), VARopt.(cnames{ii})] = VARmodel(DATA(:,1:2),nlag,const,DATA(:,3:end));
    VARopt.(cnames{ii}).vnames = VARnames; 
    VARopt.(cnames{ii}).vnames_ex = VARnames_ex;
    [table, beta] = VARprint(VAR.(cnames{ii}),VARopt.(cnames{ii}),approx);
    % Select 3nd column only (=GDP) and add country name
    table(1,3) = cnames(ii);
    TABLEy = [TABLEy table(:,3)]; % 3rd column, because first column is labels and second is VOL!
    % Save residuals from equation 2 (GDP)
    DYN.RESy(fo+max(nlag)+1:end-lo,ii) = VAR.(cnames{ii}).eq2.resid; % "eq2" is because second equation is GDP
    DYN.fo_y(ii) = fo; DYN.lo_y(ii) = lo;
    % Save Durbin-Watson test from equation 2 (GDP)
    DYN.DWy(ii,1) = VAR.(cnames{ii}).eq2.dw; % "eq2" is because second equation is GDP
    % Save VARout.(cnames{ii})
    DYN.VARy.(cnames{ii}) = VAR.(cnames{ii});
end
% Save VARX table
DYN.OLSy = [table(:,1) TABLEy];
% Save impact of factors (f only, these are 4 rows at the end of DYN.OLSy, before nobs and R2)
DYN.IMPy = [DYN.OLSy(1,:); DYN.OLSy(end-6:end-3,:)]';
% Save Pairwise Correlation of output innovations (e)
[DYN.PCe_num, DYN.PCe] = PairCorrUnbalanced(DYN.RESy,cnames_long);
% Compute unbalanced corr matrix of GDP innovations (e)
[DYN.CORe, OUTtemp] = CorrUnbalanced(DYN.RESy,cnames);
DYN.CORe_N = OUTtemp.N;
% Compute the CD statistics of GDP
DYN.CDe = CDstat(DYN.CORe,DYN.CORe_N);% Write
if savexls
    % xlswrite(fout,DYN.OLSy,'DYN_OLSy');
    xlswrite(fout,DYN.IMPy,'DYN_IMPy');
    xlswrite(fout,DYN.PCe, 'DYN_PCe');
    xlswrite(fout,TabPrint(DYN.RESy,cnames',date,approx),'DYN_RESy');
end

% 2.3 Regression of n on e
%--------------------------------------------------------------------------
for ii=1:length(cnames)
    [aux, fo, lo] = CommonSample([DYN.RESv(:,ii) DYN.RESy(:,ii)]);
    aux_norm = zscore(aux);
    OLSout = OLSmodel(aux_norm(:,1),aux_norm(:,2),0);
    OUTne(1,ii) = OLSout.beta;
    OUTne(2,ii) = OLSout.bstd;
    OUTne(3,ii) = OLSout.tstat;
    OUTne(4,ii) = OLSout.tprob;
    OUTne(5,ii) = OLSout.rbar;
    temp = corrcoef(aux(:,1),aux(:,2));
    OUTne(6,ii) = temp(1,2);
    conf(1,ii) = OLSout.beta-OLSout.bint(1);
end
[aux1, aux2] = sort(OUTne(6,:));
DYN.CORRne = TabPrint(aux1,cnames(aux2)',{'Corr'})';
DYN.corrne = [aux1; conf(aux2); aux2]';
DYN.OLSne = TabPrint(OUTne,cnames',{'beta';'std err';'tstat';'p';'adj R2';'Corr'},approx);
if savexls
    xlswrite(fout,DYN.CORRne,'DYN_CORRne');
    xlswrite(fout,DYN.OLSne,'DYN_OLSne');
end

% 2.4 Lead-lag between n on e
%--------------------------------------------------------------------------
lead_lag = 4; steps = lead_lag*2+1;
StoreCorr = nan(ncountry,steps);% CORR(u,e_t+j)
for ii=1:ncountry 
    aux = CommonSample([DYN.RESy(:,ii),DYN.RESv(:,ii)]);
    for j=1:steps
        tempnum = corrcoef(aux(1+lead_lag:end-lead_lag,1),aux(1+j-1:end+j-steps,2));
        if isnan(tempnum)==1
            StoreCorr(ii,j) = NaN;
        else
            StoreCorr(ii,j) = tempnum(1,2);
        end
    end
end
DYN.XCORRne.avg = nanmean(StoreCorr,1); 
DYN.XCORRne.med = nanmedian(StoreCorr,1); 
DYN.XCORRne.stdev = sqrt(nanvar(StoreCorr,1)/(ncountry));
DYN.XCORRne.low = prctile(StoreCorr,10,1); DYN.XCORRne.upp = prctile(StoreCorr,90,1);
DYN.XCORRne.CORR = StoreCorr;


%% 3. DYNAMIC MODEL (with f factor only, ie "gnot")
%==========================================================================
% Results are stored in VARgnot and DYNgnot
%--------------------------------------------------------------------------

% Initialize
approx = 6;
VARnames = {'v','y'}; 
% Create data for zlag
[~, aux] = VARmakexy([vwbar ywbar],nlag,0);
zlag = [nan(nlag,2*nlag); aux];
% Create labels for zlag
index = 1; clear aux;
for jj=1:nlag
    for ii=1:length(VARxnames)
        aux{index} = [VARxnames{ii} '(-' num2str(jj) ')'];
        index = index + 1;
    end
end

% 3.1 Volatility equation
%--------------------------------------------------------------------------
% Create labels for full list of exogenous [zlag f g]
VARnames_ex = [aux,'fn'];
% Initialize
DYNgnot.RESv = nan(nobs,ncountry);
TABLEv = [];
% Estimate VOL equation
for ii=1:length(cnames)
    % Common sample
    [DATA, fo, lo] = CommonSample([v(:,ii) y(:,ii) zlag fn]);
    % Estimate VARout.(cnames{ii})
    [VARgnotv.(cnames{ii}), VARgnotvopt.(cnames{ii})] = VARmodel(DATA(:,1:2),nlag,const,DATA(:,3:end));
    VARgnotvopt.(cnames{ii}).vnames = VARnames; 
    VARgnotvopt.(cnames{ii}).vnames_ex = VARnames_ex;
    [table, beta] = VARprint(VARgnotv.(cnames{ii}),VARgnotvopt.(cnames{ii}),approx);
    % Select 2nd column only (=VOL) and add country name
    table(1,2) = cnames(ii);
    TABLEv = [TABLEv table(:,2)]; % 2nd column!!
    % Save residuals from equation 1 (VOL)
    DYNgnot.RESv(fo+max(nlag)+1:end-lo,ii) = VARgnotv.(cnames{ii}).eq1.resid; % "eq1" is because first equation is VOL
    DYNgnot.fo_v(ii) = fo; DYNgnot.lo_v(ii) = lo;
    % Save VARout.(cnames{ii})
    DYNgnot.VAR.(cnames{ii}) = VARgnotv.(cnames{ii});
end
% Save VARX table
DYNgnot.OLSv = [table(:,1) TABLEv];
% Save impact of factors (f and g, these are 4 rows at the end of DYN.OLSv, before nos and R2)
DYNgnot.IMPv = [DYNgnot.OLSv(1,:); DYNgnot.OLSv(end-6:end-3,:)]';
% Save Pairwise Correlation of volatility innovations (u)
[DYNgnot.PCn_num, DYNgnot.PCn] = PairCorrUnbalanced(DYNgnot.RESv,cnames_long);
% Compute unbalanced corr matrix of VOL innovations (e)
[DYNgnot.CORu, OUTtemp] = CorrUnbalanced(DYNgnot.RESv,cnames);
DYNgnot.CORu_N = OUTtemp.N;
% Compute the CD statistics of VOL
DYNgnot.CDu = CDstat(DYNgnot.CORu,DYNgnot.CORu_N);
% Write
% Write
if savexls
    % xlswrite(fout,DYNgnot.OLSv,'DYNgnot_OLSv');
    xlswrite(fout,DYNgnot.IMPv,'DYNgnot_IMPv');
    xlswrite(fout,DYNgnot.PCn, 'DYNgnot_PCn');
    xlswrite(fout,TabPrint(DYNgnot.RESv,cnames',date,approx),'DYNgnot_RESv');
end

% 3.2 Output equation
%--------------------------------------------------------------------------
% Create labels for full list of exogenous [zlag f g]
VARnames_ex = [aux,'fn'];
% Initialize
DYNgnot.RESy = nan(nobs,ncountry);
TABLEy = [];
% Estimate VOL equation
for ii=1:length(cnames)
    % Common sample
    [DATA, fo, lo] = CommonSample([v(:,ii) y(:,ii) zlag fn]);
    % Estimate VARout.(cnames{ii})
    [VARgnoty.(cnames{ii}), VARgnotyopt.(cnames{ii})] = VARmodel(DATA(:,1:2),nlag,const,DATA(:,3:end));
    VARgnotyopt.(cnames{ii}).vnames = VARnames; 
    VARgnotyopt.(cnames{ii}).vnames_ex = VARnames_ex;
    [table, beta] = VARprint(VARgnoty.(cnames{ii}),VARgnotyopt.(cnames{ii}),approx);
    % Select 3nd column only (=GDP) and add country name
    table(1,3) = cnames(ii);
    TABLEy = [TABLEy table(:,3)]; % 3rd column!!
    % Save residuals from equation 2 (GDP)
    DYNgnot.RESy(fo+max(nlag)+1:end-lo,ii) = VARgnoty.(cnames{ii}).eq2.resid; % "eq2" is because second equation is GDP
    DYNgnot.fo_y(ii) = fo; DYNgnot.lo_y(ii) = lo;
    % Save VARout.(cnames{ii})
    DYNgnot.VAR.(cnames{ii}) = VARgnoty.(cnames{ii});
end
% Save VARX table
DYNgnot.OLSy = [table(:,1) TABLEy];
% Save impact of factors (f only, these are 4 rows at the end of DYN.OLSy, before nos and R2)
DYNgnot.IMPy = [DYNgnot.OLSy(1,:); DYNgnot.OLSy(end-6:end-3,:)]';
% Save Pairwise Correlation of output innovations (e)
[DYNgnot.PCe_num, DYNgnot.PCe] = PairCorrUnbalanced(DYNgnot.RESy,cnames_long);
% Compute unbalanced corr matrix of GDP innovations (e)
[DYNgnot.CORe, OUTtemp] = CorrUnbalanced(DYNgnot.RESy,cnames);
DYNgnot.CORe_N = OUTtemp.N;
% Compute the CD statistics of GDP
DYNgnot.CDe = CDstat(DYNgnot.CORe,DYNgnot.CORe_N);
% Write
if savexls
    % xlswrite(fout,DYNgnot.OLSy,'DYNgnot_OLSy');
    xlswrite(fout,DYNgnot.IMPy,'DYNgnot_IMPy');
    xlswrite(fout,DYNgnot.PCe, 'DYNgnot_PCe');
    xlswrite(fout,TabPrint(DYNgnot.RESy,cnames',date,approx),'DYNgnot_RESy');
end

% 3.3 Regression of n on e
%--------------------------------------------------------------------------
for ii=1:length(cnames)
    [aux, fo, lo] = CommonSample([DYNgnot.RESv(:,ii) DYNgnot.RESy(:,ii)]);
    aux_norm = zscore(aux);
    OLSout = OLSmodel(aux_norm(:,1),aux_norm(:,2),0);
    OUTne(1,ii) = OLSout.beta;
    OUTne(2,ii) = OLSout.bstd;
    OUTne(3,ii) = OLSout.tstat;
    OUTne(4,ii) = OLSout.tprob;
    OUTne(5,ii) = OLSout.rbar;
    temp = corrcoef(aux(:,1),aux(:,2));
    OUTne(6,ii) = temp(1,2);
    conf(1,ii) = OLSout.beta-OLSout.bint(1);
end
[aux1, aux2] = sort(OUTne(6,:));
DYNgnot.CORRne = TabPrint(aux1,cnames(aux2)',{'Corr'})';
DYNgnot.corrne = [aux1; conf(aux2); aux2]';
DYNgnot.OLSne = TabPrint(OUTne,cnames',{'beta';'std err';'tstat';'p';'adj R2';'Corr'},approx);
if savexls
    xlswrite(fout,DYNgnot.CORRne,'DYNgnot_CORRne');
    xlswrite(fout,DYNgnot.OLSne,'DYNgnot_OLSne');
end

% 3.4 Lead-lag between n on e
%--------------------------------------------------------------------------
lead_lag = 4; steps = lead_lag*2+1;
StoreCorr = nan(ncountry,steps);% CORR(u,e_t+j)
for ii=1:ncountry 
    aux = CommonSample([DYNgnot.RESy(:,ii),DYNgnot.RESv(:,ii)]);
    for j=1:steps
        tempnum = corrcoef(aux(1+lead_lag:end-lead_lag,1),aux(1+j-1:end+j-steps,2));
        if isnan(tempnum)==1
            StoreCorr(ii,j) = NaN;
        else
            StoreCorr(ii,j) = tempnum(1,2);
        end
    end
end
DYNgnot.XCORRne.avg = nanmean(StoreCorr,1); 
DYNgnot.XCORRne.med = nanmedian(StoreCorr,1); 
DYNgnot.XCORRne.stdev = sqrt(nanvar(StoreCorr,1)/(ncountry));
DYNgnot.XCORRne.low = prctile(StoreCorr,10,1); DYNgnot.XCORRne.upp = prctile(StoreCorr,90,1);
DYNgnot.XCORRne.CORR = StoreCorr;

disp('Done!')


