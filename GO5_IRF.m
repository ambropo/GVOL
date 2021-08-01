%--------------------------------------------------------------------------%
% This code computes IRFs to global shocks, country specific shocks, and   %
% for comparison the impulse responses froma  simple bivariate VAR in      %
% [v,y]                                                                    %
%--------------------------------------------------------------------------%
clc
disp('Compute IRFs of GVOL model...') 


%% 1. PRELIMNS
%%=========================================================================
% Create a cell array with steps for exporting
hlabel = num2cell(1:nsteps)';
% Create an index that picks only VOL or GDP
indv = 1:2:2*ncountry;
indy = 2:2:2*ncountry;

%% 2. IMPULSE RESPONSE TO "F" FACTOR (Wold)
%==========================================================================
% Impulse vector selects f from [f g]
sel_shock = [1;0];
% Compute IRs from Wold
irf_aux = zeros(2*ncountry,nsteps);
for ii=1:nsteps
    irf_aux(:,ii) = C(:,:,ii)*(GVOL.B*sel_shock);
end
% Store responses in IR
IR.IRf_VOL = irf_aux(indv,:);
IR.IRf_GDP = irf_aux(indy,:);
% Save in cell array for exporting
IRcel.IRf_VOL = TabPrint(irf_aux(indv,:)',cnames',hlabel,approx);
IRcel.IRf_GDP = TabPrint(irf_aux(indy,:)',cnames',hlabel,approx);
% Averages
ww = transpose(weights).^2;
IRw.IRf_wVOL = transpose(weights)*IR.IRf_VOL; 
IRw.IRf_wGDP = transpose(weights)*IR.IRf_GDP; 
for jj=1:nsteps
    std_VOL(jj) = sqrt(ww*(IR.IRf_VOL(:,jj)-IRw.IRf_wVOL(jj)).^2);
    std_GDP(jj) = sqrt(ww*(IR.IRf_GDP(:,jj)-IRw.IRf_wGDP(jj)).^2);
end
IRw.IRf_inf_wVOL = IRw.IRf_wVOL - 2.*std_VOL; IRw.IRf_sup_wVOL = IRw.IRf_wVOL + 2.*std_VOL; 
IRw.IRf_inf_wGDP = IRw.IRf_wGDP - 2.*std_GDP; IRw.IRf_sup_wGDP = IRw.IRf_wGDP + 2.*std_GDP; 

    
%% 3. IMPULSE RESPONSE TO "G" FACTOR (Wold)
%==========================================================================
% Impulse vector selects g from [f g]
sel_shock = [0;1];
% Compute IRs from Wold
irf_aux = zeros(2*ncountry,nsteps);
for ii=1:nsteps
    irf_aux(:,ii) = C(:,:,ii)*(GVOL.B*sel_shock);
end
% Store responses in IR
IR.IRg_VOL = irf_aux(indv,:);
IR.IRg_GDP = irf_aux(indy,:);
% Save in cell array for exporting
IRcel.IRg_VOL = TabPrint(irf_aux(indv,:)',cnames',hlabel,approx);
IRcel.IRg_GDP = TabPrint(irf_aux(indy,:)',cnames',hlabel,approx);
% Averages
ww = transpose(weights).^2;
IRw.IRg_wVOL = transpose(weights)*IR.IRg_VOL; 
IRw.IRg_wGDP = transpose(weights)*IR.IRg_GDP; 
for jj=1:nsteps
    std_VOL(jj) = sqrt(ww*(IR.IRg_VOL(:,jj)-IRw.IRg_wVOL(jj)).^2);
    std_GDP(jj) = sqrt(ww*(IR.IRg_GDP(:,jj)-IRw.IRg_wGDP(jj)).^2);
end
IRw.IRg_inf_wVOL = IRw.IRg_wVOL - 2.*std_VOL; IRw.IRg_sup_wVOL = IRw.IRg_wVOL + 2.*std_VOL; 
IRw.IRg_inf_wGDP = IRw.IRg_wGDP - 2.*std_GDP; IRw.IRg_sup_wGDP = IRw.IRg_wGDP + 2.*std_GDP; 


%% 4. IMPULSE RESPONSE TO COUNTRY-SPECIFIC VOL
%==========================================================================
% Impulse response are computed with GIRF methodology, using block diagonal
% covariance matrix
sel_shock = [1;0];
for oo=1:ncountry
    eps = CommonSample([DYN.RESv(:,oo) DYN.RESy(:,oo)]);
    sigma = (eps'*eps)/length(eps);
    IMPi = sigma*diag(1./sqrt(diag(sigma))); % GIRF as in equation 6, Warne. Equal to Cholesky: Bi = chol(sigma)'
    IMP = zeros(2*ncountry,1);
    IMP(2*oo-1:2*oo) = IMPi*sel_shock;
    % Compute IRs from Wold
    irf_aux = zeros(2*ncountry,nsteps);
    for ii=1:nsteps
        irf_aux(:,ii) = C(:,:,ii)*IMP;
    end
    % Store responses in IRn (this is country i's response to its own VOL shock
    IR.IRn_VOL(oo,:) = irf_aux(2*oo-1,:);
    IR.IRn_GDP(oo,:) = irf_aux(2*oo,:);
end
% Save in cell array for exporting
IRcel.IRn_VOL = TabPrint(IR.IRn_VOL',cnames',hlabel,approx);
IRcel.IRn_GDP = TabPrint(IR.IRn_GDP',cnames',hlabel,approx);
% Averages
ww = transpose(weights).^2;
IRw.IRn_wVOL = transpose(weights)*IR.IRn_VOL; 
IRw.IRn_wGDP = transpose(weights)*IR.IRn_GDP; 
for jj=1:nsteps
    std_VOL(jj) = sqrt(ww*(IR.IRn_VOL(:,jj)-IRw.IRn_wVOL(jj)).^2);
    std_GDP(jj) = sqrt(ww*(IR.IRn_GDP(:,jj)-IRw.IRn_wGDP(jj)).^2);
end
IRw.IRn_inf_wVOL = IRw.IRn_wVOL - 2.*std_VOL; IRw.IRn_sup_wVOL = IRw.IRn_wVOL + 2.*std_VOL; 
IRw.IRn_inf_wGDP = IRw.IRn_wGDP - 2.*std_GDP; IRw.IRn_sup_wGDP = IRw.IRn_wGDP + 2.*std_GDP; 


%% 5. IMPULSE RESPONSE TO COUNTRY-SPECIFIC GDP
%==========================================================================
% Impulse response are computed with GIRF methodology, using block diagonal
% covariance matrix
sel_shock = [0;1];
for oo=1:ncountry
    eps = CommonSample([DYN.RESv(:,oo) DYN.RESy(:,oo)]);
    sigma = (eps'*eps)/length(eps);
    IMPi = sigma*diag(1./sqrt(diag(sigma))); % GIRF as in equation 6, Warne. Equal to Cholesky: Bi = chol(sigma)'
    IMP = zeros(2*ncountry,1);
    IMP(2*oo-1:2*oo) = IMPi*sel_shock;
    % Compute IRs from Wold
    irf_aux = zeros(2*ncountry,nsteps);
    for ii=1:nsteps
        irf_aux(:,ii) = C(:,:,ii)*IMP;
    end
    % Store responses in IRe (this is country i's response to its own GDP shock
    IR.IRe_VOL(oo,:) = irf_aux(2*oo-1,:);
    IR.IRe_GDP(oo,:) = irf_aux(2*oo,:);
end
% Save in cell array for exporting
IRcel.IRe_VOL = TabPrint(IR.IRe_VOL',cnames',hlabel,approx);
IRcel.IRe_GDP = TabPrint(IR.IRe_GDP',cnames',hlabel,approx);
% Averages
ww = transpose(weights).^2;
IRw.IRe_wVOL = transpose(weights)*IR.IRe_VOL; 
IRw.IRe_wGDP = transpose(weights)*IR.IRe_GDP; 
for jj=1:nsteps
    std_VOL(jj) = sqrt(ww*(IR.IRe_VOL(:,jj)-IRw.IRe_wVOL(jj)).^2);
    std_GDP(jj) = sqrt(ww*(IR.IRe_GDP(:,jj)-IRw.IRe_wGDP(jj)).^2);
end
IRw.IRe_inf_wVOL = IRw.IRe_wVOL - 2.*std_VOL; IRw.IRe_sup_wVOL = IRw.IRe_wVOL + 2.*std_VOL; 
IRw.IRe_inf_wGDP = IRw.IRe_wGDP - 2.*std_GDP; IRw.IRe_sup_wGDP = IRw.IRe_wGDP + 2.*std_GDP; 


%% 6. BIVARIATE - VOL SHOCK
%==========================================================================
% Impulse response are computed with Cholesky within each country
for oo=1:ncountry
    aux = CommonSample([v(:,oo) y(:,oo)]);
    [BIV.cnames{oo}, BIVopt.cnames{oo}] = VARmodel(aux,nlag,1);
    BIVopt.cnames{oo}.nsteps = 20;
    [irf_aux, BIV.cnames{oo}] = VARir(BIV.cnames{oo},BIVopt.cnames{oo});
    IR.IRbiv_VOL(oo,:) = irf_aux(:,1,1)';
    IR.IRbiv_GDP(oo,:) = irf_aux(:,2,1)';
end
% Averages
ww = transpose(weights).^2;
IRw.IRbiv_wVOL = transpose(weights)*IR.IRbiv_VOL; 
IRw.IRbiv_wGDP = transpose(weights)*IR.IRbiv_GDP; 
for jj=1:nsteps
    std_VOL(jj) = sqrt(ww*(IR.IRbiv_VOL(:,jj)-IRw.IRbiv_wVOL(jj)).^2);
    std_GDP(jj) = sqrt(ww*(IR.IRbiv_GDP(:,jj)-IRw.IRbiv_wGDP(jj)).^2);
end
IRw.IRbiv_inf_wVOL = IRw.IRbiv_wVOL - 2.*std_VOL; IRw.IRbiv_sup_wVOL = IRw.IRbiv_wVOL + 2.*std_VOL; 
IRw.IRbiv_inf_wGDP = IRw.IRbiv_wGDP - 2.*std_GDP; IRw.IRbiv_sup_wGDP = IRw.IRbiv_wGDP + 2.*std_GDP; 


%% 6. EXPORT
%==========================================================================
disp('Export mean group estimates...')
if savexls
    % f factor
    xlswrite(fout, IRcel.IRf_VOL,'IRf_v'); 
    xlswrite(fout, IRcel.IRf_GDP,'IRf_y'); 
    % g factor
    xlswrite(fout, IRcel.IRg_VOL,'IRg_v'); 
    xlswrite(fout, IRcel.IRg_GDP,'IRg_y'); 
end

close all
disp('Done!')



