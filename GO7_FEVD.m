%--------------------------------------------------------------------------%
% This code computes FEVDs with a bloc-diagonal covariance matrix. Within  %
% country covariance matrices are othogonalized with a Cholesky            %
% decomposition.                                                           %
%--------------------------------------------------------------------------%
clc
disp('Compute FEVDs (block-diagonal VCV with Cholesky) of GVOL model...') 


%% 1. PRELIMNS
%%=========================================================================
% Choose the VCV matrix you want
SIGMA = GVOL.VCV_BlockDiag;
% Create a cell array with steps for exporting
hlabel = num2cell(1:nsteps)';
% Create an index that picks only VOL or GDP
indv = 1:2:2*ncountry;
indy = 2:2:2*ncountry;
% Initialize the FEVD matrices
VDf  = nan(2*ncountry,nsteps);
VDg  = nan(2*ncountry,nsteps);
VDij = nan(2*ncountry,2*ncountry,nsteps);

%% 2. FEVD
%%=========================================================================
disp('  Compute FEVDs...') 
% Initialize some stuff
sf = [1;0];
sg = [0;1];
num_f = zeros(2*ncountry,1);
num_g = zeros(2*ncountry,1);
den1 = zeros(2*ncountry,1); 
den2 = zeros(2*ncountry,1);
num_ij  = zeros(2*ncountry,2*ncountry);
% Compute FEVD
for nn=1:nsteps
    for ii=1:2*ncountry
        % Indicator for variable i (of 64 total variables)
        ei = zeros(2*ncountry,1); ei(ii) = 1;
        % Denominator
        den1(ii) = den1(ii) + ei'*C(:,:,nn)*(GVOL.B)*(GVOL.B)'*C(:,:,nn)'*ei;
        den2(ii) = den2(ii) + ei'*C(:,:,nn)*SIGMA*C(:,:,nn)'*ei;
        % Share of variable i explained by shock to f
        num_f(ii) = num_f(ii) + (ei'*C(:,:,nn)*(GVOL.B*sf))^2;
        VDf(ii,nn) = num_f(ii)./(den1(ii)+den2(ii));
        % Share of variable i explained by shock to g
        num_g(ii) = num_g(ii) + (ei'*C(:,:,nn)*(GVOL.B*sg))^2;
        VDg(ii,nn) = num_g(ii)./(den1(ii)+den2(ii));
        % Share due to other countries
        for jj=1:2*ncountry
            % Indicator for variable j (of 64 total variables)
            ej = zeros(2*ncountry,1); ej(jj) = 1;
            % Share of variable i explained by shock to variable j
            %num_ij(ii,jj) = num_ij(ii,jj) + (1/SIGMA(jj,jj))*(ei'*C(:,:,nn)*SIGMA*ej)^2; % This is standard GFEVD
            num_ij(ii,jj) = num_ij(ii,jj) + (ei'*C(:,:,nn)*chol(SIGMA,'lower')*ej)^2;     % This is GFEVD with Cholesky on block-diagonal VCV
            VDij(ii,jj,nn) = num_ij(ii,jj)./(den1(ii)+den2(ii));
        end
    end
end    

%% 3. RESHAPE (ncountry x nsteps)
%%=========================================================================
disp('  Reshape FEVDs...') 
% GFEVD VOL
v2f = VDf(indv,:); v2g = VDg(indv,:);
v2ETA = VDij(indv,indv,:); v2EPS = VDij(indv,indy,:);
v2eta_i = nan(ncountry,nsteps); v2eps_i = nan(ncountry,nsteps);
v2eta_other = nan(ncountry,nsteps); v2eps_other = nan(ncountry,nsteps);
for jj=1:ncountry
    ff = v2f(jj,:);
    gg = v2g(jj,:);
    eta = squeeze(v2ETA(jj,jj,:))';
    eps = squeeze(v2EPS(jj,jj,:))';
    eta_all = squeeze(sum(v2ETA(jj,:,:),2))';
    eps_all = squeeze(sum(v2EPS(jj,:,:),2))';
    eta_other = (eta_all' - eta')';
    eps_other = (eps_all' - eps')';
    other = eta_other + eps_other;
    % Store
    v2eta_i(jj,:) = eta;
    v2eps_i(jj,:) = eps;
    v2eta_other(jj,:) = eta_other;
    v2eps_other(jj,:) = eps_other;
end
FEVD_BDChol.v2f = v2f; 
FEVD_BDChol.v2g = v2g;
FEVD_BDChol.v2eta_i = v2eta_i; 
FEVD_BDChol.v2eps_i = v2eps_i;
FEVD_BDChol.v2eta_other = v2eta_other; 
FEVD_BDChol.v2eps_other = v2eps_other;

% Averages
v2f_w = transpose(weights)*v2f;                 FEVD_BDChol.v2f_w = v2f_w;     
v2g_w = transpose(weights)*v2g;                 FEVD_BDChol.v2g_w = v2g_w;
v2eta_i_w = transpose(weights)*v2eta_i;         FEVD_BDChol.v2eta_i_w = v2eta_i_w; 
v2eps_i_w = transpose(weights)*v2eps_i;         FEVD_BDChol.v2eps_i_w = v2eps_i_w;
v2eta_other_w = transpose(weights)*v2eta_other; FEVD_BDChol.v2eta_other_w = v2eta_other_w; 
v2eps_other_w = transpose(weights)*v2eps_other; FEVD_BDChol.v2eps_other_w = v2eps_other_w;

% GFEVD GDP
y2f = VDf(indy,:); y2g = VDg(indy,:);
y2EPS = VDij(indy,indy,:); y2ETA = VDij(indy,indv,:);
y2eta_i = nan(ncountry,nsteps); y2eps_i = nan(ncountry,nsteps);
y2eta_other = nan(ncountry,nsteps); y2eps_other = nan(ncountry,nsteps);
for jj=1:ncountry
    ff = y2f(jj,:);
    gg = y2g(jj,:);
    eta = squeeze(y2ETA(jj,jj,:))';
    eps = squeeze(y2EPS(jj,jj,:))';
    eta_all = squeeze(sum(y2ETA(jj,:,:),2))';
    eps_all = squeeze(sum(y2EPS(jj,:,:),2))';
    eta_other = (eta_all' - eta')';
    eps_other = (eps_all' - eps')';
    other = eta_other + eps_other;
    % Store
    y2eta_i(jj,:) = eta;
    y2eps_i(jj,:) = eps;
    y2eta_other(jj,:) = eta_other;
    y2eps_other(jj,:) = eps_other;
end
FEVD_BDChol.y2f = y2f; 
FEVD_BDChol.y2g = y2g;
FEVD_BDChol.y2eta_i = y2eta_i; 
FEVD_BDChol.y2eps_i = y2eps_i;
FEVD_BDChol.y2eta_other = y2eta_other; 
FEVD_BDChol.y2eps_other = y2eps_other;

% Averages
y2f_w = transpose(weights)*y2f;                 FEVD_BDChol.y2f_w = y2f_w;     
y2g_w = transpose(weights)*y2g;                 FEVD_BDChol.y2g_w = y2g_w;
y2eta_i_w = transpose(weights)*y2eta_i;         FEVD_BDChol.y2eta_i_w = y2eta_i_w; 
y2eps_i_w = transpose(weights)*y2eps_i;         FEVD_BDChol.y2eps_i_w = y2eps_i_w;
y2eta_other_w = transpose(weights)*y2eta_other; FEVD_BDChol.y2eta_other_w = y2eta_other_w; 
y2eps_other_w = transpose(weights)*y2eps_other; FEVD_BDChol.y2eps_other_w = y2eps_other_w;


disp('Export mean group estimates...')
if savexls
    xlswrite(fout, TabPrint(100.*[v2g_w' v2eta_i_w' v2eta_other_w' v2f_w' v2eps_i_w' v2eps_other_w'],...
                            {'$\hat{\xi}$','$\hat{\eta}_i$','$\sum \hat{\eta}_j$','$\hat{\zeta}$','$\hat{\varepsilon}_i$','$\sum \hat{\varepsilon}_j$'},...
                            num2cell([1:20]')),'FEVDv'); 
    xlswrite(fout, TabPrint(100.*[y2g_w' y2eta_i_w' y2eta_other_w' y2f_w' y2eps_i_w' y2eps_other_w'],...
                            {'$\hat{\xi}$','$\hat{\eta}_i$','$\sum \hat{\eta}_j$','$\hat{\zeta}$','$\hat{\varepsilon}_i$','$\sum \hat{\varepsilon}_j$'},...
                            num2cell([1:20]')),'FEVDy'); 
end

disp('Done!')








