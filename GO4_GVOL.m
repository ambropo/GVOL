%--------------------------------------------------------------------------%
% This code creates the GVOL model by combining the country-specific       %
% models in GO2_REGS.                                                           %
%                                                                          %
%--------------------------------------------------------------------------%
clc
disp('Construct GVOL model...') 

% Steps for impulse responses
nsteps = 20;

%% 1. ENDOGENOUS (PHI1, PHI2, ..., PSIp)
%==========================================================================
J=1+const; % drop the constant
for jj=1:nlag
    tempcel{1} = ['PHI' num2str(jj)];
    GVOL.(tempcel{1}) = zeros(2*ncountry,2*ncountry);
    I=1; 
    for ii=1:ncountry
        % Get coefficients from VAR for VOL
        Fv = DYN.VARv.(cnames{ii}).Ft';
        GVOL.(tempcel{1})(I,I:I+1) = Fv(1,J:J+1);
        % Get coefficients from VAR for GDP
        Fy = DYN.VARy.(cnames{ii}).Ft';
        GVOL.(tempcel{1})(I+1,I:I+1) = Fy(2,J:J+1);
        % Count
        I = I+2;
    end
    J = J+2;
end


%% 2. LAGGED EXOGENOUS (psi1, psi2, ..., psip)
%==========================================================================
J=1+const+nlag*2; % drop the constant + lags of the endogenous (2*nlag)
for jj=1:nlag
    tempcel{1} = ['psi' num2str(jj)];
    GVOL.(tempcel{1}) = zeros(2*ncountry,2);
    I=1; 
    for ii=1:ncountry
        % Get coefficients from VAR for VOL
        Fv = DYN.VARv.(cnames{ii}).Ft';
        GVOL.(tempcel{1})(I,:) = Fv(1,J:J+1);
        % Get coefficients from VAR for GDP
        Fy = DYN.VARy.(cnames{ii}).Ft';
        GVOL.(tempcel{1})(I+1,:) = Fy(2,J:J+1);
        % Count
        I = I+2;
    end
    J = J+2;
end


%% 3. Weights (W)
%==========================================================================
GVOL.W = zeros(2,2*ncountry);
I=1; 
for ii=1:ncountry
    GVOL.W(1,I) = weights(ii)'; 
    GVOL.W(2,I+1) = weights(ii)';
    GVOL.Wppp(1,I) = weights_ppp(ii)'; 
    GVOL.Wppp(2,I+1) = weights_ppp(ii)';
    I = I+2;
end


%% 4. SUM OF ENDO AND EXO (PSI1, PSI2, ..., PSIp), where PSI=PHI+psi*W
%==========================================================================
for jj=1:nlag
    tempcel{1} = ['PSI' num2str(jj)];
    aux_phi{1} = ['PHI' num2str(jj)];
    aux_psi{1} = ['psi' num2str(jj)];
    GVOL.(tempcel{1}) = GVOL.(aux_phi{1}) + GVOL.(aux_psi{1})*GVOL.W;
end
% Rewrite PSI1, PSI2, ... as PSI = [PSI1 PSI2 .. PSIp]
GVOL.PSI = [];
for jj=1:nlag
    tempcel{1} =  ['PSI' num2str(jj)];
    GVOL.PSI = [GVOL.PSI GVOL.(tempcel{1})];
end


%% 5. Factor Impact (B)
%==========================================================================
J=1+const+nlag*2+nlag*2; % drop the constant + lags of the endogenous (2*nlag) + lags of the exogenous (2*nlag)
B = zeros(2*ncountry,2);
I=1; 
for ii=1:ncountry
    % Get coefficients from VAR for VOL
    Fv = DYN.VARv.(cnames{ii}).Ft';
    GVOL.B(I,:) = Fv(1,J:J+1);
    % Get coefficients from VAR for GDP
    Fy = DYN.VARy.(cnames{ii}).Ft';
    GVOL.B(I+1,1) = Fy(2,J);
    % Count
    I = I+2;
end


%% 6. Residuals (RES) and data (ENDO)
%==========================================================================
I=1; 
for ii=1:ncountry
    GVOL.RES(I,:)= DYN.RESv(2:end,ii)';   % 2:end is because we take first diff of GDP
    GVOL.RES(I+1,:)= DYN.RESy(2:end,ii)'; % 2:end is because we take first diff of GDP

    GVOL.ENDO(I,:)= v(2:end,ii)';   % 2:end is because we take first diff of GDP
    GVOL.ENDO(I+1,:)= y(2:end,ii)'; % 2:end is because we take first diff of GDP

    I = I+2;
end


%% 7. COVARIANCE MATRIX
%==========================================================================

% 7.1 Full (balanced and unbalanced) covariance matrix
%--------------------------------------------------------------------------
GVOL.VCV = cov(CommonSample(GVOL.RES'));
GVOL.COR = corr(CommonSample(GVOL.RES'));
GVOL.EndoCORR = corr(CommonSample(GVOL.ENDO'));

% 7.2 Diagonal covariance matrix
%--------------------------------------------------------------------------
GVOL.VCV_Diag = diag(diag(GVOL.VCV));
GVOL.VCV_Diag_Chol = chol(GVOL.VCV_Diag)';

% 7.3 Block diagonal covariance matrix
%--------------------------------------------------------------------------
GVOL.VCV_BlockDiag = zeros(2*ncountry,2*ncountry);
for ii=1:2:2*ncountry
    % Create a block-diagonal VCV matrix for the GVOL
    GVOL.VCV_BlockDiag(ii:ii+1,ii:ii+1) = GVOL.VCV(ii:ii+1,ii:ii+1);
    % Create the Cholesky factor of each block
    GVOL.VCV_BlockDiag_Chol(ii:ii+1,ii:ii+1) = chol(GVOL.VCV(ii:ii+1,ii:ii+1))';
end

% 7.4 Threshold correlation and covariance matrix
%--------------------------------------------------------------------------
N = 2*ncountry;
R = GVOL.COR;
c_delta = 1; delta = 1;
c_d = 1; d = 1;
for ii=1:N
    for jj=1:N
        T = nobs;
        p = 0.05;
        fN = c_delta*N^delta;
        cpN = norminv(1-p/(2*fN));
        test = abs(R(ii,jj)) - T^(-0.5)*cpN;
        if test>0
            Rmt(ii,jj) = R(ii,jj);
        else
            Rmt(ii,jj) = 0;
        end            
    end
end
GVOL.COR_Thres = Rmt;
GVOL.VCV_Thres = GVOL.VCV_Diag_Chol*Rmt*GVOL.VCV_Diag_Chol;

% 7.5 List of the within-country correlations (with significance at 95% level)
%--------------------------------------------------------------------------
app_table(1,1:2)=NaN;
for ii=2:N
    if Rmt(ii,ii-1)==0
        app_table(ii,1)=GVOL.COR(ii,ii-1);
        app_table(ii,2)=0;
    else
        app_table(ii,1)=Rmt(ii,ii-1);
        app_table(ii,2)=1;
    end
end
app_table=app_table(2:2:end,:);
GVOL.COR_ThresWithinList = TabPrint(app_table,{'Corr','Significance'},cnames);
if savexls
    xlswrite(fout,{'List of within-country correlations (0=insign/1=signif'},'GVOL_CORR_WithinList','A1'); 
    xlswrite(fout,GVOL.COR_ThresWithinList,'GVOL_CORR_WithinList','A2'); 
end
    
% 7.6 List of all correlations that survive thresholding approach
%--------------------------------------------------------------------------
Rmt_vec = vec(tril(Rmt,-1));
index = Rmt_vec~=0;
Rmt_vec_nonzero = Rmt_vec(index);
idx = 1; 
for ii=1:N
    for jj=1:N
        % If jj is odd use VOL, if even use GDP
        if mod(jj,2)
            label_two(idx,1) = {[cnames{round(jj/2)} ' VOL']};
        else
            label_two(idx,1) = {[cnames{round(jj/2)} ' GDP']};
        end
        % If ii is odd use VOL, if even use GDP
        if mod(ii,2)
            label_one(idx,1) = {[cnames{round(ii/2)} ' VOL']};
        else
            label_one(idx,1) = {[cnames{round(ii/2)} ' GDP']};
        end
        idx = idx + 1; 
    end
end
GVOL.COR_ThresList = [label_one(index) label_two(index) num2cell(Rmt_vec_nonzero)];
if savexls
    xlswrite(fout,GVOL.COR_Thres,'GVOL_CORR'); 
    xlswrite(fout,{'List of all correlations that survive thresholding approach'},'GVOL_CORR_List','A1'); 
    xlswrite(fout,GVOL.COR_ThresList,'GVOL_CORR_List','A2'); 
end

%% CORR data and residuals
%--------------------------------------------------------------------------
N = 2*ncountry;
R = GVOL.EndoCORR;
c_delta = 1; delta = 1;
c_d = 1; d = 1;
for ii=1:N
    for jj=1:N
        T = nobs;
        p = 0.05;
        fN = c_delta*N^delta;
        cpN = norminv(1-p/(2*fN));
        test = abs(R(ii,jj)) - T^(-0.5)*cpN;
        if test>0
            Rmt(ii,jj) = R(ii,jj);
        else
            Rmt(ii,jj) = 0;
        end            
    end
end
GVOL.EndoCORR_Thres = Rmt;

% PLot
aux1 = tril(GVOL.EndoCORR,-1);
aux1(aux1==0)=NaN; aux1(aux1<0)=-1; aux1(aux1>0)=1;

aux2 = tril(GVOL.COR_Thres,-1);
aux2(aux2==0)=NaN; aux2(aux2<0)=-1; aux2(aux2>0)=1;
I = 1;
for ii=1:ncountry
    lab(I,1) =   {[cnames{ii} ' VOL']};
    lab(I+1,1) = {[cnames{ii} ' GDP']};
    I = I+2;
end

for ii=1:N
    if sum(isnan(aux2(ii,:)))<64
        lab_y(1,ii) = lab(ii);
    else
        lab_y(1,ii) = {''};
    end
end

for ii=1:N
    if sum(isnan(aux2(:,ii)))<64
        lab_x(1,ii) = lab(ii);
    else
        lab_x(1,ii) = {''};
    end
end


%% 8. COMPANION MATRIX 
%==========================================================================
GVOL.PSIcomp = [GVOL.PSI; eye(2*ncountry*(nlag-1)) zeros(2*ncountry*(nlag-1),2*ncountry)];


%% 9. WOLD REPRESENTATION
%==========================================================================
% Initialize Wold multipliers
C = zeros(2*ncountry,2*ncountry,nsteps);
% A = zeros(2*ncountry,2,nsteps);
% Re-write PSI matrix to compute multipliers, ie PSI1=(:,:,1), PSI2=(:,:,2),... 
GVOL.PSIp = zeros(2*ncountry,2*ncountry,nlag);
I = 1;
for ii=1:nsteps
    if ii<=nlag
        GVOL.PSIp(:,:,ii) = GVOL.PSI(:,I:I+2*ncountry-1);
    else
        GVOL.PSIp(:,:,ii) = zeros(2*ncountry,2*ncountry);
    end
    I = I + 2*ncountry;
end
% Compute multipliers
C(:,:,1) = eye(2*ncountry);
% A(:,:,1) = GVOL.B;
for ii=2:nsteps
    jj=1;
    tempnum = 0;
    while jj<ii
        tempnum = tempnum + C(:,:,ii-jj)*GVOL.PSIp(:,:,jj);
        jj=jj+1;
    end
    C(:,:,ii) = tempnum;
end
% Update GVOL with Wold multipliers
GVOL.C = C;


%% 10. EIGENVALUES OF COMPANION MATRIX
%==========================================================================
GVOL.EIG = sort(abs(eig(GVOL.PSIcomp)));

%% 11. NUMBER OF VARIABLES IN THE GVOL
%==========================================================================
GVOL.OBS = length(GVOL.PSI1);

disp('Done!')



