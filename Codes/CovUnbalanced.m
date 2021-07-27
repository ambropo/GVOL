function [COV, OUT] = CovUnbalanced(DATA,labels)
% =======================================================================
% Compute correlation of a panel of time series DATA (with T 
% observations and N variables) using rows with no missing values in 
% column i or j. 
% 
% note: for each pair of columns, the correlation is computed using the
% maximum amount of available observations. This is to deal with unbalanced
% panels. Therefore, NaN are accepted.
% =======================================================================
% [CORR, OUT] = CorrUnbalanced(DATA,labels)
% -----------------------------------------------------------------------
% INPUT
%	- DATA: matrix DATA T (observations) x N (variables)
%------------------------------------------------------------------------
% OPTIONAL INPUT
%   - labels: Default "Variable", names of each variable j
%------------------------------------------------------------------------
% OUPUT
%	- CORR: matrix of correlation N x N 
%	- OUT.N: number of observations used for each pair
%	- OUT.table: formatted table of correlation matrix with titles
% =======================================================================
% EXAMPLE
% DATA = rand(50,4);
% [CORR, OUT] = CorrUnbalanced(DATA)
% =======================================================================
% Ambrogio Cesa Bianchi, April 2017
% ambrogio.cesabianchi@gmail.com



%% Preliminaries: define the dimension of the matrix of interest
% =========================================================================
[nobs, nvar] = size(DATA);

% If no names are provided set it to 'Variable'
if ~exist('labels','var')
    labels(1,1:nvar) = {'Variable'};
end

% If labels are entered as jx1 vector, transpose it
if size(labels,1) > 1
    labels = labels';
end

%% Compute pairwise correlation
% =========================================================================
COV = nan(nvar,nvar);
OUT.N = nan(nvar,nvar);
% Compute the correlation matrix using the maximum amount of avaliable obs
for ii=1:nvar
    X1 = DATA(:,ii);
    for jj=1:nvar
        X2 = DATA(:,jj);
        Y = CommonSample([X1 X2]);
        if isempty(Y)
            COV(ii,jj) = NaN;
            OUT.N(ii,jj) = 0;
        else
            aux = cov(Y);
            OUT.N(ii,jj) = length(Y);
            COV(ii,jj) = aux(1,2); % take element 1,2
        end
    end
end



% %Write the table in PairCorr.xls with titles
% title = {'' , 'Level', 'First Diff.'};
% TABLE = [labels  ; num2cell(PC')];
% TABLE = [title ; TABLE'];
