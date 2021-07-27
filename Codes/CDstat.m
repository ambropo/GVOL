function OUT = CDstat(DATA,T)
% =======================================================================
% Compute CD test statistic.
% =======================================================================
% OUT = CDstat(DATA)
% -----------------------------------------------------------------------
% INPUT
%	- DATA: correlation matrix (N x N)
% -----------------------------------------------------------------------
% OPTIONAL INPUT
%	- T: number of observations per pair over which correlation has been
%        computed. If T is scalar, same number for each pair. If T is 
%        (N x N) different number for each pair
%------------------------------------------------------------------------
% OUPUT
%	- OUT: CD test statistic
% =======================================================================
% Ambrogio Cesa Bianchi, April 2017
% ambrogio.cesabianchi@gmail.com



%% Preliminaries: define the dimension of the matrix of interest
% =========================================================================
[nvar, ~] = size(DATA);

if length(T)==1
    T(1:nvar,1:nvar) = T;
end

aux = 0;
for ii=1:nvar-1
    for jj=ii+1:nvar
        aux = aux + sqrt(T(ii,jj))*(DATA(ii,jj));
    end
end
OUT = sqrt(2/(nvar*(nvar-1)))*aux;
