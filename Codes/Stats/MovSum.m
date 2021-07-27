function OUT = MovSum(DATA,window,na)
% =======================================================================
% Moving sum of the vector (or matrix) DATA (T obs x N variables). If 
% DATA is a matrix moving sum is computed down each column.
% =======================================================================
% OUT = MovAvg(DATA,window)
% -----------------------------------------------------------------------
% INPUT
%   DATA: T observations x N variables
%   window: window of the moving sum 
%------------------------------------------------------------------------
% OUPUT
%   OUT: T observations x N variables matrix (the first windows-1 
%       obseravations are NaN)
% =======================================================================
% Ambrogio Cesa Bianchi, November 2017
% ambrogio.cesabianchi@gmail.com

if ~exist('window','var')
    error('Must provide window size');
end  

if ~exist('na','var')
    na = 0;
end    
if min(size(DATA))==1
    DATA = DATA(:); % forces DATA to be a column vector
end

[nobs,nvar] = size(DATA);
if window>nobs
    error('window must not be greater than the length of DATA.')
end

temp=[];
if na
    for row=1:(nobs-window+1)
        temp = [temp; nansum(DATA(row:(row+window-1),:))];
    end
else
    for row=1:(nobs-window+1)
        temp = [temp; sum(DATA(row:(row+window-1),:))];
    end
end

OUT = temp;
OUT = [nan(window-1,nvar); OUT]; % add nans to make conformable to original 
