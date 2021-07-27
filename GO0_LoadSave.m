%--------------------------------------------------------------------------%
% This code reads from DATA.xlsx and imports the data needed to            %
% construct the GVOL model                                                 %
%                                                                          %
%--------------------------------------------------------------------------%
clc
disp('Load data and save...')
 
%% LOAD GDP DATA
% This is real GDP
[xlsdata, xlstext] = xlsread('DATA.xlsx','GDP');
date = xlstext(2:end,1);
cnames = xlstext(1,2:end)';
Y = Num2NaN(xlsdata);

%% LOAD VOL DATA
% This is the vol of equity
[xlsdata, ~] = xlsread('DATA.xlsx','VEQ');
VEQ = Num2NaN(xlsdata);

%% LOAD OTHER STUFF
[PPP, xlstext] = xlsread('DATA.xlsx','weights');
weights = nan(length(PPP),1);
for ii=1:length(PPP)
    weights(ii,1) = 1./length(PPP);
    weights_ppp(ii,1) = PPP(ii)./sum(PPP);
end
cnames_long = xlstext(2:end,1);
[xlsdata, xlstext] = xlsread('DATA.xlsx','MACRO');
macro_names = xlstext(1,2:end);
for ii=1:size(xlsdata,2)
    MACRO.(macro_names{ii}) = Num2NaN(xlsdata(:,ii));
end
 
%% PRELIMS
nobs = length(date);
ncountry = length(weights);
fout = 'Output\GVOL_RESULTS.xlsx';
v = log(VEQ);
y = 100*XoX(Y,1,'diff');

% Write
if savexls 
    xlswrite(fout, TabPrint(v,cnames',date,6),'v'); 
    xlswrite(fout, TabPrint(y,cnames',date,6),'y'); 
end

clear xlstext xlsdata
disp('Done!')

