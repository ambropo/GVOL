%--------------------------------------------------------------------------%
%                         Replication codes for                            %
%    "Uncertainty and economic activity: A multi-country perspective"      %
%            by A. Cesa-Bianchi, M.H. Pesaran, and A. Rebucci              %
%                             July, 2019                                   %
%                                                                          %
%--------------------------------------------------------------------------%
% This code replicates the results in Cesa-Bianchi, M.H. Pesaran, and A. 
% Rebucci (2019) "Uncertainty and economic activity: A multi-country 
% perspective," fortcoming in the Review of Financial Studies. 
%--------------------------------------------------------------------------%
% BEFORE RUNNING THE CODE: 
% Please make sure to add the folder "Code" and all its subfolders to your 
% Matlab path. The code has been tested with Matlab R2017b.
%--------------------------------------------------------------------------%
% WARNING FOR MAC USERS: 
% Because of compatibility issues you need to set savexls = 0 to be able to
% run the code.
%--------------------------------------------------------------------------%


%% Prelims
clearvars; close all; warning('off'); clc
mkdir('Output')
savexls = 1; % if =1 saves results to excel, =0 otherwise
saveplt = 1; % if =1 saves charts to pdf, =0 otherwise

%% Codes
GO0_LoadSave
GO1_Factors
GO2_REGS
GO3_REGSplot
GO4_GVOL
GO5_IRF
GO6_IRFplot
GO7_FEVD
GO8_FEVDplot

%% Clean and save
clc; 
save OUT


