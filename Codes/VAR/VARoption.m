function VARopt = VARoption
%==========================================================================
% Optional inputs for VAR analysis
%==========================================================================
% Ambrogio Cesa Bianchi, March 2017
% ambrogio.cesabianchi@gmail.com

VARopt.vnames    = [];    % endogenous variables' names
VARopt.vnames_ex = [];    % exogenous variables' names
VARopt.snames    = [];    % shocks' names (for sign restriction only)
VARopt.nsteps    = 40;    % number of steps for IRFs and FEVDs
VARopt.impact    = 0;     % size of the shock for IRFs: 0=1stdev, 1=unit shock
VARopt.shut      = 0;     % forces the IRF of one variable to zero
VARopt.ident     = 'oir'; % identification method for IRFs ('oir' short-run restr, 'bq' long-run restr, 'sr' sign restr)
VARopt.recurs    = 'wold';% method for computation of recursive stuff ('wold' form MA representation, 'comp' for companion form)
VARopt.ndraws    = 100;   % draws for bootstrap and sign restrictions
VARopt.pctg      = 95;    % confidence bands for bootstrap
VARopt.method    = 'bs';  % methodology for error bands, 'bs' bootstrap, 'wild' wild bootstrap
VARopt.pick      = 0;     % selects one variable for IRFs and FEVDs plots (0 => plot all)
VARopt.quality   = 0;     % quality of exported figures: 1=high (ghostscript required), 0=low
VARopt.suptitle  = 0;     % title on top of figures (Figure Toolbox required)
VARopt.firstdate = [];    % initial date of the sample in format 1999.75 => 1999Q4 (both for annual and quarterly data)
VARopt.frequency = 'q';   % frequency of the data: 'm' monthly, 'q' quarterly, 'y' yearly
VARopt.figname   = [];    % string for names of exported figures