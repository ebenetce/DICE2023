function params = LoadParams(np, params)
% LOADPARAMS Load Dice parameters. The function accepts the number of
% periods (in increments of 5 years) and it allows overriding any of the
% parameters.
%
% params = LoadParams(81, a2base = .01);

%% Parameters
arguments
    np (1,1) double {mustBeInteger} = 81

    % Population and technology
    params.gama    (1,1) double = 0.300;  % Capital elasticity in production function
    params.pop1    (1,1) double = 7752.9; % Initial world population 2020 (millions)
    params.popadj  (1,1) double = 0.145;  % Growth rate to calibrate to 2050 pop projection
    params.popasym (1,1) double = 10825;  % Asymptotic population (millions)
    params.dk      (1,1) double = 0.100;  % Depreciation rate on capital (per year)
    params.q1      (1,1) double = 135.7;  % Initial world output 2020 (trill 2019 USD)
    params.AL1     (1,1) double = 5.84;   % Initial level of total factor productivity
    params.gA1     (1,1) double = 0.066;  % Initial growth rate for TFP per 5 years
    params.delA    (1,1) double = 0.0015; % Decline rate of TFP per 5 years
    
    % Emissions parameters and Non-CO2 GHG with sigma = emissions/output
    params.gsigma1   (1,1) double = -0.015; % Initial growth of sigma (per year)
    params.delgsig   (1,1) double =  0.96;  % Decline rate of gsigma per period
    params.asymgsig  (1,1) double = -0.005; % Asymptotic gsigma
    params.e1        (1,1) double = 37.56;  % Industrial emissions 2020 (GtCO2 per year)
    params.miu1      (1,1) double =  0.05;  % Emissions control rate historical 2020
    params.fosslim   (1,1) double =  6000;  % Maximum cumulative extraction fossil fuels (GtC)
    params.CumEmiss0 (1,1) double =  633.5; % Cumulative emissions 2020 (GtC)

    % Climate damage parameters
    params.a1     (1,1) double = 0;        % Damage intercept
    params.a2base (1,1) double = 0.003467; % Damage quadratic term rev 01-13-23
    params.a3     (1,1) double = 2.00;     % Damage exponent

    % Abatement cost
    params.expcost2  (1,1) double = 2.6;    % Exponent of control cost function
    params.pback2050 (1,1) double = 515;    % Cost of backstop 2019$ per tCO2 2050
    params.gback     (1,1) double = -0.012; % Initial cost decline backstop cost per year
    params.cprice1   (1,1) double = 6;      % Carbon price 2020 2019$ per tCO2
    params.gcprice   (1,1) double = 0.025;  % Growth rate of base carbon price per year

    % Limits on emissions controls
    params.limmiu2070 (1,1) double = 1.0;  % Emission control limit from 2070
    params.limmiu2120 (1,1) double = 1.1;  % Emission control limit from 2120
    params.limmiu2200 (1,1) double = 1.05; % Emission control limit from 2220
    params.limmiu2300 (1,1) double = 1.0;  % Emission control limit from 2300
    params.delmiumax  (1,1) double = 0.12; % Emission control delta limit per period

    % Preferences, growth uncertainty, and timing
    params.betaclim  (1,1) double = 0.5;   % Climate beta
    params.elasmu    (1,1) double = 0.95;  % Elasticity of marginal utility of consumption
    params.prstp     (1,1) double = 0.001; % Pure rate of social time preference rho
    params.pi        (1,1) double = 0.05;  % Capital risk premium    
    params.k0        (1,1) double = 295;   % Initial capital stock calibrated (1012 2019 USD)

    params.siggc1    (1,1) double = 0.01;  %  Annual standard deviation of consumption growth    

    % Scaling so that MU(C(1)) = 1 and objective function = PV consumption
    params.tstep  (1,1) double = 5;           % Years per Period
    params.SRF    (1,1) double = 1000000;    % Scaling factor discounting
    params.scale1 (1,1) double = 0.00891061; % Multiplicative scaling coefficient
    params.scale2 (1,1) double = -6275.91;   % Additive scaling coefficient

    % Parameters for non-industrial emission
    % Assumes abateable share of non-CO2 GHG is 65%
    params.eland0         (1,1) double = 5.9;     % Carbon emissions from land 2015 (GtCO2 per year)
    params.deland         (1,1) double = 0.1;     % Decline rate of land emissions (per period)
    params.F_Misc2020     (1,1) double = -0.054;  % Non-abatable forcings 2020F
    params.F_Misc2100     (1,1) double = 0.265;   % Non-abatable forcings 2100
    params.F_GHGabate2020 (1,1) double = 0.518;   % Forcings of abatable nonCO2 GHG
    params.F_GHGabate2100 (1,1) double = 0.957;   % Forcings of abatabx0.CAle nonCO2 GHG

    params.ECO2eGHGB2020  (1,1) double = 9.96;    % Emis of abatable nonCO2 GHG GtCO2e  2020
    params.ECO2eGHGB2100  (1,1) double = 15.5;    % Emis of abatable nonCO2 GHG GtCO2e  2100
    params.emissrat2020   (1,1) double = 1.40;    % Ratio of CO2e to industrial CO2 2020
    params.emissrat2100   (1,1) double = 1.21;    % Ratio of CO2e to industrial CO2 2020
    params.Fcoef1         (1,1) double = 0.00955; % Coefficient of nonco2 abateable emissions
    params.Fcoef2         (1,1) double = 0.861;   % Coefficient of nonco2 abateable emissions

    % FAIR Parameters
    params.yr0 (1,1) double = 2020;         % Calendar year that corresponds to model year zero
    
    params.emshare0 (1,1) double = 0.2173;  % Carbon emissions share into Reservoir 0
    params.emshare1 (1,1) double = 0.224;   % Carbon emissions share into Reservoir 1
    params.emshare2 (1,1) double = 0.2824;  % Carbon emissions share into Reservoir 2
    params.emshare3 (1,1) double = 0.2763;  % Carbon emissions share into Reservoir 3

    params.tau0 (1,1) double =     1000000; % Decay time constant for R0  (year)
    params.tau1 (1,1) double =     394.4;   % Decay time constant for R1  (year)
    params.tau2 (1,1) double =     36.53;   % Decay time constant for R2  (year)
    params.tau3 (1,1) double =     4.304;   % Decay time constant for R3  (year)

    params.teq1 (1,1) double =    0.324;    % Thermal equilibration parameter for box 1 (m^2 per KW)
    params.teq2 (1,1) double =    0.44;     % Thermal equilibration parameter for box 2 (m^2 per KW)
    params.d1   (1,1) double =     236;     % Thermal response timescale for deep ocean (year)
    params.d2   (1,1) double =     4.07;    % Thermal response timescale for upper ocean (year)

    params.irf0   (1,1) double = 32.4;    % Pre-industrial IRF100 (year)
    params.irC    (1,1) double =  0.019;  % Increase in IRF100 with cumulative carbon uptake (years per GtC)
    params.irT    (1,1) double =  4.165;  % Increase in IRF100 with warming (years per degree K)
    params.fco22x (1,1) double =  3.93;   % Forcings of equilibrium CO2 doubling (Wm-2)

    % INITIAL CONDITIONS TO BE CALIBRATED TO HISTORY
    params.mat0   (1,1) double = 886.5128014; % Initial concentration in atmosphere in 2020 (GtC)
    params.res00  (1,1) double = 150.093;     % Initial concentration in Reservoir 0 in 2020 (GtC)
    params.res10  (1,1) double = 102.698;     % Initial concentration in Reservoir 1 in 2020 (GtC)
    params.res20  (1,1) double = 39.534;      % Initial concentration in Reservoir 2 in 2020 (GtC)
    params.res30  (1,1) double = 6.1865;      % Initial concentration in Reservoir 3 in 2020 (GtC)

    params.mateq  (1,1) double =  588;     % Equilibrium concentration atmosphere  (GtC)
    params.tbox10 (1,1) double = 0.1477;   % Initial temperature box 1 change in 2020 (C from 1765)
    params.tbox20 (1,1) double = 1.099454; % Initial temperature box 2 change in 2020 (C from 1765)
    params.tatm0  (1,1) double = 1.24715;  % Initial atmospheric temperature change in 2020

    % Variable Bounds 
    params.SLower (1,1) double = 0;
    params.SUpper (1,1) double = 0.28;
    params.AlphaUpperBound (1,1) double = 100;
    params.AlphaLowerBound (1,1) double = 0.1;
end

%% Derived PARAMETERS
params.rartp = exp( params.prstp + params.betaclim*params.pi)-1;  % Risk-adjusted rate of time preference exp(rho*(T=0))-1
params.sig1  = params.e1/(params.q1*(1-params.miu1));             % Carbon intensity 2020 kgCO2-output 2020

L = params.pop1*ones(np, 1);      % Level of population and labor
aL = params.AL1*ones(np, 1);      % Level of total factor productivity
sigma = params.sig1*ones(np, 1);  % CO2-emissions output ratio
sigmatot = zeros(np, 1);          % Emissions output ratio for CO2e
gA = zeros(np, 1);                % Growth rate of productivity from

gsig = zeros(np, 1);           % Change in sigma (rate of decarbonization)
eland          = zeros(np, 1); % Emissions from deforestation (GtCO2 per year)
cost1tot       = zeros(np, 1); % Abatement cost adjusted for backstop and sigma
PBACKTIME = zeros(np, 1);      % Backstop price 2019$ per ton CO2

cpricebase = zeros(np, 1); % Carbon price in base case
%         ppm(t)         Atmospheric concentrations parts per million
%         atfrac2020(t)  Atmospheric share since 2020
%         atfrac1765(t)  Atmospheric fraction of emissions since 1765
%         abaterat(t)    Abatement cost per net output

% Precautionary dynamic parameters
varpcc   = zeros(np, 1); % Variance of per capita consumption
rprecaut = zeros(np, 1); % Precautionary rate of return
RR1      = zeros(np, 1); % STP with precautionary factor
RR       = zeros(np, 1); % STP factor without precautionary factor;

% Parameters emissions and non-CO2
CO2E_GHGabateB = zeros(np, 1); % Abateable non-CO2 GHG emissions base
F_Misc         = zeros(np, 1); % Non-abateable forcings (GHG and other)
emissrat       = zeros(np, 1); % Ratio of CO2e to industrial emissions
%         FORC_CO2(t)               CO2 Forcings

%% Time preference for climate investments and precautionary effect

% Optimal long-run savings rate used for transversality
optlrsav        =(params.dk + .004)/(params.dk + .004*params.elasmu + params.rartp)*params.gama;
for t = 1:np
    % Precautionary dynamic parameters
    varpcc(t)     =  min(params.siggc1^2*5*(t-1),params.siggc1^2*5*47);
    rprecaut(t)   = -0.5*varpcc(t)*params.elasmu^2;
    RR1(t)        = 1/((1+params.rartp)^(params.tstep*(t-1)));
    RR(t)         = RR1(t)*(1+rprecaut(t))^(-params.tstep*(t-1));

    % Time preference for climate investments and precautionary effect
    gA(t)         = params.gA1*exp(-params.delA*5*(t-1));
    cpricebase(t) = params.cprice1*(1+params.gcprice)^(5*(t-1));
    if t <= 7
        PBACKTIME(t) = params.pback2050*exp(-0.05*(t-7));
    else
        PBACKTIME(t) = params.pback2050*exp(-0.005*(t-7));
    end
    gsig(t)       = min(params.gsigma1*params.delgsig^(t-1),params.asymgsig);
    if t < np
        L(t+1)          = L(t)*(params.popasym/L(t))^params.popadj;
        aL(t+1)         = aL(t)/(1-gA(t));
        sigma(t+1)      = sigma(t)*exp(5*gsig(t));
    end

    % Parameters emissions and non-CO2
    eland(t)          = params.eland0*(1-params.deland)^(t-1);
    if t <= 16
        CO2E_GHGabateB(t) = params.ECO2eGHGB2020+((params.ECO2eGHGB2100-params.ECO2eGHGB2020)/16)*(t-1);
        F_Misc(t)         = params.F_Misc2020 +((params.F_Misc2100-params.F_Misc2020)/16)*(t-1);
        emissrat(t) = params.emissrat2020 +((params.emissrat2100-params.emissrat2020)/16)*(t-1);
    else
        CO2E_GHGabateB(t) = params.ECO2eGHGB2100;
        F_Misc(t)         = params.F_Misc2100;
        emissrat(t)       = params.emissrat2100;
    end
    sigmatot(t) = sigma(t)*emissrat(t);
    cost1tot(t) = PBACKTIME(t)*sigmatot(t)/params.expcost2/1000;
end

%% Emissions limits
% Upper bound on miu
miuup = zeros(np,1);
miuup(1) = 0.05;
miuup(2) = 0.10;

idx = 1:np;
miuup(idx > 2)  = params.delmiumax*(idx(3:end)-1);
miuup(idx > 8)  = 0.85+.05*(idx(9:end)-8);
miuup(idx > 11) = params.limmiu2070;
miuup(idx > 20) = params.limmiu2120;
miuup(idx > 37) = params.limmiu2200;
miuup(idx > 57) = params.limmiu2300;

%% Capital Limits
% sLBounds = zeros(np,1);
sLBounds = params.SLower*ones(np,1);
sLBounds(38:end) = 0.28;
sUBounds = inf(np,1);
sUBounds(38:end) = 0.28;

%% Collect parameters
params.L = L;
params.F_Misc = F_Misc;
params.RR = RR;
params.aL = aL;
params.gA = gA;
params.CO2E_GHGabateB = CO2E_GHGabateB;
params.sigmatot = sigmatot;
params.cost1tot = cost1tot;
params.sigma = sigma;
params.eland = eland;
params.miuup = miuup;
params.sLBounds = sLBounds;
params.sUBounds = sUBounds;
params.eco2Param = aL.*(L/1000).^(1-params.gama);
params.pbacktime = PBACKTIME;
params.optlrsav = optlrsav;

%% Solve for Alpha0
proba0 = eqnproblem();
opts = optimset('fzero');
opts.Display = 'off';
a0 = optimvar('a0',1);
proba0.Equations = params.irf0 + params.irC*(params.CumEmiss0-(params.mat0-params.mateq))+params.irT*params.tatm0 == ...
    (a0*params.emshare0*params.tau0*(1-exp(-100/(a0*params.tau0)))) + ...
    (a0*params.emshare1*params.tau1*(1-exp(-100/(a0*params.tau1))))+ ...
    (a0*params.emshare2*params.tau2*(1-exp(-100/(a0*params.tau2))))+ ...
    (a0*params.emshare3*params.tau3*(1-exp(-100/(a0*params.tau3))));
a0 = solve(proba0, struct('a0',0.5), Options = opts);
a0 = a0.a0;
params.a0 = a0;