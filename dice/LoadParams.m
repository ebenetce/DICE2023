function params = LoadParams(np)
% Population and technology

gama    = 0.300;  % Capital elasticity in production function
pop1    = 7752.9; % Initial world population 2020 (millions)
popadj  = 0.145;  % Growth rate to calibrate to 2050 pop projection
popasym = 10825;  % Asymptotic population (millions)
dk      = 0.100;  % Depreciation rate on capital (per year)
q1      = 135.7;  % Initial world output 2020 (trill 2019 USD)
AL1     = 5.84;   % Initial level of total factor productivity
gA1     = 0.066;  % Initial growth rate for TFP per 5 years
delA    = 0.0015; % Decline rate of TFP per 5 yearst
% Emissions parameters and Non-CO2 GHG with sigma = emissions/output

gsigma1   = -0.015; % Initial growth of sigma (per year)
delgsig   =  0.96;  % Decline rate of gsigma per period
asymgsig  = -0.005; % Asympototic gsigma
e1        = 37.56;  % Industrial emissions 2020 (GtCO2 per year)
miu1      =  0.05;  % Emissions control rate historical 2020
fosslim   =  6000;  % Maximum cumulative extraction fossil fuels (GtC)
CumEmiss0 =  633.5; % Cumulative emissions 2020 (GtC)
% Climate damage parameters

a1     = 0;        % Damage intercepte
a2base = 0.003467; % Damage quadratic term rev 01-13-23
a3     = 2.00;     % Damage exponent
% Abatement cost

expcost2  = 2.6;    % Exponent of control cost function
pback2050 = 515;    % Cost of backstop 2019$ per tCO2 2050
gback     = -0.012; % Initial cost decline backstop cost per year
cprice1   = 6;      % Carbon price 2020 2019$ per tCO2
gcprice   = 0.025;  % Growth rate of base carbon price per year
% Limits on emissions controls

limmiu2070 = 1.0;  % Emission control limit from 2070
limmiu2120 = 1.1;  % Emission control limit from 2120
limmiu2200 = 1.05; % Emission control limit from 2220
limmiu2300 = 1.0;  % Emission control limit from 2300
delmiumax  = 0.12; % Emission control delta limit per period
% Preferences, growth uncertainty, and timing

betaclim  = 0.5;   % Climate beta
elasmu    = 0.95;  % Elasticity of marginal utility of consumption
prstp     = 0.001; % Pure rate of social time preference rho
pi        = 0.05;  % Capital risk premium
rartp     = exp( prstp + betaclim*pi)-1;  % Risk-adjusted rate of time preference exp(rho*(T=0))-1
k0        = 295;   % Initial capital stock calibrated (1012 2019 USD)
siggc1    = 0.01;  %  Annual standard deviation of consumption growth
sig1      = e1/(q1*(1-miu1)); % Carbon intensity 2020 kgCO2-output 2020
% Scaling so that MU(C(1)) = 1 and objective function = PV consumption

params.tstep = 5;           % Years per Period
params.SRF    = 1000000;    % Scaling factor discounting
params.scale1 = 0.00891061; % Multiplicative scaling coefficient
params.scale2 = -6275.91;   % Additive scaling coefficient
% Parameters for non-industrial emission
% Assumes abateable share of non-CO2 GHG is 65%

eland0         = 5.9;     % Carbon emissions from land 2015 (GtCO2 per year)
deland         = 0.1;     % Decline rate of land emissions (per period)
F_Misc2020     = -0.054;  % Non-abatable forcings 2020F
F_Misc2100     = 0.265;   % Non-abatable forcings 2100
F_GHGabate2020 = 0.518;   % Forcings of abatable nonCO2 GHG
F_GHGabate2100 = 0.957;   % Forcings of abatabx0.CAle nonCO2 GHG

ECO2eGHGB2020  = 9.96;    % Emis of abatable nonCO2 GHG GtCO2e  2020
ECO2eGHGB2100  = 15.5;    % Emis of abatable nonCO2 GHG GtCO2e  2100
emissrat2020   = 1.40;    % Ratio of CO2e to industrial CO2 2020
emissrat2100   = 1.21;    % Ratio of CO2e to industrial CO2 2020
Fcoef1         = 0.00955; % Coefficient of nonco2 abateable emissions
Fcoef2         = 0.861;   % Coefficient of nonco2 abateable emissions
% FAIR Parameters

yr0 = 2020;         % Calendar year that corresponds to model year zero
emshare0 = 0.2173;  % Carbon emissions share into Reservoir 0
emshare1 = 0.224;   % Carbon emissions share into Reservoir 1
emshare2 = 0.2824;  % Carbon emissions share into Reservoir 2
emshare3 = 0.2763;  % Carbon emissions share into Reservoir 3
tau0 =     1000000; % Decay time constant for R0  (year)
tau1 =     394.4;   % Decay time constant for R1  (year)
tau2 =     36.53;   % Decay time constant for R2  (year)
tau3 =     4.304;   % Decay time constant for R3  (year)

teq1 =    0.324;    % Thermal equilibration parameter for box 1 (m^2 per KW)
teq2 =    0.44;     % Thermal equilibration parameter for box 2 (m^2 per KW)
d1  =     236;      % Thermal response timescale for deep ocean (year)
d2  =     4.07;     % Thermal response timescale for upper ocean (year)

irf0  =   32.4;   % Pre-industrial IRF100 (year)
irC  =    0.019;  % Increase in IRF100 with cumulative carbon uptake (years per GtC)
irT   =   4.165;  % Increase in IRF100 with warming (years per degree K)
fco22x  = 3.93;   % Forcings of equilibrium CO2 doubling (Wm-2)
% PARAMETERS

L = pop1*ones(np, 1);      % Level of population and labor
aL = AL1*ones(np, 1);      % Level of total factor productivity
sigma = sig1*ones(np, 1);  % CO2-emissions output ratio
sigmatot = zeros(np, 1);   % Emissions output ratio for CO2e
gA = zeros(np, 1);         % Growth rate of productivity from
%         gL(t)          Growth rate of labor and population
%         gcost1         Growth of cost factor
gsig = zeros(np, 1);       % Change in sigma (rate of decarbonization)
eland          = zeros(np, 1);  % Emissions from deforestation (GtCO2 per year)
cost1tot       = zeros(np, 1); % Abatement cost adjusted for backstop and sigma
PBACKTIME = zeros(np, 1);  % Backstop price 2019$ per ton CO2
%         optlrsav       Optimal long-run savings rate used for transversality
%         scc(t)         Social cost of carbon
cpricebase = zeros(np, 1); % Carbon price in base case
%         ppm(t)         Atmospheric concentrations parts per million
%         atfrac2020(t)  Atmospheric share since 2020
%         atfrac1765(t)  Atmospheric fraction of emissions since 1765
%         abaterat(t)    Abatement cost per net output
%         miuup(t)       Upper bound on miu
%         gbacktime(t)   Decline rate of backstop price
% Precautionary dynamic parameters

varpcc   = zeros(np, 1); % Variance of per capita consumption
rprecaut = zeros(np, 1); % Precautionary rate of return
RR1      = zeros(np, 1); % STP with precautionary factor
RR       = zeros(np, 1); % STP factor without precautionary factor;
% Parameters emissions and non-CO2 

CO2E_GHGabateB = zeros(np, 1); % Abateable non-CO2 GHG emissions base
%         CO2E_GHGabateact(t)       Abateable non-CO2 GHG emissions base (actual)
F_Misc         = zeros(np, 1); % Non-abateable forcings (GHG and other)
emissrat       = zeros(np, 1); % Ratio of CO2e to industrial emissions
%         FORC_CO2(t)               CO2 Forcings
% ;
% Time preference for climate investments and precautionary effect

optlrsav        =(dk + .004)/(dk + .004*elasmu + rartp)*gama;
for t = 1:np
    % Precautionary dynamic parameters
    varpcc(t)     =  min(siggc1^2*5*(t-1),siggc1^2*5*47);
    rprecaut(t)   = -0.5*varpcc(t)*elasmu^2;
    RR1(t)        = 1/((1+rartp)^(params.tstep*(t-1)));
    RR(t)         = RR1(t)*(1+rprecaut(t))^(-params.tstep*(t-1));

    % Time preference for climate investments and precautionary effect
    gA(t)         = gA1*exp(-delA*5*(t-1));
    cpricebase(t) = cprice1*(1+gcprice)^(5*(t-1));
    if t <= 7
        PBACKTIME(t) = pback2050*exp(-0.05*(t-7));
    else
        PBACKTIME(t) = pback2050*exp(-0.005*(t-7));
    end
    gsig(t)       = min(gsigma1*delgsig^(t-1),asymgsig);
    if t < np
        L(t+1)          = L(t)*(popasym/L(t))^popadj;
        aL(t+1)         = aL(t)/(1-gA(t));
        sigma(t+1)      = sigma(t)*exp(5*gsig(t));
    end

    % Parameters emissions and non-CO2
    eland(t)          = eland0*(1-deland)^(t-1);
    if t <= 16
        CO2E_GHGabateB(t) = ECO2eGHGB2020+((ECO2eGHGB2100-ECO2eGHGB2020)/16)*(t-1);
        F_Misc(t)         = F_Misc2020 +((F_Misc2100-F_Misc2020)/16)*(t-1);
        emissrat(t) = emissrat2020 +((emissrat2100-emissrat2020)/16)*(t-1);
    else
        CO2E_GHGabateB(t) = ECO2eGHGB2100;
        F_Misc(t)         = F_Misc2100;
        emissrat(t)       = emissrat2100;
    end
    sigmatot(t) = sigma(t)*emissrat(t);
    cost1tot(t) = PBACKTIME(t)*sigmatot(t)/expcost2/1000;
end
% Emissions limits

miuup = zeros(np,1);
miuup(1) = 0.05;
miuup(2) = 0.10;

idx = 1:np;
miuup(idx > 2)  = delmiumax*(idx(3:end)-1);
miuup(idx > 8)  = 0.85+.05*(idx(9:end)-8);
% miuup(idx > 2)  = delmiumax*(np-1);
% miuup(idx > 8)  = 0.85+.05*(np-8);
miuup(idx > 11) = limmiu2070;
miuup(idx > 20) = limmiu2120;
miuup(idx > 37) = limmiu2200;
miuup(idx > 57) = limmiu2300;
% Capital Limits

sLBounds = zeros(np,1);
sLBounds(38:end) = 0.28;
sUBounds = inf(np,1);
sUBounds(38:end) = 0.28;
% INITIAL CONDITIONS TO BE CALIBRATED TO HISTORY

mat0   = 886.5128014; %  Initial concentration in atmosphere in 2020 (GtC)
res00  = 150.093;    % Initial concentration in Reservoir 0 in 2020 (GtC)
res10  = 102.698; % Initial concentration in Reservoir 1 in 2020 (GtC)
res20  = 39.534; % Initial concentration in Reservoir 2 in 2020 (GtC)
res30  = 6.1865; % Initial concentration in Reservoir 3 in 2020 (GtC)

mateq  =  588; %    Equilibrium concentration atmosphere  (GtC)
tbox10 = 0.1477;  % Initial temperature box 1 change in 2020 (C from 1765)
tbox20 = 1.099454;  % Initial temperature box 2 change in 2020 (C from 1765)
tatm0  = 1.24715;   % Initial atmospheric temperature change in 2020
%% 
% Collect parameters

params.L = L;
params.F_Misc = F_Misc;
params.RR = RR;
params.aL = aL;
params.gA = gA;
params.gama = gama;
params.dk = dk;
params.a1 = a1;
params.a2base = a2base;
params.a3 = a3;
params.k0 = k0;
params.fco22x = fco22x;
params.Fcoef2 = Fcoef2;
params.Fcoef1 = Fcoef1;
params.CO2E_GHGabateB = CO2E_GHGabateB;
params.sigmatot = sigmatot;
params.cost1tot = cost1tot;
params.expcost2 = expcost2;
params.sigma = sigma;
params.eland = eland;
params.F_GHGabate2020 = F_GHGabate2020;
params.CumEmiss0 = CumEmiss0;
params.d1 = d1;
params.d2 = d2;
params.teq1 = teq1;
params.teq2 = teq2;
params.emshare0 = emshare0;
params.emshare1 = emshare1;
params.emshare2 = emshare2;
params.emshare3 = emshare3;

params.tau0 = tau0;
params.tau1 = tau1;
params.tau2 = tau2;
params.tau3 = tau3;

params.irf0 = irf0;
params.irC = irC;
params.irT = irT;

params.elasmu = elasmu;

params.mat0 = mat0;
params.res00 = res00;
params.res10 = res10;
params.res20 = res20;
params.res30 = res30;

params.mateq = mateq;
params.tbox10 = tbox10;
params.tbox20 = tbox20;
params.tatm0 = tatm0;

params.miuup = miuup;
params.sLBounds = sLBounds;
params.sUBounds = sUBounds;
params.eco2Param = aL.*(L/1000).^(1-gama);

params.pbacktime = PBACKTIME;


%Solve for Alpha0
proba0 = eqnproblem();
opts = optimset('fzero');
opts.Display = 'off';
a0 = optimvar('a0',1);
proba0.Equations = irf0+irC*(CumEmiss0-(mat0-mateq))+irT*tatm0 == ((a0*emshare0*tau0*(1-exp(-100/(a0*tau0))))+(a0*emshare1*tau1*(1-exp(-100/(a0*tau1))))+(a0*emshare2*tau2*(1-exp(-100/(a0*tau2))))+(a0*emshare3*tau3*(1-exp(-100/(a0*tau3)))));
a0 = solve(proba0, struct('a0',0.5), Options = opts);
a0 = a0.a0;
params.a0 = a0;