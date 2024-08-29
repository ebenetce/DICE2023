function allVars = recoverAllVars(sol, params)
np = numel(sol)/3;

MIU = sol(1:np);
S = sol(np+1:2*np);
alpha = sol(2*np+1:3*np);

% Unpack params
sigma = params.sigma;
eco2Param = params.eco2Param;
eland = params.eland;
mateq = params.mateq;
fco22x = params.fco22x;
F_Misc = params.F_Misc;
irf0 = params.irf0;
irC = params.irC;
irT = params.irT;
Fcoef1 = params.Fcoef1;
Fcoef2 = params.Fcoef2;
CO2E_GHGabateB = params.CO2E_GHGabateB;

d1 = params.d1;
d2 = params.d2;
dk = params.dk;
aL = params.aL;

L = params.L;
RR = params.RR;
a1 = params.a1;
a2base = params.a2base;
a3 = params.a3;
cost1tot = params.cost1tot;
expcost2 = params.expcost2;
% scale1 = params.scale1;
% scale2 = params.scale2;
elasmu = params.elasmu;

teq1 = params.teq1;
teq2 = params.teq2;

tstep = params.tstep;

emshare0 = params.emshare0;
emshare1 = params.emshare1;
emshare2 = params.emshare2;
emshare3 = params.emshare3;

tau0 = params.tau0;
tau1 = params.tau1;
tau2 = params.tau2;
tau3 = params.tau3;

gama = params.gama;
pbacktime = params.pbacktime;
SRF = params.SRF;

years = 2020 + (params.tstep:params.tstep:params.tstep*np)';

CCATOT = zeros(np,1); CCATOT(1) = params.CumEmiss0;
K = zeros(np,1); K(1) = params.k0;
I = zeros(np,1);
Y = zeros(np,1);
YNET = zeros(np,1);
ABATECOST = zeros(np,1);
YGROSS = zeros(np,1);
DAMFRAC = zeros(np,1);
F_GHGabate = zeros(np,1); F_GHGabate(1) = params.F_GHGabate2020;
RES0 = zeros(np,1); RES0(1) = params.res00;
RES1 = zeros(np,1); RES1(1) = params.res10;
RES2 = zeros(np,1); RES2(1) = params.res20;
RES3 = zeros(np,1); RES3(1) = params.res30;
MAT = zeros(np,1); MAT(1) = params.mat0;
TBOX1 = zeros(np,1); TBOX1(1) = params.tbox10;
TBOX2 = zeros(np,1); TBOX2(1) = params.tbox20;
TATM = zeros(np,1); TATM(1) = params.tatm0;
IRFt = zeros(np,1); IRFt(1) = irf0+irC*(CCATOT(1)-(MAT(1)-mateq))+irT*TATM(1);

DAMFRAC(1) = a1*TATM(1)+a2base*TATM(1).^a3;
YGROSS(1) = eco2Param(1)*(K(1)^gama);
YNET(1) = YGROSS(1)*(1-DAMFRAC(1));
ABATECOST(1) = YGROSS(1)*cost1tot(1)*(MIU(1)^expcost2);
Y(1) = YNET(1) - ABATECOST(1);
I(1) = S(1)*Y(1);

for i = 2 : np
    CCATOT(i) = CCATOT(i-1) + ((sigma(i-1)*(eco2Param(i-1)*(K(i-1)^gama)) + eland(i-1))*(1-MIU(i-1)))*5/3.666;

    K(i) = (1-dk)^tstep*K(i-1) + tstep*I(i-1);

    F_GHGabate(i) = Fcoef2*F_GHGabate(i-1) + Fcoef1*CO2E_GHGabateB(i-1)*(1-(MIU(i-1)));

    RES0(i) =  (emshare0*tau0*alpha(i)*( (sigma(i)*eco2Param(i)*(K(i)^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau0*alpha(i))))+RES0(i-1)*exp(-tstep/(tau0*alpha(i))); % Reservoir 0 law of motion
    RES1(i) =  (emshare1*tau1*alpha(i)*( (sigma(i)*eco2Param(i)*(K(i)^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau1*alpha(i))))+RES1(i-1)*exp(-tstep/(tau1*alpha(i))); % Reservoir 1 law of motion
    RES2(i) =  (emshare2*tau2*alpha(i)*( (sigma(i)*eco2Param(i)*(K(i)^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau2*alpha(i))))+RES2(i-1)*exp(-tstep/(tau2*alpha(i))); % Reservoir 2 law of motion
    RES3(i) =  (emshare3*tau3*alpha(i)*( (sigma(i)*eco2Param(i)*(K(i)^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau3*alpha(i))))+RES3(i-1)*exp(-tstep/(tau3*alpha(i))); % Reservoir 3 law of motion

    MAT(i) = mateq+RES0(i)+RES1(i)+RES2(i)+RES3(i);

    TBOX1(i) =  TBOX1(i-1)*exp(-tstep/d1) + teq1*(fco22x*log(MAT(i)/mateq)/log(2)+F_Misc(i)+F_GHGabate(i))*(1-exp(-tstep/d1));
    TBOX2(i) =  TBOX2(i-1)*exp(-tstep/d2) + teq2*(fco22x*log(MAT(i)/mateq)/log(2)+F_Misc(i)+F_GHGabate(i))*(1-exp(-tstep/d2));

    TATM(i)  = TBOX1(i)+TBOX2(i);

    DAMFRAC(i) = a1*TATM(i)+a2base*TATM(i)^a3;

    YGROSS(i) = eco2Param(i)*(K(i)^gama);
    ABATECOST(i) = YGROSS(i)*cost1tot(i)*(MIU(i)^expcost2);

    YNET(i) = YGROSS(i)*(1-DAMFRAC(i));
    Y(i) = YNET(i) - ABATECOST(i);
    I(i) = S(i)*Y(i);
    C(i) = Y(i) - I(i);

    IRFt(i) = irf0+irC*(CCATOT(i)-(MAT(i)-mateq))+irT*TATM(i);
end

CACC = CCATOT-(MAT-mateq);                    % Accumulated carbon in sinks equation
DAMFRAC = (a1*TATM)+(a2base*TATM.^a3) ;       % Equation for damage fraction
YGROSS = aL.*(L/1000).^(1-gama).*(K.^gama);       % Output gross equation
ABATECOST = YGROSS.*(cost1tot.*(MIU.^expcost2));  % Cost of emissions reductions equation
DAMAGES = YGROSS.*DAMFRAC;                            % Damage equation
YNET = YGROSS.*(1-DAMFRAC);                           % Output net of damages equation
Y = YNET - ABATECOST;                                 % Output net equation
I = S.*Y;                                         % Savings rate equation
C = Y - I;                                            % Consumption equation
CPC =  1000*C./L;                                     % Per capita consumption definition
MCABATE = pbacktime.*(MIU.^(expcost2-1));         % Equation for MC abatement
CPRICE  = pbacktime.*(MIU.^(expcost2-1));         % Carbon price equation from abatement

RFACTLONG = SRF.*ones(np,1);
RFACTLONG(2:end) = SRF*(CPC(2:end)/CPC(1)).^(-elasmu).*RR(2:end); % Long interest factor
RLONG = zeros(np,1);
RLONG(2:end) = -log(RFACTLONG(2:end)/SRF)./(5*(1:np-1)');           % Long-run interest rate equation
RSHORT = zeros(np,1);
RSHORT(2:end) = -log(RFACTLONG(2:end)./RFACTLONG(1:end-1))/5;     % Short-run interest rate equation

FORC =  fco22x*log(MAT/mateq)/log(2)+F_Misc+F_GHGabate;   % Radiative forcing equation
ECO2 = (sigma.*(eco2Param.*(K.^gama)) + eland).*(1-MIU);  % CO2 Emissions equation
FORC_CO2 = fco22x*(log((MAT/mateq))/log(2));
PERIODU = ((C*1000./L).^(1-elasmu)-1)/(1-elasmu)-1; % Instantaneous utility function equation
TOTPERIODU  = PERIODU.*L.* RR;                      % Period utility

SCC = zeros(np,1);
W0 = sum(TOTPERIODU);

dWdC = (1000*C./L).^(-elasmu).*RR;

for i = 2 : np
    params2 = params;
    params2.eland(i) = params2.eland(i)+0.1;

    [~,Cd,Kd] = diceTrajectory(params2, np, MIU, S, alpha);

    PERIODUd = ((Cd*1000./params2.L).^(1-elasmu)-1)/(1-elasmu)-1; % Instantaneous utility function equation
    TOTPERIODUd  = PERIODUd.*params2.L.* params2.RR;                      % Period utility

    W1 = sum(TOTPERIODUd);

    YGROSSd = (aL(i)*(L(i)/1000)^(1-gama))*(Kd(i)^gama);
    ECO2d = (sigma(i)*YGROSSd + params2.eland(i))*(1-(MIU(i)));
    SCC(i) = -(W1-W0)/(ECO2d - ECO2(i))/dWdC(i);
end
SCC(1) = SCC(2)*0.85;

ABATERAT = ABATECOST./Y;

allVars = table(years, MAT, TATM, FORC, F_Misc, ...
    F_GHGabate, CPRICE, MIU, SCC, Y, RSHORT, L, aL, YGROSS,  ...
    params.gA, K, S, I, YNET, CPC, C, DAMFRAC, DAMAGES, ABATECOST, ...
    ABATERAT, sigma, params.sigmatot, pbacktime, ECO2, ...
    params.CO2E_GHGabateB, eland, CCATOT, ...
    FORC_CO2, RES0, RES1,RES2,RES3, TBOX1, TBOX2, alpha, IRFt, CACC, ...
    RLONG, ...
    VariableNames= ...
    ["Year", ...
    "Atmospheric concentrations GtC", ...
    "Atmospheric temperature (deg c above preind)", ...
    "Total forcings w/m2", ...
    "Forcings, exogenous w/m2", ...
    "Actual other abatable GHG forcings w/m2", ...
    "Carbon price (2019 $ per t CO2)" , ...
    "Emissions control rate", ...
    "Social cost of carbon $/tCO2", ...
    "Output, net net trill 2019$", ...
    "Short Interest rate, %/yr", ...
    "Population", ...
    "TFP", ...
    "Output, gross-gross, 2019$", ...
    "Change TFP, %/year", ...
    "Capital stock, 2019$" , ...
    "Savings rate, fraction gross output", ...
    "Gross investment, 2019$", ...
    "Y gross-net, 2019$", ...
    "Consumption per capita, 2019$" , ...
    "Consumption", ...
    "Climate damages, fraction of output", ...
    "Damages, 2019$", ...
    "Abatement, 2019$", ...
    "Abatement/Output", ...
    "Sigmabase (CO2/output, no controls, industrial CO2)", ...
    "Sigmatot,(CO2/output, no controls, all CO2)", ...
    "Cost, backstop technology ($/tCO2)", ...
    "Total CO2 Emissions, GTCO2/year", ...
    "Base abateable non-CO2 emission, GTCO2-E/year", ...
    "Land emissions, GtCO2/year", ...
    "Cumulative CO2 emissions, GtC", ...
    "CO2 forcings w/m2", ...
    "Permanent C box", ...
    "Slow C box", ...
    "Medium C box", ...
    "Fast C box", ...
    "Temp Box 1", ...
    "Temp Box 2", ...
    "Alpha", ...
    "IFR", ...
    "cacc", ...
    "Long interest rate(%)"]);

[~,idx]=sort(upper(allVars.Properties.VariableNames));
allVars = allVars(:,idx);

end