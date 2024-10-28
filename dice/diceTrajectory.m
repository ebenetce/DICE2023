function [UTILITY, C, K] = diceTrajectory(params, np, MIU, S, alpha)

% Unpack params
sigma = params.sigma;
eco2Param = params.eco2Param;
eland = params.eland;
mateq = params.mateq;
fco22x = params.fco22x;
F_Misc = params.F_Misc;

Fcoef1 = params.Fcoef1;
Fcoef2 = params.Fcoef2;
CO2E_GHGabateB = params.CO2E_GHGabateB;

d1 = params.d1;
d2 = params.d2;
dk = params.dk;

a1 = params.a1;
a2base = params.a2base;
a3 = params.a3;
cost1tot = params.cost1tot;
expcost2 = params.expcost2;

teq1 = params.teq1;
teq2 = params.teq2;

emshare0 = params.emshare0;
emshare1 = params.emshare1;
emshare2 = params.emshare2;
emshare3 = params.emshare3;

tau0 = params.tau0;
tau1 = params.tau1;
tau2 = params.tau2;
tau3 = params.tau3;

gama = params.gama;

L = params.L;
RR = params.RR;

scale1 = params.scale1;
scale2 = params.scale2;
elasmu = params.elasmu;

tstep = params.tstep;

% Define equations
C = zeros(np,1);
K = zeros(np,1);

[RES0, RES1, RES2, RES3, TBOX1, TBOX2, CCATOT, F_GHGabate, K(1,1), I, C(1,1)] = getInitialState(params, S(1), MIU(1));
alpha(1) = params.a0;

for i = 2 : np
    [C(i,1), CCATOT, K(i,1), I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2] = diceForward(i, MIU, S, alpha, CCATOT, K(i-1,1), I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2, ...
        sigma, eco2Param, eland, mateq, fco22x, F_Misc, Fcoef1, Fcoef2, CO2E_GHGabateB, d1, d2, dk, a1, a2base, a3, cost1tot, expcost2, teq1, teq2, tstep, emshare0, emshare1, emshare2, emshare3, ...
    tau0, tau1, tau2, tau3, gama);
end

PERIODU = ((C*1000./L).^(1-elasmu)-1)./(1-elasmu)-1;
TOTPERIODU = PERIODU.*L.*RR;
UTILITY = tstep * scale1 * sum( TOTPERIODU ) + scale2;

end