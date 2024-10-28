function prob = diceFunc(np, params, nvp)

arguments
    np (1,1) double {mustBeInteger, mustBePositive}
    params (1,1) struct
    nvp.MIU (:,1) double = [];
    nvp.S (:,1) double = [];
    nvp.Alpha (:,1) double = [];
    nvp.TempUpperConstraint (1,1) double = 20;
    nvp.TempLowerConstraint (1,1) double = 0.5;
end

%% VARIABLES
% Control variables
% Emission control rate GHGs
if isempty (nvp.MIU)
    MIU = optimvar('MIU',np, 'LowerBound', params.MIULowerBound, 'UpperBound', params.miuup(1:np)); % Emission control rate GHGs
else
    MIU = nvp.MIU;
end

% Gross savings rate as fraction of gross world product
if isempty(nvp.S)
    S   = optimvar('S',np, 'LowerBound', params.sLBounds, 'UpperBound', params.sUBounds);
else
    S = nvp.S;
end

% Carbon decay time scaling factor
if isempty(nvp.Alpha)
    alpha = optimvar('alpha',np, 'LowerBound', params.AlphaLowerBound, 'UpperBound', params.AlphaUpperBound);
else
    alpha = nvp.Alpha;
end

L = params.L;
RR = params.RR;

scale1 = params.scale1;
scale2 = params.scale2;
elasmu = params.elasmu;

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
mat0 = params.mat0;
gama = params.gama;

tstep = params.tstep;

% Define equations
C = optimexpr(np, 1);
TATM = params.tatm0*ones(np,1);

[RES0, RES1, RES2, RES3, TBOX1, TBOX2, CCATOT, F_GHGabate, K, I, C(1,1)] = getInitialState(params, S(1), MIU(1));

[C, CCATOT,MAT,TATM] = fcn2optimexpr(@diceLoop, np,CCATOT,sigma,eco2Param, ...
    K,gama,eland,MIU,dk,tstep,I,F_GHGabate,Fcoef2,Fcoef1,CO2E_GHGabateB, ...
    RES0,emshare0,tau0,alpha,RES1,emshare1,tau1,RES2,emshare2,tau2,RES3, ...
    emshare3,tau3,mateq,TBOX1,d1,teq1,fco22x,F_Misc,TBOX2,d2,teq2,TATM, ...
    a1,a2base,a3,cost1tot,expcost2,S,C, mat0);

PERIODU = ((C*1000./L).^(1-elasmu)-1)./(1-elasmu)-1;
TOTPERIODU = PERIODU.*L.*RR;
UTILITY = tstep * scale1 * sum( TOTPERIODU ) + scale2;

%% Define problem
prob = optimproblem('Objective', UTILITY, 'ObjectiveSense', 'max');
prob.Constraints.IRFeqLHS = params.irf0+params.irC*(CCATOT-(MAT-mateq))+params.irT*TATM  ==  ...
    alpha*emshare0*tau0.*(1-exp(-100./(alpha*tau0))) + ...
    alpha*emshare1*tau1.*(1-exp(-100./(alpha*tau1))) + ...
    alpha*emshare2*tau2.*(1-exp(-100./(alpha*tau2))) + ...
    alpha*emshare3*tau3.*(1-exp(-100./(alpha*tau3)));
prob.Constraints.TATMup = TATM <= nvp.TempUpperConstraint;
prob.Constraints.TATMlow = TATM >= nvp.TempLowerConstraint;

end

function [C, CCATOTv,MATv,TATMv] = diceLoop(np,CCATOT,sigma,eco2Param,K,gama,eland,MIU,dk,tstep,I,F_GHGabate,Fcoef2,Fcoef1,CO2E_GHGabateB,RES0,emshare0,tau0,alpha,RES1,emshare1,tau1,RES2,emshare2,tau2,RES3,emshare3,tau3,mateq,TBOX1,d1,teq1,fco22x,F_Misc,TBOX2,d2,teq2,TATM,a1,a2base,a3,cost1tot,expcost2,S,C, mat0)

CCATOTv = CCATOT*ones(np,1);
MATv = mat0*ones(np,1);
TATMv = TATM;
for i = 2 : np
    [C(i,1), CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2, MAT, TATM] = diceForward(i, MIU, S, alpha, CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2, ...
    sigma, eco2Param, eland, mateq, fco22x, F_Misc, Fcoef1, Fcoef2, CO2E_GHGabateB, d1, d2, dk, a1, a2base, a3, cost1tot, expcost2, teq1, teq2, tstep, emshare0, emshare1, emshare2, emshare3, ...
    tau0, tau1, tau2, tau3, gama);
    
    CCATOTv(i,1) = CCATOT;
    MATv(i,1) = MAT;
    TATMv(i,1) = TATM;

end

end