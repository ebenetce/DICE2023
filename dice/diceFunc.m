function [UTILITY, IRFtl, ITFtr] = diceFunc(x, params)

np = numel(x)/3;
MIU = x(1:np);
S = x(np+1:np*2);
alpha = x(2*np+1: end);

L = params.L;
RR = params.RR;

scale1 = params.scale1;
scale2 = params.scale2;
elasmu = params.elasmu;

mateq = params.mateq;
irf0 = params.irf0;
irC = params.irC;
irT = params.irT;

emshare0 = params.emshare0;
emshare1 = params.emshare1;
emshare2 = params.emshare2;
emshare3 = params.emshare3;

tau0 = params.tau0;
tau1 = params.tau1;
tau2 = params.tau2;
tau3 = params.tau3;


tstep = params.tstep;

% Define equations
C = zeros(np,1);
IRFtl = zeros(np,1);
ITFtr = zeros(np,1);
for i = 1 : np
    if i == 1
        [RES0, RES1, RES2, RES3, TBOX1, TBOX2, CCATOT, F_GHGabate, K, I, C(1,1)] = getInitialState(params, S(1), MIU(1));
        alpha(1) = params.a0;
    else
        [C(i,1), CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2, MAT, TATM] = diceForward(i, params, MIU, S, alpha, CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2);
        IRFtl(i) =  irf0+irC*(CCATOT-(MAT-mateq))+irT*TATM;
        ITFtr(i) =  alpha(i)*emshare0*tau0.*(1-exp(-100./(alpha(i)*tau0)))+ alpha(i)*emshare1*tau1.*(1-exp(-100./(alpha(i)*tau1))) + alpha(i)*emshare2*tau2.*(1-exp(-100./(alpha(i)*tau2))) + alpha(i)*emshare3*tau3.*(1-exp(-100./(alpha(i)*tau3)));
    end
end

PERIODU = ((C*1000./L).^(1-elasmu)-1)./(1-elasmu)-1;
TOTPERIODU = PERIODU.*L.*RR;
UTILITY = tstep * scale1 * sum( TOTPERIODU ) + scale2;