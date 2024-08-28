function UTILITY = diceFunc(x, params)

np = numel(x)/3;
MIU = x(1:np);
S = x(np+1:np*2);
alpha = x(2*np+1: end);

L = params.L;
RR = params.RR;

scale1 = params.scale1;
scale2 = params.scale2;
elasmu = params.elasmu;

tstep = params.tstep;

% Define equations
C = zeros(np,1);
for i = 1 : np
    if i == 1
        [RES0, RES1, RES2, RES3, TBOX1, TBOX2, CCATOT, F_GHGabate, K, I, C(1,1)] = getInitialState(params, S(1), MIU(1));
        alpha(1) = params.a0;
    else
        [C(i,1), CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2] = diceForward(i, params, MIU, S, alpha, CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2);
    end
end

PERIODU = ((C*1000./L).^(1-elasmu)-1)./(1-elasmu)-1;
TOTPERIODU = PERIODU.*L.*RR;
UTILITY = tstep * scale1 * sum( TOTPERIODU ) + scale2;