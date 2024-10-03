function [UTILITY, C, K] = diceTrajectory(params, np, MIU, S, alpha)

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
    [C(i,1), CCATOT, K(i,1), I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2] = diceForward(i, params, MIU, S, alpha, CCATOT, K(i-1,1), I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2);
end

PERIODU = ((C*1000./L).^(1-elasmu)-1)./(1-elasmu)-1;
TOTPERIODU = PERIODU.*L.*RR;
UTILITY = tstep * scale1 * sum( TOTPERIODU ) + scale2;

end