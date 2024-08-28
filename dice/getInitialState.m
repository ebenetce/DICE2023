function [RES0, RES1, RES2, RES3, TBOX1, TBOX2, CCATOT, F_GHGabate, K, I, C] = getInitialState(params, S, MIU)

RES0 = params.res00;
RES1 = params.res10;
RES2 = params.res20;
RES3 = params.res30;

TBOX1 = params.tbox10;
TBOX2 = params.tbox20;

CCATOT = params.CumEmiss0;

F_GHGabate = params.F_GHGabate2020;

K = params.k0;

TATM = params.tatm0;

DAMFRAC = params.a1*TATM+params.a2base*TATM.^params.a3;

YGROSS = params.eco2Param(1)*(K^params.gama);
YNET = YGROSS*(1-DAMFRAC);
ABATECOST = YGROSS*params.cost1tot(1)*(MIU^params.expcost2);
Y = YNET - ABATECOST;
I = S*Y;
C = Y - I;

end