function [C, CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2, MAT, TATM] = diceForward(i, MIU, S, alpha, CCATOT, K, I, F_GHGabate, RES0, RES1, RES2, RES3, TBOX1, TBOX2, ...
    sigma, eco2Param, eland, mateq, fco22x, F_Misc, Fcoef1, Fcoef2, CO2E_GHGabateB, d1, d2, dk, a1, a2base, a3, cost1tot, expcost2, teq1, teq2, tstep, emshare0, emshare1, emshare2, emshare3, ...
    tau0, tau1, tau2, tau3, gama)

% Define equations
CCATOT = CCATOT + ((sigma(i-1)*(eco2Param(i-1)*(K^gama)) + eland(i-1))*(1-MIU(i-1)))*5/3.666;

K = (1-dk)^tstep*K + (tstep)*I;

F_GHGabate = Fcoef2*F_GHGabate + Fcoef1*CO2E_GHGabateB(i-1)*(1-(MIU(i-1)));

RES0 =  (emshare0*tau0*alpha(i)*( (sigma(i)*eco2Param(i)*(K^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau0*alpha(i))))+RES0*exp(-tstep/(tau0*alpha(i))); % Reservoir 0 law of motion
RES1 =  (emshare1*tau1*alpha(i)*( (sigma(i)*eco2Param(i)*(K^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau1*alpha(i))))+RES1*exp(-tstep/(tau1*alpha(i))); % Reservoir 1 law of motion
RES2 =  (emshare2*tau2*alpha(i)*( (sigma(i)*eco2Param(i)*(K^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau2*alpha(i))))+RES2*exp(-tstep/(tau2*alpha(i))); % Reservoir 2 law of motion
RES3 =  (emshare3*tau3*alpha(i)*( (sigma(i)*eco2Param(i)*(K^gama) + eland(i))*(1-MIU(i))/3.667) )*(1-exp(-tstep/(tau3*alpha(i))))+RES3*exp(-tstep/(tau3*alpha(i))); % Reservoir 3 law of motion

MAT = mateq+RES0+RES1+RES2+RES3;

TBOX1 =  TBOX1*exp(-tstep/d1) + teq1*(fco22x*log(MAT/mateq)/log(2)+F_Misc(i)+F_GHGabate)*(1-exp(-tstep/d1));
TBOX2 =  TBOX2*exp(-tstep/d2) + teq2*(fco22x*log(MAT/mateq)/log(2)+F_Misc(i)+F_GHGabate)*(1-exp(-tstep/d2));

TATM  = TBOX1+TBOX2;

DAMFRAC = a1*TATM+a2base*TATM^a3;

YGROSS = eco2Param(i)*(K^gama);
ABATECOST = YGROSS*cost1tot(i)*(MIU(i)^expcost2);

YNET = YGROSS*(1-DAMFRAC);
Y = YNET - ABATECOST;
I = S(i)*Y;

C = Y - I;
end