function [C, CCATOTv,MATv,TATMv] = diceLoop(np,CCATOT,sigma,eco2Param,K,gama,eland,MIU,dk,tstep,I,F_GHGabate,Fcoef2,Fcoef1,CO2E_GHGabateB,RES0,emshare0,tau0,alpha,RES1,emshare1,tau1,RES2,emshare2,tau2,RES3,emshare3,tau3,mateq,TBOX1,d1,teq1,fco22x,F_Misc,TBOX2,d2,teq2,TATM,a1,a2base,a3,cost1tot,expcost2,S,C, mat0)

CCATOTv = CCATOT*ones(np,1);
MATv = mat0*ones(np,1);
TATMv = TATM;
for i = 2 : np
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
    C(i,1) = Y - I;

    CCATOTv(i,1) = CCATOT;
    MATv(i,1) = MAT;
    TATMv(i,1) = TATM;

end

end