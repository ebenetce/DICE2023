function [c,ceq] = IRFtcon(x, params)

np = numel(x)/3;
MIU = x(1:np);
S = x(np+1:np*2);
alpha = x(2*np+1:end);

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

a1 = params.a1;
a2base = params.a2base;
a3 = params.a3;
cost1tot = params.cost1tot;
expcost2 = params.expcost2;
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
%% Check eqs
np = length(MIU);

c = [];
ceq = zeros(np,1);

C = zeros(np,1);
for i = 1 : np
    if i == 1
        RES0 = params.res00;
        RES1 = params.res10;
        RES2 = params.res20;
        RES3 = params.res30;

        TBOX1 = params.tbox10;
        TBOX2 = params.tbox20;

        CCATOT = params.CumEmiss0;

        F_GHGabate = params.F_GHGabate2020;

        % MAT = params.mat0;
        
        K = params.k0;

        TATM = params.tatm0;
        % IRFt =  irf0+irC*(CCATOT-(MAT-mateq))+irT*TATM;

        DAMFRAC = a1*TATM+a2base*TATM.^a3;

        YGROSS = eco2Param(1)*(K^gama);
        YNET = YGROSS*(1-DAMFRAC);
        ABATECOST = YGROSS*cost1tot(1)*(MIU(1)^expcost2);
        Y = YNET - ABATECOST;
        I = S(1)*Y;
        C(1,1) = Y - I;
     
        %Solve for Alpha0
        alpha(1) = params.a0;

    else
        CCATOT = CCATOT + ((sigma(i-1)*(eco2Param(i-1)*(K^gama)) + eland(i-1))*(1-MIU(i-1)))*5/3.666;

        K = (1-dk)^tstep*K + tstep*I;
        
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
        ceq(i) =  irf0+irC*(CCATOT-(MAT-mateq))+irT*TATM - ( alpha(i)*emshare0*tau0.*(1-exp(-100./(alpha(i)*tau0)))+ alpha(i)*emshare1*tau1.*(1-exp(-100./(alpha(i)*tau1))) + alpha(i)*emshare2*tau2.*(1-exp(-100./(alpha(i)*tau2))) + alpha(i)*emshare3*tau3.*(1-exp(-100./(alpha(i)*tau3))) );
    end
end