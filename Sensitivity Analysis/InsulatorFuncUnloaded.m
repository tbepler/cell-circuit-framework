function dx = InsulatorFuncUnloaded(t,x,u,Kdox,del,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2)
 


%% Inputs
dox = u(1,1);
pt = u(1,2); 

Xt = 0.0712; %2570/(6.02e23)*1/(60e-15);
Wt = 0.1752; %6330/(6.02e23)*1/(60e-15);
%% States
Z = x(1);
Zp = x(2);
Wp = x(3);
Xp = x(4);
Xpp = x(5);
Cp = x(6);
Cpp = x(7);
GFP = x(8);
X = Xt - Xp - Xpp - Cp - Cpp;
p = pt - Cp - Cpp;
%% ODEs 
dx = zeros(8,1);
dx(1) = km*dox^n1/(Kdox + dox^n1) - del*Z - k2*Wp*Z + k1*Zp*(Wt-Wp) - kp*Z + kpp*Zp;
dx(2) = -k1*Zp*(Wt-Wp) + k2*Wp*Z + kp*Z - kpp*Zp - del*Zp;
dx(3) = k1*Zp*(Wt-Wp) - k2*Wp*Z - k3*X*Wp + k4*Xp*(Wt-Wp) - k3*Xp*Wp + k4*Xpp*(Wt-Wp) -k7*Wp ...
        +-k9*Cp*Wp + k10*Cpp*(Wt - Wp);
dx(4) = k3*X*Wp - k4*Xp*(Wt-Wp) - k3*Xp*Wp +k4*Xpp*(Wt-Wp) - k5*Xp + k6*Xpp ...
        -kon2*Xp*p + koff*Cp;
dx(5) = k3*Xp*Wp - k4*Xpp*(Wt - Wp) - k6*Xpp...
        -(kon*Xpp*p + koff*Cpp)*0;
dx(6) = kon2*Xp*p - koff*Cp - k9*Cp*Wp + k10*Cpp*(Wt - Wp)+ k8*Cpp;
dx(7) = kon*Xpp*p - koff*Cpp + k9*Cp*Wp - k10*Cpp*(Wt - Wp)- k8*Cpp;
dx(8) = ksgfp + kg*Xpp^n2/(Kdgfp + Xpp^n2) - kdgfp*GFP;

end