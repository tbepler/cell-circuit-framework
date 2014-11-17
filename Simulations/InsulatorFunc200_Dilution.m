function dx = InsulatorFunc200_Dilution(t,x,pt)



%% Input Data
if t<=150 || t >=250 && t<300 || t>=450 && t<500 || t>= 650 && t<700 || t>= 850 && t<900 || t >=1050 && t<1100 || t>=1250 && t<1300 || t>= 1450 && t<1500
       
    dox = 20;
else 
    dox = 0;
end

global Kdox delz delx delg delw delc kp kpp k1 k2 k3 k4 k5 k6 kon koff ksgfp Kdgfp k7 km kg n1 n2 
 
X0 = 0.0712; %2570/(6.02e23)*1/(60e-15);
W0 = 0.1752; %6330/(6.02e23)*1/(60e-15);
kx = X0*delx;
kw = W0*delw;
%% States
Zt = x(1);
Zp = x(2);
 
Wt = x(3);
Wp = x(4);
 
Xt = x(5);
Xp = x(6);
Xpp = x(7);
 
Cp = x(8);
Cpp = x(9);
 
GFP = x(10);
 
Z = Zt - Zp;
X = Xt - Xp - Xpp - Cp - Cpp;
p = pt - Cp - Cpp;
%%
dx = zeros(10,1);
dx(1) = km*dox^n1/(Kdox + dox^n1) - delz*Zt;
dx(2) = -k1*Zp*(Wt-Wp) + k2*Wp*Z + kp*Z - kpp*Zp - delz*Zp;
dx(3) = kw - delw*Wt;
dx(4) = k1*Zp*(Wt-Wp) - k2*Wp*Z - k3*X*Wp + k4*Xp*(Wt-Wp) - k3*Xp*Wp ...
        +k4*Xpp*(Wt-Wp) -k7*Wp - delw*Wp;
dx(5) = kx - delx*Xt;
dx(6) = k3*X*Wp - k4*Xp*(Wt-Wp) - k3*Xp*Wp +k4*Xpp*(Wt-Wp) - k5*Xp + k6*Xpp -delx*Xp ...
        -kon*Xp*p + koff*Cp + delc*Cp;
dx(7) = k3*Xp*Wp - k4*Xpp*(Wt - Wp) - k6*Xpp - delx*Xpp...
        -kon*Xpp*p + koff*Cpp + delc*Cpp;
dx(8) = kon*Xp*p - koff*Cp - delc*Cp;
dx(9) = kon*Xpp*p - koff*Cpp - delc*Cpp;
dx(10) = ksgfp + kg*Xpp^n2/(Kdgfp + Xpp^n2) - delg*GFP;

end