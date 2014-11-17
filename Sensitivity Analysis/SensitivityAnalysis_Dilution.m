%% Sensitivity Analysis
function dx = SensitivityAnalysis_Dilution(t,x,pt)

global Kdox delz delx delg delw delc kp kpp k1 k2 k3 k4 k5 k6 kon koff ksgfp Kgfp k7 km kg n1 n2 k8 k9 k10 k11 k12
Xt = 0.0712;
Wt = 0.1752;
%%Input

if t<=58 || t >=500 && t<550 || t>=1000 && t<1050 || t>= 1500 && t<1550 
       
    dox = 20;
else 
    dox = 0;
end
%% States

Z = x(1);
Zp = x(2);
Wp = x(3);
Xp = x(4);
Xpp = x(5);
Cp = x(6);
Cpp = x(7);
GFP = x(8);

X = Xt - Xp - Xpp  - Cp - Cpp;
W = Wt - Wp;
p = pt - Cp - Cpp;

Sx(1,1:23) = x(10:32);
Sx(2,1:23) = x(33:55);
Sx(3,1:23) = x(56:78);
Sx(4,1:23) = x(79:101);
Sx(5,1:23) = x(102:124);
Sx(6,1:23) = x(125:147);
Sx(7,1:23) = x(148:170);
%% ODEs
dx = zeros(170,1);

dx(1) = km*dox^n1/(Kdox + dox^n1) - delz*Z - k2*Wp*Z + k1*Zp*(Wt-Wp) - kp*Z + kpp*Zp;
dx(2) = -k1*Zp*(Wt-Wp) + k2*Wp*Z + kp*Z - kpp*Zp - delz*Zp;
dx(3) = k1*Zp*(Wt-Wp) - k2*Wp*Z - k3*X*Wp + k4*Xp*(Wt-Wp) - k3*Xp*Wp + k4*Xpp*(Wt-Wp) -k7*Wp -delw*Wp...
        + k12*Cp*(Wt - Wp) -k9*Cp*Wp + k10*Cpp*(Wt - Wp);
dx(4) = k3*X*Wp - k4*Xp*(Wt-Wp) - k3*Xp*Wp +k4*Xpp*(Wt-Wp) - k5*Xp -delx*Xp + k6*Xpp ...
        -kon*Xp*p + koff*Cp + delc*Cp;
dx(5) = k3*Xp*Wp - k4*Xpp*(Wt - Wp) - k6*Xpp -delx*Xpp...
        -kon*Xpp*p + koff*Cpp + delc*Cpp;
dx(6) = kon*Xp*p - koff*Cp - delc*Cp - k9*Cp*Wp + k10*Cpp*(Wt - Wp) + k8*Cpp -k11*Cp - k12*Cp*(Wt-Wp);
dx(7) = kon*Xpp*p - koff*Cpp - delc*Cpp + k9*Cp*Wp - k10*Cpp*(Wt - Wp)- k8*Cpp;
dx(8) = ksgfp + kg*Xpp^n2/(Kgfp + Xpp^n2) - delg*GFP;

dfdZ = [-delz-k2*Wp-kp;k2*Wp+kp;-k2*Wp;0;0;0;0];
dfdZp = [k1*W+kpp;-k1*W-kpp-delz;k1*W;0;0;0;0];
dfdWp = [-k2*Z-k1*Zp;k1*Zp+k2*Z;-k1*Zp-k2*Z-k3*X-k4*Xp-k4*Xpp-k3*Xp-(k7+delw)-k9*Cp-k10*Cpp-Cp*k12;k3*X+k4*Xp-k3*Xp-k4*Xpp;k3*Xp+k4*Xpp;-k9*Cp-k10*Cpp+k12*Cp;k9*Cp+k10*Cpp];
dfdXp = [0;0;k4*W;-2*k3*Wp-k4*W-(k5+delx)-kon*p;k3*Wp;kon*p;0];
dfdXpp = [0;0;k4*W+k3*Wp;-k3*Wp+k4*W+k6;-k4*W-(k6+delx)-kon*p;0;kon*p];
dfdCp = [0;0;k3*Wp-k9*Wp+k12*W;-k3*Wp+kon*Xp+(koff+delc);kon*Xpp;-kon*Xp-(koff+delc)-k9*Wp-k11-k12*W;-kon*Xpp+k9*Wp];
dfdCpp = [0;0;k3*Wp+k10*W;-k3*Wp+kon*Xp;kon*Xpp+(koff+delc);-kon*Xp+k8+k10*W;-kon*Xpp-(koff+delc)-k8-k10*W];
dfdxL = [dfdZ dfdZp dfdWp dfdXp dfdXpp dfdCp dfdCpp];

dfdp1 = [-Z;-Zp;0;0;0;0;0];
dfdp2 = [0;0;-Wp;0;0;0;0];
dfdp3 = [0;0;0;-Xp;-Xpp;0;0];
dfdp4 = [-Z;Z;0;0;0;0;0];
dfdp5 = [Zp;-Zp;0;0;0;0;0];
dfdp6 = [Zp*W;-Zp*W;Zp*W;0;0;0;0];
dfdp7 = [-Wp*Z;Wp*Z;-Wp*Z;0;0;0;0];
dfdp8 = [0;0;-X*Wp-Xp*Wp;X*Wp-Xp*Wp;Xp*Wp;0;0];
dfdp9 = [0;0;Xp*W+Xpp*W;-Xp*W+Xpp*W;-Xpp*W;0;0];
dfdp10 = [0;0;0;-Xp;0;0;0];
dfdp11= [0;0;0;Xpp;-Xpp;0;0];
dfdp12 = [0;0;-Wp;0;0;0;0];
dfdp13 = [0;0;-k3*Wp;k3*Wp;0;0;0];
dfdp14 = [k1*Zp;-k1*Zp;k1*Zp+k4*Xp+k4*Xpp+k10*Cpp+k12*Cp;-k4*Xp+k4*Xpp;-k4*Xpp;k10*Cp-k12*Cp;-k10*Cpp];
dfdp15 = [0;0;0;0;0;Cpp;-Cpp];
dfdp16 = [0;0;-Cp*Wp;0;0;-Cp*Wp;Cp*Wp];
dfdp17 = [0;0;Cpp*W;0;0;Cpp*W;-Cpp*W];
dfdp18 = [0;0;0;0;0;-Cp;0];
dfdp19 = [0;0;Cp*W;0;0;-Cp*W;0];
dfdp20 = [0;0;0;Cp;Cpp;-Cp;-Cpp];
dfdp21 = [0;0;0;-Xp*p;-Xpp*p;Xp*p;Xpp*p];
dfdp22 = [0;0;0;Cp;Cpp;-Cp;-Cpp];
dfdp23 = [0;0;0;-kon*Xp;-kon*Xpp;kon*Xp;kon*Xpp];



dfdpL = [dfdp1 dfdp2 dfdp3 dfdp4 dfdp5 dfdp6 dfdp7 dfdp8 dfdp9 dfdp10 dfdp11 dfdp12 dfdp13 dfdp14 dfdp15 dfdp16 dfdp17 dfdp18 dfdp19 dfdp20 dfdp21 dfdp22 dfdp23];

dSdt = dfdxL*Sx + dfdpL;

dx(10:32) = dSdt(1,1:23);
dx(33:55) = dSdt(2,1:23);
dx(56:78) = dSdt(3,1:23);
dx(79:101) = dSdt(4,1:23);
dx(102:124) = dSdt(5,1:23);
dx(125:147) = dSdt(6,1:23);
dx(148:170) = dSdt(7,1:23);

end
