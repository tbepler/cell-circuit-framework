function dx = InsulatorFunc500G1G2(t,x,G1,G2,u,pt,Kdox,del,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2,k11,k12)

 

%% Input Data
% if t<=58 || t >=500 && t<550 || t>=1000 && t<1050 || t>= 1500 && t<1550 
%        
%     dox = 20;
% else 
%     dox = 0;
% end
time = 0:1:6000;
dox = interp1(time,u,t);

Xt = G1*0.0712; %G1*0.0712;
pt = G2*pt; 
Wt = G1*0.1752; %G1*0.1752; 

%% States Unloaded
Z = x(1);
Zp = x(2);
Wp = x(3);
Xp = x(4);
Xpp = x(5);
Cp = x(6);
Cpp = x(7);
GFP = x(8);
X = Xt - Xp - Xpp - Cp - Cpp;
p = 0 - Cp - Cpp;
%% ODEs 
dx = zeros(16,1);

dx(1) = km*dox^n1/(Kdox + dox^n1) - del*Z +1*(- k2*Wp*Z + k1*Zp*(Wt-Wp) - kp*Z + kpp*Zp);
dx(2) = 1*(-k1*Zp*(Wt-Wp) + k2*Wp*Z + kp*Z - kpp*Zp) - del*Zp;
dx(3) = 1*(k1*Zp*(Wt-Wp) - k2*Wp*Z - k3*X*Wp + k4*Xp*(Wt-Wp) - k3*Xp*Wp + k4*Xpp*(Wt-Wp) -k7*Wp);
dx(4) = 1*(k3*X*Wp - k4*Xp*(Wt-Wp) - k3*Xp*Wp +k4*Xpp*(Wt-Wp) - k5*Xp + k6*Xpp) ...
        -kon2*Xp*p + koff*Cp;
dx(5) = 1*(k3*Xp*Wp - k4*Xpp*(Wt - Wp) - k6*Xpp)...
        -kon*Xpp*p + koff*Cpp;
dx(6) = kon2*Xp*p - koff*Cp;
dx(7) = kon*Xpp*p - koff*Cpp;
dx(8) = ksgfp + kg*Xpp^n2/(Kdgfp + Xpp^n2) - kdgfp*GFP;

% zT = x(1);
% zp = x(2);
% wp = x(3);
% xp = x(4);
% xpp = x(5);
% cp = x(6);
% cpp = x(7);
% g = x(8);
% 
% dx(1) = km*dox/(Kdox + dox) - del*zT;
% dx(2) = G1*(-c1*del*zp*(1 - wp) + c2*del*wp*(zT - zp) + kappap*del*(zT - zp) )-kapap*del*zp - del*zp;
% dx(3) = G1*((Z0/Wt)*zp*(1-wp) - c2*del*(Z0/Wt)*wp*(zT - zp) - c3*alpha*del*(1 - xp - xpp - (pt/Xt)*cp - (pt/Xt)*cpp)*wp ) ...
%        + G1*( alpha*del*xp*(1 - wp) - c3*alpha*del*xp*wp + alpha*del*xpp*(1 - wp) -    )
%% States Loaded
ZL = x(9);
ZpL = x(10);
WpL = x(11);
XpL = x(12);
XppL = x(13);
CpL = x(14);
CppL = x(15);
GFPL = x(16);
XL = Xt - XpL - XppL - CpL - CppL;
pL = pt - CpL - CppL;
%% ODEs 

dx(9) = km*dox^n1/(Kdox + dox^n1) - del*ZL +1*(- k2*WpL*ZL + k1*ZpL*(Wt-WpL) - kp*ZL + kpp*ZpL);
dx(10) = 1*(-k1*ZpL*(Wt-WpL) + k2*WpL*ZL + kp*ZL - kpp*ZpL) - del*ZpL;
dx(11) = 1*(k1*ZpL*(Wt-WpL) - k2*WpL*ZL - k3*XL*WpL + k4*XpL*(Wt-WpL) - k3*XpL*WpL + k4*XppL*(Wt-WpL) -k7*WpL);
dx(12) = 1*(k3*XL*WpL - k4*XpL*(Wt-WpL) - k3*XpL*WpL +k4*XppL*(Wt-WpL) - k5*XpL + k6*XppL) ...
        -kon2*XpL*pL + koff*CpL;
dx(13) = 1*(k3*XpL*WpL - k4*XppL*(Wt - WpL) - k6*XppL)...
        -kon*XppL*pL + koff*CppL;
dx(14) = kon2*XpL*pL - koff*CpL;
dx(15) = kon*XppL*pL - koff*CppL;
dx(16) = ksgfp + kg*XppL^n2/(Kdgfp + XppL^n2) - kdgfp*GFPL;




end