function dx = ControlFunc200_Dilution(t,x,pt)


%% Input Data
if t<=135 || t >=250 && t<300 || t>=450 && t<500 || t>= 650 && t<700 || t>= 850 && t<900 || t >=1050 && t<1100 || t>=1250 && t<1300 || t>= 1450 && t<1500
       
    u = 20;
else 
    u = 0;
end

global km Kdox delc delx delg kon koff ksgfp kg Kgfp

%% States 
   X = x(1);  % Skn7m
   C = x(2);  % Load Complex    
   G = x(3);  % GFP
   dx = zeros(3,1);
 
%% ODEs
   dx(1) = km*u/(Kdox + u)- delx*X - kon*X*(pt - C) + koff*C + delc*C;
    dx(2) = kon*X*(pt - C) - koff*C - delc*C; 
    dx(3) = ksgfp+kg*X/(Kgfp+X) -delg*G;
end