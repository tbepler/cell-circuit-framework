function dx = ControlFunc150_Dilution(t,x,pt)


%% Input Data
if t<=200 || t >=350 && t<400 || t>=500 && t<550 || t>= 650 && t<700 || t>= 800 && t<850 || t >=950 && t<1000 || t>=1100 && t<1150 || t>=1250 && t<1300 || t>= 1400 && t<1450 || t>=1550 && t<1600 || t>=1700 && t<1750 || t>=1850 && t<1900
       
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