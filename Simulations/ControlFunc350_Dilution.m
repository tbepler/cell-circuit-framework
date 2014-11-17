function dx = ControlFunc350_Dilution(t,x,pt)


%% Input Data
if t<=60 || t >=350 && t<400 || t>=700 && t<750 || t>= 1050 && t<1100 || t>= 1400 && t<1450 || t>= 1750 && t<1800       
    u = 20;
else 
    u = 0;
end


global km Kdox delx delg kon koff ksgfp kg Kgfp delc

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