function dx = Control_Dilution(t,x,u,pt)

%% Parameters
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