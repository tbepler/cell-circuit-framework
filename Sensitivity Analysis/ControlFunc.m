function dx = ControlFunc(t,x,u,del,Kdox,kon,koff,ksgfp,kdgfp,Kgfp,n2,n1,km,kg)

%% States
 
%% Input Data
dox1 = u(1,1); 
pt = u(1,2);
    
% %% States 
%    X = x(1);  % Skn7m
%    C = x(2);  % Load Complex    
%    G = x(3);  % GFP
%    dx = zeros(3,1);
% 
% %% ODEs
%     dx(1) = km*dox1^n1/(Kdox + dox1^n1) - del*X - kon*X*(pt1 - C) + koff*C;
%     dx(2) = kon*X*(pt1 - C) - koff*C - del*C; 
%     dx(3) = ksgfp+kg*X^n2/(Kgfp+X^n2) - del*G;

%% States 
   X = x(1);  % Skn7m
   C = x(2);  % Load Complex    
   G = x(3);  % GFP
   dx = zeros(3,1);

%% ODEs
    dx(1) = km*dox1^n1/(Kdox + dox1^n1) - del*X - kon*X*(pt - C) + koff*C;
    dx(2) = kon*X*(pt - C) - koff*C; 
    dx(3) = ksgfp+kg*X^n2/(Kgfp+X^n2) - kdgfp*G;
end