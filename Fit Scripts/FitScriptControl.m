%% GreyBoxFitScript Control
clear all
clc
close all
ContDataRevised; % Loading Data File


%% Parameters Initial Estimate
global delx kon Kgfp koff km kg Kdox del ksgfp delg

%% New parameter fit
    del = 0;
    kon = 6; 
    km = 6e-5;
    kg = 0.0088;
    Kdox = 3.0240;
    Kgfp = 0.001;
    ksgfp = 0;
    delx = 0.0130;
    delg = 0.004;
    koff  = kon*Kgfp;
%% Estimation GrayBox Model Control

order  = [1 2 3]; % Model order [# outputs, # inputs, # states ]

% InitialStates = {0 0 0};
InitialStates = struct('Name',    {'s1' 's2' 's3'}, ...
                       'Unit',    {'u1' 'u2' 'u3'},     ...
                       'Value',   {0 0 0},                                         ...
                       'Minimum', {0 0 0},                                             ...
                       'Maximum', {Inf Inf Inf},                                         ...
                       'Fixed',   {false false false}); %Fixing initial states to zero   
                   
Ts = 0; % Sampling time for continuous data

parameters    = {del,kon,km,kg,Kdox,Kgfp,ksgfp,delx,delg};

Control_model= idnlgrey('ControlGrey_c',order,parameters,InitialStates,Ts);
%%%                             del,    kon,    km,         kg,     Kdox,   Kgfp,  ksgfp,   delx,   delg
% setpar(Control_model,'Minimum',{.0017,  0.006,  3.47e-5,    1e-9,   1e-3,   1e-4,  0,       0.0017,  0.0017});
% setpar(Control_model,'Maximum',{.0347,  6,      1e3,        1e3,    20,     10e-2, 0,       0.0347, 0.0347});
setpar(Control_model,'Fixed',{  true,  false,  false,      false,  false,  false, true,    false,  false});
 
Control_model.Algorithm.SearchMethod = 'lsqnonlin';
 %%% Trust Region based algorithms, lsqnonlin = Newton Method
 %%%                                 lm = Levenberg-Marquardt 
Control_model.Algorithm.SimulationOptions.Solver = 'ode15s';
Control_model.Algorithm.SimulationOptions.RelTol = 1e-6;


%%
Data1 = merge(ContSq3,ContSq3L,ContSq4,ContSq4L,ContSq5,ContSq5L); %,Cont20W,Cont20WL
Data2 = merge(Cont20,Cont20L); 
Control_model = pem(Control_model,ContDataPt,'Display','Full','MaxIter',200,'Tolerance',1e-6);  
%Control_model = pem(ContDataPt,Control_model,'Display','Full','MaxIter',200,'Tolerance',1e-6);              

%%


%%% Printing Parameters
for a = 1:size(parameters,2)
par(a) = Control_model.Parameters(a).Value;
end


    del = par(1)
    kon = par(2) 
    km = par(3)
    kg = par(4)   
    Kdox = par(5)
    Kgfp = par(6)
    ksgfp = par(7)
    delx = par(8)
    delg = par(9)
    koff  = kon*Kgfp
    ksgfp = 0

values = [del kon koff km kg Kdox Kgfp ksgfp delx delg];
save('Parameters','values'); %% Saving parameter values
