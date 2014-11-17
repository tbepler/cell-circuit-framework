%% GreyBoxFitScriptIns Insulator
clear all
clc
close all

InsData;

%% Parameters Original
%%%Parameter%%%Range                 %%%Reaction
G = .7;
Gw = .4;
Gx = .3;
km = .6e-4;         % Fitted                 %Dox Maximum Velocity

Kdox = 3.024;       % Fitted                 %Dox Km 
n1 = 1.0;           % Fitted                 %Dox Hill coefficient
 
kg =  14.0743;      % Fitted                 %Gfp Maximum Velocity
n2 = 1.0;           % Fitted                 %Gfp Hill coefficient
delg = 0.004;       % 0.0039< del <0.0069    %Gfp Decay
ksgfp = 0.0012;     % Fited                 %Bassal Gfp expression

kon = 6;            % kon < 6                %Skn7**/DNA association rate
koff = 0.006;       % oder of min/hr         %Skn7**/DNA dissociation rate
Kdgfp = koff/kon; 
 
delz = 0.013;       % From Control           %Snl1 Decay
kp = 0.6323;        % Fitted                 %Snl1 phosphorylation
kpp = 0.0373;       % Fitted                 %Snl1 dephosphorylation
k1 = 350.1322;      % Fitted                 %(Snl1*/Ypd1) Ypd1 phosphorylation
k2 = 881.2300;      % Fitted                 %(Snl1/Ypd1*) Ypd1 dephosphorylation

delx = 0;           % From Control           %SKN7 decay
k3 = 335.1981;      % Fitted                 %(Skn7/Ypd1*) Skn7 phosphorylation     
k4 = 42.0015;       % Fitted                 %(Skn7*/Ypd1) Skn7 dephosphorylation
k5 =   0.0024;      % Fitted                %Skn7* spontaneous dephosphorylation
k6 =  0.4818;       % Fitted                 %Skny** spontaneous dephosphorylation
 
k7 = 0.0024;        % Fitted                %Ypd1 spontaneous dephosphorylation
delw = 0;           % Ref                    %Ypd1 decay

delc = 0;       % Added                  %Complex TF clearance
%%% No laod phosphorylation
k8 = 0; k9 = 0; k10 = 0; k11 = 0; k12 = 0;
kon2 = kon;

pt = 0.0220;
Xt = 0.0712;
Wt = 0.1752;
%% Creating GrayBox Model Insulator
order  = [1 4 8];
parameters    = {Kdox,delz,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,delg,k7,km,kg,n1,n2,k8,k9,k10,k11,k12,kon2};           
% initial_states = [0 0 0 0 0 0 0 0]';             
Ts = 0;
Insulator_model= idnlgrey('InsulatorGrey',order,parameters);
                                % Kdox delz    kp    kpp   k1    k2    k3    k4    k5    k6    kon  koff ksgfp Kdgfp  delg   k7    km    kg    n1    n2    k8    k9    k10   k11   k12 kon2
setpar(Insulator_model,'Minimum',{.001 0.004 0.000 0.000 500.0 500.0 60.00 60.00 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 0.000 1.000 1.000 0.004 0.004 0.004 0.004 0.004 1e-3});
                                % Kdox   delz    kp    kpp   k1   k2    k3  k4    k5    k6    kon  koff ksgfp Kdgfp  delg   k7    km    kg    n1    n2     k8    k9    k10   k11   k12  kon2
% setpar(Insulator_model,'Maximum',{10.00 0.0161 1.000 1.000 3000 3000 1000 500.0 1.000 1.000   inf   inf inf    inf    inf  1.000  inf   inf   2.000 2.000 0.100 0.100  1.000 0.100 20.00 6});
                               %Kdox delz  kp    kpp    k1   k2    k3    k4      k5    k6   kon  koff  ksgfp  Kdgfp    delg    k7    km   kg   n1   n2    k8    k9    k10   k11  k12  kon2
setpar(Insulator_model,'Fixed',{true,true,false,  false,true,true,  true,  true, false,false, true,  true,true,  true, false, false,  true,false,true,true,true,true,true,true,true,true});

Insulator_model.Algorithm.SearchMethod = 'lsqnonlin'; %lsqnonlin, lm
 %%% Trust Region based algorithms, lsqnonlin = Newton Method
 %%%                                 lm = Levenberg-Marquardt 
Insulator_model.Algorithm.SimulationOptions.Solver = 'ode15s';
Insulator_model.Algorithm.SimulationOptions.RelTol = 1e-6;

Insulator_model = pem(merge(InsDataRed,InsDataPt),Insulator_model,'Display','Full','MaxIter',100,'Tolerance',1e-6);

Kdox = Insulator_model.Parameters(1).Value
delz = Insulator_model.Parameters(2).Value
kp =  Insulator_model.Parameters(3).Value
kpp =  Insulator_model.Parameters(4).Value
k1 =  Insulator_model.Parameters(5).Value
k2 = Insulator_model.Parameters(6).Value
k3 = Insulator_model.Parameters(7).Value
k4 = Insulator_model.Parameters(8).Value
k5 = Insulator_model.Parameters(9).Value
k6 = Insulator_model.Parameters(10).Value
kon = Insulator_model.Parameters(11).Value
koff = Insulator_model.Parameters(12).Value
ksgfp = Insulator_model.Parameters(13).Value
Kdgfp = Insulator_model.Parameters(14).Value
delg = Insulator_model.Parameters(15).Value
k7 = Insulator_model.Parameters(16).Value
km = Insulator_model.Parameters(17).Value
kg = Insulator_model.Parameters(18).Value
n1 = Insulator_model.Parameters(19).Value
n2 = Insulator_model.Parameters(20).Value
k8 = Insulator_model.Parameters(21).Value
k9 = Insulator_model.Parameters(22).Value
k10 = Insulator_model.Parameters(23).Value
k11 = Insulator_model.Parameters(24).Value
k12 = Insulator_model.Parameters(25).Value
kon2 = Insulator_model.Parameters(26).Value
