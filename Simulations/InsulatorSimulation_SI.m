%% GreyBoxFitScriptIns Insulator
clear all
clc
close all
%% Parameters
global Kdox delz delx delg delw delc kp kpp k1 k2 k3 k4 k5 k6 kon koff ksgfp Kdgfp k7 km kg n1 n2

%%%Parameter        %%%Range                 %%%Reaction
km = 6e-5;          % Fitted                 %Dox Maximum Velocity

Kdox = 3.024;       % Fitted                 %Dox Km 
n1 = 1.0;           % Fitted                 %Dox Hill coefficient
 
Kdgfp = 0.001;      % Fitted                 %Gfp Km
kg =  8.8;          % Fitted                 %Gfp Maximum Velocity
n2 = 1.0;           % Fitted                 %Gfp Hill coefficient
delg = 0.004;       % 0.0039< del <0.0069    %Gfp Decay
ksgfp = 0.0012;     % Fited                 %Bassal Gfp expression

kon = 6;            % kon < 6                %Skn7**/DNA association rate
koff = 0.006;       % oder of min/hr         %Skn7**/DNA dissociation rate
 
delz = 0.013;       % From Control           %Snl1 Decay
kp = 0.6323;        % Fitted                 %Snl1 phosphorylation
kpp = 0.0373;       % Fitted                 %Snl1 dephosphorylation
k1 = 350.1322;      % Fitted                 %(Snl1*/Ypd1) Ypd1 phosphorylation
k2 = 881.2300;      % Fitted                 %(Snl1/Ypd1*) Ypd1 dephosphorylation

delx = 0.004;       % From Control           %SKN7 decay
k3 = 335.1981;      % Fitted                 %(Skn7/Ypd1*) Skn7 phosphorylation     
k4 = 42.0015;       % Fitted                 %(Skn7*/Ypd1) Skn7 dephosphorylation
k5 =   0.001;       % By Hand                %Skn7* spontaneous dephosphorylation
k6 =  0.4818;       % Fitted                 %Skny** spontaneous dephosphorylation
 
k7 = 0.001;         % By Hand                %Ypd1 spontaneous dephosphorylation
delw = 0.004;       % Ref                    %Ypd1 decay

delc = 0.004;       % Added                  %Complex TF clearance

pt = 0.0220;
Xt = 0.0712;
Wt = 0.1752;
%%
DynamicAnalysisInsulator_SI
