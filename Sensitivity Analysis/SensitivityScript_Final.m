%% Sensitivity Script
clear all
close all
clc
set(0,'DefaultAxesFontname', 'Times New Roman')
set(0,'DefaultAxesFontSize', 28)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 28)
set(0,'DefaultLineLinewidth',1)
set(0,'DefaultFigureColor','w')
set(0,'defaulttextinterpreter','latex')
set(0,'defaulttextinterpreter','latex')

global Kdox delz delx delg delw delc kp kpp k1 k2 k3 k4 k5 k6 kon koff ksgfp Kgfp k7 km kg n1 n2 k8 k9 k10 k11 k12
%%%Parameter        %%%Range                 %%%Reaction
km = 6e-5;          % Fitted                 %Dox Maximum Velocity

Kdox = 3.024;       % Fitted                 %Dox Km 
n1 = 1.0;           % Fitted                 %Dox Hill coefficient
 
Kgfp = 0.001;      % Fitted                 %Gfp Km
kg =  0.0141;      % Fitted                 %Gfp Maximum Velocity
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

%%% No complex phosphorylation
k8 = 0; k9 = 0; k10 = 0; k11 = 0; k12 = 0;

%%% Species concentrations
pt = 0.0220; % TR-SSRE Binding Sites
Xt = 0.0712; %SKN7 total
Wt = 0.1752; %YPD1 total
%% Unloaded Dynamic Sensitivity
lambda = diag([delz delw delx kp kpp k1 k2 k3 k4 k5 k6 k7 Xt Wt k8 k9 k10 k11 k12 delc kon koff pt]);

time = 0:1:1600;
[t,y] = ode23s(@SensitivityAnalysis_Dilution,time,zeros(170,1),[],0);
Xpp = y(:,5);
SXpp = y(:,102:124)*lambda;
scaleUnloaded = 6e-3;

figure
plot(t,SXpp(:,1),t,SXpp(:,2),t,SXpp(:,3),t,SXpp(:,4),t,SXpp(:,5))
legend('$\delta_Z$','$\delta_W$','$\delta_X$','$k_p$','$k_p$''')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X^{**}(t,\lambda)$')
ylim([-1 1]*scaleUnloaded)
xlim([0 1600])

figure
plot(t,SXpp(:,6),t,SXpp(:,7),t,SXpp(:,8),t,SXpp(:,9),t,SXpp(:,10))
legend('$k_1$','$k_2$','$k_3$','$k_4$','$k_5$')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X^{**}(t,\lambda)$')
ylim([-1 1]*scaleUnloaded)
xlim([0 1600])

figure
plot(t,SXpp(:,11),t,SXpp(:,12),t,SXpp(:,13),t,SXpp(:,14))
legend('$k_6$','$k_7$','$X_T$','$W_T$')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X^{**}(t,\lambda)$')
ylim([-1 1]*scaleUnloaded)
xlim([0 1600])

%% Loaded Dynamic Sensitivity

[tL,yL] = ode23s(@SensitivityAnalysis_Dilution,time,zeros(170,1),[],pt);

XppL = yL(:,5);
SXppL = yL(:,102:124)*lambda;

scaleLoad = 6e-3;

figure
plot(t,SXppL(:,1),t,SXppL(:,2),t,SXppL(:,3),t,SXppL(:,4),t,SXppL(:,5))
legend('$\delta_Z$','$\delta_W$','$\delta_X$','$k_p$','$k_p$''')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X_L^{**}(t,\lambda)$')
ylim([-1 1]*scaleLoad)
xlim([0 1600])

figure
plot(t,SXppL(:,6),t,SXppL(:,7),t,SXppL(:,8),t,SXppL(:,9),t,SXppL(:,10))
legend('$k_1$','$k_2$','$k_3$','$k_4$','$k_5$')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X_L^{**}(t,\lambda)$')
ylim([-1 1]*scaleLoad)
xlim([0 1600])

figure
plot(t,SXppL(:,11),t,SXppL(:,12),t,SXppL(:,13),t,SXppL(:,14))
legend('$k_6$','$k_7$','$X_T$','$W_T$')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X_L^{**}(t,\lambda)$')
ylim([-1 1]*scaleLoad)
xlim([0 1600])

figure
plot(t,SXppL(:,20),t,SXppL(:,21),t,SXppL(:,22),t,SXppL(:,23))
legend('$\delta_C$','$k_{on}$','$k_{off}$','$p_T$')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Sensitivity $\bar{S}_{\lambda}$ of $X_L^{**}(t,\lambda)$')
ylim([-1 1]*scaleLoad)
xlim([0 1600])

%% Output Sensitivity
 
SE = 2*diag(Xpp-XppL)/max(Xpp)*(SXpp(:,1:14) - SXppL(:,1:14));

scaleOut = 2.5e-4;
figure
plot(t,SE(:,1),t,SE(:,2),t,SE(:,3),t,SE(:,4),t,SE(:,5))
legend('$\delta_Z$','$\delta_W$','$\delta_X$','$k_p$','$k_p$''')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Error Sensitivity')
ylim([-6e-4 4e-4])
xlim([0 1600])

figure
plot(t,SE(:,6),t,SE(:,7),t,SE(:,8),t,SE(:,9),t,SE(:,10))
legend('$k_1$','$k_2$','$k_3$','$k_4$','$k_5$')
ylabel('$[\mu M]$')
xlabel('time [min]')
title('Error Sensitivity')
ylim([-6e-4 4e-4])
xlim([0 1600])

figure
plot(t,SE(:,11),t,SE(:,12),t,SE(:,13),t,SE(:,14))
legend('$k_6$','$k_7$','$X_T$','$W_T$')
ylabel('$[\mu M]$')
xlabel('Time [min]')
title('Error Sensitivuty')
ylim([-6e-4 4e-4])
xlim([0 1600])

%% Output Error Timescale and Absolute Concentration change

time = 0:1:6000;
%%% Creating DOX input
dox = zeros(max(time)+1,1);
input = 20*ones(50+1,1);
            for c = 1:ceil(6000/500)
                dox(1+(c-1)*500:1+(c-1)*500+50,1) = input;
            end
doxTrimed = dox(1:6001,1);    
%%% Choice of Gains
G1 = .1:.02:1; %1:.5:10; %= 1:.2:3; %1:.5:10; %.1:.02:1; %logspace(-1.3,0,20);
G2 = 1:.5:10;

%%% Main Loop
 matlabpool open
for a = 1:size(G1,2)
     parfor b = 1:size(G2,2)
 [tE,yE] = ode23s(@InsulatorFunc500G1G2,time,zeros(16,1),[],G1(a),G2(b),doxTrimed,pt,Kdox,delz,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kgfp,delg,k7,km,kg,n1,n2,k8,k9,k10,0,k11,k12);  
    X = yE(:,5);
    XL = yE(:,13);
    e = ((X - XL)/max(X(5000:6001))).^2;
    E(a,b) = trapz(time,e)/6000;
    disp(b)
    end
disp(a)    
end
 matlabpool close

%%% Subplot Choice of Gains 
G1_t = [1 5]; G2_t = [9 4];
% G1_t = [1e-2 1e-1];  G2_t = [9 3]
% G1_t = [.1 .4]; G2_t = [9 3];
% G1_t = [1 2.5]; G2_t = [9 6];
for a = 1:2
[tE,yE] = ode23s(@InsulatorFunc500G1G2,time,zeros(16,1),[],G1_t(a),G2_t(a),doxTrimed,pt,Kdox,delz,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kgfp,delg,k7,km,kg,n1,n2,k8,k9,k10,0,k11,k12);  
    X(:,a) = yE(:,5);
    XL(:,a) = yE(:,13);
end

figure
subplot(1,1,1)
caxis([0 0.3])
contour(G2,G1,E,20)
colorbar
% colorbar('YTicklabel',{'1E-1','1E+0','1E+1','1E+2','1E+3','1E+4','1E+5','1E+6'})
% set(gca,'yscale','log')
xlabel('\gamma [Fold increase in promoter sites]')
% ylabel('G [Fold increase in timescale]')
ylabel('\rho [Fold decrease in X_T and W_T]')
title('Error Level Map')

axes('position',[0.5 0.4 0.2 0.15]) ; % inset
plot(time,X(:,1),'k',time,XL(:,1),'r')
xlim([0 6000])
ylim([-max(X(:,1))*.1 max(X(:,1))*1.1])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

axes('position',[0.2 0.7 0.2 0.15]) ; % inset
plot(time,X(:,2),'k',time,XL(:,2),'r')
xlim([0 6000])
ylim([-max(X(:,2))*.1 max(X(:,2))*1.1])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])

% 
figure
surf(G2,G1,E)
% ylim([10^(-1.3) 10^0])
% set(gca,'yScale','log')
view([0 90])
colorbar
% colorbar('YTicklabel',{'1E-1','1E+0','1E+1','1E+2','1E+3','1E+4','1E+5','1E+6'})
% colorbar('YScale','log')
xlabel('\gamma [Fold increase in promoter sites]')
% ylabel('G [Fold increase in timescale]')
 ylabel('\rho [Fold decrease in X_T and W_T]')
title('Error Heat Map')
   
    



return