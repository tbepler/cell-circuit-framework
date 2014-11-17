%% Steady state sensitivity
clear all
close all
clc
set(0,'DefaultAxesFontname', 'Times New Roman')
set(0,'DefaultAxesFontSize', 28)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 28)
set(0,'DefaultLineLinewidth',1)
set(0,'DefaultFigureColor','w')
%% Parameters
Kdox = 3.024;       % Fitted                 %Dox Km 
km = 1.2e-04;       % Fitted                 %Dox Maximum Velocity
n1 = 1.0;           % Fitted                 %Dox Hill coefficient
 
Kdgfp = 0.0023;     % Fitted                 %Gfp Km 
kg =  0.0034;       % Fitted                 %Gfp Maximum Velocity
n2 = 1.0;           % Fitted                 %Gfp Hill coefficient
kdgfp = 0.004;      % 0.0039< del <0.0069    %Gfp Decay
ksgfp = 0.0012;     % Fitted                 %Bassal Gfp expression
 
kon = 6;            % kon < 6                %Skn7**/DNA association rate
koff = 0.0138;      % oder of min/hr         %Skn7**/DNA dissociation rate
 
 
del = 0.0065;       % From Control           %Snl1 Decay
kp =  0.9033;       % Fitted                 %Snl1 phosphorylation
kpp = 0.0533;       % Fitted                 %Snl1 dephosphorylation
k1 = 500.1888;      % Fitted                 %(Snl1*/Ypd1) Ypd1 phosphorylation
k2 = 1.2589e+03;    % Fitted                 %(Snl1/Ypd1*) Ypd1 dephosphorylation
 
k3 = 478.8544;      % Fitted                 %(Skn7/Ypd1*) Skn7 phosphorylation     
k4 = 60.0021;       % Fitted                 %(Skn7*/Ypd1) Skn7 dephosphorylation
k5 =  0.01;         % By Hand                %Skn7* spontaneous dephosphorylation
k6 = 0.6948;        % Fitted                 %Skny** spontaneous dephosphorylation
 
k7 = 0.01;          % By Hand                %Ypd1 spontaneous dephosphorylation
 
%%% No complex phosphorylation
k8 = 0; k9 = 0; k10 = 0; k11 = 0; k12 = 0;
kon2 = kon;         % Fitted                 %Skn7*/DNA association rate
%%
dox = logspace(-3,2);
for a = 1:size(dox,2)
    
[t,ydelp] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del*1.5,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
[t,ydelm] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del*.5,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
Gdelp(a) = ydelp(end,8);
Gdelm(a) = ydelm(end,8);

[t,ykpp] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del,kp*1.5,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
[t,ykpm] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del,kp*.5,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
Gkpp(a) = ykpp(end,8);
Zkpp(a) = ykpp(end,1);
Zpkpp(a) = ykpp(end,2);

Gkpm(a) = ykpm(end,8);
Zkpm(a) = ykpm(end,1);
Zpkpm(a) = ykpm(end,2);

[t,yk6p] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del,kp,kpp,k1,k2,k3,k4,k5,k6*1.5,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
[t,yk6m] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del,kp,kpp,k1,k2,k3,k4,k5,k6*.5,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
Gk6p(a) = yk6p(end,8);
Gk6m(a) = yk6m(end,8);

[t,ydelxp] = ode23s(@ControlFunc,[0 5000],[0 0 0]',[],[dox(a) 0],del*1.5,Kdox,kon,koff,ksgfp,kdgfp,Kdgfp,n2,n1,km,kg);
[t,ydelxm] = ode23s(@ControlFunc,[0 5000],[0 0 0]',[],[dox(a) 0],del*.5,Kdox,kon,koff,ksgfp,kdgfp,Kdgfp,n2,n1,km,kg);
Gdelxp(a) = ydelxp(end,3);
Gdelxm(a) = ydelxm(end,3);

[t,y] = ode23s(@InsulatorFuncUnloaded,[0 5000],[0 0 0 0 0 0 0 0]',[],[dox(a),0],Kdox,del,kp,kpp,k1,k2,k3,k4,k5,k6,kon,koff,ksgfp,Kdgfp,kdgfp,k7,km,kg,n1,n2,k8,k9,k10,kon2);
[tc,yc] = ode23s(@ControlFunc,[0 5000],[0 0 0]',[],[dox(a) 0],del,Kdox,kon,koff,ksgfp,kdgfp,Kdgfp,n2,n1,km,kg);
G(a) = y(end,8);
Z(a) = y(end,1);
Zp(a) = y(end,2);
Gm(a) = yc(end,3);
end

figure
semilogx(dox,Gdelp,'--',dox,Gdelm,'--',dox,G)
xlabel('Dox [\mu M]')
ylabel('GFP [\mu M]')
legend('\delta +50%','\delta -50%','\delta')
title('Buffered system steady state change with \delta')
xlim([1e-3 1e2])

figure
semilogx(dox,Gkpp,'--',dox,Gkpm,'--',dox,G)
xlabel('Dox [\mu M]')
ylabel('GFP [\mu M]')
title('Steady State Change with k_p')
legend('k_p +50%','k_p -50%','k_p')
xlim([1e-3 1e2])

figure
semilogx(dox,Gk6p,'--',dox,Gk6m,'--',dox,G)
xlabel('Dox [\mu M]')
ylabel('GFP [\mu M]')
title('Buffered system steady state change with k_6')
legend('k_6 +50%','k_6 -50%','k_6')
xlim([1e-3 1e2])

figure
semilogx(dox,Gdelxp,'--',dox,Gdelxm,'--',dox,Gm)
xlabel('Dox [\mu M]')
ylabel('GFP [\mu M]')
title('Unbuffered system steady state change with \delta_M')
legend('\delta_M +50%','\delta_M -50%','\delta_M')
xlim([1e-3 1e2])
 