ContDataRevised; %% Loading data
Inputs;
set(0,'DefaultAxesFontname', 'Times New Roman')
set(0,'DefaultAxesFontSize', 16)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 16)
set(0,'DefaultLineLinewidth',1)
set(0,'DefaultFigureColor','w')
%% Concentration of Skn7, Ypd1 and Load 
Xt = 0.0712; %2570/(6.02e23)*1/(60e-15);
pt = 0.022; %200/(6.02e23)*1/(60e-15); 800 site * 1/(6.02e23) (moles/sites) * 1/(60 e-15) (1/L)
Wt = 0.1752; %6330/(6.02e23)*1/(60e-15);
%% Steady State 
doxl = logspace(-3,2);
time = 1:5000;
for a = 1:size(doxl,2)
[t,y] = ode23s(@Control_Dilution,time,[0 0 ksgfp/(delg)]',[],doxl(a),0);
[tL1,yL1] = ode23s(@Control_Dilution,time,[0 0 ksgfp/(delg)]',[],doxl(a),pt/2);
[tL2,yL2] = ode23s(@Control_Dilution,time,[0 0 ksgfp/(delg)]',[],doxl(a),pt);
Gfp(a) = y(end,3);
GfpL1(a) = yL1(end,3);
GfpL2(a) = yL2(end,3);
end
figure
semilogx(doxl,Gfp,'k',doxl,GfpL2,'r')
xlabel('DOX')
ylabel('Fluoresence (A.U.)')
title('Unbuffered Steady State Fit')
legend('Model','Model+2x')

%% Step Up
time = 0:2000;
[tInd,yInd] = ode23s(@Control_Dilution,time,[0 0 ksgfp/(delg)],[],20,0);
[tIndL1,yIndL1] = ode23s(@Control_Dilution,time,[0 0 ksgfp/(delg)],[],20,pt/2);
[tIndL2,yIndL2] = ode23s(@Control_Dilution,time,[0 0 ksgfp/(delg)],[],20,pt);
GfpInd = yInd(:,3);
GfpIndL1 = yIndL1(:,3);
GfpIndL2 = yIndL2(:,3);

figure
plot(tInd,GfpInd,'k',tInd,GfpIndL1,'b',tInd,GfpIndL2,'r')
legend('Model C+0x','Model C+1x','Model C+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
title('Control 20 uM DOX Induction Model')


%% Step Down
time = 0:10000;
[tWash,yWash] = ode23s(@Control_Dilution,time,yInd(end,:),[],0,0);
[tWashL1,yWashL1] = ode23s(@Control_Dilution,time,yIndL1(end,:),[],0,pt/2);
[tWashL2,yWashL2] = ode23s(@Control_Dilution,time,yIndL2(end,:),[],0,pt);
GfpWash = yWash(:,3);
GfpWashL1 = yWashL1(:,3);
GfpWashL2 = yWashL2(:,3);


figure
plot(tWash,GfpWash,'k',tWash,GfpWashL1,'b',tWash,GfpWashL2,'r')
legend('Model C+0x','Model C+1x','Model C+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
title('Control 20 uM DOX Induction Model')
xlim([0 9000])

%% Dynamic Experiments
Experiment = [3 1 2 4 5];
data = {ContSq3,ContSq1,ContSq2,ContSq4,ContSq5,ContSq3L,ContSq1L,ContSq2L,ContSq4L,ContSq5L};
points = {s1'*1e-3  s2'*1e-3   s3'*1e-3   s4'*1e-3   s5'*1e-3   s1L'*1e-3   s2L'*1e-3   s3L'*1e-3   s4L'*1e-3   s5L'*1e-3  }; %% Unloaded System
tSq = {t1' t2' t3' t4' t5'};
tDox = {t1Dox' t2Dox' t3Dox' t4Dox' t5Dox'};
sDox = {DOX1' DOX2' DOX3' DOX4' DOX5'};

%% F = 350
[tSim,ySim] = ode23s(@ControlFunc350_Dilution,[0 1800],[0 0 ksgfp/delg],[],0);
Gfp350= ySim(:,3);
[tSimL,ySimL] = ode23(@ControlFunc350_Dilution,[0 1800],[0 0 ksgfp/delg],[],pt);
Gfp350L= ySimL(:,3);

figure
plot(tSim,Gfp350,'k',tSimL,Gfp350L,'r',tDox{4},sDox{4}*max(Gfp350*1.1))
title('Control F = 350 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')

save('350.mat')
%% F = 200
[tSim,ySim] = ode23s(@ControlFunc200_Dilution,[0 1600],[0 0 ksgfp/delg],[],0);
Gfp200= ySim(:,3);
[tSimL,ySimL] = ode23s(@ControlFunc200_Dilution,[0 1600],[0 0 ksgfp/delg],[],pt);
Gfp200L= ySimL(:,3);

figure
plot(tSim,Gfp200,'k',tSimL,Gfp200L,'r',tDox{3},sDox{3}*max(Gfp200*1.1))
title('Control F = 200 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])

save('200.mat')
%% F = 150
[tSim,ySim] = ode23s(@ControlFunc150_Dilution,[0 1600],[0 0 ksgfp/delg],[],0);
Gfp150= ySim(:,3);
[tSimL,ySimL] = ode23s(@ControlFunc150_Dilution,[0 1600],[0 0 ksgfp/delg],[],pt);
Gfp150L= ySimL(:,3);

figure
plot(tSim,Gfp150,'k',tSimL,Gfp150L,'r',tDox{2},sDox{2}*max(Gfp150*1.1))
title('Control F = 150 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])
 
save('150.mat')
%% F = 250
[tSim,ySim] = ode23s(@ControlFunc250_Dilution,[0 1600],[0 0 ksgfp/delg],[],0);
Gfp250= ySim(:,3);
[tSimL,ySimL] = ode23s(@ControlFunc250_Dilution,[0 1600],[0 0 ksgfp/delg],[],pt);
Gfp250L= ySimL(:,3);

figure
plot(tSim,Gfp250,'k',tSimL,Gfp250L,'r',tDox{1},sDox{1}*max(Gfp250*1.1))
title('Control F = 250 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])

save('250.mat')

%% F = 500
[tSim,ySim] = ode23s(@ControlFunc500_Dilution,[0 1600],[0 0 ksgfp/delg]',[],0);
Gfp500= ySim(:,3);
[tSimL,ySimL] = ode23(@ControlFunc500_Dilution,[0 1600],[0 0 ksgfp/delg]',[],pt);
Gfp500L= ySimL(:,3);

figure
plot(tSim,Gfp500,'k',tSimL,Gfp500L,'r',tDox{5},sDox{5}*max(Gfp500*1.1))
title('Control F = 500 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])