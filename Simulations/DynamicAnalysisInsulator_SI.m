%% Dynamic Analysis SI
set(0,'DefaultAxesFontname', 'Times New Roman')
set(0,'DefaultAxesFontSize', 20)
set(0,'DefaultTextFontname', 'Times New Roman')
set(0,'DefaultTextFontSize', 20)
set(0,'DefaultLineLinewidth',1)
set(0,'DefaultFigureColor','w')

InsData; Inputs;

Xt = 0.0712; %2570/(6.02e23)*1/(60e-15);
pt = 0.0220; %200/(6.02e23)*1/(60e-15); 200 site * 1/(6.02e23) (moles/sites) * 1/(60 e-15) (1/L)
Wt = 0.1752; %6330/(6.02e23)*1/(60e-15);
%% Steady State

doxl = logspace(-3,2);
time = 1:5000;
for a = 1:size(doxl,2)
[tL0,yL0] = ode23s(@InsulatorFunc_Dilution,time,[0 0 0.4*Wt 0 0.3*Xt 0 0 0 0 0.3063],[],doxl(a),0,Xt*.3,Wt*.4);
[tL1,yL1] = ode23s(@InsulatorFunc_Dilution,time,[0 0 0.4*Wt 0 0.3*Xt 0 0 0 0 0.3063],[],doxl(a),pt/2,Xt*.3,Wt*.4);
[tL2,yL2] = ode23s(@InsulatorFunc_Dilution,time,[0 0 0.4*Wt 0 0.3*Xt 0 0 0 0 0.3063],[],doxl(a),pt,Xt*.3,Wt*.4);
GfpL0(a) =  yL0(end,10);
GfpL1(a) =  yL1(end,10);
GfpL2(a) =  yL2(end,10);


Zm(a) = yL0(end,1) - yL0(end,2);
Zpm(a)= yL0(end,2);
Wpm(a) = yL0(end,4);
Xpm(a) = yL0(end,6);
Xppm(a) = yL0(end,7);

end

figure
semilogx(doxl,GfpL0,'k',doxl,GfpL2,'r')
xlabel('DOX')
ylabel('GFP [\muM]')
title('Insulattor  Steady State Fit')
legend('Model I+0x','Model I+2x')

%% Step Up
[tInd,yInd] = ode23s(@InsulatorFunc_Dilution,0:1000,[0 0 Wt 0 Xt 0 0 0 0 0.3063]',[],20,0,Xt,Wt);
[tIndL1,yIndL1] = ode23s(@InsulatorFunc_Dilution,0:1000,[0 0 Wt 0 Xt 0 0 0 0 0.3063]',[],20,pt/2,Xt,Wt);
[tIndL2,yIndL2] = ode23s(@InsulatorFunc_Dilution,0:1000,[0 0 Wt 0 Xt 0 0 0 0 0.3063]',[],20,pt,Xt,Wt);

figure
plot(tInd,yInd(:,10),'k',tInd,yIndL1(:,10),'b',tInd,yIndL2(:,10),'r')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
legend('Model I+0x','Model I+1x','Model I+2x')
title('Insulator 20 uM DOX Induction Model')

%%% Reduced Xt,Wt Step up
[tu0x,yu0x] = ode23s(@InsulatorFunc_Dilution,0:1000,[0 0 Wt*.4 0 Xt*.3 0 0 0 0 0.3063]',[],20,0,Xt*.3,Wt*.4);
[tu1x,yu1x] = ode23s(@InsulatorFunc_Dilution,0:1000,[0 0 Wt*.4 0 Xt*.3 0 0 0 0 0.3063]',[],20,pt/2,Xt*.3,Wt*.4);
[tu2x,yu2x] = ode23s(@InsulatorFunc_Dilution,0:1000,[0 0 Wt*.4 0 Xt*.3 0 0 0 0 0.3063]',[],20,pt,Xt*.3,Wt*.4);

figure
plot(tu0x,yu0x(:,10),'k',tu1x,yu1x(:,10),'b',tu2x,yu2x(:,10),'r')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
legend('Model I+0x','Model I+1x','Model I+2x')
title('30% SKN7 and 40% YPD1')

figure
plot(tu0x,(yu0x(:,10) - min(yu0x(:,10)))/max(yu0x(:,10) - min(yu0x(:,10))),'k',...
    tu1x,(yu1x(:,10) - min(yu1x(:,10)))/max(yu1x(:,10) - min(yu1x(:,10))),'b',...
    tu2x,(yu2x(:,10) - min(yu2x(:,10)))/max(yu2x(:,10) - min(yu2x(:,10))),'r')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
legend('Model I+0x','Model I+1x','Model I+2x')
title('30% SKN7 and 40% YPD1 - Normalized')

%% Step Down
[tWash,yWash] = ode23s(@InsulatorFunc_Dilution,0:2100,yInd(end,:),[],0,0,Xt,Wt);
[tWashL1,yWashL1] = ode23s(@InsulatorFunc_Dilution,0:2100,yIndL1(end,:),[],0,pt/2,Xt,Wt);
[tWashL2,yWashL2] = ode23s(@InsulatorFunc_Dilution,0:2100,yIndL2(end,:),[],0,pt,Xt,Wt);

GfpWash= yWash(:,10);
GfpWashL1 = yWashL1(:,10);
GfpWashL2 = yWashL2(:,10);

figure
plot(tWash,GfpWash,'k',tWash,GfpWashL1,'b',tWash,GfpWashL2,'r')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
legend('Model I+0x','Model I+1x','Model I+2x')
title('Step Down')

%%% Reduced Xt,Wt wash

[td0x,yd0x] = ode23s(@InsulatorFunc_Dilution,0:2100,yu0x(end,:),[],0,0,Xt*.3,Wt*.4);
[td1x,yd1x] = ode23s(@InsulatorFunc_Dilution,0:2100,yu1x(end,:),[],0,pt/2,Xt*.3,Wt*.4);
[td2x,yd2x] = ode23s(@InsulatorFunc_Dilution,0:2100,yu2x(end,:),[],0,pt,Xt*.3,Wt*.4);

figure
plot(td0x,yd0x(:,10),'k',td1x,yd1x(:,10),'b',td2x,yd2x(:,10),'r')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
legend('Model I+0x','Model I+1x','Model I+2x')
title('30% SKN7 and 40% YPD1')

figure
plot(td0x,(yd0x(:,10) - min(yd0x(:,10)))/max(yd0x(:,10) - min(yd0x(:,10))),'k',...
    td1x,(yd1x(:,10) - min(yd1x(:,10)))/max(yd1x(:,10) - min(yd1x(:,10))),'b',...
    td2x,(yd2x(:,10) - min(yd2x(:,10)))/max(yd2x(:,10) - min(yd2x(:,10))),'r')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
legend('Model I+0x','Model I+1x','Model I+2x')
title('30% SKN7 and 40% YPD1 - Normalized')


%% Square Wave input

Experiment = [3 1 2 4 5];
data = {InsSq3,InsSq1,InsSq2,InsSq4,InsSq5,InsSq3L,InsSq1L,InsSq2L,InsSq4L,InsSq5L};
points = {s1' s2' s3' s4' s5' s1L' s2L' s3L' s4L' s5L'}; 
tSq = {t1' t2' t3' t4' t5'};
tDox = {t1Dox' t2Dox' t3Dox' t4Dox' t5Dox'};
sDox = {DOX1' DOX2' DOX3' DOX4' DOX5'};

%% 350
[tSim,ySim] = ode23s(@InsulatorFunc350_Dilution,0:.1:1800,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],0);
Gfp350= ySim(:,10);
[tSimL,ySimL] = ode23s(@InsulatorFunc350_Dilution,0:.1:1800,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],pt);
Gfp350L= ySimL(:,10);

figure
plot(tSim,Gfp350,'k',tSimL,Gfp350L,'r',tDox{4},sDox{4}*max(Gfp350*1.1))
title('Insulattor  F = 350 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1400])

save('350.mat');
%% 200
[tSim,ySim] = ode23s(@InsulatorFunc200_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],0);
Gfp200= ySim(:,10);
[tSimL,ySimL] = ode23s(@InsulatorFunc200_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],pt);
Gfp200L= ySimL(:,10);

figure
plot(tSim,Gfp200,'k',tSimL,Gfp200L,'r',tDox{3},sDox{3}*max(Gfp200*1.1))
title('Insulattor  F = 200 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])

save('200.mat');
%% 150
[tSim,ySim] = ode23s(@InsulatorFunc150_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],0);
Gfp150= ySim(:,10);
[tSimL,ySimL] = ode23s(@InsulatorFunc150_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],pt);
Gfp150L= ySimL(:,10);

figure
plot(tSim,Gfp150,'k',tSimL,Gfp150L,'r',tDox{2},sDox{2}*max(Gfp150*1.1))
title('Insulattor  F = 150 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])

save('150.mat');
%% 250
[tSim,ySim] = ode23s(@InsulatorFunc250_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],0);
Gfp250= ySim(:,10);
[tSimL,ySimL] = ode23s(@InsulatorFunc250_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],pt);
Gfp250L= ySimL(:,10);

figure
plot(tSim,Gfp250,'k',tSimL,Gfp250L,'r',tDox{1},sDox{1}*max(Gfp250*1.1))
title('Insulattor  F = 250 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])

save('250.mat');
%% 500
time = 0:.1:1600;
[tSim,ySim] = ode23s(@InsulatorFunc500_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],0);
Gfp500= ySim(:,10);
[tSimL,ySimL] = ode23s(@InsulatorFunc500_Dilution,0:.1:1600,[0 0 Wt 0 Xt 0 0 0 0 0.3063],[],pt);
Gfp500L= ySimL(:,10);


figure
plot(tSim,Gfp500,'k',tSimL,Gfp500L,'r',tDox{5},sDox{5}*max(Gfp500*1.1))
title('Insulattor  F = 500 min')
legend('Model I+0x','Model I+2x')
ylabel('Fluoresence (A.U.)')
xlabel('Time (min)')
xlim([0 1600])

save('500.mat');
return