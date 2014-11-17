%% GreyBoxFitScript Control
clear all
clc
close all

%% Parameters
global km Kdox delc delx delg kon koff ksgfp kg Kgfp

%% Old parameters
%%Parameter    %%%Range                %%%%Reaction
Kdox = 3.024;   % Fitted                %Dox Km
km = 6e-5;      % Fitted                %Dox maximum velocity


kg = 8.8;       % Fitted                %Gfp maximum velocity

kon = 6;        % kon < 6               %Skn7/DNA association rate
koff = 0.006; 
Kgfp = koff/kon;% Fitted                %Gfp Km

ksgfp = 0.0012; % Fitted                %Bassal Gfp expression

delx = 0.013;                          %SKN7m degradation
delg = 0.004;                          %GFP degradation
delc = 0.004;   % Added                %Complex TF clearance

%% Calling Simulation Script
DynamicAnalysisControl_SI