%% Generating Tube Shape Comparison
clc; clear; close all;

% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
% Choose Property Set 
parameters = PropertySets(196,:);
% Pick desired theta readings and force readings

%% For Tip First Image: DSC 2320
n = 4;
theta_des = deg2rad([12.95479079; 21.02054754; 35.44172543; 37.65617083]); %[deg]
force_input = 1.342364; %[N]
% Make Wrist
% wrist = MakeWrist('TipFirstTube');

%% For 150 Tube
% n = 5;
% theta_des = deg2rad([	19.81294513	17.02661081	15.41854107	13.90919924	14.35484336]'); %[deg]
% force_input = 1.198876; %[N]
% % Make Wrist

%% For 90 Tube
n = 5;
theta_des = deg2rad([17.02091355	15.54938205	14.6406139	13.33544722	13.30206208]'); %[deg]
force_input = 1.362315; %[N]
% Make Wrist
wrist = MakeWrist('90Tube');


%%
% Getting S and kappa based on theta
wrist.E_lin = parameters.E_lin;
wrist.E_se = parameters.E_se;
wrist.strain_lower = parameters.Strain_Lower;
wrist.mu = parameters.Mu;
s_des = wrist.h - (wrist.ybar).*abs(theta_des);  %[m]
kappa_des = theta_des./s_des; % [1/m]

%% Plot Experimental Data

wrist.theta = theta_des;
wrist.s = s_des;
wrist.F = force_input;
wrist.kappa = kappa_des;

[~,T] = wrist.robot_kin();
wrist.plot_stick_model;

%% Linear and NO friction
% wrist.E_lin = (parameters.E_lin+parameters.E_se)/2;
wrist.use_friction = false;
wrist.use_non_linear = false;
[~,T1] = wrist.fwkin([force_input,0,0]);
gca;
hold on;
wrist.plot_stick_model;

%% Linear and friction 
% wrist.E_lin = (parameters.E_lin+parameters.E_se)/2;
wrist.use_friction = true;
wrist.use_non_linear = false;
[~,T2] = wrist.fwkin([force_input,0,0]);
gca;
hold on;
wrist.plot_stick_model;

%% Non Linear and Friction
wrist.E_lin = parameters.E_lin;
wrist.use_friction = true;
wrist.use_non_linear = true;
[~,T3] = wrist.fwkin([force_input,0,0]);
gca;
hold on;
wrist.plot_stick_model;

view(0,0);
title('90 Tube Comparisons of different kinematic algorithms with E_{lin} = 10GPa','FontSize',16)'
legend('Experiment','No Friction All Linear', 'Friction All Linear','Full Model','FontSize',14)