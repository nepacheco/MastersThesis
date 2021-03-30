%% Generating Tube Shape Comparison
clc; clear; close all;

% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
% Choose Property Set 
parameters = PropertySets(196,:);
% Pick desired theta readings and force readings

%% For Tip First Image: 2324
n = 4;
theta_des = deg2rad([31.5282 38.0678 37.6685 37.6272]'); %[deg]
force_input = 1.8852; %[N]
% Make Wrist
wrist = MakeWrist('TipFirstTube',true);

%% For 150 Tube
% n = 5;
% theta_des = deg2rad([31.01822282	30.25792933	25.50073754	19.95987288	19.03623321]'); %[deg]
% force_input = 1.55515; %[N]
% % Make Wrist
% wrist = MakeWrist('150Tube',true);

%% For 90 Tube
% n = 5;
% theta_des = deg2rad([17.02091355	15.54938205	14.6406139	13.33544722	13.30206208]'); %[deg]
% force_input = 1.362315; %[N]
% % Make Wrist
% wrist = MakeWrist('90Tube');


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
wrist.E_lin = 7.51E9;
wrist.mu = 0.05;
wrist.use_friction = false;
wrist.use_non_linear = false;
[~,T1] = wrist.fwkin([force_input,0,0]);
gca;
hold on;
wrist.plot_stick_model;

%% Linear and friction 
wrist.E_lin = 7.51E9;
wrist.mu = 0.05;
wrist.use_friction = true;
wrist.use_non_linear = false;
[~,T2] = wrist.fwkin([force_input,0,0]);
gca;
hold on;
wrist.plot_stick_model;

%% Non Linear and Friction
wrist.E_lin = parameters.E_lin;
wrist.mu = parameters.Mu;
wrist.use_friction = true;
wrist.use_non_linear = true;
[~,T3] = wrist.fwkin([force_input,0,0]);
gca;
hold on;
wrist.plot_stick_model;

view(0,0);
title('Tip First Tube Comparisons of different kinematic algorithms with E_{lin} = 7.51GPa','FontSize',16)'
legend('Experiment','No Friction All Linear', 'Friction All Linear','Full Model','FontSize',14)