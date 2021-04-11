%% Generating Tube Shape Comparison
clc; clear; close all;

% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
% Choose Property Set
parameters = PropertySets(196,:);
% Pick desired theta readings and force readings

for i = 1:3
    figure()
    %% For Tip First Image: 2324
    if i == 1
        n = 4;
        theta_des = deg2rad([31.5282 38.0678 37.6685 37.6272]'); %[deg]
        force_input = 1.8852; %[N]
        % Make Wrist
        wrist = MakeWrist('TipFirstTube',true);
        wrist.name = 'Wrist C';
    end
    %% For 150 Tube
    if i == 2
        n = 5;
        theta_des = deg2rad([31.01822282	30.25792933	25.50073754	19.95987288	19.03623321]'); %[deg]
        force_input = 1.55515; %[N]
%         theta_des = deg2rad([6.7270	7.6068	6.7920	6.4308	7.9189]');
%         force_input = 0.4942;
        % Make Wrist
        wrist = MakeWrist('150Tube',true);
        wrist.name = 'Wrist A';
    end
    %% For 90 Tube : 2799
    if i == 3
        n = 5;
        theta_des = deg2rad([18.9249 19.6790 18.2729 16.1427 16.3920]'); %[deg]
        force_input = 1.6002; %[N]
        % Make Wrist
        wrist = MakeWrist('90Tube',true);
        wrist.name = 'Wrist B';
    end
    
    
    
    %%
    % Getting S and kappa based on theta
    wrist.E_lin = parameters.E_lin;
    wrist.E_se = parameters.E_se;
    wrist.strain_lower = parameters.Strain_Lower;
    wrist.mu = parameters.Mu;
    s_des = wrist.h - (wrist.ybar).*abs(theta_des);  %[m]
    kappa_des = theta_des./s_des; % [1/m]
    
    %% Plot Experimental Data
    wrist.E_lin = parameters.E_lin;
    wrist.mu = parameters.Mu;
    wrist.use_friction = true;
    wrist.use_non_linear = true;
    wrist.theta = theta_des;
    wrist.s = s_des;
    wrist.F = force_input;
    wrist.kappa = kappa_des;
    color = [0.4940 0.1840 0.5560];
    [~,T] = wrist.robot_kin();
    wrist.plot_stick_model('Color',color);
    
    %% Non Linear and Friction
    wrist.E_lin = parameters.E_lin;
    wrist.mu = parameters.Mu;
    wrist.use_friction = true;
    wrist.use_non_linear = true;
    [~,T3] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = [0 0.4470 0.7410];
    wrist.plot_stick_model('Marker','none','Color',color);
    %% Linear and NO friction
    wrist.E_lin = 7.51E9;
    wrist.mu = 0.05;
    wrist.use_friction = false;
    wrist.use_non_linear = false;
    [~,T1] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = [0.8500 0.3250 0.0980];
    wrist.plot_stick_model('Marker','none','Color',color);
    
    %% Linear and friction
    wrist.E_lin = 7.51E9;
    wrist.mu = 0.05;
    wrist.use_friction = true;
    wrist.use_non_linear = false;
    [~,T2] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = [0.9290 0.6940 0.1250];
    wrist.plot_stick_model('Marker','none','Color',color);
    
    %%
    set(gcf,'Position',[100 100 500 500])
    set(gca, 'FontSize',16,'FontName','CMU Serif');
    view(0,0);
    xlim([0,7]);
    zlim([0,6]);
    axis equal
    title(sprintf('%s Shape Analysis',wrist.name),'FontSize',20,'FontName','CMU Serif');
%     legend('Experiment','Non-Linear with Friction','Linear without Friction', 'Linear with Friction','FontSize',16,'Location','south','FontName','CMU Serif')
end

%% Test
% wrist.E_lin = parameters.E_lin;
% wrist.mu = parameters.Mu;
% wrist.use_friction = true;
% wrist.use_non_linear = true;
% [~,T3] = wrist.fwkin([2.05,0,0]);
% gca;
% hold on;
% wrist.plot_stick_model;
% view(0,0);