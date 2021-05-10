%% Generating Tube Shape Comparison
clc; clear; close all;

% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
SaveDestination = sprintf("RA-L_Figures/TipAnalysis");

% Choose Property Set
parameters = PropertySets(202,:);
% Pick desired theta readings and force readings
c = distinguishable_colors(10);
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
    wristType = wrist.name;
    
    
    %%
    % Getting S and kappa based on theta
    s_des = wrist.h - (wrist.ybar).*abs(theta_des);  %[m]
    kappa_des = theta_des./s_des; % [1/m]
    
    %% Plot Experimental Data
    wrist.theta = theta_des;
    wrist.s = s_des;
    wrist.F = force_input;
    wrist.kappa = kappa_des;
    color = c(5,:);
    [~,T] = wrist.robot_kin();
    wrist.plot_stick_model('Color',color);
    
    %% Non Linear and Friction
    wrist.E_lin = parameters.E_lin;
    wrist.E_se = parameters.E_se;
    wrist.strain_lower = parameters.Strain_Lower;
    wrist.mu = parameters.Mu;
    
    wrist.use_friction = true;
    wrist.use_non_linear = true;
    
    [~,T3] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = c(1,:);
    wrist.plot_stick_model('Marker','none','Color',color);
    %% Linear and NO friction
    wrist.E_lin = PropertySets.E_lin(198);
    wrist.mu = PropertySets.Mu(198);
    wrist.use_friction = false;
    wrist.use_non_linear = false;
    [~,T1] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = c(3,:);
    wrist.plot_stick_model('Marker','none','Color',color);
    
    %% Linear and friction
    wrist.E_lin = PropertySets.E_lin(198);
    wrist.mu = PropertySets.Mu(198);
    wrist.use_friction = true;
    wrist.use_non_linear = false;
    [~,T2] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = c(2,:);
    wrist.plot_stick_model('Marker','none','Color',color);
    
    %% Nonlinear without friction
    wrist.E_lin = PropertySets.E_lin(201);
    wrist.E_se = PropertySets.E_se(201);
    wrist.strain_lower = PropertySets.Strain_Lower(201);
    wrist.mu = PropertySets.Mu(201);
    wrist.use_friction = false;
    wrist.use_non_linear = true;
    [~,T4] = wrist.fwkin([force_input,0,0]);
    gca;
    hold on;
    color = c(4,:);
    wrist.plot_stick_model('Marker','none','Color',color);
    %%
    set(gcf,'Position',[100 100 500 500])
    set(gca, 'FontSize',16,'FontName','CMU Serif');
    view(0,0);
    xlim([0,7]);
    zlim([0,6]);
    axis equal
    title(sprintf('%s Shape Analysis',wrist.name),'FontSize',20,'FontName','CMU Serif');
    legend('Experiment','Non-Linear with Friction','Linear without Friction', 'Linear with Friction','NonLinear without Friction','FontSize',16,'Location','south','FontName','CMU Serif')
    
    
    destdirectory = sprintf("%s/",SaveDestination);
    if ~exist(destdirectory, 'dir')
        mkdir(destdirectory);
    end
    saveas(gcf,sprintf("%s/%s_TubeAnalysis.png",SaveDestination,wristType));
    saveas(gcf,sprintf("%s/%s_TubeAnalysis.fig",SaveDestination,wristType));
    saveas(gcf,sprintf("%s/%s_TubeAnalysis.svg",SaveDestination,wristType));
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