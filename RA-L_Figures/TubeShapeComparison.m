%% Generating Tube Shape Comparison
clc; clear; close all;

% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
SaveDestination = sprintf("RA-L_Figures/TipAnalysis");

% Choose Property Set
parameters = PropertySets(196,:);
% Pick desired theta readings and force readings
c = distinguishable_colors(10);
for i = 1:3
    figure()
    %% For Tip First Image
    if i == 1
        n = 4;
        % Image 2324
%         theta_des = deg2rad([31.5282 38.0678 37.6685 37.6272]'); %[deg]
%         force_input = 1.8852; %[N]
%         tendon_disp = 5E-3; % [m]
        % Image 2870
        theta_des = deg2rad([8.0606	10.9139	17.3038	28.8193]');
        force_input = 0.9138;
        tendon_disp = 2.25E-3; % [m]
        % Make Wrist
        wrist = MakeWrist('TipFirstTube',true);
        wrist.name = 'Wrist C';
    end
    %% For 150 Tube
    if i == 2
        n = 5;
        % Image 2258
%         theta_des = deg2rad([31.01822282	30.25792933	25.50073754	19.95987288	19.03623321]'); %[deg]
%         force_input = 1.55515; %[N]
%         tendon_disp = 4.5E-3; % [m]
        % Image 2871
        theta_des = deg2rad([12.9007	11.6051	11.4855	10.4380	11.5129]');
        force_input = 0.8840;
        tendon_disp = 2.25E-3; % [m]
        % Make Wrist
        wrist = MakeWrist('150Tube',true);
        wrist.name = 'Wrist A';
    end
    %% For 90 Tube 
    if i == 3
        n = 5;
        % Image 2799
%         theta_des = deg2rad([18.9249 19.6790 18.2729 16.1427 16.3920]'); %[deg]
%         force_input = 1.6002; %[N]
%         tendon_disp = 3.25E-3; % [m] 
        % Image 2808
        theta_des = deg2rad([10.2108	10.1908	9.6146	9.2993	9.7316]'); %[deg]
        force_input = 1.0876; %[N]
        tendon_disp = 2E-3; % [m] 
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
    color = c(4,:);
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
%     wrist.E_lin = PropertySets.E_lin(201);
%     wrist.E_se = PropertySets.E_se(201);
%     wrist.strain_lower = PropertySets.Strain_Lower(201);
%     wrist.mu = PropertySets.Mu(201);
%     wrist.use_friction = false;
%     wrist.use_non_linear = true;
%     [~,T4] = wrist.fwkin([force_input,0,0]);
%     gca;
%     hold on;
%     color = c(4,:);
%     wrist.plot_stick_model('Marker','none','Color',color);
    
    %% Geometry
    wrist.E_lin = parameters.E_lin;
    wrist.E_se = parameters.E_se;
    wrist.strain_lower = parameters.Strain_Lower;
    wrist.mu = parameters.Mu;
    
    wrist.use_friction = true;
    wrist.use_non_linear = true;
    
    [~,T3] = wrist.fwkin([tendon_disp,0,0],'Type','geometry');
    gca;
    hold on;
    color = c(5,:);
    wrist.plot_stick_model('Marker','none','Color',color);
    
    %%
%     set(gcf,'Position',[100 100 500 500])
    set(gca, 'FontSize',16,'FontName','CMU Serif');
    view(0,0)
    zlim([0,10])
    xlim([0,7.5])
    axis equal
    title(sprintf('%s Shape Analysis',wrist.name),'FontSize',20,'FontName','CMU Serif');
    legend('Experiment','Our Model',...
        'Linear without Friction', 'Linear with Friction',...
        'York et al. Model',...
        'FontSize',16,'Location','south','FontName','CMU Serif')
    
    
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