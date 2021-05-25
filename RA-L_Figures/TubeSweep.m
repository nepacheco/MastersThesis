clc; clear; close all;

% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
SaveDestination = sprintf("RA-L_Figures/TipAnalysis");

% Choose Property Set
parameters = PropertySets(196,:);
for w = 1:3
    figure()
    %% For Tip First Image: 2324
    if w == 1
        legendlocation = 'northeast';
        n = 4;
        % Make Wrist
        wrist = MakeWrist('TipFirstTube',true);
        file_name = "RA-L_Figures/02-03-2021_Results.xlsx";
        wrist.name = 'Wrist C';
    end
    %% For 150 Tube
    if w == 2
        legendlocation = 'northeast';
        n = 5;
        wrist = MakeWrist('150Tube',true);
        file_name = "RA-L_Figures/27-6-9_Trial1_Results.xlsx";
        wrist.name = 'Wrist A';
    end
    %% For 90 Tube : 2799
    if w == 3
        legendlocation = 'south';
        n = 5;
        % Make Wrist
        wrist = MakeWrist('90Tube',true);
        file_name = "RA-L_Figures/02-17-2021_Results.xlsx";
        wrist.name = 'Wrist B';
    end
    wristType = wrist.name;
    opts = detectImportOptions(file_name);
    opts.Sheet = 'AvgMeasurements';
    file = readcell(file_name,opts);
    [force_vec,notch_mat,tendon_disp] = ParseExperimentFile(file,n);
    num_elem = 2;
    for i = 2:size(force_vec,1)
        if force_vec(i) < force_vec(i-1)
            num_elem = i-1;
            break;
        end
    end
    c = distinguishable_colors(10);
    num_examples = 9;
    counter = 0;
    for i = num_elem:-ceil(num_elem/num_examples):1
        counter = counter + 1;
%         fprintf("Notch Angles %f ", notch_mat(i,1:end-1))
%         fprintf("\nForce Value %f \n",force_vec(i));
        theta_des = deg2rad(notch_mat(i,1:end-1)');
        s_des = wrist.h - (wrist.ybar).*abs(theta_des);  %[m]
        kappa_des = theta_des./s_des; % [1/m]

        %% Plot Experimental Data
        wrist.theta = theta_des;
        wrist.s = s_des;
        wrist.F = force_vec(i);
        wrist.kappa = kappa_des;
        color = c(4,:);
        color = color + (ones(1,3)-color)/(num_examples*1.05)*(counter); % Lets me apply a gradient
        [~,T] = wrist.robot_kin();
        hold on
        wrist.plot_stick_model('Color',color);


        %% Non Linear and Friction
        wrist.E_lin = parameters.E_lin;
        wrist.E_se = parameters.E_se;
        wrist.strain_lower = parameters.Strain_Lower;
        wrist.mu = parameters.Mu;

        wrist.use_friction = true;
        wrist.use_non_linear = true;

        [~,T3] = wrist.fwkin([force_vec(i),0,0]);
        gca;
        hold on;
        color = c(1,:);
        color = color + (ones(1,3)-color)/(num_examples*1.05)*(counter); % Lets me apply a gradient
        wrist.plot_stick_model('Marker','none','Color',color);

        %% Geometry
%         wrist.E_lin = parameters.E_lin;
%         wrist.E_se = parameters.E_se;
%         wrist.strain_lower = parameters.Strain_Lower;
%         wrist.mu = parameters.Mu;
% 
%         wrist.use_friction = true;
%         wrist.use_non_linear = true;
% 
%         [~,T3] = wrist.fwkin([tendon_disp(i)*1E-3,0,0],'Type','geometry');
%         gca;
%         hold on;
%         color = c(5,:);
%         color = color + (ones(1,3)-color)/(num_examples*1.05)*(counter); % Lets me apply a gradient
%         wrist.plot_stick_model('Marker','none','Color',color);

    end
    view(0,0)
    zlim([0,10])
    xlim([0,7.5])
    set(gca, 'FontSize',16,'FontName','CMU Serif');
    title(sprintf('%s Bending',wrist.name),'FontSize',20,'FontName','CMU Serif');
    legend('Experiment','Our Model','FontSize',16,'Location',legendlocation,'FontName','CMU Serif')
    
    
    SaveDestination = "RA-L_Figures/TotalBend";
    destdirectory = sprintf("%s/",SaveDestination);
    if ~exist(destdirectory, 'dir')
        mkdir(destdirectory);
    end
    saveas(gcf,sprintf("%s/%s_TotalBend.png",SaveDestination,wristType));
    saveas(gcf,sprintf("%s/%s_TotalBend.fig",SaveDestination,wristType));
    saveas(gcf,sprintf("%s/%s_TotalBend.svg",SaveDestination,wristType));
end
