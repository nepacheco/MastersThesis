%% Setup
clc; clear; close all;
load('PropertySets.mat')
parameters = PropertySets(198,:);
wristTypes = {'150Tube','90Tube','TipFirstTube'};
wristNames = {'A', 'B', 'C'};
expFiles = ["\27-6-9_Trial1_Results.xlsx","\02-17-2021_Results.xlsx","\02-03-2021_Results.xlsx"];
fontsize = 20;
markerSize = 500;
Force = 2.5;

az = -32.7;
el = 33.8;

c = distinguishable_colors(50);
%%
for w = 1 : 3
    wristType = wristTypes{w};
    SaveDestination = sprintf("RA-L_Figures/%s",wristType);
    experimentFiles = expFiles(w);
    wrist = MakeWrist(wristType,true);
    wrist.use_non_linear = false;
    
    
    %% Parse Experiment Files
    numFiles = size(experimentFiles,1);
    force_cell = cell(1,numFiles);
    notch_cell = cell(1,numFiles);
    experimentStr = "";
    average = zeros(wrist.n+1,numFiles);
    for i = 1:numFiles
        % Parsing File
        backslash_indicies = strfind(experimentFiles(i),"\");
        period_indicies = strfind(experimentFiles(i),".");
        experimentStr = experimentStr + extractBetween(experimentFiles(i),backslash_indicies(end)+1,period_indicies(end)-1);
        opts = detectImportOptions(experimentFiles(i));
        opts.Sheet = 'AvgMeasurements';
        file = readcell(experimentFiles(i),opts);
        [force_vec, notch_data] = ParseExperimentFile(file,wrist.n);
        
        force_cell(:,i) = {force_vec};
        notch_cell(:,i) = {notch_data};
        average(:,i) = mean(notch_data)';
    end
    %% Plotting is done for each set of material properties
    for m = 1:size(parameters,1)
        %% Defining the material properties of the wrist for this experiment
        wrist.E_lin = table2array(parameters(m,'E_lin'));
        wrist.E_se = table2array(parameters(m,'E_se'));
        wrist.strain_lower = table2array(parameters(m,'Strain_Lower'));
        wrist.mu = table2array(parameters(m,'Mu'));
        
        %% Plot The Experiments along with the desired parameters
        points = 100;
        theta_mat_force = zeros(wrist.n,points);
        F_vec = linspace(0,Force,points); % Getting list of forces (x axis values)
        for k = 1:points
            wrist.fwkin([F_vec(k),0,0],'Type','force');
            theta_mat_force(:,k) = wrist.theta;
        end
        
        plot_options = {'k.','bx','r*','mo'};
        
        figure, hold on
        %        l = {};
        
        for p = 1:wrist.n
            %title(sprintf("Notch %d",p),'FontSize',fontsize);
            %hold on
            for i = 1:numFiles
                scatter3(p*ones(1,length(force_cell{1,i})),force_cell{1,i},notch_cell{1,i}(:,p),50,.2*ones(1,3),'filled');
                hold on
            end
            plot3(p*ones(1,length(F_vec)),F_vec,rad2deg(theta_mat_force(p,:)),'Color',c(p+5,:),'Linewidth',3);
            hold on
            xlabel('Notch number')
            ylabel("Force (N)",'FontSize',fontsize);
            zlabel("Deflection (deg)",'FontSize',fontsize)
            ax = gca;
            set(ax,'FontSize',fontsize)
            axis tight
            
            ylim([0,Force]);
            zlim([0 inf]);
            
            if w == 2
               ylim([0,2]); 
            end
            
            if w == 3
                xlim([1 4]);
                zlim([0 40]);
                xticks(1:4);
                zticks([0 10 20 30 40]);
            elseif w == 2
                xlim([1 5]);
                zlim([0 20]);
                xticks(1:5);
                zticks([0 10 20 30]);
            else 
                xlim([1 5]);
                zlim([0 35]);
                xticks(1:5);
                zticks([0 10 20 30 35]);
                zticklabels({'0', '10', '20', '30', ''});
            end
            
            view(az,el);
            %ylim([0 inf]);
            %xlim([0,Force]);
            %ax.PlotBoxAspectRatio= [1,0.5,1];
            %axis tight
            ax = gca;
            %ax.PlotBoxAspectRatio= [0.5,1,0.5];
            set(gca, 'FontName', 'CMU Serif');
            grid on
        end
        
        title(['Exp. Data vs. Model - Wrist ' wristNames{w}]);
        hold off
        
        % Saving the figure
        destdirectory = sprintf("%s/",SaveDestination);
        if ~exist(destdirectory, 'dir')
           mkdir(destdirectory);
        end
    saveas(gcf,sprintf("%s/%s_3D%s.png",SaveDestination,wristType,experimentStr));
    saveas(gcf,sprintf("%s/%s_3D%s.fig",SaveDestination,wristType,experimentStr));
    saveas(gcf,sprintf("%s/%s_3D%s.svg",SaveDestination,wristType,experimentStr));
    end
end