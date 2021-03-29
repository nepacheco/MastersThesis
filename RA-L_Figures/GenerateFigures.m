%% Setup
load('ExperimentFiles.mat');
load('PropertySets.mat');
wristTypes = {'150Tube','90Tube','TipFirstTube'};
expFiles = [experimentFiles150(end),experimentFiles90(end),experimentFilesTip(end)];
parameters = PropertySets(196,:);
fontsize = 22;
markerSize = 100;
Force = 2.5;
%%
for w = 1:3
    wristType = wristTypes{w};
    SaveDestination = sprintf("RA-L_Figures/%s",wristType);
    experimentFiles = expFiles(w);
    wrist = MakeWrist(wristType,true);
    
    
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
        figure('WindowState','Maximize');
        for p = 1:wrist.n+1
            subplot(2,3,p);
            if p < wrist.n+1
                title(sprintf("Notch %d",p),'FontSize',fontsize);
                hold on
                for i = 1:numFiles
                    scatter(force_cell{1,i},notch_cell{1,i}(:,p),markerSize,plot_options{i});
                end
                plot(F_vec,rad2deg(theta_mat_force(p,:)),'Color','#EDB120','Linewidth',3);
                xlabel("Force (N)",'FontSize',fontsize);
                ylabel("Deflection (deg)",'FontSize',fontsize)
                ax = gca;
                set(ax,'FontSize',fontsize)
                axis tight
                ylim([0 inf]);
                xlim([0,Force]);
                ax.PlotBoxAspectRatio= [1,0.5,1];
                grid on
            else
                title(sprintf("Total Deflection"),'FontSize',fontsize);
                hold on
                legend_entries = cell(1,numFiles+1);
                for i = 1:numFiles
                    scatter(force_cell{1,i},notch_cell{1,i}(:,p),markerSize,plot_options{i});
                end
                plot(F_vec,rad2deg(sum(theta_mat_force(:,:))),'Color','#EDB120','Linewidth',3);
                xlabel("Force (N)",'FontSize',fontsize);
                ylabel("Deflection (deg)",'FontSize',fontsize)
                set(gca,'FontSize',fontsize)
                axis tight
                ylim([0 inf]);
                xlim([0,Force]);
                ax = gca;
                ax.PlotBoxAspectRatio= [1,0.5,1];
                grid on
            end
            hold off
        end
        
        
        %% Saving the figure
        destdirectory = sprintf("%s/",SaveDestination);
        if ~exist(destdirectory, 'dir')
            mkdir(destdirectory);
        end
        saveas(gcf,sprintf("%s/%s_%s.png",SaveDestination,wristType,experimentStr));
        saveas(gcf,sprintf("%s/%s_%s.fig",SaveDestination,wristType,experimentStr));
    end
end