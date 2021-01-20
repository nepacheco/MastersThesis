function [rmse_total] = CompareModel(wristType,experimentFiles,usePrecurve,parameters)
%COMPAREMODEL Compares the model to a set of experimental data using the
%given parameters
%   Wrist is the wrist class to use
%   experimentFiles is a vector containing the absolute path to each
%   experiment excel file
arguments
    wristType char {mustBeMember(wristType,{'90Tube','150Tube','TipFirstTube'})}
    experimentFiles (:,1) string
    usePrecurve logical
    parameters table
end

wrist = MakeWrist(wristType,usePrecurve);
% Parse Experiment Files
numFiles = size(experimentFiles,1);
for i = 1:numFiles
    % Parsing File
    opts = detectImportOptions(experimentFiles(i));
    opts.Sheet = 'AvgMeasurements';
    file = readcell(experimentFiles(i),opts);
    [force_vec notch_data] = ParseExperimentFile(file,wrist.n);
    
    force_mat(:,i) = force_vec;
    notch_mat(:,:,i) = notch_data;
end

rmse_total = zeros(wrist.n+1,numFiles,size(parameters,1));
% Calculating RMSE and Plotting is done for each set of material properties
for m = 1:size(parameters,1)
    % Defining the material properties of the wrist for this experiment
    wrist.E_lin = table2array(parameters(m,'E_lin'));
    wrist.E_se = table2array(parameters(m,'E_se'));
    wrist.strain_lower = table2array(parameters(m,'Strain_Lower'));
    wrist.mu = table2array(parameters(m,'Mu'));
    % Generating RMSE
    for i = 1:numFiles
        diff = zeros(wrist.n+1,length(force_mat(:,i)));
        for j = 1:length(force_mat(:,i))
            input = force_mat(j,i);
            wrist.fwkin([input,0,0],'Type','force');
            diff(:,j) = notch_mat(j,:,i)' - rad2deg([wrist.theta; sum(wrist.theta)]);
        end
        se = diff.^2;
        mse = mean(se,2);
        rmse_total(:,i,m) = sqrt(mse);
    end
    
    % Plot The Experiments along with the desired parameters
    points = 100;
    theta_mat_force = zeros(wrist.n,points);
    F_vec = linspace(0,5.5,points); % Getting list of forces (x axis values)
    for k = 1:points
        wrist.fwkin([F_vec(k),0,0],'Type','force');
        theta_mat_force(:,k) = wrist.theta;
    end
    
    plot_options = {'b.','rx','k*','mo'};
    figure('WindowState','Maximize');
    for p = 1:wrist.n+1
        subplot(2,3,p);
        if p < wrist.n+1
            title(sprintf("Notch %d Experimental Results Tip First",i),'FontSize',16);
            hold on
            legend_entries = cell(1,numFiles+1);
            for i = 1:numFiles
                scatter(force_mat(:,i),notch_mat(:,p,i),plot_options{i});
                legend_entries(i) = {sprintf("Experiment %d",i)};
            end
            plot(F_vec,rad2deg(theta_mat_force(p,:)),'g','Linewidth',2);
            legend_entries(end) = {"Model"};
            legend(legend_entries,'Location','southeast','FontSize',14)
            xlabel("Force (N)",'FontSize',14);
            ylabel("Notch Deflection (deg)",'FontSize',14)
            set(gca,'FontSize',12)
            grid on
        else
            title(sprintf("Total Deflection Experimental Results"),'FontSize',18);
            hold on
            legend_entries = cell(1,numFiles+1);
            for i = 1:numFiles
                scatter(force_mat(:,i),notch_mat(:,p,i),plot_options{i});
                legend_entries(i) = {sprintf("Experiment %d",i)};
            end
            plot(F_vec,rad2deg(sum(theta_mat_force(:,:))),'g','Linewidth',2);
            legend_entries(end) = {"Model"};
            legend(legend_entries,'Location','southeast','FontSize',14)
            xlabel("Force (N)",'FontSize',16);
            ylabel("Tip Deflection (deg)",'FontSize',16)
            set(gca,'FontSize',14)
            grid on
        end
        hold off
    end
    
    destdirectory = sprintf("ComparisonImages/PropertySets/PropertySet%d/",table2array(parameters(m,'ID')));
    mkdir(destdirectory);
    saveas(gcf,sprintf("ComparisonImages/PropertySets/PropertySet%d/%s.png",table2array(parameters(m,'ID')),wristType));
    saveas(gcf,sprintf("ComparisonImages/PropertySets/PropertySet%d/%sfig",table2array(parameters(m,'ID')),wristType));
end


end

