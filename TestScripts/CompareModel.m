function [rmse_total, r2_total] = CompareModel(wrist,experimentData,parameters,SaveDestination,options)
%COMPAREMODEL Compares the model to a set of experimental data using the
%given parameters
%   Wrist is the wrist class to use
%   experimentFiles is a vector containing the absolute path to each
%   experiment excel file
arguments
    wrist Wrist
    experimentData (1,:) cell
    parameters table
    SaveDestination string = "ComparisonImages/PropertySets"
    options.Force (1,1) double = 3;
    options.Plot logical = true;
    options.numFiles = 1;
end
numFiles = experimentData{1,end};
force_cell = experimentData(1,1:numFiles);
notch_cell = experimentData(1,numFiles+1:numFiles*2);
experimentStr = convertCharsToStrings(experimentData{1,end-1});
expAverage = experimentData{1,end-2};

%% Calculating RMSE and Plotting is done for each set of material properties
rmse_total = zeros(wrist.n+1,numFiles,size(parameters,1));
r2_total = zeros(wrist.n+1,numFiles,size(parameters,1));
for m = 1:size(parameters,1)
    %% Defining the material properties of the wrist for this experiment
    wrist.E_lin = table2array(parameters(m,'E_lin'));
    wrist.E_se = table2array(parameters(m,'E_se'));
    wrist.strain_lower = table2array(parameters(m,'Strain_Lower'));
    wrist.mu = table2array(parameters(m,'Mu'));
    %% Generating RMSE
    for i = 1:numFiles
        diff = zeros(wrist.n+1,length(force_cell{1,i}));
        
        for j = 1:length(force_cell{1,i})
            input = force_cell{1,i}(j);
            wrist.fwkin([input,0,0],'Type','force');
            diff(:,j) = (notch_cell{1,i}(j,:))' - rad2deg([wrist.theta; sum(wrist.theta)]);
        end

        se = diff.^2;
        mse = mean(se,2);
        rmse_total(:,i,m) = sqrt(mse);
        
        ssTot = sum((notch_cell{1,i}' - expAverage(:,i)).^2,2);
        ssRes = sum(se,2);
        r2_total(:,i,m) = ones(wrist.n+1,1) - (ssRes./ssTot);
    end
    
    %% Plot The Experiments along with the desired parameters
    if options.Plot
        points = 100;
        theta_mat_force = zeros(wrist.n,points);
        F_vec = linspace(0,options.Force,points); % Getting list of forces (x axis values)
        for k = 1:points
            wrist.fwkin([F_vec(k),0,0],'Type','force');
            theta_mat_force(:,k) = wrist.theta;
        end

        plot_options = {'rx','b.','k*','mo','rx','rx','rx'};
        figure('WindowState','Maximize');
        for p = 1:wrist.n+1
            subplot(2,3,p);
            if p < wrist.n+1
                title(sprintf("Notch %d Experimental Results",p),'FontSize',16);
                hold on
                legend_entries = cell(1,numFiles+1);
                for i = 1:numFiles
                    scatter(force_cell{1,i},notch_cell{1,i}(:,p),plot_options{i});
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
                    scatter(force_cell{1,i},notch_cell{1,i}(:,p),plot_options{i});
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


        %% Saving the figure
        wristType = wrist.name;
        if ~wrist.use_friction
            experimentStr = experimentStr + "_noFriction";
        end
        if ~wrist.use_non_linear
            experimentStr = experimentStr + "_linear";
        end
        destdirectory = sprintf("%s/PropertySet%d/",SaveDestination,table2array(parameters(m,'ID')));
        if ~exist(destdirectory, 'dir')
            mkdir(destdirectory);
        end
        saveas(gcf,sprintf("%s/PropertySet%d/%s_%s.png",SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));
        saveas(gcf,sprintf("%s/PropertySet%d/%s_%sfig",SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));

        writematrix(rmse_total(:,:,m),...
            sprintf("%s/PropertySet%d/RMSE_Values_%s_%s.xlsx",SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));
        writematrix(r2_total(:,:,m),...
            sprintf("%s/PropertySet%d/R2_Values_%s_%s.xlsx",SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));
        close gcf
    end
end

end

