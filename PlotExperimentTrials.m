function p = PlotExperimentTrials(experimentFile,n)
opts = detectImportOptions(experimentFile);
opts.Sheet = 'AvgMeasurements';
file = readcell(experimentFile,opts);
[force_vec notch_data] = ParseExperimentFile(file,n);

trialForceList = cell(1,1);
trialNotchDataList = cell(1,1);
last = 1;
for i = 1:length(force_vec)-1
    if force_vec(i) > force_vec(i+1)
        if last == 1
            trialForceList = {force_vec(last:i)};
            trialNotchDataList = {notch_data(last:i,:)};
            last = i+1;
        else
            trialForceList = [trialForceList, {force_vec(last:i)}];
            trialNotchDataList = [trialNotchDataList {notch_data(last:i,:)}];
            last = i+1;
        end
    elseif i == length(force_vec) - 1
            trialForceList = [trialForceList, {force_vec(last:i+1)}];
            trialNotchDataList = [trialNotchDataList {notch_data(last:i+1,:)}];
    end
end



figure('Name','Experiment Trials Plotted')
hold on
for i = 1:n+1
    subplot(3,2,i)
    hold on
    title_str = sprintf("Notch %d Deflection",i);
    if i == n+1
        title_str = sprintf("Total tip Deflection");
    end
    title(title_str,'FontSize',16);
    legend_entries = cell(1,size(trialForceList,2));
    for t = 1:size(trialForceList,2)
        plot(trialForceList{1,t},trialNotchDataList{1,t}(:,i),...
            '-','LineWidth',2,'Marker','.');
        legend_entries(1,t) = {sprintf('Trial %d',t)};
    end
    xlabel('Force (N)','FontSize',12);
    ylabel('Deflection (deg)','FontSize',12);
    legend(legend_entries,'Location','southeast','FontSize',12);
    hold off
end

SaveDestination = "ComparisonImages/TrialImages"
destdirectory = sprintf("%s/",SaveDestination);
if ~exist(destdirectory, 'dir')
    mkdir(destdirectory);
end
saveas(gcf,sprintf("%s/",SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));
saveas(gcf,sprintf("%s/",SaveDestination,SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));

writematrix(rmse_total(:,:,m),...
    sprintf("%s/PropertySet%d/RMSE_Values_%s_%s.xlsx",SaveDestination,table2array(parameters(m,'ID')),wristType,experimentStr));
end