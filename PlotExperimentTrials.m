function plot_ = PlotExperimentTrials(experimentFiles,n)

experimentStr = "";
line_type = {'-','--',':','-.'};
figure('Name','Experiment Trials Plotted','WindowState','Maximize');
for i = 1:length(experimentFiles)
    force_vec = [];
    notch_data = [];
    backslash_indicies = strfind(experimentFiles(i),"\");
    period_indicies = strfind(experimentFiles(i),".");
    experimentStr = experimentStr + extractBetween(experimentFiles(i),backslash_indicies(end)+1,period_indicies(end)-1);
    
    opts = detectImportOptions(experimentFiles(i));
    opts.Sheet = 'AvgMeasurements';
    file = readcell(experimentFiles(i),opts);
    [force_vec_i notch_data_i] = ParseExperimentFile(file,n);
    force_vec = [force_vec; force_vec_i];
    notch_data = [notch_data; notch_data_i];
    
    
    trialForceList = cell(1,1);
    trialNotchDataList = cell(1,1);
    last = 1;
    for f = 1:length(force_vec)-1
        if force_vec(f) > force_vec(f+1)
            if last == 1
                trialForceList = {force_vec(last:f)};
                trialNotchDataList = {notch_data(last:f,:)};
                last = f+1;
            else
                trialForceList = [trialForceList, {force_vec(last:f)}];
                trialNotchDataList = [trialNotchDataList {notch_data(last:f,:)}];
                last = f+1;
            end
        elseif f == length(force_vec) - 1
            trialForceList = [trialForceList, {force_vec(last:f+1)}];
            trialNotchDataList = [trialNotchDataList {notch_data(last:f+1,:)}];
        end
    end
    
    
    hold on
    for p = 1:n+1
        subplot(3,2,p)
        hold on
        title_str = sprintf("Notch %d Deflection",p);
        if p == n+1
            title_str = sprintf("Total tip Deflection");
        end
        title(title_str,'FontSize',16);
        for t = 1:size(trialForceList,2)
            disp_name = sprintf('Exp%d Trial %d',i,t);
            plot(trialForceList{1,t},trialNotchDataList{1,t}(:,p),...
                line_type{i},'LineWidth',2,'Marker','.',...
                'DisplayName',disp_name);
        end
        xlabel('Force (N)','FontSize',12);
        ylabel('Deflection (deg)','FontSize',12);
        legend;
        hold off
    end
end
SaveDestination = "ComparisonImages/TrialImages";
destdirectory = sprintf("%s/",SaveDestination);
if ~exist(destdirectory, 'dir')
    mkdir(destdirectory);
end
saveas(gcf,sprintf("%s/%s_Trials.png",SaveDestination,experimentStr));
saveas(gcf,sprintf("%s/%s_Trials.fig",SaveDestination,experimentStr));

end