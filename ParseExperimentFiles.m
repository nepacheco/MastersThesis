function experimentData = ParseExperimentFiles(experimentFiles,n,options)
%ParseExperimentFiles - takes in a list of paths to experiment files and
%parses all of them into a single cell array
arguments
    experimentFiles
    n = 5 % number of notches
    options.force_index = 4 % the force index
    options.tendon_index = 3 
    options.notch1_index = 5
    options.Sheet = 'AvgMeasurements'
end

% Parse Experiment Files
numFiles = size(experimentFiles,1);
force_cell = cell(1,numFiles);
notch_cell = cell(1,numFiles);
tendon_cell = cell(1,numFiles);
experimentStr = "";
average = zeros(n+1,numFiles);
for i = 1:numFiles
    % Parsing File
    backslash_indicies = strfind(experimentFiles(i),"\");
    period_indicies = strfind(experimentFiles(i),".");
    experimentStr = experimentStr + extractBetween(experimentFiles(i),backslash_indicies(end)+1,period_indicies(end)-1);
    opts = detectImportOptions(experimentFiles(i));
    opts.Sheet = options.Sheet;
    file = readcell(experimentFiles(i),opts);
    [force_vec, notch_data, tendon_data] = ParseExperimentFile(file,n,...
        'force_index', options.force_index,'tendon_index',...
        options.tendon_index, 'notch1_index', options.notch1_index);
    
    force_cell(1,i) = {force_vec};
    notch_cell(1,i) = {notch_data};
    tendon_cell(1,i) = {tendon_data};
    average(:,i) = mean(notch_data)';
end

experimentData = struct('force_data',{force_cell},...
    'notch_data',{notch_cell},'notch_averages',{average},...
    'tendon_data',{tendon_cell},...
    'experiment_str',convertStringsToChars(experimentStr),'num_files',numFiles);

end

