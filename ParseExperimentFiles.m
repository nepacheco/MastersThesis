function experimentData = ParseExperimentFiles(experimentFiles,n)
%ParseExperimentFile - parses the NxM cell passed into a force vector and notch value
%matrix
%   force vector output is Nx1 and notch_mat is Nxn where n is the number
%   of notches in the tube.

% Parse Experiment Files
numFiles = size(experimentFiles,1);
force_cell = cell(1,numFiles);
notch_cell = cell(1,numFiles);
experimentStr = "";
average = zeros(n+1,numFiles);
for i = 1:numFiles
    % Parsing File
    backslash_indicies = strfind(experimentFiles(i),"\");
    period_indicies = strfind(experimentFiles(i),".");
    experimentStr = experimentStr + extractBetween(experimentFiles(i),backslash_indicies(end)+1,period_indicies(end)-1);
    opts = detectImportOptions(experimentFiles(i));
    opts.Sheet = 'AvgMeasurements';
    file = readcell(experimentFiles(i),opts);
    [force_vec, notch_data] = ParseExperimentFile(file,n);
    
    force_cell(1,i) = {force_vec};
    notch_cell(1,i) = {notch_data};
    average(:,i) = mean(notch_data)';
end

experimentData = [force_cell, notch_cell,average,convertStringsToChars(experimentStr),numFiles];

end

