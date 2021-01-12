function CreateExcelFile(path,varargin)
%CREATEEXCELFILE - creates the excel file for storing the results of the
%force trials.
%   It will go through and average all the force measurements and add them
%   as well as the image numbers to the excel sheet

%****** INPUT PARSING *********************
% default values
isRelative = false;
saveLocation = "TestResults.xlsx";
numNotches = 5;
imgNum = 0;

p = inputParser();
addRequired(p,'path',@isstring);
addOptional(p, 'isRelative', isRelative, @islogical);
addParameter(p,'SaveLocation',saveLocation,@isstring);
addParameter(p,'NumNotches',numNotches,@isnumeric);
addParameter(p,'ImageNumber',imgNum,@isnumeric);

parse(p,path,varargin{:});

isRelative = p.Results.isRelative;
saveLocation = p.Results.SaveLocation;
numNotches = p.Results.NumNotches;
imgNum = p.Results.ImageNumber;
%*********************************************

if isRelative
    path = pwd + "\" + path;
    saveLocation = pwd + "\" + saveLocation;
end

% Analyzing multiple files in the directory
filesAndFolders = dir(path);
foldersInDir = filesAndFolders(([filesAndFolders.isdir]));
numOfFolders = length(foldersInDir);


writecell({},saveLocation);
row_header = {'','Picture','Tendon Displacement','Force'};
for i = 1:numNotches+1
    row_header(end+1) = {sprintf('Notch %d',i)};
    if i > numNotches
        row_header(end+1) = {'Tip Displacement'};
    end
end

for i = 1:numOfFolders
    % For loop through the files in the directory and analyze each file
    folder_path = path+foldersInDir(i).name+'\';
        
    row_header(1,1) = {foldersInDir(i).name};
    writecell(row_header,saveLocation,'WriteMode','append');
    

    
    F_vec = GetForceReadings(folder_path,false);
    data = cell(length(F_vec),5+numNotches);
    for k = 1:length(F_vec)
        data(k,1:4) = {'','','',F_vec(k)};
        data(k,2) = {sprintf('DSC_%d',imgNum)};
        imgNum = imgNum + 1;
        data(k,5:end) = {''};
    end
    
    writecell(data,saveLocation,'WriteMode','append');
end

end

