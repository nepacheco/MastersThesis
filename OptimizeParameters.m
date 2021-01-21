function parameters = OptimizeParameters(wristType,experimentFiles,usePrecurve,numSets,paramRange)
%OPTIMIZEPARAMETERS - optimizes the material property of a tube
%   wristType is the type of wrist you want to use (90Tube, 150Tube,
%   TipFirstTube)
%   experimentFiles is the list of experiment files to use for optimization
%   usePrecurve is a boolean that determines whether we want to add
%   precurvature to the tube
%   numSets it the number of sets we want to save

arguments
    wristType char {mustBeMember(wristType,{'90Tube','150Tube','TipFirstTube'})}
    experimentFiles (:,1) string
    usePrecurve (1,1) logical
    numSets double
    paramRange.E_lin (1,:) double = 15E9;
    paramRange.scale (1,:) double = 0.2;
    paramRange.strain_lower (1,:) double = 0.02;
    paramRange.mu (1,:) double = 0.4;
    
end
% construct wrist
wrist = MakeWrist(wristType,usePrecurve);

%% read experiment files
numFiles = size(experimentFiles,1);
force_readings = [];
notch_mat = [];
for i = 1:numFiles
    % Parsing File
    opts = detectImportOptions(experimentFiles(i));
    opts.Sheet = 'AvgMeasurements';
    file = readcell(experimentFiles(i),opts);
    [force_vec notch_data] = ParseExperimentFile(file,wrist.n);
    
    force_readings = [force_readings; force_vec];
    notch_mat = [notch_mat; notch_data];
end

%% Optimize Parameters
tic
min_norm_rmse = 100*ones(numSets,5);
for E_lin = paramRange.E_lin
    for scale = paramRange.scale
        for strain_lower = paramRange.strain_lower
            for mu = paramRange.mu
                wrist.E_lin = E_lin;
                wrist.E_se = scale*E_lin;
                wrist.strain_lower = strain_lower;
                wrist.mu = mu;
                diff = zeros(wrist.n+1,length(force_readings));
                for i = 1:length(force_readings)
                    input = force_readings(i);
                    wrist.fwkin([input,0,0],'Type','force');
                    diff(:,i) = notch_mat(i,:)' -...
                        rad2deg([wrist.theta; sum(wrist.theta)]);
                end
                se = diff.^2;
                mse = mean(se,2);
                rmse = sqrt(mse);
                % Replace the largest value in the min_norm_rmse vector
                if (norm(rmse) < max(min_norm_rmse(:,1)))
                    index = find(min_norm_rmse(:,1) == max(min_norm_rmse(:,1)));
                    min_norm_rmse(index(1),1) = norm(rmse);
                    min_norm_rmse(index(1),2) = E_lin/(1E9);
                    min_norm_rmse(index(1),3) = scale;
                    min_norm_rmse(index(1),4) = strain_lower;
                    min_norm_rmse(index(1),5) = mu;
                end
            end
        end
    end
end
fprintf("Optimization duration: %f seconds \n",toc);

%% Add the new unique parameters to the Property Sets table
tic
load('PropertySets.mat','PropertySets');
parameter_time = string(datetime(now,'ConvertFrom','datenum'));
parameters = table();
init_row = size(PropertySets,1);
for m = 1:numSets
    set_num = PropertySets.ID(end) + 1;
    E_lin = min_norm_rmse(m,2)*1E9;
    E_se = min_norm_rmse(m,3)*E_lin;
    strain_lower = min_norm_rmse(m,4);
    mu = min_norm_rmse(m,5);
    
    new_row = {set_num,E_lin,E_se,strain_lower,mu,wristType,...
        usePrecurve,experimentFiles',parameter_time};
    % Don't add if duplicate
    
    PropertySets = [PropertySets; new_row];
    [~,ia] = unique(PropertySets(:,{'E_lin','E_se','Strain_Lower','Mu','Tube','Precurvature'}),'stable');
    PropertySets = PropertySets(ia,:);
end
% Save the property set and add it to the workspace 
save('PropertySets.mat','PropertySets');
assignin ('base', 'PropertySets', PropertySets);
fprintf("Table Update duration: %f seconds \n",toc);

% output
if (size(PropertySets,1) > init_row) % we added to the table
    parameters = PropertySets(init_row+1:end,:); % return the new parameters we added
end
end