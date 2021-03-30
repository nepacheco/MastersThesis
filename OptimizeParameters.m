function [parameters,min_norm_rmse] = OptimizeParameters(experimentDataTip,experimentData90,experimentData150,wristTipFirst,wrist90,wrist150,useTipDeflection,numSets,paramRange)
%OPTIMIZEPARAMETERS - optimizes the material property of a tube
%   wristType is the type of wrist you want to use (90Tube, 150Tube,
%   TipFirstTube)
%   experimentFiles is the list of experiment files to use for optimization
%   usePrecurve is a boolean that determines whether we want to add
%   precurvature to the tube
%   numSets it the number of sets we want to save

arguments
    experimentDataTip (1,:) cell
    experimentData90 (1,:) cell
    experimentData150 (1,:) cell
    wristTipFirst Wrist
    wrist90 Wrist
    wrist150 Wrist
    useTipDeflection (1,1) logical
    numSets double = 1
    paramRange.E_lin (1,:) double = 15E9;
    paramRange.E_se (1,:) double = 3E9; % typical range is 3E9-5E9
    paramRange.strain_lower (1,:) double = 0.02;
    paramRange.mu (1,:) double = 0.4;
    
end
% construct wrist


%% Optimize Parameters
tic
min_norm_rmse = 100*ones(numSets,5);
for E_lin = paramRange.E_lin
    for E_se = paramRange.E_se
        for strain_lower = paramRange.strain_lower
            for mu = paramRange.mu
                set = [E_lin,E_se,strain_lower,mu];
                output = ObjectiveFunction(set,experimentDataTip,experimentData90,experimentData150, wristTipFirst,wrist90,wrist150,useTipDeflection);
                if (output) < max(min_norm_rmse(:,1))
                    index = find(min_norm_rmse(:,1) == max(min_norm_rmse(:,1)));
                    min_norm_rmse(index(1),1) = output;
                    min_norm_rmse(index(1),2) = E_lin/(1E9);
                    min_norm_rmse(index(1),3) = E_se/(1E9);
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
    try
        set_num = PropertySets.ID(end) + 1;
    catch e
        set_num = 1;
    end
    E_lin = min_norm_rmse(m,2)*1E9;
    E_se = min_norm_rmse(m,3)*1E9;
    strain_lower = min_norm_rmse(m,4);
    mu = min_norm_rmse(m,5);
    wristType = 'all';
    new_row = {set_num,E_lin,E_se,strain_lower,mu,wristType,...
        usePrecurve,experimentFiles',parameter_time};
    % Don't add if duplicate
    
    PropertySets = [PropertySets; new_row];
    [~,ia] = unique(PropertySets(:,{'E_lin','E_se','Strain_Lower','Mu','Precurvature'}),'stable');
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