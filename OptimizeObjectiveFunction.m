%% Using MATLAB To optimize parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb = [5E9,1.75E9,0.01,0.05];
ub = [30E9,4E9,0.04,0.35];
x0 = [10E9,2.25E9,0.0272,0.2389];

problem = createOptimProblem('fmincon',...
    'objective',@(x)ObjectiveFunction(x),...
    'x0',x0,'lb', lb, 'ub', ub,'options',...
    optimoptions(@fmincon,'Algorithm','sqp','Display','off'));
tic
gs = GlobalSearch('Display','iter','NumTrialPoints',3000,'MaxTime',36000);
rng(14,'twister') % for reproducibility
[x,fval] = run(gs,problem);
toc

load('PropertySets.mat')
load('ExperimentFiles.mat')
set_num = PropertySets.ID(end) + 1;
E_lin = x(1);
E_se = x(2);
strain_lower = x(3);
mu = x(4);
parameter_time = string(datetime(now,'ConvertFrom','datenum'));
Precurve = true;
expFiles = "all";
Tube = "multiple";

new_row = {set_num,E_lin,E_se,strain_lower,mu,Tube,...
        Precurve,expFiles',parameter_time};
% Don't add if duplicate
PropertySets = [PropertySets; new_row];
valsTip = CompareModel('TipFirstTube',experimentFilesTip(3),true, PropertySets(end,:),'Force',3);
vals150 = CompareModel('150Tube',experimentFiles150(6),true, PropertySets(end,:),'Force',3);
vals90 = CompareModel('90Tube',experimentFiles90(2),true,PropertySets(end,:),'Force',3);
% save('PropertySets.mat','PropertySets')

%% Getting Objective Function data

numPoints = 50;
E_lin = 10E9;
E_se = linspace(1.0E9,3.5E9,numPoints);
strain_lower = linspace(0.02,0.03,numPoints)';
mu = linspace(0.1,0.3,numPoints)';
tic
function_mat = zeros(numPoints,numPoints);
parfor i = 1:numPoints
    for j = 1:numPoints
        input = [E_lin,E_se(i),strain_lower,mu(j)];
        function_mat(i,j) = ObjectiveFunction(input);
    end
end
toc

%% Graphing Objective funciton with 3D Plot
[E,m] = meshgrid(E_se,mu);
surf(E,m,function_mat)
xlabel('E_se')
ylabel('mu')
zlabel('RMSE')  

%% Graphing Objective function with a 4D plot
cla                           
numPoints = 10;
E_lin = 10E9;
E_se = linspace(1.0E9,3.5E9,numPoints);
strain_lower = linspace(0.02,0.03,numPoints)';
mu = linspace(0.1,0.3,numPoints)';
tic
function_mat = zeros(numPoints,numPoints,numPoints);
parfor i = 1:numPoints
    for j = 1:numPoints
        for k = 1:numPoints
            input = [E_lin,E_se(i),strain_lower(j),mu(k)];
            function_mat(i,j,k) = ObjectiveFunction(input);
        end
    end
end
[E,strain,m] = meshgrid(E_se,strain_lower,mu);
toc
%%
E = reshape(E,[numPoints^3,1]);
strain = reshape(strain,[numPoints^3,1]);
m = reshape(m,[numPoints^3,1]);
function_mat = reshape(function_mat,[numPoints^3,1]);
scatter3(E,strain,m,40,function_mat,'filled')    % draw the scatter plot
xlabel('E_se')
ylabel('Strain')
zlabel('mu')

cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'RMSE';