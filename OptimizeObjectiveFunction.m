%% Using MATLAB To optimize parameters to reduce angle error
load('ExperimentFiles.mat');
tic
% Getting file data and creating wrists
experimentData90 = ParseExperimentFiles(experimentFiles90(2),5);
wrist90 = MakeWrist('90Tube',true);
wrist90.use_non_linear = true;
wrist90.use_friction = true;

experimentData150 = ParseExperimentFiles(experimentFiles150(6),5);
wrist150 = MakeWrist('150Tube',true);
wrist150.use_non_linear = true;
wrist150.use_friction = true;

experimentDataTipFirst = ParseExperimentFiles(experimentFilesTip(3),4);
wristTipFirst = MakeWrist('TipFirstTube',true);
wristTipFirst.use_non_linear = true;
wristTipFirst.use_friction = true;

%%
save
A = [];
b = [];
Aeq = [];
beq = [];
% lb = [5E9,1.75E9,0.01,0.05];
lb = [7E9,2.25E9,0.01 0.05];
% ub = [10E9,4E9,0.04,0.35];
ub = [60E9,10E9,0.05,0.5];
x0 = [7E9,2.25E9,0.0272,0.2389];

problem = createOptimProblem('fmincon',...
    'objective',@(x)ObjectiveFunction(x,experimentDataTipFirst,experimentData90,experimentData150,...
    wristTipFirst,wrist90,wrist150),...
    'x0',x0,'lb', lb, 'ub', ub,'options',...
    optimoptions(@fmincon,'Algorithm','sqp','Display','off'));
tic
gs = GlobalSearch('Display','iter','NumTrialPoints',1000,'MaxTime',360);
rng(14,'twister') % for reproducibility
[x,fval] = run(gs,problem);
toc
save
load('PropertySets.mat')
load('ExperimentFiles.mat')

set_num = PropertySets.ID(end) + 1;
E_lin = x(1);
E_se = x(2);
strain_lower = x(3);
mu = x(4);
parameter_time = string(datetime(now,'ConvertFrom','datenum'));
Precurve = true;
expFiles = "Just90";
Tube = "Just 90 with relaxed assumptions";

new_row = {set_num,E_lin,E_se,strain_lower,mu,Tube,...
        Precurve,expFiles',parameter_time};
% Don't add if duplicate
PropertySets = [PropertySets; new_row];

valsTip = CompareModel(wristTipFirst,experimentDataTipFirst,PropertySets(end,:));
vals150 = CompareModel(wrist150,experimentData150, PropertySets(end,:));
vals90 = CompareModel(wrist90,experimentData90, PropertySets(end,:));
save('PropertySets.mat','PropertySets')
%% Optimizing for Tip Error
clc; clear; close all;
load('ExperimentFiles.mat');
load('PropertySets.mat')
tic
% Getting file data and creating wrists
experimentData90 = ParseExperimentFiles(experimentFiles90(2),5);
wrist90 = MakeWrist('90Tube',true);

experimentData150 = ParseExperimentFiles(experimentFiles150(6),5);
wrist150 = MakeWrist('150Tube',true);

experimentDataTipFirst = ParseExperimentFiles(experimentFilesTip(3),4);
wristTipFirst = MakeWrist('TipFirstTube',true);
toc

CompareTipModel(wrist150,experimentData150,PropertySets(196,:))
CompareTipModel(wrist90,experimentData90,PropertySets(196,:))
CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(196,:))

%% 
A = [];
b = [];
Aeq = [];
beq = [];
lb = [5E9,1.75E9,0.01,0.05];
ub = [20E9,5E9,0.04,0.35];
x0 = [10E9,3E9,0.0272,0.1275];
tic
problem = createOptimProblem('fmincon',...
    'objective',@(x)ObjectiveFunction2(x,experimentDataTipFirst,experimentData90,experimentData150,...
    wristTipFirst,wrist90,wrist150),...
    'x0',x0,'lb', lb, 'ub', ub,'options',...
    optimoptions(@fmincon,'Algorithm','sqp','Display','off'));
toc
tic
gs = GlobalSearch('Display','iter','NumTrialPoints',5000,'MaxTime',36000);
rng(14,'twister') % for reproducibility
[x,fval] = run(gs,problem);
toc
save

load('PropertySets.mat')
load('ExperimentFiles.mat')

set_num = PropertySets.ID(end) + 1;
E_lin = x(1);
E_se = x(2);
strain_lower = x(3);
mu = x(4);
parameter_time = string(datetime(now,'ConvertFrom','datenum'));
Precurve = true;
expFiles = "TipRMSE";
Tube = "multiple";

new_row = {set_num,E_lin,E_se,strain_lower,mu,Tube,...
        Precurve,expFiles',parameter_time};
% Don't add if duplicate
PropertySets = [PropertySets; new_row];

valsTip = CompareModel(wristTipFirst,experimentDataTipFirst,PropertySets(end,:));
vals150 = CompareModel(wrist150,experimentData150, PropertySets(end,:));
vals90 = CompareModel(wrist90,experimentData90, PropertySets(end,:));
save('PropertySets.mat','PropertySets')
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
numPoints = 12;
E_lin = linspace(5E9,60E9,numPoints);
E_se = linspace(1.5E9,18E9,numPoints);
% E_se = 3E9;
strain_lower = linspace(0.010,0.032,numPoints)';
% mu = linspace(0.10,0.21,numPoints)';
function_mats = cell(4,1);
for index = 1:4
mu = 0.13*index/2 + 0.13*1/2;
tic
function_mat = zeros(numPoints,numPoints,numPoints);
E = zeros(numPoints,numPoints,numPoints);
strain = E;
E_s = E;
parfor i = 1:numPoints
    for j = 1:numPoints
        for k = 1:numPoints
            input = [E_lin(i),E_se(k),strain_lower(j),mu];
            E(i,j,k) = E_lin(i);
            strain(i,j,k) = strain_lower(j);
            E_s(i,j,k) = E_se(k);
            function_mat(i,j,k) = ObjectiveFunction(input,experimentDataTipFirst,experimentData90,experimentData150,...
    wristTipFirst,wrist90,wrist150);
        end
    end
end
function_mats{index,1} = function_mat;
fprintf("Done running points\n");
toc

E_reshaped = reshape(E,[numPoints^3,1]);
strain_reshaped = reshape(strain,[numPoints^3,1]);
E_s_reshaped = reshape(E_s,[numPoints^3,1]);
function_mat_reshaped = reshape(function_mat,[numPoints^3,1]);
subplot(2,2,index);
scatter3(E_reshaped,strain_reshaped,E_s_reshaped,40,function_mat_reshaped,'filled')    % draw the scatter plot
xlabel('E_{lin}','FontSize',14,'FontName','CMU Serif')
ylabel('Strain','FontSize',14,'FontName','CMU Serif')
zlabel('E_{se}','FontSize',14,'FontName','CMU Serif')



cb = colorbar;                                     % create and label the colorbar
cb.Label.String = 'RMSE';
cb.Label.FontName = 'CMU Serif';
end

%% Getting Objective Function data using 2D plot

numPoints = 5000;
E_lin = linspace(1E9,60E9,numPoints);
E_se = 3E9;
strain_lower = 0.028;
mu = 0.13;
tic
function_mat = zeros(numPoints,1);
parfor i = 1:numPoints
        input = [E_lin(i),E_se,strain_lower,mu];
        function_mat(i) = ObjectiveFunction(input,experimentDataTipFirst,experimentData90,experimentData150,...
    wristTipFirst,wrist90,wrist150);
end
toc
%%
plot(E_lin./1E9,function_mat,'LineWidth',2);
xlabel('E_{lin} (GPa)','FontSize',14,'FontName','CMU Serif');
ylabel('RMSE (degrees)','FontSize',14,'FontName','CMU Serif'); 
title('RMSE of Model based on Linear Elastic Modulus','FontSize',18,'FontName','CMU Serif')