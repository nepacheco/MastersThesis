clc;clear;close all;
c = distinguishable_colors(50);
load('ExperimentFiles.mat');
load('PropertySets.mat');
tic
% Getting file data and creating wrists
experimentData90 = ParseExperimentFiles(experimentFiles90(2),5);
wrist90 = MakeWrist('90Tube',true);

experimentData150 = ParseExperimentFiles(experimentFiles150(6),5);
wrist150 = MakeWrist('150Tube',true);

experimentDataTipFirst = ParseExperimentFiles(experimentFilesTip(3),4);
wristTipFirst = MakeWrist('TipFirstTube',true);


% Getting Data for FRICTION and NONLINEAR
nonlin150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(202,:),'Plot',false);
nonlin90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(202,:),'Plot',false);
nonlinTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(202,:),'Plot',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Data for FRICTION and LINEAR
wrist90.use_non_linear = false;
wrist150.use_non_linear = false;
wristTipFirst.use_non_linear = false;


lin150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(198,:),'Plot',false);
lin90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(198,:),'Plot',false);
linTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(198,:),'Plot',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting Data for NO FRICTION and LINEAR
wrist90.use_non_linear = false;
wrist150.use_non_linear = false;
wristTipFirst.use_non_linear = false;
wrist90.use_friction = false;
wrist150.use_friction = false;
wristTipFirst.use_friction = false;


nof150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(198,:),'Plot',false);
nof90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(198,:),'Plot',false);
nofTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(198,:),'Plot',false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Getting data for NO FRICTION and NONLINEAR
wrist90.use_non_linear = true;
wrist150.use_non_linear = true;
wristTipFirst.use_non_linear = true;
wrist90.use_friction = false;
wrist150.use_friction = false;
wristTipFirst.use_friction = false;

nofnonlin150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(201,:),'Plot',false);
nofnonlin90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(201,:),'Plot',false);
nofnonlinTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(201,:),'Plot',false);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting
x = categorical({'Wrist A' ,'Wrist B', 'Wrist C'});
y = 1000*[nonlin150_rmse nof150_rmse lin150_rmse nofnonlin150_rmse ;nonlin90_rmse  nof90_rmse lin90_rmse nofnonlin90_rmse; nonlinTip_rmse  nofTip_rmse linTip_rmse nofnonlinTip_rmse];
b = bar(x,y);


b(1).FaceColor = c(1,:);
b(2).FaceColor = c(3,:);
b(3).FaceColor = c(2,:);
b(4).FaceColor = c(4,:);


set(gca,'FontSize',12,'FontName','CMU Serif')
ylabel('RMSE (mm)','FontName','CMU Serif','FontSize',14)
title('Tip Tracking Accuracy','FontName','CMU Serif','FontSize',16)
legend('Nonlinear with Friction', 'Linear without Friction', 'Linear with Friction', 'Nonlinear without Friction','FontName','CMU Serif','FontSize',14)

%% Calculate rmse of rmse

numInputs90 = size(experimentData90{1,1},1);
numInputs150 = size(experimentData150{1,1},1);
numInputsTip = size(experimentDataTipFirst{1,1},1);


% non linear rmse across tubes
error_vecNonLin = [nonlin150_rmse*ones(numInputs150,1); nonlin90_rmse*ones(numInputs90,1); nonlinTip_rmse*ones(numInputsTip,1)];
rmse_NonLin = sqrt(mean(error_vecNonLin.^2));


% Frictionless linear rmse across tubes
error_veclin = [lin150_rmse*ones(numInputs150,1); lin90_rmse*ones(numInputs90,1); linTip_rmse*ones(numInputsTip,1)];
rmse_lin = sqrt(mean(error_veclin.^2));

% non linear rmse across tubes
error_vecNoF = [nof150_rmse*ones(numInputs150,1); nof90_rmse*ones(numInputs90,1); nofTip_rmse*ones(numInputsTip,1)];
rmse_NoF = sqrt(mean(error_vecNoF.^2));

fprintf("%0.5f, %0.5f, %0.5f\n",rmse_NonLin,rmse_lin,rmse_NoF)