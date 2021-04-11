clc;clear;close all;

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

nonlin150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(196,:),'Plot',false);
nonlin90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(196,:),'Plot',false);
nonlinTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(196,:),'Plot',false);

wrist90.use_non_linear = false;
wrist150.use_non_linear = false;
wristTipFirst.use_non_linear = false;


lin150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(198,:),'Plot',false);
lin90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(198,:),'Plot',false);
linTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(198,:),'Plot',false);


wrist90.use_friction = false;
wrist150.use_friction = false;
wristTipFirst.use_friction = false;

    
nof150_rmse = CompareTipModel(wrist150,experimentData150,PropertySets(198,:),'Plot',false);
nof90_rmse = CompareTipModel(wrist90,experimentData90,PropertySets(198,:),'Plot',false);
nofTip_rmse = CompareTipModel(wristTipFirst,experimentDataTipFirst,PropertySets(198,:),'Plot',false);


x = categorical({'Wrist A' ,'Wrist B', 'Wrist C'});
y = 1000*[nonlin150_rmse nof150_rmse lin150_rmse  ;nonlin90_rmse  nof90_rmse lin90_rmse ; nonlinTip_rmse  nofTip_rmse linTip_rmse];
bar(x,y)
set(gca,'FontSize',12,'FontName','CMU Serif')
ylabel('RMSE (mm)','FontName','CMU Serif','FontSize',14)
title('Tip Tracking Accuracy','FontName','CMU Serif','FontSize',16)
legend('Non-Linear with Friction', 'Linear without Friction', 'Linear with Friction','FontName','CMU Serif','FontSize',14)

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