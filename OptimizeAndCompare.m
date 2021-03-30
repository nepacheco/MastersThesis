clc;clear;close all;

%%
load('ExperimentFiles.mat');
% Getting file data and creating wrists
experimentData90 = ParseExperimentFiles(experimentFiles90(2),5);
wrist90 = MakeWrist('90Tube',true);
wrist90.use_non_linear = false;

experimentData150 = ParseExperimentFiles(experimentFiles150(6),5);
wrist150 = MakeWrist('150Tube',true);
wrist150.use_non_linear = false;

experimentDataTipFirst = ParseExperimentFiles(experimentFilesTip(3),4);
wristTipFirst = MakeWrist('TipFirstTube',true);
wristTipFirst.use_non_linear = false;
E_linRange = linspace(6E9,7E9,50);
% % E_seRange = linspace(2.00E9,4.25E9,10);
% % strain_lowerRange = linspace(0.005,0.03,10);
muRange = linspace(0.05,0.35,10);
numSets = 70;

OptimizeParameters(experimentDataTipFirst,experimentData90,experimentData150,...
    wristTipFirst,wrist90,wrist150,...
    false,numSets,...
    'E_lin',E_linRange,'mu',muRange);
%%
load('PropertySets.mat');
valsTip = CompareModel('TipFirstTube',experimentFilesTip(3),true, PropertySets,'Force',3,'Plot',false);
vals150 = CompareModel('150Tube2',experimentFiles150_2(1),true, PropertySets,'Force',3,'Plot',false);
vals90 = CompareModel('90Tube',experimentFiles90(2),true,PropertySets,'Force',3,'Plot',false);

%%
min_norm = [100,0];
min_norm_noTip = [100, 0];
for i = 1:size(valsTip,3)
    if (norm(vals150(:,1,i))) < min_norm(1)
        min_norm(1) = norm(vals150(:,1,i));
        min_norm(2) = i;
    end
%     if(norm(valsTip(1:end-1,1,i)) + norm(vals150(1:end-1,1,i))+ norm(vals90(1:end-1,1,i))) < min_norm_noTip(1)
%         min_norm_noTip(1) = norm(valsTip(1:end-1,1,i)) + norm(vals150(1:end-1,1,i)) + norm(vals90(1:end-1,1,i));
%         min_norm_noTip(2) = i;
%     end
end
%%
min_norm = [100,0];
min_norm_noTip = [100, 0];
for i = 1:size(valsTip,3)
    if (norm(valsTip(:,1,i)) + norm(vals90(:,1,i))) < min_norm(1)
        min_norm(1) = norm(valsTip(:,1,i))  + norm(vals90(:,1,i));
        min_norm(2) = i;
    end
    if(norm(valsTip(1:end-1,1,i)) + norm(vals90(1:end-1,1,i))) < min_norm_noTip(1)
        min_norm_noTip(1) = norm(valsTip(1:end-1,1,i))  + norm(vals90(1:end-1,1,i));
        min_norm_noTip(2) = i;
    end
end
%%
load('PropertySets.mat')
load('ExperimentFiles.mat')
vals150_1 = CompareModel('150Tube',experimentFiles150(end),true,PropertySets([116,187],:));

