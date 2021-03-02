clc;clear;close all;

%%
load('ExperimentFiles.mat')

E_linRange = linspace(10E9,40E9,10);
E_seRange = linspace(2.00E9,4.25E9,10);
strain_lowerRange = linspace(0.005,0.03,10);
muRange = linspace(0.1,0.35,10);
numSets = 70;

% Optimizing using Tip First Tube
% OptimizeParameters('TipFirstTube',experimentFilesTip(3),true,false,numSets,...
%     'E_lin',E_linRange,'E_se',E_seRange,...
%     'strain_lower',strain_lowerRange,'mu',muRange);
% 
% % Optimizing using 150 Tube 2
% OptimizeParameters('150Tube2',experimentFiles150_2(1),true,false,numSets,...
%     'E_lin',E_linRange,'E_se',E_seRange,...
%     'strain_lower',strain_lowerRange,'mu',muRange);

% Optimizing using 150 Tube 2
OptimizeParameters('90Tube',experimentFiles90(2),true,false,numSets,...
    'E_lin',E_linRange,'E_se',E_seRange,...
    'strain_lower',strain_lowerRange,'mu',muRange);
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

