%% Plot nitinol stress strain curve
clc; clear; close all;
% Load experiment files and property set data;
load('ExperimentFiles.mat');
load('PropertySets.mat');
% Choose Property Set
parameters = PropertySets(196,:);

wrist90 = MakeWrist('90Tube');
wrist90.E_lin = parameters.E_lin;
wrist90.E_se = parameters.E_se;
wrist90.strain_lower = parameters.Strain_Lower;
wrist90.mu = parameters.Mu;

points = 1000;
strain = linspace(0,0.12,points);
stress = zeros(points,1);
for i = 1:points
stress(i,1) = wrist90.get_stress(strain(i));
end
    
plot(strain,stress/1e6,'r','Linewidth',2)
% title('Nitinol Stress Strain Curve','FontSize',20,'FontName','CMU Serif')
xlabel('Strain ','FontSize',16,'FontName','CMU Serif')
ylabel('Stress [MPa]','FontSize',16,'FontName','CMU Serif')
ax = gca;
ax.FontSize = 14;
ax.FontName = 'CMU Serif';
ax.LineWidth = 1.5;
grid on;
