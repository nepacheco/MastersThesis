%% Comparison Script
% This script is meant to compare the outputs from force measurements of
% the nitinol notched wrist bending to the output predicted by Josh's model
clc; clear; close all;
%% Tube Creation
% This is where we put all the tube parameters for the tube that bends up
% to 150° with uniform notches
od = 1.62E-3; % [m] - outer diameter of tube
id = 1.4E-3; % [m] - inner diameter of tube
n = 5; % number of notches
phi = zeros(n,1);
% tip first bending tube
% g = [1.36,1.39,1.42,1.45].*1E-3; 
% h = 1.0E-3*ones(n,1);
% c = 1.0E-3*ones(n,1);
% 150 bend tube and 90 bend tube
g = 1.4E-3*ones(n,1); 
h = 0.8E-3*ones(n,1);
c = 1.2E-3*ones(n,1);

E_lin = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
E_se = 0.08*E_lin; % [N/m^2] - Slope of Super Elastic Region for Nitinol
mu = 0.2; % coefficient of friction for capstan
T_0 = eye(4,4);
b = c(1);    


%% Comparing my Wrist Class to Josh's
wrist = Wrist(od,id,n,h,phi,c,g,'CutType','off-axis'); % Nick's wrist class
sheath = JoshWrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h,'h_c',h,...
    'r_o',od/2,'r_i',id/2,'n',n,'material','nitinol','plotkin',false,'verbose',false);
% **** PLOTTING NOTCH ANGLE WRT FORCE **********
theta_last = zeros(n,1);
points = 100; % how many points to plot
diff_values = zeros(n+1,points);
theta_mat_sheath = zeros(n,points); % Initializing empty array to store values
theta_mat_wrist = zeros(n,points);
F_vec = linspace(0,5.5,points); % Getting list of forces (x axis values)
for i = 1:points
    sheath.GetKinematicsForce(F_vec(i)); % Updating tube position
    theta_mat_sheath(:,i) = sheath.theta; % Getting tube position
    wrist.fwkin([F_vec(i),0,0],'Type','force');
    theta_mat_wrist(:,i) = wrist.theta;
    if (norm(wrist.theta-sheath.theta) > 1E-4)
        diff_values(:,i) = [wrist.theta-sheath.theta; F_vec(i)];
    end
end

for i = 1:n
    subplot(2,3,i);
    title(sprintf("Notch %d comparison",i));
    hold on
    plot(F_vec,theta_mat_sheath(i,:));
    plot(F_vec,theta_mat_wrist(i,:));
    legend('Josh Model','Nick Model');
    grid on
    hold off
end

fprintf("Notch difference between josh and my wrist class: %g\n",sqrt(mean((diff_values(1:5,:).^2)')));

points = 10000;
strain = linspace(0,.12,points);
my_func = zeros(points,2);
josh_func = zeros(points,2);
for i = 1:points
    [my_func(i,1),my_func(i,2)] = wrist.get_stress(strain(i));
    [josh_func(i,1), josh_func(i,2)] = sheath.SuperElastic(strain(i));
    if isnan(my_func(i,1)) || isnan(josh_func(i,1))
        my_func(i,1) = my_func(i-1,1);
        josh_func(i,1) = josh_func(i-1,1);
    end
    
end
diff_mat = my_func - josh_func;
fprintf("Difference in stress: %f, and eta: %f\n",mean(diff_mat));

%% Comparing Test Results to Wrist Model
close all
cutType = 'on-axis';
wrist = Wrist(od,id,n,h,phi,c,g,'CutType',cutType); % Nick's wrist class

% tip first bending01-12-20201 Experiment
% precurvature = deg2rad([2.191;2.264;2.534;3.062]); 

% this is the aversage initial reading from the 12-12 experiment 
precurvature = deg2rad([2.39580099;2.268315378;2.433246067;1.724263869;2.334074353]); 

% Precurvature for experiment on 12-19-2020
% precurvature = deg2rad([1.521853776;1.320520452;1.255100512;1.149834263;1.336681514]);

wrist.precurve_theta = precurvature;
wrist.E_lin = 10E9;
wrist.E_se = 0.35*wrist.E_lin;
wrist.strain_lower = 0.03;
wrist.mu = 0.4;
file_path2 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-12-2020_Experiment\12-12-2020_Results.xlsx";
file_path1 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-7-2020_Experiment\12-7-2020_Results.xlsx";
% file_path1 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-19-2020_Experiment\12-19-2020_Results.xlsx";
% file_path1 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-29-2020_Experiment\12-29-2020_Results.xlsx";
% file_path2 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\01-12-2021_Experiment\01-12-2021_Results.xlsx";
opts = detectImportOptions(file_path1);
opts.Sheet = 'AvgMeasurements';

file1 = readcell(file_path1,opts);
file2 = readcell(file_path2,opts);

fprintf("Cut depth :%f, Cut Height: %f, Cut Type: %s\n",g(1),h(1),cutType);
[force_vec1,notch_mat1] = parseFile(file1,n);
[force_vec2,notch_mat2] = parseFile(file2,n);

points = 100; % how many points to plot
diff_values = zeros(n+1,points);
theta_mat_force = zeros(n,points);
F_vec = linspace(0,6,points); % Getting list of forces (x axis values)
for i = 1:points
    wrist.fwkin([F_vec(i),0,0],'Type','force');
    theta_mat_force(:,i) = wrist.theta;
end

diff = zeros(n+1,length(force_vec1));
for i = 1:length(force_vec1)
    input = force_vec1(i);
    wrist.fwkin([input,0,0],'Type','force');
    diff(:,i) = notch_mat1(i,:)' -...
        rad2deg([wrist.theta; sum(wrist.theta)]);
end
se = diff.^2;
mse = mean(se,2);
disp("Experiment 1");
rmse = sqrt(mse)

% Generating RMSE
diff = zeros(n+1,length(force_vec2));
for i = 1:length(force_vec2)
    input = force_vec2(i);
    wrist.fwkin([input,0,0],'Type','force');
    diff(:,i) = notch_mat2(i,:)' -...
        rad2deg([wrist.theta; sum(wrist.theta)]);
end
se = diff.^2;
mse = mean(se,2);
disp("Experiment 2");
rmse = sqrt(mse)

% Plotting Force 
figure(2)
for i = 1:n+1
    subplot(2,3,i);
    if i < n+1
        title(sprintf("Notch %d Experimental Results Tip First",i),'FontSize',16);
        hold on
        scatter(force_vec1,notch_mat1(:,i),'b.');
        scatter(force_vec2,notch_mat2(:,i),'rx');
        plot(F_vec,rad2deg(theta_mat_force(i,:)),'g','Linewidth',2);
        legend("Experiment1","Experiment2","Model",'Location','southeast','FontSize',12)
        xlabel("Force (N)",'FontSize',14);
        ylabel("Notch Deflection (deg)",'FontSize',14)
        set(gca,'FontSize',12)
        grid on
    else
        title(sprintf("Total Deflection Experimental Results"),'FontSize',18);
        hold on
        scatter(force_vec1,notch_mat1(:,i),'b.');
        scatter(force_vec2,notch_mat2(:,i),'rx');
        plot(F_vec,rad2deg(sum(theta_mat_force(:,:))),'g','Linewidth',2);
        legend("Experiment1","Experiment2","Model",'Location','southeast','FontSize',14)
        xlabel("Force (N)",'FontSize',16);
        ylabel("Tip Deflection (deg)",'FontSize',16)
        set(gca,'FontSize',14)
        grid on
    end
    hold off
end

% DATA = [notch_mat1(:,6),force_vec1(:)];%notch_mat2(:,6), force_vec2(:)];
% Wrist_statics
%% Functions
function [force_vec, notch_mat] = parseFile(file,n)
%PARSEFILE - parses the NxM cell passed into a force vector and notch value
%matrix
force_index = 4;
force_vec = [];
notch1_index = 5;
notch_mat = [];
[N,M] = size(file);
for i = 1:N
    if (isnumeric(file{i,force_index}) && isnumeric(file{i,notch1_index}))
        force_vec = [force_vec; file{i,force_index}];
        notch_mat = [notch_mat; file{i,notch1_index:notch1_index+n}];
    end
end


end