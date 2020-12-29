%% Solving for E_lin, E_se, mu, and strain lower
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
g = 1.4E-3*ones(n,1); % Diagrma sent by Pulse might be 1.45 [m] - Max depth which we assign to notch n.
h = 0.8E-3*ones(n,1);
c = 1.2E-3*ones(n,1);
E_lin = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
E_se = 0.08*E_lin; % [N/m^2] - Slope of Super Elastic Region for Nitinol
mu = 0.2; % coefficient of friction for capstan
T_0 = eye(4,4);
b = c(1);


wrist = Wrist(od,id,n,h,phi,c,g,'CutType','on-axis'); % Nick's wrist class
file_path2 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-12-2020_Experiment\12-12-2020_Results.xlsx";
file_path1 = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-7-2020_Experiment\12-7-2020_Results.xlsx";
opts = detectImportOptions(file_path2);
opts.Sheet = 'AvgMeasurements';
file2 = readcell(file_path2,opts);
file1 = readcell(file_path1,opts);

[force_vec2,notch_mat2] = parseFile(file2,n);
% [force_vec1,notch_mat1] = parseFile(file1,n);
tic
min_norm_rmse = [100 0 0 0 0];
for E_lin = linspace(10E9,40E9,5)
    for strain_lower = linspace(0.010,0.030,5)
        for scale = linspace(0.1,0.35,5)
            for mu = linspace(0.1,0.5,5)
                wrist.E_lin = E_lin;
                wrist.E_se = scale*E_lin;
                wrist.strain_lower = strain_lower;
                wrist.mu = mu;
                diff = zeros(n+1,length(force_vec2));
                for i = 1:length(force_vec2)
                    input = force_vec2(i);
                    wrist.fwkin([input,0,0],'Type','force');
                    diff(:,i) = notch_mat2(i,:)' -...
                        rad2deg([wrist.theta; sum(wrist.theta)]);
                end
                se = diff.^2;
                mse = mean(se,2);
                rmse = sqrt(mse);
                if (norm(rmse) < min_norm_rmse(1))
                    min_norm_rmse(1) = norm(rmse);
                    min_norm_rmse(2) = E_lin/(1E9);
                    min_norm_rmse(3) = strain_lower;
                    min_norm_rmse(4) = scale;
                    min_norm_rmse(5) = mu;
                end
            end
        end
    end
end
toc
points = 100;
wrist.E_lin = min_norm_rmse(2)*1E9;
wrist.E_se = min_norm_rmse(4)*min_norm_rmse(2)*1E9;
wrist.strain_lower = min_norm_rmse(3);
wrist.mu = min_norm_rmse(5);
theta_mat_force = zeros(n,points);
F_vec = linspace(0,6,points); % Getting list of forces (x axis values)
for i = 1:points
    wrist.fwkin([F_vec(i),0,0],'Type','force');
    theta_mat_force(:,i) = wrist.theta;
end
figure(2)
for i = 1:n+1
    subplot(3,2,i);
    if i < n+1
        title(sprintf("Notch %d Experimental Results",i),'FontSize',16);
        hold on
        %         scatter(force_vec1,notch_mat1(:,i),'bo');
        scatter(force_vec2,notch_mat2(:,i),'rx');
        plot(F_vec,rad2deg(theta_mat_force(i,:)),'g','Linewidth',2);
        legend("Experiment2","Model",'Location','southeast','FontSize',12)
        xlabel("Force (N)",'FontSize',14);
        ylabel("Notch Deflection (deg)",'FontSize',14)
        set(gca,'FontSize',12)
        grid on
    else
        title(sprintf("Total Deflection Experimental Results"),'FontSize',18);
        hold on
        %         scatter(force_vec1,notch_mat1(:,i),'bo');
        scatter(force_vec2,notch_mat2(:,i),'rx');
        plot(F_vec,rad2deg(sum(theta_mat_force(:,:))),'g','Linewidth',2);
        legend("Experiment2","Model",'Location','southeast','FontSize',14)
        xlabel("Force (N)",'FontSize',16);
        ylabel("Tip Deflection (deg)",'FontSize',16)
        set(gca,'FontSize',14)
        grid on
    end
    hold off
end
%%
diff = zeros(n+1,length(force_vec2));
for i = 1:length(force_vec2)
    input = force_vec2(i);
    wrist.fwkin([input,0,0],'Type','force');
    diff(:,i) = notch_mat2(i,:)' -...
        rad2deg([wrist.theta; sum(wrist.theta)]);
end
se = diff.^2;
mse = mean(se,2);
rmse = sqrt(mse);
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