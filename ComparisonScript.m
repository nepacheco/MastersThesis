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
g = 1.46E-3*ones(n,1); % Diagrma sent by Pulse might be 1.45 [m] - Max depth which we assign to notch n.
h = 0.8E-3*ones(n,1);
c = 2.2E-3*ones(n,1);
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
sheath.FindMaxForce(1,5);
theta_last = zeros(n,1);
points = 100; % how many points to plot
diff_values = zeros(n+1,points);
theta_mat_sheath = zeros(n,points); % Initializing empty array to store values
theta_mat_wrist = zeros(n,points);
F = linspace(0,sheath.F_max,points); % Getting list of forces (x axis values)
for i = 1:points
    sheath.GetKinematicsForce(F(i)); % Updating tube position
    theta_mat_sheath(:,i) = sheath.theta; % Getting tube position
    wrist.fwkin([F(i),0,0],'Type','force');
    theta_mat_wrist(:,i) = wrist.theta;
    if (norm(wrist.theta-sheath.theta) > 1E-4)
        diff_values(:,i) = [wrist.theta-sheath.theta; F(i)];
    end
end

for i = 1:n
    subplot(2,3,i);
    title(sprintf("Notch %d comparison",i));
    hold on
    plot(F,theta_mat_sheath(i,:));
    plot(F,theta_mat_wrist(i,:));
    legend('Josh Model','Nick Model');
    grid on
    hold off
end

disp(diff_values)

%% Comparison of single force spot with both models
F = 1.15;
wrist = Wrist(od,id,n,h,phi,c,g,'CutType','off-axis'); % Nick's wrist class
sheath = JoshWrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h,'h_c',h,...
    'r_o',od/2,'r_i',id/2,'n',n,'material','nitinol','plotkin',false,'verbose',false);
sheath.GetKinematicsForce(F); % Updating tube position
disp("MY MODEL")
wrist.fwkin([F,0,0],'Type','force');

diff = wrist.theta - sheath.theta

wrist_fwkin = readmatrix("wrist_fwkin.csv");
josh_fwkin = readmatrix("josh_fwkin.csv");

% If results look good move on to next section

%% Comparing Test Results to Wrist Model
wrist = Wrist(od,id,n,h,phi,c,g,'CutType','on-axis'); % Nick's wrist class
file_path = "C:\Users\nickp\OneDrive - Worcester Polytechnic Institute (wpi.edu)\School Files\Thesis\ForceTest\12-7-2020_Experiment\12-7-2020_Results.xlsx";
opts = detectImportOptions(file_path);
opts.Sheet = 'AvgMeasurements';
file = readmatrix(file_path);

force = rmmissing([file(1:47,3); file(63:75,3); file(91:106,3); file(124:135,3)]);
notch_matrix = rmmissing(file(:,4:8))';


points = 100; % how many points to plot
diff_values = zeros(n+1,points);
theta_mat_force = zeros(n,points);
F = linspace(0,6,points); % Getting list of forces (x axis values)
for i = 1:points
    wrist.fwkin([F(i),0,0],'Type','force');
    theta_mat_force(:,i) = wrist.theta;
    sheath.GetKinematicsForce(F(i)); % Updating tube position
    theta_mat_sheath(:,i) = sheath.theta; % Getting tube position
end

% Generating RMSE
diff = zeros(n+1,length(force));
for i = 1:length(force)
    input = force(i);
    wrist.fwkin([F(i),0,0],'Type','force');
    diff(:,i) = [notch_matrix(:,i); sum(notch_matrix(:,i))] -...
        [wrist.theta; sum(wrist.theta)];
end
se = diff.^2;
mse = mean(diff')';
rmse = sqrt(mse);

% Plotting
figure(2)
for i = 1:n+1
    subplot(2,3,i);
    if i < n+1
        title(sprintf("Notch %d Experimental Results",i),'FontSize',16);
        hold on
        scatter(force,notch_matrix(i,:));
        plot(F,rad2deg(theta_mat_force(i,:)));
%         plot(F,rad2deg(theta_mat_sheath(i,:)));
        legend("Experiment","Model",'Location','southeast','FontSize',12)
        xlabel("Force (N)",'FontSize',14);
        ylabel("Notch Deflection (deg)",'FontSize',14)
        grid on
    else
        title(sprintf("Total Deflection Experimental Results"),'FontSize',16);
        hold on
        scatter(force,sum(notch_matrix(:,:)));
        plot(F,rad2deg(sum(theta_mat_force(:,:))));
        legend("Experiment","Model",'Location','southeast','FontSize',12)
        xlabel("Force (N)",'FontSize',14);
        ylabel("Tip Deflection (deg)",'FontSize',14)
        grid on
    end
    hold off
end