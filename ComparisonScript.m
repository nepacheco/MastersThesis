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

wrist = Wrist(od,id,n,h,phi,c,g,'CutType','off-axis'); % Nick's wrist class

T_0 = eye(4,4);
b = c(1);    
sheath = JoshWrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h,'h_c',h,...
    'r_o',od/2,'r_i',id/2,'n',n,'material','nitinol','plotkin',false,'verbose',false);

%% Comparing my Wrist Class to Josh's


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

% If results look good move on to next section

%% Comparing Test Results to Wrist Model
wrist = Wrist(od,id,n,h,phi,c,g,'CutType','on-axis'); % Nick's wrist class
