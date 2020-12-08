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

wrist = Wrist(od,id,n,h,phi,c,g,'CutType','on-axis'); % Nick's wrist class

T_0 = eye(4,4);
b = c(1);    
sheath = JoshWrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h,'h_c',h,...
    'r_o',od/2,'r_i',id/2,'n',n,'material','nitinol','plotkin',false,'verbose',false);

%% Comparing my Wrist Class to Josh's


% **** PLOTTING NOTCH ANGLE WRT FORCE **********
sheath.FindMaxForce(1,5);
des_points = zeros(n,2);
theta_last = zeros(n,1);
points = 100; % how many points to plot
theta_mat = zeros(n,points); % Initializing empty array to store values
F = linspace(0,sheath.F_max,points); % Getting list of forces (x axis values)
for i = 1:points
    sheath.GetKinematicsForce(F(i)); % Updating tube position
    theta_mat(:,i) = sheath.theta; % Getting tube position
    
    % *** DATA COLLECTION TO SEE IF WE CROSSED DESIRED ANGLE VALUE ***
    diff = (theta_last-theta_des') .* (sheath.theta-theta_des'); 
    for k = 1:n % check each notch individually
        if diff(k) <= 0 % We have crosed over the desired angle
            % Grab the point that was closest to the desired angle
            if abs(theta_last(k)-theta_des(k)) < abs(sheath.theta(k) - theta_des(k))
                des_points(k,1) = F(i-1);
                des_points(k,2) = rad2deg(theta_last(k));
            else
                des_points(k,1) = F(i);
                des_points(k,2) = rad2deg(sheath.theta(k));
            end
        end
    end
    theta_last = sheath.theta;
    % ****************************************************************
end
