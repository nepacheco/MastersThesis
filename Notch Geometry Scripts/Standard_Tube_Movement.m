%% Initial Conditions
clc; clear; close all;
maxBendPerNotch = 30; % [deg] - angle per notch
od = 1.62E-3; % [m] - outer diameter of tube
id = 1.4E-3; % [m] - inner diameter of tube
n = 5; % number of notches
maxG = 0.85*od;
g_vec = maxG.*ones(n,1); % [m] - Max depth which we assign to notch n.
[~,h,u] = GetNotchSynthesis(maxBendPerNotch,maxG,od,id,n);
E_lin = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
E_se = 0.08*E_lin; % [N/m^2] - Slope of Super Elastic Region for Nitinol
mu = 0.2; % coefficient of friction for capstan

T_0 = eye(4);
N = n;                 % Number of notches
r_o = od/2;         % Tube Outer Radius [m]
r_i = id/2;         % Tube Inner Radius [m]
h_w = linspace(h,h,N);    % Notch Width Vector [m]
h_c = linspace(h,h,N);    % Notch Collision Width Vector (<h) [m]
g = g_vec;     % Notch depth vector which we just determined [m]
c = (u).*ones(N,1);             % Notch spacing vector [m]
b = u;             % Distal offset [m]
FOS = 1;                % Factor of Safety

% Create instance of wrist object
sheath = Wrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h_w,'h_c',h_c,...
    'r_o',r_o,'r_i',r_i,'n',N,'material','nitinol','plotkin',false,'verbose',false);

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
end
theta_mat = theta_mat.*(180/pi); % Convert to deg
close all;
figure();
plot(F,theta_mat);
title("Notch angles with respect to force applied at tendon",'FontSize', 16)
xlabel("Force (N)",'FontSize', 12)
ylabel("Angle (deg)", 'FontSize', 12)
% *** LABELING OUR DESIRED POINTS ***
labels = {};
legendLables = {};
for (i = 1:n)
    % Creating Legend Labels as well
    legendLabels(i) = cellstr(...
        sprintf("theta %u",i));
end
legend(legendLabels,'Location','northwest');
