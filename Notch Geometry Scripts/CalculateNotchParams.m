%% Initial Conditions
clc; clear; close all;

n = 5; % number of notches
totalWristBend = 90*pi/180; % [rad] - total wrist deflection
maxBendPerNotch = totalWristBend*2/(n*(n+1)).*(1:n); % [rad] - maximum bending angle per notch. Tip Dominant
% maxBendPerNotch = totalWristBend./(n*ones(1,n));  % [rad] - maximum bending angle per notch. Constant Curvature

% Desired theta setpoint for each notch.
% This creates a distance during tendon displacement where the notches in
% the wrist will take on these values. 
% theta_des = maxBendPerNotch/2;  % [rad] - Uniform Bending
theta_des = maxBendPerNotch; % [rad] - close together (tip first)

wristLength = 10e-3; %[m] Length of wrist excluding any base.
od = 1.1E-3; % [m] - outer diameter of tube
id = 0.9E-3; % [m] - inner diameter of tube
phi = zeros(n,1);
maxG =  0.875*od; % [m] - Max depth which we assign to notch n.
cutType = 'on-axis';
[~,h,u] = GetNotchSynthesis(maxBendPerNotch,maxG,od,id,'CutType',cutType,...
    'L',wristLength);
maxForce = 5.5; %[N] for plotting

% Based on Pacheco et al. JMRR. 2021.
E_lin = 10E9; % [N/m^2] - Elastic Modulus of Nitinol
E_se = 3e9; % [N/m^2] - Slope of Super Elastic Region for Nitinol
mu = 0.13; % coefficient of friction for capstan
strain_lower = 0.028;
% Default Values
% E_lin = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
% E_se = 0.08*E_lin; % [N/m^2] - Slope of Super Elastic Region for Nitinol
% mu = 0.2; % coefficient of friction for capstan
% strain_lower = 0.02;


%% Gradient Descent to determine the cut depth per notch

% Get the resulting ybar and I from our chosen cut depth
[maxYbar, minI] = GetNeutralAxis(od/2, id/2, maxG,'CutType',cutType);
% Determin the approximated linear elastic modulus
[strainFin, stressFin, E] = GetStrainInformation(theta_des(n), h(n), od/2, maxYbar,...
    'E_lin',E_lin, 'E_se',E_se,'strainLower',strain_lower);
% Determine force necessary to achieve desired bend at specifically notch n
% Inverting equation 12 in Pacheco et al. JMRR. 2021
Fdesired = theta_des(n)*E*minI/(h(n)*(id/2 + maxYbar)*exp(-mu*sum(theta_des)));


% *** FOR DATA COLLECTION ***
E_exp = zeros(1,n);
ybar_exp = zeros(1,n);
I_exp = zeros(1,n);
stress_exp = zeros(1,n);
strain_exp = zeros(1,n);
% initialize gradient descent To determine cut depth
g_vec = zeros(1,n);
for i = 1:n
    g = 0.9*od;
    stepsize = 0.00001; % Step size needs to be small
    error = 10;
    steps = 0;
    while(abs(error) > 1E-9)
        % Same process as up above to determine how much force would be
        % necessary to pull notch (i) given its current cut depth, theta
        % target, and other wrist parameters.
        [ybar, I] = GetNeutralAxis(od/2, id/2, g,'CutType',cutType);
        [strain, stress, E] = GetStrainInformation(theta_des(i), h(i), od/2, ybar,...
                'E_lin',E_lin, 'E_se',E_se,'strainLower',strain_lower);
        Fpw = theta_des(i)*E*I/(h(i)*(id/2 + ybar)*exp(-mu*sum(theta_des(1:i))));
        
        % Need our pullwire force to get close to our desired force
        error = Fpw - Fdesired;
        g = g + stepsize*error;
        steps = steps + 1;
        % *** DATA COLLECTION ***
        E_exp(i) = E;
        ybar_exp(i) = ybar;
        I_exp(i) = I;
        stress_exp(i) = stress;
        strain_exp(i) = strain;
        % **************************
        if (steps >= 800)
            disp("max steps reached")
            disp(error);
            break;
        end
    end
    g_vec(i) = g;
end

%% *** TESTING THE ABOVE RESULTS ***

% Create instance of wrist object
sheath = Wrist(od,id,n,h,phi,u,g_vec,'CutType',cutType,'Name','TestWrist');
sheath.E_lin = E_lin;
sheath.E_se = E_se;
sheath.mu = mu;
sheath.strain_lower = strain_lower;


% **** PLOTTING NOTCH ANGLE WRT FORCE **********
% sheath.FindMaxForce(1,5);
des_points = zeros(n,2);
theta_last = zeros(n,1);
points = 200; % how many points to plot
theta_mat = zeros(n,points); % Initializing empty array to store values
F = linspace(0,maxForce,points); % Getting list of forces (x axis values)
for i = 1:points
%     sheath.GetKinematicsForce(F(i)); % Updating tube position
    sheath.fwkin([F(i);0;0]);
    theta_mat(:,i) = sheath.theta; % Getting tube position
    
    % *** DATA COLLECTION TO SEE IF WE CROSSED DESIRED ANGLE VALUE ***
    diff = (theta_last-theta_des') .* (sheath.theta-theta_des'); 
    for k = 1:n % check each notch individually
        if diff(k) <= 0% We have crosed over the desired angle
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
theta_mat = theta_mat.*(180/pi); % Convert to deg
close all;
figure();
plot(F,theta_mat,'LineWidth',2);
title(sprintf("Notch angles with respect to force applied at tendon\ndesigned to close at the same time"),'FontSize',16)
xlabel("Force (N)",'FontSize',14)
ylabel("Angle (deg)",'FontSize',14)
% *** LABELING OUR DESIRED POINTS ***
labels = cell(n,1);
legendLabels = cell(n,1);
for i = 1:n
    labels(i) = cellstr(...
        sprintf("Angle: %.1f",theta_des(i)*180/pi));
    % Creating Legend Labels as well
    legendLabels(i) = cellstr(...
        sprintf("theta %u",i));
end
hold on
% Plot our desired force and the desired theta for each wrist
stem(Fdesired.*ones(1,n),theta_des.*180/pi,'ok','MarkerSize',10)
% Plot the locations where our wrist actually hit the desired theta.
stem(des_points(:,1),des_points(:,2),'.r','MarkerSize',10);
text(Fdesired.*ones(1,n)-0.01,theta_des.*180/pi,labels,'VerticalAlignment','bottom',...
    'HorizontalAlignment','right');
hold off
legend(legendLabels,'Location','northwest','FontSize',12);

% *** MEASURING THETA DIFFERENCE ***
% Commented out because it messes with previous plot
% Apply our desired force to the wrist tendon
 sheath.fwkin([Fdesired;0;0]);

% Print out percentage error for each notch
percentage_error = 100*(theta_des' - sheath.theta)./(theta_des');
s = sprintf("");
for i = 1:n
    s = s + sprintf("Error for notch %u: %0.2f %%\n",i,percentage_error(i));
end
disp(s)

disp(['Notch Width (w): ' num2str(g_vec*1e3) ' mm']);
disp(['Notch Height (h) : ' num2str(h'*1e3) ' mm']);
disp(['Notch Height (u) : ' num2str(u'*1e3) ' mm']);
