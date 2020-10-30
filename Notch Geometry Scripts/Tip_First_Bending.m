%% Initial Conditions
clc; clear; close all;
maxBendPerNotch = 30; % [deg] - angle per notch
od = 1.62E-3; % [m] - outer diameter of tube
id = 1.4E-3; % [m] - inner diameter of tube
n = 5; % number of notches
maxG = 0.9*od; % [m] - Max depth which we assign to notch n.
[~,h,u] = GetNotchSynthesis(maxBendPerNotch,maxG,od,id,n);
E_lin = 60E9; % [N/m^2] - Elastic Modulus of Nitinol
E_se = 0.08*E_lin; % [N/m^2] - Slope of Super Elastic Region for Nitinol
mu = 0.2; % coefficient of friction for capstan


%% Tip First Bending but incorporating different elastic modulus per notch

% Get the resulting ybar and I from our chosen cut depth
[maxYbar, minI] = GetNeutralAxis(od/2, id/2, maxG);
% Desired theta for each notch. If n changes change the size of this vector
theta_des = [15 16 17 18 19].*pi/180;
% Determin the approximated linear elastic modulus
[stress, strain, E] = GetStrainInformation(theta_des(5), h, od/2, maxYbar,E_lin,E_se);
% Determine force necessary to achieve desired bend at specifically notch n
Fdesired = theta_des(n)*E*minI/(h*(id/2 + maxYbar)*exp(-mu*sum(theta_des)));


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
        % Same process as up above
        [ybar, I] = GetNeutralAxis(od/2, id/2, g);
        [stress, strain, E] = GetStrainInformation(theta_des(i), h, od/2, ybar,E_lin,E_se);
        Fpw = theta_des(i)*E*I/(h*(id/2 + ybar)*exp(-mu*sum(theta_des(1:i))));
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
g_vec

% *** TESTING THE ABOVE RESULTS ***
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
    % *** DATA COLLECTION TO SEE IF WE CROSSED DESIRED ANGLE VALUE ***
    diff = (theta_last-theta_des') .* (sheath.theta-theta_des'); 
    for k = 1:n % check each notch individually
        if diff(k) <= 0 % We have crosed over the desired angle
            disp(theta_last)
            % Grab the point that was closest to the desired point
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
end
theta_mat = theta_mat.*(180/pi); % Convert to deg
close all;
figure();
plot(F,theta_mat);
title("Notch angles with respect to force applied at tendon Using Josh's Model")
xlabel("Force (N)")
ylabel("Angle (rad)")
% *** LABELING OUR DESIRED POINTS ***
labels = {};
legendLables = {};
for (i = 1:n)
    labels(i) = cellstr(...
        sprintf("Angle: %.1f",theta_des(i)*180/pi));
    % Creating Legend Labels as well
    legendLabels(i) = cellstr(...
        sprintf("theta %u",i));
end
hold on
stem(Fdesired.*ones(1,n),theta_des.*180/pi,'ok','MarkerSize',10)
stem(des_points(:,1),des_points(:,2),'r','MarkerSize',5);
text(Fdesired.*ones(1,n)-0.01,theta_des.*180/pi,labels,'VerticalAlignment','bottom',...
    'HorizontalAlignment','right');
hold off
legend(legendLabels,'Location','northwest');

% % *** MEASURING THETA DIFFERENCE ***
% % Commented out because it messes with previous plot
% % Apply our desired force to the wrist tendon
% sheath.GetKinematicsForce(Fdesired);
% 
% % Print out percentage error for each notch
% percentage_error = sheath.theta./(theta_des');
% s = sprintf("");
% for i = 1:n
%     s = s + sprintf("Percentage error for notch %u, %f\n",i,percentage_error(i));
% end
% disp(s)

disp(['Notch Width (w): ' num2str(g_vec) ' m']);
disp(['Notch Height (h) : ' num2str(h) ' m']);
disp(['Notch Height (u) : ' num2str(u) ' m']);
