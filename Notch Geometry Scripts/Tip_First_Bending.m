%% Initial Conditions
clc; clear; close all;
maxBendPerNotch = 10; % [deg] - angle per notch
od = 1.62E-3; % [m] - outer diameter of tube
id = 1.4E-3; % [m] - inner diameter of tube
n = 9; % number of notches
maxG = 0.9*od; % [m] - Max depth per notch
[~,h,u] = GetNotchSynthesis(maxBendPerNotch,maxG,od,id,n);
E = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
mu = 0.2; % coefficient of friction for capstan


%% Tip First Bending but incorporating different elastic modulus per notch
% Determine force necessary to achieve desired bend
[maxYbar, minI] = GetNeutralAxis(od/2, id/2, maxG);
theta_des = [0*ones(1,n)].*pi/180;
[stress, strain, E] = GetStrainInformation(theta_des(5), h, od/2, maxYbar);
Fclosed = theta_des(5)*E*minI/(h*(id/2 + maxYbar)*exp(-mu*sum(theta_des)));

% initialize gradient descent To determine cut depth
g_vec = zeros(1,n);
E_exp = zeros(1,n);
ybar_exp = zeros(1,n);
I_exp = zeros(1,n);
stress_exp = zeros(1,n);
strain_exp = zeros(1,n);
for i = 1:n
    g = 0.9*od;
    stepsize = 0.00001;
    error = 10;
    steps = 0;
    while(abs(error) > 1E-9)
        [ybar, I] = GetNeutralAxis(od/2, id/2, g);
        [stress, strain, E] = GetStrainInformation(theta_des(i), h, od/2, ybar);
        Fpw = theta_des(i)*E*I/(h*(id/2 + ybar)*exp(-mu*sum(theta_des(1:i))));
        error = Fpw - Fclosed;
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
g_vec;

% *** TESTING THE ABOVE RESULTS ***
T_0 = eye(4);
N = n;                 % Number of notches
r_o = od/2;         % Tube Outer Radius [m]
r_i = id/2;         % Tube Inner Radius [m]
h_w = linspace(h,h,N);    % Notch Width Vector [m]
h_c = linspace(h,h,N);    % Notch Collision Width Vector (<h) [m]
g = g_vec;     % Notch depth vector [m]
c = (u).*ones(N,1);             % Notch spacing vector [m]
b = u;             % Distal offset [m]
FOS = 1;                % Factor of Safety

% Create instance of wrist object
sheath = Wrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h_w,'h_c',h_c,'r_o',r_o,'r_i',r_i,'n',N,'material','nitinol','plotkin',false);
sheath.GetKinematicsForce(Fclosed);

% Print out percentage error for each notch
percentage_error = sheath.theta./(theta_des');
s = sprintf("");
for i = 1:n
    s = s + sprintf("Percentage error for notch %u, %f\n",i,percentage_error(i));
end
disp(s)

% **** PLOTTING NOTCH ANGLE WRT FORCE **********
sheath.FindMaxForce(1,5);
points = 100;
theta_mat = zeros(n,points);
F = linspace(0,sheath.F_max,points);
for i = 1:points
    sheath.GetKinematicsForce(F(i));
    theta_mat(:,i) = sheath.theta;
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
        sprintf("Angle: %.2f\nForce: %.2f",theta_des(i)*180/pi,Fclosed));
    legendLabels(i) = cellstr(...
        sprintf("theta %u",i));
end
hold on
plot(Fclosed.*ones(1,n),theta_des.*180/pi,'ok','MarkerSize',12)
text(Fclosed.*ones(1,n),theta_des.*180/pi,labels,'VerticalAlignment','bottom',...
    'HorizontalAlignment','right');
hold off
legend(legendLabels,'Location','northwest');


disp(['Notch Width (w): ' num2str(g_vec) ' m']);
disp(['Notch Height (h) : ' num2str(h) ' m']);
disp(['Notch Height (u) : ' num2str(u) ' m']);