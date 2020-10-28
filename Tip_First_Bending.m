%% Initial Conditions
clc; clear; close all;
od = 1.1E-3; % [m] - outer diameter of tube
id = 0.9E-3; % [m] - inner diameter of tube
h = 0.53096E-3; % [m] - notch height
E = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
n = 5; % number of notches
mu = 0.2; % coefficient of friction for capstan

% General theta equation is
% $$theta = \frac{F(ro + ybar)h}{(EI)}$$
%% Determining necessary g for same angle at proximal and distal notches
% Determining the necessary pullwire force to have all notches closed at
% 30° with a 5 notch tube. Assume that the depth for the most distal notch
% is 90% of the outer diameter.
anglePerNotch = 20*pi/180; % [rad] - maximum angle per notch
maxG = 0.9*od; %[mm]
[maxYbar, minI] = GetNeutralAxis(od/2, id/2, maxG);
Fclosed = anglePerNotch*E*minI/(h*(id/2 + maxYbar)*exp(-mu*n*anglePerNotch));

% initialize gradient descent
g_vec = zeros(1,5);
for i = 1:5
    g = 0.9*od;
    stepsize = 0.00001;
    error = 10;
    steps = 0;
    Fpw = 0;
    while(abs(error) > 1E-9)
        [ybar, I] = GetNeutralAxis(od/2, id/2, g);
        Fpw = anglePerNotch*E*I/(h*(id/2 + ybar)*exp(-mu*i*anglePerNotch));
        error = Fpw - Fclosed;
        g = g + stepsize*error;
        steps = steps + 1;
        if (steps >= 800)
            disp("max steps reached")
            break;
        end
    end
    g_vec(i) = g;
end
g_vec

% *** TESTING THE ABOVE RESULTS ***
T_0 = eye(4);
N = 5;                 % Number of notches
r_o = 1.1E-3/2;         % Tube Outer Radius [m]
r_i = 0.9E-3/2;         % Tube Inner Radius [m]
h_w = linspace(0.53096E-3,0.53096E-3,N);    % Notch Width Vector [m]
h_c = linspace(0.53096E-3,0.53096E-3,N);    % Notch Collision Width Vector (<h) [m]
g = g_vec;     % Notch depth vector [m]
c = (1E-3).*ones(N,1);             % Notch spacing vector [m]
b = 1.5E-3;             % Distal offset [m]
FOS = 1;                % Factor of Safety

% Create instance of wrist object
sheath = Wrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h_w,'h_c',h_c,'r_o',r_o,'r_i',r_i,'n',N,'material','nitinol','plotkin',false);
sheath.GetKinematicsForce(Fclosed);
theta = (sheath.theta.*180/pi)'


%% Tip First Bending but incorporating different elastic modulus per notch
clc;
anglePerNotch = 15*pi/180; % [rad] - maximum angle per notch
maxG = 9.9E-4; %[mm]
[maxYbar, minI] = GetNeutralAxis(od/2, id/2, maxG);
[stress, strain, E] = GetStrainInformation(anglePerNotch, h, od/2, maxG, maxYbar);
theta_des = anglePerNotch.*ones(1,5);
Fclosed = anglePerNotch*E*minI/(h*(id/2 + maxYbar)*exp(-mu*sum(theta_des)));

% initialize gradient descent To determine cut depth
g_vec = zeros(1,5);
E_exp = zeros(1,5);
ybar_exp = zeros(1,5);
I_exp = zeros(1,5);
stress_exp = zeros(1,5);
strain_exp = zeros(1,5);
for i = 1:5
    g = 0.9*od;
    stepsize = 0.00001;
    error = 10;
    steps = 0;
    while(abs(error) > 1E-9)
        [ybar, I] = GetNeutralAxis(od/2, id/2, g);
        [stress, strain, E] = GetStrainInformation(anglePerNotch, h, od/2, g, ybar);
        Fpw = anglePerNotch*E*I/(h*(id/2 + ybar)*exp(-mu*sum(theta_des(1:i))));
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
N = 5;                 % Number of notches
r_o = 1.1E-3/2;         % Tube Outer Radius [m]
r_i = 0.9E-3/2;         % Tube Inner Radius [m]
h_w = linspace(0.53096E-3,0.53096E-3,N);    % Notch Width Vector [m]
h_c = linspace(0.53096E-3,0.53096E-3,N);    % Notch Collision Width Vector (<h) [m]
g = g_vec;     % Notch depth vector [m]
c = (.51E-3).*ones(N,1);             % Notch spacing vector [m]
b = 1.5E-3;             % Distal offset [m]
FOS = 1;                % Factor of Safety

% Create instance of wrist object
sheath = Wrist();
sheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h_w,'h_c',h_c,'r_o',r_o,'r_i',r_i,'n',N,'material','nitinol','plotkin',false);
sheath.GetKinematicsForce(Fclosed);

% Print out percentage error for each notch
theta_exp = anglePerNotch*ones(1,5);
percentage_error = sheath.theta./(theta_exp');
s = sprintf("");
for i = 1:5
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
figure();
plot(F,theta_mat);
title("Notch angles with respect to force applied at Tendon Using Josh's Model")
xlabel("Force (N)")
ylabel("Angle (rad)")
label = sprintf("Angle: %.2f\nForce: %.2f",anglePerNotch*180/pi,Fclosed);
hold on
plot(Fclosed,anglePerNotch*180/pi,'ok','MarkerSize',12)
text(Fclosed,anglePerNotch*180/pi,label,'VerticalAlignment','bottom',...
    'HorizontalAlignment','right');
legend("theta1 (most proximal)","theta2","theta3","theta4","theta5 (most distal)")
%% Graphing Fclosed
anglePerNotch = 0.4;
maxG = 9.9E-4; %[mm]
[ybar, I] = GetNeutralAxis(od/2, id/2, maxG);
[stress, strain, E] = GetStrainInformation(anglePerNotch, h, od/2, maxG, ybar);
ratio_init = anglePerNotch*exp(mu*n*anglePerNotch)*E;


points = 200;
x = zeros(1,points^2);
theta_vec= [];
E_vec = [];
y = zeros(1,points^2);
z = zeros(1,points^2);
count = 1;
for anglePerNotch = linspace(0.001,30*pi/180,points)
    [stress, strain, E] = GetStrainInformation(anglePerNotch, h, od/2, maxG, ybar);
    for E_eff = linspace(E,40E9,points)
        Fclosed = anglePerNotch*E*I/(h*(od/2 + ybar)*exp(-mu*n*anglePerNotch));
        
        x(count) = anglePerNotch*exp(mu*n*anglePerNotch);
        y(count) = E_eff;
        z(count) = x(count)*y(count);
        if abs(ratio_init-z(count)) < 1
            theta_vec = [theta_vec anglePerNotch];
            E_vec = [E_vec E_eff];
        end
        count = count + 1;
    end
end

xv = linspace(min(x), max(x), 50);
yv = linspace(min(y), max(y), 50);
[X,Y] = meshgrid(xv, yv);
Z = griddata(x,y,z,X,Y);
surf(X,Y,Z)
xlabel("Angle per Notch")
ylabel("Effective Modulus Range")
zlabel("Force on pullwire")
