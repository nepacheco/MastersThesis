%% Test script for Wrist() class

clc;clear all;
T_0 = eye(4);   % Arbitrary initial transform

% NOTE: the values below are from the York et al. design for comparison to
% experimental results. If you change the values, be sure to comment out
% the plotting of the experimental results (as they won't be equivalent)

N = 5;                 % Number of notches
r_o = 1.16E-3/2;         % Tube Outer Radius [m]
r_i = 0.86E-3/2;         % Tube Inner Radius [m]
h = linspace(0.51E-3,0.51E-3,N);    % Notch Width Vector [m]
h_c = linspace(.51E-3,.51E-3,N);    % Notch Collision Width Vector (<h) [m]
g = linspace(0.97E-3,0.97E-3,N);     % Notch depth vector [m]
c = (.51E-3).*ones(N,1);             % Notch spacing vector [m]
b = 1.5E-3;             % Distal offset [m]
FOS = 1;                % Factor of Safety

% Create instance of wrist object
newSheath = JoshWrist();
newSheath.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h,'h_c',h_c,'r_o',r_o,'r_i',r_i,'n',N,'material','nitinol','plotkin',true);

% Find maximum force for desired factor-of-safety
newSheath.FindMaxForce(FOS,5);
F_max = newSheath.F_max;
newSheath.PrintResults();

%% Create solution vector for force sweep
% Here we are just linearly sweeping up to the maximum force for parametric
% analysis

% Create vector up to maximum force
n_pts = 40;
Force = linspace(0,F_max,n_pts);

% Allocate space for solution vectors
stress_vec = [];
strain_vec = [];
theta_vec = [];

fprintf('\nCreating force sweep\n');
for i = 1:n_pts
    if i == n_pts
        % Only plot kinematics on last iteration
        figure(1);
        newSheath.plotkin = true;
        fprintf('%i\n',i);
    else
        newSheath.plotkin = false;
        fprintf('%i,',i);
    end
    
    % Run mechanics/kinematics model and extract results
    newSheath.GetKinematicsForce(Force(i));
    stress_vec(i,:) = newSheath.sig;
    strain_vec(i,:) = newSheath.strain;
    theta_vec(i) = sum(newSheath.theta);
end


%%
% Data from York et al (comment out if deviating from York design)
DATA = [5.9	0.004154105
11	0.357560665
15.4	0.776502946
22.9	1.20881778
28.8	1.622529453
35.3	1.947777032
45	2.220738834
55	2.409689019
65.6	2.541027722
75.1	2.66422041
85.1	2.816352783
95.1	2.949542383
105.6	3.118125643
114.9	3.247102412
122.9	3.405919401
131.6	3.632075577
137.7	3.880093348
143.7	4.202531803];
% 

theta_exp = DATA(:,1)-6;
Force_exp = DATA(:,2);

figure(2)
subplot(1,3,1)
plot(theta_vec.*180./pi,Force,'-');hold on;
plot(theta_exp,Force_exp,'ko','MarkerFaceColor','w');
legend('Model','Experiment (York et al.)','Location','SouthEast');
xlabel('Tip Deflection [deg]');
ylabel('Tendon Force [N]');
grid on;

subplot(1,3,2)
surf([1:1:N],theta_vec.*180./pi,stress_vec./1000000,'EdgeAlpha',.2)
xlabel('Joint');
ylabel('Deflection [deg]');
zlabel('Joint Stress [MPa]');
subplot(1,3,3)
surf([1:1:N],theta_vec.*180./pi,strain_vec,'EdgeAlpha',.2);
xlabel('Joint');
ylabel('Deflection [deg]');
zlabel('Joint Strain');
set(gcf,'Color','w');
set(gcf,'Position',[100 100 1200 250]);

% % Uncomment to save as .pdf
% set(gcf,'renderer','painters');
% set(gcf,'Units','Inches');
% pos = get(gcf,'Position');
% set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(gcf,'model_nitinol_tapered','-dpdf')

