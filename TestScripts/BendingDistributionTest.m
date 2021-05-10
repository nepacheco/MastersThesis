clc; clear; close all

% TUBE PARAMS %
n = 5; % number of notches
ro = (1.62E-3)/2; % [mm] outer radius
ri = (1.4E-3)/2; % [mm] inner radius
h = 0.7915E-3; % [mm] notch height
g = 1.4E-3; % [mm] notch depth
c = 2.2085E-3; % [mm] uncut section height
E = 40E9; % [Pa] linear elastic modulus
mu = 0.2; % coefficient of friction
delta_l = 0.5E-3; % [mm] tendon displacement
theta_vec_base = deg2rad(2)*ones(n,1); % this represents the precurvature in each notch

wrist = Wrist(ro*2,ri*2,n,h*ones(n,1),zeros(n,1),c*ones(n,1),g*ones(n,1));
ybar = wrist.get_neutral_axis();
p = zeros(3,2);

points = 3;
theta_max = deg2rad(90);
theta = zeros(n,points);
theta(:,1) = theta_max/n*ones(n,1);
theta(:,2) = [0.27*theta_max 0.22*theta_max 0.19*theta_max 0.172*theta_max 0.148*theta_max]';
theta(:,3) = [0.1,0.14,0.19,0.25,0.32].*theta_max;

colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250];

for j = 1:points
    if(sum(theta(:,j)) - theta_max > 1E-6)
        disp("error")
    end
    kappa = theta(:,j)./(h-ybar(1).*theta(:,j));
    s = zeros(n,1);
    for i = 1:n
        if(kappa(i)>0)
            s(i) = theta(i,j)./kappa(i);
        else
            s(i) = h;
        end
    end
    wrist.kappa = kappa;
    wrist.s = s;
    wrist.alpha = 0;
    wrist.tau = 0;
    [~,T] = wrist.robot_kin();
    p(:,j) = T(1:3,4)*1000;
    wrist.plot_stick_model('LineWidth',3,'Marker','none','Color',colors(j,:));
    hold on
end
view(0,0)
legend('Equal Bending','Base First Bending','Tip First Bending','FontSize',16,'FontName','CMU Serif','Location','southeast');
title(sprintf('Tube position based on \ndistribution of wrist deflection'),'FontSize',20,'FontName','CMU Serif');
ax = gca;
ax.FontSize = 16;
ax.FontName = 'CMU Serif';

fprintf("Difference in tip positions: %f\n",norm(p(:,1)-p(:,2)))
fprintf("Difference in tip positions: %f\n",norm(p(:,1)-p(:,3)))