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
theta_max = deg2rad(60);
theta = zeros(n,points);
theta(:,1) = theta_max/n*ones(n,1);
theta(:,2) = [0.45*theta_max 0.25*theta_max 0.15*theta_max 0.1*theta_max 0.05*theta_max]';
theta(:,3) = flip(theta(:,2));

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
    wrist.plot_stick_model();
    xlabel('Position (m)');
    zlabel('Position (m)');
    xlim([0,15e-3]);
    ylim([0,1e-3]);
    zlim([0,15e-3]);
    hold on
end
legend('Equal Bending','Base First Bending','Tip First Bending');
fprintf("Difference in tip positions: %f\n",norm(p(:,1)-p(:,2)))
fprintf("Difference in tip positions: %f\n",norm(p(:,1)-p(:,3)))