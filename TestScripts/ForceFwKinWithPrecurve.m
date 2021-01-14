clear; close all;
% TUBE PARAMETERS %
n = 5; % number of notches
ro = (1.62E-3)/2; % [mm] outer radius
ri = (1.4E-3)/2; % [mm] inner radius
h = 0.8E-3; % [mm] notch height
g = 1.46E-3; % [mm] notch depth
c = 2.2E-3; % [mm] uncut section height
E = 40E9; % [Pa] linear elastic modulus
mu = 0.2; % coefficient of friction
F = 1.15; % [N] input Force
theta_vec_base = deg2rad(0)*ones(n,1); % this represents the precurvature in each notch

% INITIALIZATION OF VECTORS %
theta_vec = theta_vec_base; % vector of each notch angle which should start at precurve value
F_vec = zeros(n,1); % vector of force experienced by each notch
M_vec = zeros(n,1); % vector of moment experienced by each notch
E_vec = E*ones(n,1); % effective elastic modulus for each notch
[ybar,I] = GetNeutralAxis(ro,ri,g); % neutral axis information (uniform notches)

% INTIALIZATION OF GRADIENT DESCENT %
theta_diff = 100; % start the change in theta high
trials = 0; % trials we have completed
maxTrials = 10; % maximum number of trials to perform
theta_last = theta_vec;
while theta_diff > 1E-6 && trials < maxTrials
    for ii = 1:n
        % Compute distal force and moment
        F_vec(ii) = F*exp(-mu*sum(theta_vec(1:ii))); % get force experienced by notch
        M_vec(ii) = F_vec(ii)*(ybar +ri); % get moment experienced by notch
        
        pct = 100;
        k = 1;
        % Gradient descent for non-linear modulus
        while k<100 && pct>1E-4
            
            % Compute angular deflection of current segment
            theta_vec(ii) = M_vec(ii)*h/(E_vec(ii)*I);
            
            if (h/(ro+ybar)<=theta_vec(ii))
                theta_vec(ii) = h/(ro+ybar);
            end
            % Compute section arc length, curvature and
            % strain
            s1 = h-(ybar)*abs(theta_vec(ii));
            kappa = theta_vec(ii)/s1;
            epsilon = (kappa.*(ro-ybar))./(1+ybar*kappa);
            
            % Update modulus via gradient descent
            [stress_eff,eta] = GetStress(abs(epsilon),E,0.08*E);
            new_E = E_vec(ii)-eta*(E_vec(ii)-stress_eff/abs(epsilon));
            
            % Percent change in modulus (for convergence
            % checking)
            pct = abs((new_E-E_vec(ii))/E_vec(ii));
            
            % Update modulus guess
            E_vec(ii) = new_E;
            
            % Increment
            k = k+1;
        end
        
    end
    theta_vec = theta_vec + theta_vec_base; % add offset for precurvature
    theta_delta = sum(theta_vec)-sum(theta_last);
    trials = trials + 1; % increment
end
rad2deg(theta_vec)

sheath = JoshWrist();
sheath.ConstructWrist('T_0',eye(4,4),'g',g*ones(n,1),'c',c*ones(n,1),'b',c,'h',h*ones(n,1),'h_c',h*ones(n,1),...
    'r_o',ro,'r_i',ri,'n',n,'material','nitinol','plotkin',false,'verbose',false);
sheath.GetKinematicsForce(F);
rad2deg(sheath.theta)

%% Showing effect of adding precurvature
clc;close all;
% TUBE PARAMETERS %
n = 5; % number of notches
ro = (1.62E-3)/2; % [mm] outer radius
ri = (1.4E-3)/2; % [mm] inner radius
h = 0.8E-3*ones(n,1); % [mm] notch height
g = 1.4E-3*ones(n,1); % [mm] notch depth
c = 1.2E-3*ones(n,1); % [mm] uncut section height
E = 10E9; % [Pa] linear elastic modulus
mu = 0.4; % coefficient of friction
F = 1.15; % [N] input Force
theta_vec_base = deg2rad(0)*ones(n,1); % this represents the precurvature in each notch

wrist = Wrist(ro*2,ri*2,n,h,zeros(n,1),c,g,'CutType','on-axis');

points = 151;
F_max = 6;
F = linspace(0,F_max,points);

notch_angles = zeros(n,points);
for i = 1:points
    wrist.fwkin([F(i),0,0]);
    notch_angles(:,i) = wrist.theta;
end


notch_angles_precurve = zeros(n,points);
wrist.fwkin([F(10),0,0]);
precurve_theta = wrist.theta;
wrist.precurve_theta = precurve_theta;
for i = 1:points
    wrist.fwkin([F(i),0,0]);
    notch_angles_precurve(:,i) = wrist.theta;
end

figure(2)
for i = 1:n+1
    subplot(2,3,i);
    if i < n+1
        title(sprintf("Notch %d Experimental Results Tip First",i),'FontSize',16);
        hold on
        plot(F,rad2deg(notch_angles(i,:)),'r','Linewidth',2);
        plot(F,rad2deg(notch_angles_precurve(i,:)),'g','Linewidth',2);
        legend("Base Model","Precurvature Model",'Location','southeast','FontSize',12)
        xlabel("Force (N)",'FontSize',14);
        ylabel("Notch Deflection (deg)",'FontSize',14)
        set(gca,'FontSize',12)
        grid on
    else
        title(sprintf("Total Deflection Experimental Results"),'FontSize',18);
        hold on
        plot(F,rad2deg(sum(notch_angles)),'r','Linewidth',2);
        plot(F,rad2deg(sum(notch_angles_precurve)),'g','Linewidth',2);
        legend("Base Model","Precurvature Model",'Location','southeast','FontSize',14)
        xlabel("Force (N)",'FontSize',16);
        ylabel("Tip Deflection (deg)",'FontSize',16)
        set(gca,'FontSize',14)
        grid on
    end
    hold off
end

