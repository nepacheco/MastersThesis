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

%% Checking to see if my GetStress funtion is the same as Josh's
points = 100000;
strain = linspace(0,.12,points);
my_func = zeros(points,2);
josh_func = zeros(points,2);
for i = 1:points
    [my_func(i,1),my_func(i,2)] = GetStress(strain(i),E,0.08*E);
    [josh_func(i,1), josh_func(i,2)] = sheath.SuperElastic(strain(i));
    if isnan(my_func(i,1)) && isnan(josh_func(i,1))
        my_func(i,1) = my_func(i-1,1);
        josh_func(i,1) = josh_func(i-1,1);
    end
    
end

diff_mat = my_func - josh_func;
mean(diff_mat)
%% Comparing Josh Model vs my model
clear; close all; clc;
% TUBE PARAMETERS %
n = 5; % number of notches
ro = (1.62E-3)/2; % [mm] outer radius
ri = (1.4E-3)/2; % [mm] inner radius
h = 0.8E-3; % [mm] notch height
g = 1.46E-3; % [mm] notch depth
c = 2.2E-3; % [mm] uncut section height
E = 40E9; % [Pa] linear elastic modulus
mu = 0.2; % coefficient of friction
F = 1.144595415; % [N] input Force
theta_vec_base = deg2rad(0)*ones(n,1); % this represents the precurvature in each notch


sheath = JoshWrist();
sheath.ConstructWrist('T_0',eye(4,4),'g',g*ones(n,1),'c',c*ones(n,1),'b',c,'h',h*ones(n,1),'h_c',h*ones(n,1),...
    'r_o',ro,'r_i',ri,'n',n,'material','nitinol','plotkin',false,'verbose',false);
wrist = Wrist(ro*2,ri*2,n,h*ones(n,1),zeros(n,1),c*ones(n,1),g*ones(n,1),'CutType','off-axis'); % Nick's wrist class
wrist.theta = zeros(n,1);
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
wrist.ybar-sheath.ybar'
wrist.I-sheath.J'
while theta_diff > 1E-6 && trials < maxTrials
    for ii = 1:n
        idx1 = ii;
        
        E_vec(ii) = wrist.E_lin;
        % Compute distal force and moment
        F_vec(ii) = F*exp(-mu*sum(wrist.theta(1:ii))); % get force experienced by notch
        M_vec(ii) = F_vec(ii)*(wrist.ybar(ii) + wrist.ID/2); % get moment experienced by notch
        
        sheath.E_eff(idx1)=40e9;
        sheath.F(idx1) = F*exp(-sheath.mu*sum(sheath.theta(1:idx1)));
        sheath.M(idx1) = sheath.F(idx1)*(sheath.ybar(idx1)+sheath.r_i);
        
        pct = 100;
        k = 1;
        % Gradient descent for non-linear modulus
        while k<100 && pct>1E-4
            % Compute angular deflection of current segment
            wrist.theta(ii) = M_vec(ii)*wrist.h(ii)/(E_vec(ii)*wrist.I(ii));
            sheath.theta(idx1) = sheath.M(idx1)*sheath.h(idx1)/(sheath.E_eff(idx1)*sheath.J(idx1));
            if norm(wrist.theta-sheath.theta) > eps
                error("Error with theta")
            end

            wrist.check_notch_limits(ii);
            sheath.CheckLimits(idx1);
            
            
            % Compute section arc length, curvature and
            % strain
            s1 = wrist.h(ii)-(wrist.ybar(ii))*abs(wrist.theta(ii));
            wrist.kappa(ii) = wrist.theta(ii)/s1;
            epsilon = (wrist.kappa(ii).*(wrist.OD/2-wrist.ybar(ii)))./(1+wrist.ybar(ii)*wrist.kappa(ii));
            
            s1_j = sheath.h(idx1)-(sheath.ybar(idx1))*abs(sheath.theta(idx1));
            kappa_j = sheath.theta(idx1)/s1_j;
            epsilon_j = (kappa_j.*(sheath.r_o-sheath.ybar(idx1)))./(1+sheath.ybar(idx1)*kappa_j);
            if norm(epsilon-epsilon_j) > eps
                error("Error with epsilon")
            end
            
            % Update modulus via gradient descent
            [stress_eff,eta] = wrist.get_stress(abs(epsilon));
            new_E = E_vec(ii)-eta*(E_vec(ii)-stress_eff/abs(epsilon));
            
            [stress_eff_j,eta_j] = sheath.SuperElastic(abs(epsilon_j));
            new_E_j = sheath.E_eff(idx1)-eta_j*(sheath.E_eff(idx1)-stress_eff_j/abs(epsilon_j));
            if norm(new_E-new_E_j) > eps
                error("Error with E")
            end
            % Percent change in modulus (for convergence
            % checking)
            pct = abs((new_E-E_vec(ii))/E_vec(ii));
            pct = abs((new_E_j-sheath.E_eff(idx1))/sheath.E_eff(idx1));
            % Update modulus guess
            E_vec(ii) = new_E;
            sheath.E_eff(idx1) = new_E_j;
            
            % Increment
            k = k+1;
        end
        
    end
    theta_vec = theta_vec + theta_vec_base; % add offset for precurvature
    theta_delta = sum(theta_vec)-sum(theta_last);
    trials = trials + 1; % increment
end
