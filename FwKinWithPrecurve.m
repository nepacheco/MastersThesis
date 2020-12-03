clear; close all;
% TUBE PARAMETERS %
n = 5; % number of notches
ro = (1.62E-3)/2; % [mm] outer radius
ri = (1.4E-3)/2; % [mm] inner radius
h = 0.7915E-3; % [mm] notch height
g = 1.4E-3; % [mm] notch depth
c = 2.2085E-3; % [mm] uncut section height
E = 40E9; % [Pa] linear elastic modulus
mu = 0.2; % coefficient of friction
F = 1; % [N] input Force
theta_vec_base = deg2rad(2)*ones(n,1); % this represents the precurvature in each notch

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
while theta_diff > 1E-6 && trials < maxTrials
    for ii = 1:n    
        % Compute distal force and moment
        F_vec(ii) = F*exp(-mu*sum(theta_vec(1:ii))); % get force experienced by notch
        M_vec(ii) = F_vec(ii)*(ybar +ri); % get moment experienced by notch
        theta_vec(ii) = M_vec(ii)*h/(E_vec(ii)*I); % update theta
        
        pct = 100;
        k = 1;
        % Gradient descent for non-linear modulus
        while k<100 && pct>1E-4
            
            % Compute angular deflection of current segment
            theta_vec(ii) = M_vec(ii)*h/(E_vec(ii)*I);
            
            % Compute section arc length, curvature and
            % strain
            s1 = (h-ybar)*abs(theta_vec(ii));
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
    trials = trials + 1; % increment
end
rad2deg(theta_vec)