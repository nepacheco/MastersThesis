clear; close all;
% TUBE PARAMETERS %
n = 5;
ro = (1.62E-3)/2; % [mm]
ri = (1.4E-3)/2; % [mm]
h = 0.7915E-3; % [mm]
g = 1.4E-3; % [mm]
c = 2.2085E-3; % [mm]
E = 40E9; % [Pa]
mu = 0.2;
F = 1; % [N] input Force
theta_vec_base = deg2rad(2)*ones(n,1);

% INITIALIZATION OF VECTORS %
theta_vec = zeros(n,1);
F_vec = zeros(n,1);
M_vec = zeros(n,1);
E_vec = E*ones(n,1);
[ybar,I] = GetNeutralAxis(ro,ri,g);

% INTIALIZATION OF GRADIENT DESCENT %
theta_diff = 100;
trials = 0;
maxTrials = 10;
while theta_diff > 1E-6 && trials <= maxTrials
    for ii = 1:n
        E_vec(ii)=E;
        
        % Compute distal force and moment
        F_vec(ii) = F*exp(-mu*sum(theta_vec(1:ii)));
        M_vec(ii) = F_vec(ii)*(ybar +ri);
        theta_vec(ii) = M_vec(ii)*h/(E_vec(ii)*I);
        
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
    theta_vec = theta_vec + theta_vec_base;
    trials = trials + 1;
end

rad2deg(theta_vec)