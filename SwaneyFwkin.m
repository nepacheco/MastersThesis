% TUBE PARAMETERS %
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

[ybar,I] = GetNeutralAxis(ro,ri,g); % neutral axis information (uniform notches)

kappa = delta_l/(h*(ri + ybar) - delta_l*ybar);
s = h/(1 + ybar*kappa);

theta = s*kappa;
rad2deg(theta)

theta_max = h/(ro + ybar);
rad2deg(theta_max);
