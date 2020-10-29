
function [strain, stress, E] = GetStrainInformation(theta, h, ro, ybar)
% GETSTRAININFORMATION - returns the stress, strain, and effective youngs
% modulus for a notch
kappa = theta/(h-theta*ybar);
y = max([ro-ybar]);
strain = abs(kappa*y/(1 + ybar*kappa));
stress = GetStress(strain, 40E9);
E = stress/strain; % Assume we are in linear elastic range
end