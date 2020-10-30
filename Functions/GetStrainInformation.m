
function [strain, stress, E] = GetStrainInformation(theta, h, ro, ybar,E_lin)
% GETSTRAININFORMATION - returns the stress, strain, and effective youngs
% modulus for a notch
kappa = theta/(h-theta*ybar);
y = max([ro-ybar]);
strain = abs(kappa*y/(1 + ybar*kappa));
stress = GetStress(strain, E_lin);
if (strain > 0)
    E = stress/strain; % Assume we are in linear elastic range
else
    E = E_lin; % This is linear modulus of nitinol
end
end