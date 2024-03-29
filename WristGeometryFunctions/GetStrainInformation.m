
function [strain, stress, E] = GetStrainInformation(theta, h, ro, ybar,opts)
% GETSTRAININFORMATION - returns the stress, strain, and effective youngs
% modulus for a notch
arguments
    theta (1,1) double  % [rad] bending angle of notch
    h (1,1) double  % [m] notch height
    ro (1,1) double  % [m] outer radius of tube
    ybar (1,1) double  % [m] neutral bending plane of tube
    opts.E_lin (1,1) double = 10e9 % [Pa] Linear Elastic Modulus
    opts.E_se (1,1) double  = 3e9% [Pa]  Super elastic modulus
    opts.strainLower (1,1) double = 0.028 % lower strain threshold
end

kappa = theta/(h-theta*ybar);
y = max([ro-ybar]);
strain = abs(kappa*y/(1 + ybar*kappa));
stress = GetStress(strain, opts.E_lin, opts.E_se, 'strainLower', opts.strainLower);
if (strain > 0)
    E = stress/strain; % Assume we are in linear elastic range
else
    E = E_lin; % This is linear modulus of nitinol
end
end