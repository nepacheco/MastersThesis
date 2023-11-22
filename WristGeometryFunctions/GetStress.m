function [stress,eta] = GetStress(strain, E_lin, E_se, opts)
% GETSTRESS - returns the stress of a single notch based on stress strain
% curve of nitinol
arguments
    strain (1,1) double  % approximate strain experienced by notch
    E_lin  (1,1) double  % the estimated linear elastic modulus of the material
    E_se   (1,1) double  % The estimated super elastic modulus of the material
    opts.strainLower  (1,1) double = 0.028  % The strain threshold where material becomes super elastic
    opts.strainUpper (1,1) double = 0.1  % the strain threshold where material becomes plastic. 
end

strain_lower = opts.strainLower;
strain_upper = opts.strainUpper;

sigma = @(e) (e<strain_lower).*e*E_lin+...
    (e >= strain_lower && e < strain_upper)*((e-strain_lower)*E_se+strain_lower*E_lin)+...
    (e >= strain_upper)*((1.0E9)*exp(-.01/(e-strain_upper))+strain_lower*1*E_lin+(strain_upper-strain_lower)*E_se);
stress = abs(sigma(strain));
eta_fun = @(e) (e<strain_lower).*.5+...
    (e >= strain_lower && e < strain_upper)*1+...
    (e >= strain_upper)*0.3;
eta = eta_fun(strain);
end