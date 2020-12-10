function [stress,eta] = GetStress(strain, E_lin, E_se)
% GETSTRESS - returns the stress of a single notch based on stress strain
% curve of nitinol
strain_lower = 0.02;
strain_upper = 0.1;

sigma = @(e) (e<strain_lower).*e*E_lin+...
    (e >= strain_lower && e < strain_upper)*((e-strain_lower)*.08*E_lin+strain_lower*E_lin)+...
    (e >= strain_upper)*((1.0E9)*exp(-.01/(e-strain_upper))+strain_lower*1*E_lin+(strain_upper-strain_lower)*.08*E_lin);
stress = abs(sigma(strain));
eta_fun = @(e) (e<strain_lower).*.5+...
    (e >= strain_lower && e < strain_upper)*1+...
    (e >= strain_upper)*0.3;
eta = eta_fun(strain);
end