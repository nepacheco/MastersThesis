function stress = GetStress(strain, E_lin)
% GETSTRESS - returns the stress of a single notch based on stress strain
% curve of nitinol
strain_lower = 0.02;
if (strain > 0.1) 
    disp("very high");
end
sigma = @(e) (e < strain_lower).*E_lin*e + ...
    (strain_lower <= e).*(E_lin*(0.08)*(e - strain_lower) + strain_lower*E_lin);
stress = abs(sigma(strain));
end