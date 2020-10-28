
function [ybar, I] = GetNeutralAxis(ro, ri, g)
%GETNEUTRALAXIS - returns the distance to neutral bending axis and area
%moment of inertia about bending axis

% Assumes we are cutting between center of tube and inner radius
phi_o = 2*acos((g - ro)/ro);
phi_i = 2*acos((g - ro)/ri);

% Equations from Swaney and York
A_o = (ro^2)*(phi_o-sin(phi_o))/2;
A_i = (ri^2)*(phi_i-sin(phi_i))/2;
ybar_o = 4*ro*sin(1/2*phi_o)^3/(3*(phi_o - sin(phi_o)));
ybar_i = 4*ri*sin(1/2*phi_i)^3/(3*(phi_i - sin(phi_i)));
ybar = (ybar_o*A_o - ybar_i*A_i)/(A_o - A_i);

% using Circular segment for outer and inner regions of tube
I_o = (phi_o - sin(phi_o) + 2*sin(phi_o)*sin(phi_o/2)^2)*ro^4/8;
I_i = (phi_i - sin(phi_i) + 2*sin(phi_i)*sin(phi_i/2)^2)*ri^4/8;
% subtract inner from outer to determine area moment of inertia for cut
% section only
Io = I_o - I_i;
% Use parallel axis theorem to shift the area moment of inertia to centroid
% of the notch
I = Io - (A_o - A_i)*ybar^2;
end