
function [ybar, I] = GetNeutralAxis(ro, ri, g, varargin)
%GETNEUTRALAXIS - returns the distance to neutral bending axis and area
%moment of inertia about bending axis

%****** INPUT PARSING *********************
% default values
cutType = 'off-axis';
cutOptions = {'off-axis','on-axis'};

p = inputParser();
addRequired(p,'outerRadius',@isnumeric);
addRequired(p,'innerRadius',@isnumeric);
addRequired(p, 'cutDepth',@isnumeric);
addOptional(p,'CutType',cutType,@(x) any(validatestring(x,cutOptions)));
parse(p,ro,ri,g,varargin{:});

cutType = p.Results.CutType;
%*********************************************

switch(cutType)
    case 'off-axis'
        % Assumes we are cutting between center of tube and inner radius
        phi_o = 2*acos((g - ro)/ro);
        phi_i = 2*acos((g - ro)/ri);
        
        % Equations from Swaney and York
        A_o = (ro^2)*(phi_o-sin(phi_o))/2;
        A_i = (ri^2)*(phi_i-sin(phi_i))/2;
        ybar_o = 4*ro*sin(1/2*phi_o)^3/(3*(phi_o - sin(phi_o)));
        ybar_i = 4*ri*sin(1/2*phi_i)^3/(3*(phi_i - sin(phi_i)));
        ybar = (ybar_o*A_o - ybar_i*A_i)/(A_o - A_i);
        (A_o - A_i);
        
        % using Circular segment for outer and inner regions of tube
        I_o = (phi_o - sin(phi_o) + 2*sin(phi_o)*sin(phi_o/2)^2)*ro^4/8;
        I_i = (phi_i - sin(phi_i) + 2*sin(phi_i)*sin(phi_i/2)^2)*ri^4/8;
        % subtract inner from outer to determine area moment of inertia for cut
        % section only
        Io = I_o - I_i;
        % Use parallel axis theorem to shift the area moment of inertia to centroid
        % of the notch
        I = Io - (A_o - A_i)*ybar^2;
    case 'on-axis'
        phi = 2*acos((g - ro)/ro);

        ybar_o = 2*ro*sin(phi/2)/(3*phi/2);
        ybar_i = 2*ri*sin(phi/2)/(3*phi/2);
        Ao = phi/2*ro^2;
        Ai = phi/2*ri^2;
        A = Ao - Ai;
        ybar = (Ao*ybar_o - Ai*ybar_i)/A;
        
        % Using Circular Sector for outer and inner regions of tube
        I_o = ro^4/4 * (phi/2 + 1/2*sin(phi));
        I_i = ri^4/4*(phi/2 + 1/2*sin(phi));
        % Subtract inner from outer to determine area moment of inertia for
        % the remaining backbone
        Iprime = I_o - I_i;
        % Parallel axis theorem to move the area moment of inertia to
        % centroid of notch
        I = Iprime - A*ybar^2;
end
end