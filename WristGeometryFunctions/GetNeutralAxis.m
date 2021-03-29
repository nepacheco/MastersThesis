
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
        alpha_in = phi_i/2;
        alpha_out = phi_o/2;
        
        % Equations from Swaney and York
        A_o = (ro^2)*(phi_o-sin(phi_o))/2;
        A_seg_out=(ro^2)*(alpha_out-sin(alpha_out)*cos(alpha_out));
        A_o_diff = A_o - A_seg_out;
        
        A_i = (ri^2)*(phi_i-sin(phi_i))/2;
        A_seg_in=(ri^2)*(alpha_in-sin(alpha_in)*cos(alpha_in));
        A_idiff = A_i - A_seg_in;
        
        ybar_o = 4*ro*sin(1/2*phi_o)^3/(3*(phi_o - sin(phi_o)));
        ybar_out=(2*ro/3)*(sin(alpha_out)^3)/(alpha_out-sin(alpha_out)*cos(alpha_out));
        ybar_o_diff = ybar_o - ybar_out;
        
        ybar_i = 4*ri*sin(1/2*phi_i)^3/(3*(phi_i - sin(phi_i)));
        ybar_in=(2*ri/3)*(sin(alpha_in)^3)/(alpha_in-sin(alpha_in)*cos(alpha_in));
        ybar_i_diff = ybar_i - ybar_in;
        
%         ybar = (ybar_o*A_o - ybar_i*A_i)/(A_o - A_i);
        A_seg=A_seg_out-A_seg_in;
        
        
        % using Circular segment for outer and inner regions of tube
        I_o = (phi_o - sin(phi_o) + 2*sin(phi_o)*sin(phi_o/2)^2)*ro^4/8;
        Jz_oo=(ro^4)/4*(alpha_out-sin(alpha_out)*cos(alpha_out)+2*sin(alpha_out)^3*cos(alpha_out));
        I_o_diff = I_o - Jz_oo;
        
        I_i = (phi_i - sin(phi_i) + 2*sin(phi_i)*sin(phi_i/2)^2)*ri^4/8;
        Jz_oi=(ri^4)/4*(alpha_in-sin(alpha_in)*cos(alpha_in)+2*sin(alpha_in)^3*cos(alpha_in));
        I_i_diff = I_i - Jz_oi;
        % subtract inner from outer to determine area moment of inertia for cut
        % section only
        Io = I_o - I_i;
        Io_diff = Io - (Jz_oo - Jz_oi);
        Jz_o=Jz_oo-Jz_oi;
        % Use parallel axis theorem to shift the area moment of inertia to centroid
        % of the notch
        ybar = (ybar_out*A_seg_out-ybar_in*A_seg_in)/A_seg;
%         I = Io - (A_o - A_i)*ybar^2;
        I = Jz_o-A_seg*ybar^2;
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