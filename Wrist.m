classdef Wrist
    %WRIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        OD = 1.62E-3;               % outer diameter of tube
        ID = 1.4E-3;               % inner diameter of tube
        n = 6;                      % number of notches
        h = 0.9E-3.*ones(6,1);        % vector containing notch base width parameters
        phi = 0.*ones(6,1);
        c = 1.0E-3.*ones(6,1);        % vector containing notch spacing parameters
        g = 1.45E-3.*ones(6,1);         % vector containing tube notch height
        E_lin = 40E9;
        E_se = 0.08*40E9;
        transformations = [];
        T_tip = zeros(4,4);
    end
    
    methods
        function obj = Wrist(OD,ID, n, h, phi, x, w)
            %WRIST Construct an instance of this class
            %   Detailed explanation goes here
            obj.OD = OD;
            obj.ID = ID;
            obj.n = n;
            obj.h = h;
            obj.phi = phi;
            obj.c = x;
            obj.g = w;
        end
        
        function [path,T_tip] = fwkin(obj,q,varargin)
            %FWKIN - returns the forward kinematics of a notched wrist
            %depending on the input paramters q and whether we are using
            %force or geometry model to determine arc paratmers
            %   q - a (3x1) vector containing tendon displace or force
            %   (depending on type), tube rotation, and tube translation.
            %   type - an name-value pair for the type of model we want to
            %   use (force, geometry)
            
            %****** INPUT PARSING *********************
            % default values
            type = 'geometry';
            typeOptions = {'force','geometry'};
            
            p = inputParser();
            addRequired(p,'obj');
            addRequired(p,'q',@(x) isnumeric(x));
            addParameter(p,'Type',type,@(x) any(validatestring(x,typeOptions)));
            parse(p,obj,q,varargin{:});
            
            type = p.Results.Type;
            %*********************************************
            
            switch(type)
                case 'force'
                case 'geometry'
                    [k,s] = obj.get_geom_arc_params(q(1));
            end
            [path,T_tip] = obj.robot_kin([k,s],q(2),q(3));
        end
        
        
        function [path,T_tip] = robot_kin(obj, arc_parameters,alpha,d)
            %ROBOT_KIN - returns the forward kinematics of a notched wrist
            %   arc parameters is a (n,2) matrix containing the curvature
            %       and arc length for each notched section
            %   alpha - represents base rotation
            %   d - base translation
            k = arc_parameters(:,1);
            s = arc_parameters(:,2);
            
            % Initial frame and translation section
            transformation = zeros(4,4,2*(obj.n+1));
            % Plotting multiple frames from multiple sections of tube
            obj.T_tip = obj.get_arc_fwdkin(0,0,0);
            obj.transformation(:,:,1) = obj.get_arc_fwdkin(0,0,0);
            % Add section of un-notched tube
            T_straight = obj.get_arc_fwdkin(0, alpha, d);
            obj.transformation(:,:,2) = T_straight;
            obj.T_tip = obj.T_tip*T_straight;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % getting position of notches
            % For loop to go through the c, h, and g, matrices to create
            for i = 1:obj.n
                T_arc = obj.get_arc_fwdkin(k(i), 0, s(i));
                transformation(:,:,2*i + 1) = T_arc;
                obj.T_tip = obj.T_tip*T_arc;
                T_straight = obj.get_arc_fwdkin(0,obj.phi(i),obj.c(i));
                transformation(:,:,2*i+2) = T_straight;
                obj.T_tip = T_tip*T_straight;
            end
            path = transformation;
        end
        
        function transform = get_arc_fwdkin(obj,kappa,phi,arc_length)
            %GET_ARC_FWDKIN - Define transformation matrices based on arc parameters
            % theta - rotation about the base
            % k - 1/r where r is the radius of curvature
            % s - arc length of the constant curvature section
            
            % rot is the twist rotation matrix as part of
            % creating the transformation matrix between
            % sections of the tube
            % phi - (radians)the base angle change around the z-axis between sections
            rot = @(theta) theta*[0 -1 0 0;...
                1  0 0 0;...
                0  0 0 0;...
                0  0 0 0];
            
            % inp is the linear translation described between
            % the two points of a tube
            % k - 1/r (mm) which is the radius of curvatature of the section
            % l - (mm) the arclength of the section of tube
            inp = @(k,s) s*[0 0 k 0;...
                0 0 0 0;...
                -k 0 0 1;...
                0 0 0 0];  % Webster says this matrix should be transposed
            % but it seems to produce the wrong matrix
            % when transposed
            
            T = @(k,theta,s) expm(rot(theta))*expm(inp(k,s));
            transform = T(kappa, phi,arc_length);
        end
        
        function [ybar, I] = get_neutral_axis(obj,index,varargin)
            %GETNEUTRALAXIS - returns the distance to neutral bending axis and area
            %moment of inertia about bending axis
            
            %****** INPUT PARSING *********************
            % default values
            cutType = 'off-axis';
            cutOptions = {'off-axis','on-axis'};
            
            p = inputParser();
            addRequired(p,'obj');
            addRequired(p,'index',@(x) isnumeric(x));
            addParameter(p,'CutType',cutType,@(x) any(validatestring(x,cutOptions)));
            parse(p,obj,index,varargin{:});
            
            cutType = p.Results.CutType;
            %*********************************************
            
            switch(cutType)
                case 'off-axis'
                    % Assumes we are cutting between center of tube and inner radius
                    phi_o = 2*acos((obj.g(index) - obj.OD/2)/(obj.OD/2));
                    phi_i = 2*acos((obj.g(index) - obj.OD/2)/(obj.ID/2));
                    
                    % Equations from Swaney and York
                    A_o = ((obj.OD/2)^2)*(phi_o-sin(phi_o))/2;
                    A_i = ((obj.ID/2)^2)*(phi_i-sin(phi_i))/2;
                    ybar_o = 4*(obj.OD/2)*sin(1/2*phi_o)^3/(3*(phi_o - sin(phi_o)));
                    ybar_i = 4*(obj.ID/2)*sin(1/2*phi_i)^3/(3*(phi_i - sin(phi_i)));
                    ybar = (ybar_o*A_o - ybar_i*A_i)/(A_o - A_i);
                    (A_o - A_i);
                    
                    % using Circular segment for outer and inner regions of tube
                    I_o = (phi_o - sin(phi_o) + 2*sin(phi_o)*sin(phi_o/2)^2)*(obj.OD/2)^4/8;
                    I_i = (phi_i - sin(phi_i) + 2*sin(phi_i)*sin(phi_i/2)^2)*(obj.ID/2)^4/8;
                    % subtract inner from outer to determine area moment of inertia for cut
                    % section only
                    Io = I_o - I_i;
                    % Use parallel axis theorem to shift the area moment of inertia to centroid
                    % of the notch
                    I = Io - (A_o - A_i)*ybar^2;
                case 'on-axis'
                    phi_o = 2*acos((obj.g(index) - (obj.OD/2))/ri);
                    
                    ybar_o = 2*(obj.OD/2)*sin(phi_o/2)/(3*phi_o/2);
                    ybar_i = 2*(obj.ID/2)*sin(phi_o/2)/(3*phi_o/2);
                    Ao = phi_o/2*(obj.OD/2)^2;
                    Ai = phi_o/2*(obj.ID/2)^2;
                    A = Ao - Ai;
                    ybar = (Ao*ybar_o - Ai*ybar_i)/A;
                    
                    % Using Circular Sector for outer and inner regions of tube
                    I_o = (obj.OD/2)^4/4 * (phi_o/2 + 1/2*sin(phi_o));
                    I_i = ri^4/4*(phi_o/2 + 1/2*sin(phi_o));
                    % Subtract inner from outer to determine area moment of inertia for
                    % the remaining backbone
                    Iprime = I_o - I_i;
                    % Parallel axis theorem to move the area moment of inertia to
                    % centroid of notch
                    I = Iprime - A*ybar^2;
            end
        end
        
        function [k, s] = get_geom_arc_params(obj, l)
            %GET_ARC_PARAMS - gets the curvature and arc length of notched
            %wrist
            %   Given a specific tendon displacement this function
            %   will determine the curvature and arc length of each notch
            %   notched wrist.
            % assumes the l given is the total tendon displacement
            k = zeros(1,obj.n);
            s = zeros(1,obj.n);
            l = l/obj.n;
            for i = 1:obj.n
                y_bar = obj.get_neutral_axis(i);
                k(i) = l/(obj.h(i)*(obj.ID/2 + y_bar) - l*y_bar);
                s(i) = obj.h(i)/(1 + y_bar*k(i));
            end
        end
        
        function maxStrain = calcMaxStrain(obj, q)
            %CALCMAXSTRAIN - takes as input a vector q of actuator variables
            %and calculates the maximum strain underwent by the material
            delta_l = q(1);
            ybar = obj.get_neutral_axis(1);
            y = max([obj.OD/2 - ybar,obj.w(1)-obj.OD/2-ybar]);
            [k, s] = obj.get_geom_arc_params(delta_l);
            maxStrain = abs(k(1)*y/(1 + ybar*k(1)));
            if maxStrain > 0.08
                warning("Exceeding maximum recoverable strain: %f",...
                    maxStrain);
            end
        end
        
        function [stress,eta] = get_stress(strain)
            % GETSTRESS - returns the stress of a single notch based on stress strain
            % curve of nitinol
            strain_lower = 0.02;
            strain_upper = 0.1;
            
            sigma = @(e) (e < strain_lower).*obj.E_lin*e + ...
                (e >= strain_lower && e < strain_upper).*(obj.E_se*(e - strain_lower) + strain_lower*obj.E_lin)+...
                (e >= strain_upper)*((1.0E9)*exp(-.01/(e-strain_upper))+strain_upper*1*obj.E_lin+(strain_upper-strain_lower)*obj.E_se);
            stress = abs(sigma(strain));
            eta_fun = @(e) (e<strain_lower).*.5+...
                (e >= strain_lower && e < strain_upper)*1+...
                (e >= strain_upper)*0.3;
            eta = eta_fun(strain);
        end
        
        function maxL = calcMaxL(obj)
            %CALCMAXL - determines the maximum tendon displacement before
            %8 percent strain is experienced by the notch
            ybar = obj.get_neutral_axis(1);
            maxStrain = 0.08;
            kappa = maxStrain/(obj.OD/2 - ybar*(1 + maxStrain));
            maxL = (kappa*obj.h(1)*(obj.ID/2 + ybar)/(1 + ybar*kappa));
        end
    end
end

