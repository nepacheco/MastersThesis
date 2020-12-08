classdef Wrist < handle
    %WRIST Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        % Tube Geometry properties
        OD = 1.62E-3;               % outer diameter of tube
        ID = 1.4E-3;               % inner diameter of tube
        n = 6;                      % number of notches
        h = 0.9E-3.*ones(6,1);        % vector containing notch base width parameters
        phi = 0.*ones(6,1);
        c = 1.0E-3.*ones(6,1);        % vector containing notch spacing parameters
        g = 1.45E-3.*ones(6,1);         % vector containing tube notch height
        E_lin = 40E9;
        E_se = 0.08*40E9;
        
        % Kinematic properties
        transformation = []; % contains every intermediate transform
        T_tip = zeros(4,4); % Transform from base to tip
        kappa = zeros(6,1); % [m^-1] Curvature for each notch
        s = zeros(6,1); % [m] Arc length for each notch
        theta = zeros(6,1); % [rad] bending angle for each notch
        alpha = 0 % [rad] Base rotation
        tau = 0; % [m] Base translation
        force = 0; % [N] Force on the wire
        delta_l = 0; % [m] Tendon displacement
        ybar = zeros(6,1);
        I = zeros(6,1);
        
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
        
        function [path,T_tip,obj] = fwkin(obj,q,varargin)
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
            
            obj.alpha = q(2);
            obj.tau = q(3);
            switch(type)
                case 'force'
                    obj.force = q(1);
                    obj.get_force_arc_params();
                case 'geometry'
                    disp("geometry")
                    obj.delta_l = q(1);
                    obj.get_geom_arc_params(q(1));
            end
            [path,T_tip] = obj.robot_kin();
        end
        
        
        function [path,T,obj] = robot_kin(obj)
            %ROBOT_KIN - returns the forward kinematics of a notched wrist
            disp("Running robot Kin")
            % Initial frame and translation section
            obj.transformation = zeros(4,4,2*(obj.n+1));
            % Plotting multiple frames from multiple sections of tube
            obj.T_tip = obj.get_arc_fwdkin(0,0,0);
            obj.transformation(:,:,1) = obj.get_arc_fwdkin(0,0,0);
            % Add section of un-notched tube
            T_straight = obj.get_arc_fwdkin(0, obj.alpha, obj.tau);
            obj.transformation(:,:,2) = T_straight;
            obj.T_tip = obj.T_tip*T_straight;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % getting position of notches
            % For loop to go through the c, h, and g, matrices to create
            for i = 1:obj.n
                T_arc = obj.get_arc_fwdkin(obj.kappa(i), 0, obj.s(i));
                obj.transformation(:,:,2*i + 1) = T_arc;
                obj.T_tip = obj.T_tip*T_arc;
                T_straight = obj.get_arc_fwdkin(0,obj.phi(i),obj.c(i));
                obj.transformation(:,:,2*i+2) = T_straight;
                obj.T_tip = obj.T_tip*T_straight;
            end
            path = obj.transformation;
            T = obj.T_tip;
        end
        
        function [transform,obj] = get_arc_fwdkin(obj,kappa,phi,arc_length)
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
        
        function [ybar_, I_,obj] = get_neutral_axis(obj,varargin)
            %GETNEUTRALAXIS - returns the distance to neutral bending axis and area
            %moment of inertia about bending axis
            
            %****** INPUT PARSING *********************
            % default values
            cutType = 'off-axis';
            cutOptions = {'off-axis','on-axis'};
            
            p = inputParser();
            addRequired(p,'obj');
            addParameter(p,'CutType',cutType,@(x) any(validatestring(x,cutOptions)));
            parse(p,obj,varargin{:});
            
            cutType = p.Results.CutType;
            %*********************************************
            obj.ybar = zeros(obj.n,1);
            obj.I = zeros(obj.n,1);
            for index = 1:obj.n
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
                        obj.ybar(index) = (ybar_o*A_o - ybar_i*A_i)/(A_o - A_i);
                        
                        % using Circular segment for outer and inner regions of tube
                        I_o = (phi_o - sin(phi_o) + 2*sin(phi_o)*sin(phi_o/2)^2)*(obj.OD/2)^4/8;
                        I_i = (phi_i - sin(phi_i) + 2*sin(phi_i)*sin(phi_i/2)^2)*(obj.ID/2)^4/8;
                        % subtract inner from outer to determine area moment of inertia for cut
                        % section only
                        Io = I_o - I_i;
                        % Use parallel axis theorem to shift the area moment of inertia to centroid
                        % of the notch
                        obj.I(index) = Io - (A_o - A_i)*obj.ybar(index)^2;
                    case 'on-axis'
                        phi_o = 2*acos((obj.g(index) - (obj.OD/2))/(obj.OD/2));
                        
                        ybar_o = 2*(obj.OD/2)*sin(phi_o/2)/(3*phi_o/2);
                        ybar_i = 2*(obj.ID/2)*sin(phi_o/2)/(3*phi_o/2);
                        Ao = phi_o/2*(obj.OD/2)^2;
                        Ai = phi_o/2*(obj.ID/2)^2;
                        A = Ao - Ai;
                        obj.ybar(index) = (Ao*ybar_o - Ai*ybar_i)/A;
                        
                        % Using Circular Sector for outer and inner regions of tube
                        I_o = (obj.OD/2)^4/4 * (phi_o/2 + 1/2*sin(phi_o));
                        I_i = ri^4/4*(phi_o/2 + 1/2*sin(phi_o));
                        % Subtract inner from outer to determine area moment of inertia for
                        % the remaining backbone
                        Iprime = I_o - I_i;
                        % Parallel axis theorem to move the area moment of inertia to
                        % centroid of notch
                        obj.I(index) = Iprime - A*obj.ybar(index)^2;
                end
            end
            ybar_ = obj.ybar;
            I_ = obj.I;
        end
        
        function [k, s,obj] = get_geom_arc_params(obj, l)
            %GET_ARC_PARAMS - gets the curvature and arc length of notched
            %wrist
            %   Given a specific tendon displacement this function
            %   will determine the curvature and arc length of each notch
            %   notched wrist.
            % assumes the l given is the total tendon displacement
            k = zeros(obj.n,1);
            s = zeros(obj.n,1);
            l = l/obj.n;
            obj.get_neutral_axis();
            for i = 1:obj.n
                k(i) = l/(obj.h(i)*(obj.ID/2 + obj.ybar(i)) - l*obj.ybar(i));
                s(i) = obj.h(i)/(1 + obj.ybar(i)*k(i));
            end
            obj.kappa = k;
            obj.s = s;
        end
        
        function obj = get_force_arc_params(obj)
            theta_vec_base = deg2rad(0)*ones(obj.n,1); % this represents the precurvature in each notch
            
            % INITIALIZATION OF VECTORS %
            obj.theta = theta_vec_base; % vector of each notch angle which should start at precurve value
            F_vec = zeros(obj.n,1); % vector of force experienced by each notch
            M_vec = zeros(obj.n,1); % vector of moment experienced by each notch
            E_vec = obj.E_lin*ones(obj.n,1); % effective elastic modulus for each notch
            mu = 0.2;
            
            theta_diff = 100; % start the change in theta high
            trials = 0; % trials we have completed
            maxTrials = 10; % maximum number of trials to perform
            while theta_diff > 1E-6 && trials < maxTrials
                for ii = 1:obj.n
                    % Compute distal force and moment
                    F_vec(ii) = F*exp(-mu*sum(obj.theta(1:ii))); % get force experienced by notch
                    M_vec(ii) = F_vec(ii)*(obj.ybar(ii) + obj.ID/2); % get moment experienced by notch
                    obj.theta(ii) = M_vec(ii)*obj.h(ii)/(E_vec(ii)*obj.I(ii)); % update theta
                    
                    pct = 100;
                    k = 1;
                    % Gradient descent for non-linear modulus
                    while k<100 && pct>1E-4
                        
                        % Compute angular deflection of current segment
                        obj.theta(ii) = M_vec(ii)*obj.h(ii)/(E_vec(ii)*obj.I(ii));
                        
                        % Compute section arc length, curvature and
                        % strain
                        s1 = (obj.h(ii)-obj.ybar(ii))*abs(obj.theta(ii));
                        obj.kappa(ii) = obj.theta(ii)/s1;
                        epsilon = (obj.kappa(ii).*(obj.OD/2-obj.ybar(ii)))./(1+obj.ybar(ii)*obj.kappa(ii));
                        
                        % Update modulus via gradient descent
                        [stress_eff,eta] = obj.get_stress(abs(epsilon));
                        new_E = E_vec(ii)-eta*(E_vec(ii)-stress_eff/abs(epsilon));
                        
                        % Percent change in modulus (for convergence
                        % checking)
                        pct = abs((new_E-E_vec(ii))/E_vec(ii));
                        
                        % Update modulus guess
                        E_vec(ii) = new_E;
                        
                        % Increment
                        k = k+1;
                    end
                    
                end
                obj.theta = obj.theta + theta_vec_base; % add offset for precurvature
                trials = trials + 1; % increment
            end
        end
        
        function maxStrain = calcMaxStrain(obj, q)
            %CALCMAXSTRAIN - takes as input a vector q of actuator variables
            %and calculates the maximum strain underwent by the material
            obj.delta_l = q(1);
            obj.get_neutral_axis();
            y = max([obj.OD/2 - obj.ybar(1),obj.w(1)-obj.OD/2-obj.ybar(1)]);
            obj.get_geom_arc_params(obj.delta_l);
            maxStrain = abs(obj.k(1)*y/(1 + obj.ybar(1)*obj.k(1)));
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
        
        function obj = plot_stick_model(obj)
            x_vec = [];
            y_vec = [];
            z_vec = [];
            T = eye(4,4);
            for i = 1:size(obj.transformation,3)
                T = T*obj.transformation(:,:,i);
                x_vec = [x_vec T(1,4)];
                y_vec = [y_vec T(2,4)];
                z_vec = [z_vec T(3,4)];
            end
            plot3(x_vec,y_vec,z_vec,'Linewidth',2.5);
            axis equal
            grid on;
        end
    end
end

