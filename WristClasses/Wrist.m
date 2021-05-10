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
        E_se = 0.08*(40E9);
        strain_lower = 0.02;
        strain_upper = 0.1;
        mu = 0.2;
        
        % Kinematic properties
        transformation = []; % contains every intermediate transform
        T_tip = zeros(4,4); % Transform from base to tip
        kappa = zeros(6,1); % [m^-1] Curvature for each notch
        s = zeros(6,1); % [m] Arc length for each notch
        theta = zeros(6,1); % [rad] bending angle for each notch
        precurve_theta = zeros(6,1);
        alpha = 0 % [rad] Base rotation
        tau = 0; % [m] Base translation
        F = 0; % [N] Force on the wire
        delta_l = 0; % [m] Tendon displacement
        ybar = zeros(6,1);
        I = zeros(6,1);
        use_friction = true;
        use_non_linear = true;
        use_precurve = false;
        name = "wrist";
        
        DEBUG = false;
        cutType = 'on-axis'
        
    end
    
    methods
        function obj = Wrist(OD,ID, n, h, phi, c, g, varargin)
            %WRIST Construct an instance of this class
            %   Detailed explanation goes here
            
            %****** INPUT PARSING *********************
            % default values
            cutType = 'off-axis';
            cutOptions = {'off-axis','on-axis'};
            
            p = inputParser();
            addRequired(p,'OD');
            addRequired(p,'ID');
            addRequired(p,'n');
            addRequired(p,'h');
            addRequired(p,'phi');
            addRequired(p,'c');
            addRequired(p,'g');
            addParameter(p,'CutType',cutType,@(x) any(validatestring(x,cutOptions)));
            addParameter(p,'Name',"wrist");
            parse(p,OD,ID,n,h,phi,c,g,varargin{:});
            
            obj.cutType = p.Results.CutType;
            obj.name = p.Results.Name;
            %*********************************************
            
            obj.OD = OD;
            obj.ID = ID;
            obj.n = n;
            obj.h = h;
            obj.phi = phi;
            obj.c = c;
            obj.g = g;
            obj.theta = zeros(n,1);
            obj.s = zeros(n,1);
            obj.kappa = zeros(n,1);
            obj.precurve_theta = zeros(n,1);
            
            obj.get_neutral_axis(); % Define geometry parameters
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
            type = 'force';
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
                    obj.F = q(1);
                    obj.get_force_arc_params();
                case 'geometry'
                    obj.delta_l = q(1);
                    obj.get_geom_arc_params(q(1));
            end
            [path,T_tip] = obj.robot_kin();
        end
        
        
        function [path,T,obj] = robot_kin(obj)
            %ROBOT_KIN - returns the forward kinematics of a notched wrist
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
        
        function [ybar_, I_,obj] = get_neutral_axis(obj)
            %GETNEUTRALAXIS - returns the distance to neutral bending axis and area
            %moment of inertia about bending axis
            ro = obj.OD/2;
            ri = obj.ID/2;
            obj.ybar = zeros(obj.n,1);
            obj.I = zeros(obj.n,1);
            for index = 1:obj.n
                switch(obj.cutType)
                    case 'off-axis'
                        % Assumes we are cutting between center of tube and inner radius
                        alpha_out=acos((obj.g(index)-ro)/ro);
                        ybar_out=(2*ro/3)*(sin(alpha_out)^3)/(alpha_out-sin(alpha_out)*cos(alpha_out));
                        A_seg_out=(ro^2)*(alpha_out-sin(alpha_out)*cos(alpha_out));
                        Jz_oo=(ro^4)/4*(alpha_out-sin(alpha_out)*cos(alpha_out)+2*sin(alpha_out)^3*cos(alpha_out));
                        
                        % Inner Circular Segment
                        alpha_in=acos((obj.g(index)-ro)/ri);
                        ybar_in=(2*ri/3)*(sin(alpha_in)^3)/(alpha_in-sin(alpha_in)*cos(alpha_in));
                        A_seg_in=(ri^2)*(alpha_in-sin(alpha_in)*cos(alpha_in));
                        Jz_oi=(ri^4)/4*(alpha_in-sin(alpha_in)*cos(alpha_in)+2*sin(alpha_in)^3*cos(alpha_in));
                        
                        A_seg=A_seg_out-A_seg_in;
                        Jz_o=Jz_oo-Jz_oi;
                        
                        % ybar and Jz
                        obj.ybar(index)=(ybar_out*A_seg_out-ybar_in*A_seg_in)/A_seg;
                        obj.I(index)=Jz_o-A_seg*obj.ybar(index)^2;
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
                        I_i = (obj.ID/2)^4/4*(phi_o/2 + 1/2*sin(phi_o));
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
            obj.theta = obj.kappa.*obj.s;
            for(i = 1:obj.n)
                obj.check_notch_limits(i);
            end
        end
        
        function obj = get_force_arc_params(obj)
            % INITIALIZATION OF VECTORS %
            obj.theta = obj.precurve_theta; % vector of each notch angle which should start at precurve value
            F_vec = zeros(obj.n,1); % vector of force experienced by each notch
            M_vec = zeros(obj.n,1); % vector of moment experienced by each notch
            E_vec = obj.E_lin*ones(obj.n,1); % effective elastic modulus for each notch
            
            if (obj.DEBUG)
                debug_mat = [];
                debug_vector = [];
            end
            
            theta_delta = 100; % start the change in theta high
            theta_last = obj.theta;
            trials = 0; % trials we have completed
            maxTrials = 10; % maximum number of trials to perform
            while theta_delta > 1E-6 && trials < maxTrials
                for ii = 1:obj.n
                    %                     obj.check_notch_limits();
                    E_vec(ii) = obj.E_lin;
                    % Compute distal force and moment
                    if obj.use_friction
                        F_vec(ii) = obj.F*exp(-obj.mu*sum(obj.theta(1:ii))); % get force experienced by notch
                    else
                        F_vec(ii) = obj.F;
                    end
                    M_vec(ii) = F_vec(ii)*(obj.ybar(ii) + obj.ID/2); % get moment experienced by notch
                    
                    pct = 100;
                    k = 1;
                    % Gradient descent for non-linear modulus
                    if obj.use_non_linear
                        while k<100 && pct>1E-4
                            % Compute angular deflection of current segment
                            obj.theta(ii) = M_vec(ii)*obj.h(ii)/(E_vec(ii)*obj.I(ii));
                            
                            % Check to see if notch angle closes the notch, and
                            % if so, limit the notch angle
                            obj.check_notch_limits(ii);
                            
                            % Compute section arc length, curvature and
                            % strain
                            obj.s(ii)  = obj.h(ii)-(obj.ybar(ii))*abs(obj.theta(ii));
                            obj.kappa(ii) = obj.theta(ii)/obj.s(ii);
                            epsilon = (obj.kappa(ii).*(obj.OD/2-obj.ybar(ii)))./(1+obj.ybar(ii)*obj.kappa(ii));
                            
                            % Update modulus via gradient descent
                            [stress_eff,eta] = obj.get_stress(abs(epsilon));
                            new_E = E_vec(ii)-eta*(E_vec(ii)-stress_eff/abs(epsilon));
                            if isnan(new_E)
                                new_E = E_vec(ii);
                            end
                            
                            % Percent change in modulus (for convergence
                            % checking)
                            pct = abs((new_E-E_vec(ii))/E_vec(ii));
                            
                            % Update modulus guess
                            E_vec(ii) = new_E;
                            
                            if(obj.DEBUG)
                                debug_vec = [trials, k, s1, obj.kappa(ii), epsilon, stress_eff, eta, new_E];
                                debug_mat = [debug_mat; debug_vec];
                            end
                            % Increment
                            k = k+1;
                        end
                    else
                        obj.theta(ii) = M_vec(ii)*obj.h(ii)/(E_vec(ii)*obj.I(ii));
                        
                        % Check to see if notch angle closes the notch, and
                        % if so, limit the notch angle
                        obj.check_notch_limits(ii);
                        obj.s(ii) = obj.h(ii)-(obj.ybar(ii))*abs(obj.theta(ii));
                        obj.kappa(ii) = obj.theta(ii)/obj.s(ii);
                    end
                    obj.theta(ii) = obj.theta(ii) + obj.precurve_theta(ii);
                    obj.check_notch_limits(ii);
                    obj.s(ii)  = obj.h(ii)-(obj.ybar(ii))*abs(obj.theta(ii));
                    obj.kappa(ii) = obj.theta(ii)/obj.s(ii);
                end
                theta_delta = sum(obj.theta)-sum(theta_last);
                %theta_last = obj.theta;
                trials = trials + 1; % increment
            end
            if (obj.DEBUG)
                writematrix(debug_mat,"wrist_fwkin.csv")
            end
        end
        
        function obj = check_notch_limits(obj,index)
            %CHECK_NOTCH_LIMITS - makes sure the notches aren't closing past
            %what they should do based on geometry.
            if obj.theta(index) >= 0 && (obj.h(index)/(obj.OD/2 + obj.ybar(index)) <= obj.theta(index))
                obj.theta(index) = obj.h(index)/(obj.OD/2 + obj.ybar(index));
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
        
        function [stress,eta,obj] = get_stress(obj,strain)
            % GETSTRESS - returns the stress of a single notch based on stress strain
            % curve of nitinol
            
            sigma = @(e) (e<obj.strain_lower).*e*obj.E_lin+...
                (e >= obj.strain_lower && e < obj.strain_upper)*((e-obj.strain_lower)*obj.E_se+obj.strain_lower*obj.E_lin)+...
                (e >= obj.strain_upper)*((1.0E9)*exp(-.01/(e-obj.strain_upper))+obj.strain_lower*1*obj.E_lin+(obj.strain_upper-obj.strain_lower)*obj.E_se);
            stress = abs(sigma(strain));
            eta_fun = @(e) (e<obj.strain_lower).*.5+...
                (e >= obj.strain_lower && e < obj.strain_upper).*1+...
                (e >= obj.strain_upper).*0.3;
            
            eta = eta_fun(strain);
        end
        
        function maxL = calcMaxL(obj)
            %CALCMAXL - determines the maximum tendon displacement before
            %8 percent strain is experienced by the notch
            obj.get_neutral_axis();
            maxStrain = 0.08;
            k = maxStrain/(obj.OD/2 - obj.ybar(1)*(1 + maxStrain));
            maxL = (k*obj.h(1)*(obj.ID/2 + obj.ybar(1))/(1 + obj.ybar(1)*k));
        end
        
        function obj = plot_stick_model(obj,options)
            arguments
                obj Wrist
                options.Linewidth (1,1) double = 2.5
                options.Marker (1,:) char = '.'
                options.MarkerSize (1,1) double = 25
                options.Color (1,3) double = [0 0.4470 0.7410];
            end
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
            plot3(1000*x_vec,1000*y_vec,1000*z_vec,'Linewidth',...
                options.Linewidth,'Marker',options.Marker,'MarkerSize',...
                options.MarkerSize,'Color',options.Color);
            axis equal
            xlabel('X axis (mm)','FontSize',18,'FontName','CMU Serif');
            zlabel('Z axis (mm)','FontSize',18,'FontName','CMU Serif');
            %             xlim([min(x_vec),max(x_vec)+eps])
            %             ylim([min(y_vec),max(y_vec)+eps])
            %             zlim([min(z_vec),max(z_vec)+eps])
            grid on;
        end
        
        function p = plot_notch_angles(obj,F_max,points)
            %PLOT_NOTCH_ANGLES plots the notch angles with respect to force
            %for the given wrist
            notch_angles = zeros(obj.n,points);
            F_vec = linspace(0,F_max,points);
            for i = 1:points
                obj.fwkin([F_vec(i),0,0]);
                notch_angles(:,i) = obj.theta;
            end
            figure()
            hold on
            title('Model output of notch angles with respect to force input','FontSize',16);
            p = plot(F_vec,rad2deg(notch_angles),'LineWidth',2);
            xlabel("Force Input (N)",'FontSize',12);
            ylabel("Notch Angle (deg)",'FontSize',12);
            legend_entries = cell(1,obj.n);
            for e = 1:obj.n
                legend_entries(1,e) = {sprintf('Notch %d',e)};
            end
            legend(legend_entries,'location','southeast');
            hold off
        end
    end
end

