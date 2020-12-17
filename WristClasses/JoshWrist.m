classdef JoshWrist < handle
% Description: This class defines a framework for constructing slotted 
% wrists with variable notch geometries and simulating their kinematics and
% mechanics. Kinematics/Mechanics modeling is a Castigliano-based
% mechanics model where each manipulator segment is treated as a beam with
% a fixed support and a distal moment of F*r, and the notch length (h),
% notch depth (g), and subsequently, the moment of inertia (J) are allowed
% to vary between notches, enabling deviations from constant-curvature
% kinematics. If tendon friction is specified, the model is modified to
% include the effects of friction loss due to the curvature of the wrist.
%
% The class implements a nonlinear stress-strain model for Nitinol and an
% iterative algorithm for estimating effective modulus based on the strain
% generated by the kinematics. Under certain circumstances, this algorithm
% may not converge especially for high strains so check the command prompt
% to see which joints may not have converged.
%
% This class can be coupled with g-code generating scripts to automatically
% generate g-code based on the inputted notch parameters through the 
% GCode() method. To enable this functionality, ensure that 
% niti_wrist_gcode_generator.m and alu_fix_gcode_generator.m are in the 
% same parent directory.
%
% EXAMPLE USAGE
%
% % Define Wrist Parameters
% T_0 = HTM(0,0,0,0,0,2E-3);    % Base Transformation
% N = 20;                       % Number of notches
% c = (.4E-3).*ones(N,1);       % Notch spacing vector
% r_o = 1.6E-3/2;               % tube outer radius
% r_i = 1.4E-3/2;               % tube inner radius
% h = linspace(0.35E-3,0.35E-3,N);  % notch base width vector
% h_c = linspace(.15E-3,.15E-3,N);  % notch collision width vector
% g = linspace(1.35E-3,1.35E-3,N);  % notch depth vector
% b = 1.5E-3;                       % distal offset
%
% % Construct Wrist object
% newWrist = Wrist();
% newWrist.ConstructWrist('T_0',T_0,'g',g,'c',c,'b',b,'h',h,'h_c',h_c,'r_o',r_o,'r_i',r_i,'n',N,'material','nitinol','plotkin',1);
%
% % Find maximum tendon force and print results
% newWrist.FindMaxForce(1,5);
% newWrist.PrintResults();
%
%
% LAST UPDATED: 01/11/2019

    % Define class parameters
    properties (Access = public)
        
        % Default wrist geometric parameters
        r_o = .55E-3;               % outer radius of tube
        r_i = .45E-3;               % inner radius of tube
        n = 6;                      % number of notches
        h = 1.5E-3.*ones(6,1);        % vector containing notch base width parameters
        h_c = 1.5E-3.*ones(6,1)       % vector containing notch top width parameters
        c = 1.0E-3.*ones(6,1);        % vector containing notch spacing parameters
        g = .8E-3.*ones(6,1);         % vector containing tube notch height 
        b = 1E-3;           % Length from last notch to tip of tube
        
        % Default Wrist Kinematics
        T_0 = eye(4,4);     % transformation matrix of initial sheath orientation in global frame
        plotkin = true;        % 1 to plot kinematics, 0 to suppress
        X = [];             % vector containing sheath joint coordinates
        T_ef = [];          % End effector transformation matrix
        d_ef = [];          % End effector direction vector
        p_ef = [];          % End effector coordinates
        theta = [];         % angular deflection vector
        delta_L = [];       % tendon displacement vector
        
        % Default Wrist Mechanics
        J = [];            % second moment of area vector
        M = [];             % Vector containing tip moments for each segment
        F = [];             % Vector containing tip forces for each segment
        ybar = [];          % neutral axis vector of tube
        A_seg = [];         % cross-sectional area of notch segment
        strain = [];        % tube strain vector
        sig = [];           % tube stress vector
        F_max = [];         % Maximum applied force
        FOS = [];           % tube strain factor-of-safety
        FOS_b = [];         % tube buckling factor-of-safety
        material = 'nitinol';               % Sheath material
        material_model = 'superelastic';    % material model
        verbose = true;        % Verbose output for force iterations
        E_eff = [];         % Vector containing effective moduli for each segment
        mu = 0.2;           % Tendon friction coefficient
    end
    
    properties (Constant)
        epsilon_up = 0.1;    % Upper plateau strain (for nitinol mechanics)
        epsilon_lp = 0.02;   % Lower plateau strain (for nitinol mechanics)
    end
    
    properties (Access = private)
        E = 40e9;           % Young's Modulus 
        strain_max = 0.08;  % Maximum recoverable strain
        K_sig = 1;          % Stress concentration factor (worst case K_sig=2)
        s_o = [];           % outer tube cross sectional arc length [m]
        s_i = [];           % inner tube cross sectional arc length [m]
        alpha_o = [];       % outer tube intermediate angle [rad]    
        alpha_i = [];       % inner tube intermediate angle [rad]   
        limit = [];         % boolean vector for joint limit (=1 if joint closed)
        yield_flag = 0;     % Yield flag
        buckle_flag = 0;    % Buckling flag
    end
        
    methods
          
        function obj = ConstructWrist(obj,varargin)
            % DESCRIPTION: Parse arguments and construct wrist object
            %
            % INPUT:
            %   varargin (various arguments)
            % 
            % OUTPUT - none (implicit class update)
            
            % Input parsing
            valid_trans = @(x) isnumeric(x) && ismatrix(x)  && numel(x)==16;
            p = inputParser;
            addRequired(p,'obj');
            addOptional(p, 'T_0', obj.T_0, valid_trans);
            addOptional(p, 'r_o', obj.r_o, @(x) isnumeric(x) && (x>0));
            addOptional(p, 'r_i', obj.r_i, @(x) isnumeric(x) && (x>0));
            addOptional(p, 'n', obj.n, @(x) isnumeric(x) && (x>0));
            addOptional(p, 'h', obj.h, @(x) isnumeric(x) && min(x)>0);
            addOptional(p, 'h_c', obj.h_c, @(x) isnumeric(x) && min(x)>0);
            addOptional(p, 'c', obj.c, @(x) isnumeric(x) && min(x)>0);
            addOptional(p, 'b', obj.b, @(x) isnumeric(x) && isscalar(x) && (x>0));
            addOptional(p, 'g', obj.g, @(x) isnumeric(x) && min(x)>0);
            addOptional(p, 'plotkin', obj.plotkin, @(x) islogical(x));
            addOptional(p, 'verbose', obj.verbose, @(x) islogical(x));
            addOptional(p, 'material', obj.material,...
                @(s) any(validatestring(s, {'nitinol','PEEK','steel','PEBAX'})));
            parse(p,varargin{:});
            
            % parse the inputparser object results
            obj.T_0 = p.Results.T_0;
            obj.r_o = p.Results.r_o;
            obj.r_i = p.Results.r_i;
            obj.n = p.Results.n;
            obj.h = p.Results.h;
            obj.h_c = p.Results.h_c;
            obj.c = p.Results.c;
            obj.b = p.Results.b;
            obj.g = p.Results.g;
            obj.plotkin = p.Results.plotkin;
            obj.verbose = p.Results.verbose;
            obj.material = p.Results.material;
            
            % Material properties
            switch obj.material
                case 'nitinol'
                    obj.E = 40e9;
                    obj.strain_max = 0.12;
                    obj.material_model = 'superelastic';
                case 'PEEK'
                    obj.E = 3.75e9;
                    obj.strain_max = 0.043;
                    obj.material_model = 'linear';
                case 'PEBAX'
                    obj.E = .51E9;
                    obj.strain_max = .18;
                    obj.material_model = 'linear';
                case 'steel'
                    obj.E = 200e9;
                    obj.strain_max = 0.002;
                    obj.material_model = 'linear';
            end
            
            % Error Checking for invalid parameter dimensions
            if length(obj.h)==obj.n && length(obj.c)==obj.n && length(obj.g)==obj.n
            elseif isscalar(obj.h) && isscalar(obj.c) && isscalar(obj.g)
                obj.h = obj.h.*ones(obj.n,1);
                obj.c = obj.c.*ones(obj.n,1);
                obj.g = obj.g.*ones(obj.n,1);
            else
                error_str = 'Check dimension of g, h, and c, must either be scalar or equal to n';
                error(error_str);
            end
            
            % Checking for feasible sheath design
            if max(obj.g)>(obj.r_o*2)
                error('Notch depth of tube (g) cannot be greater than tube diameter!');
            end
            
            % Initialize joint angle vector
            obj.theta = zeros(obj.n,1);
            
            % If sheath design is feasible, compute second moments and
            % neutral axes of flexible sections, set notch limit flags to zero
            obj.GetInertialProperties();
            obj.limit = zeros(obj.n,1);
            
        end
        
        function obj = PrintResults(obj)
            % DESCRIPTION: Print results after study
            %
            % INPUT - none (implicit class update)
            % 
            % OUTPUT - none (implicit class update)
            
            if isempty(obj.T_ef)
                error('Cannot print results without running GetKinematicsForce() or MaxForce() first!');
            else
                fprintf('TUBE:\n');
                fprintf('Maximum Strain: %f\n',max(abs(obj.strain)));
                fprintf('Maximum Stress: %f\n',max(abs(obj.sig)));
                fprintf('Minimum Buckling FOS: %f\n',min(obj.FOS_b));
                fprintf('Overall Deflection Angle: %f deg\n',sum(obj.theta).*180/pi);
                fprintf('Overall Pullwire Displacement: %f [mm]\n',sum(obj.delta_L)*1000);
                fprintf('OVERALL RESULTS:\n');
                fprintf('Yield: %i %i\n',obj.yield_flag);
                fprintf('Buckling: %i\n',obj.buckle_flag);
                
                % If there are any notch closures, list them
                if sum(obj.limit)>0
                    fprintf('Notch closure on joint ');
                    lim_idx=find(obj.limit);
                    for i = 1:length(lim_idx)
                        fprintf('%i, ',lim_idx(i));
                    end
                    fprintf('\n');
                end
            end
            
        end
        
        function obj = GetKinematicsForce(obj,F)
            % DESCRIPTION: Generate full kinematics and mechanics of the
            % steerable wrist subject to tendon force F
            %
            % INPUT
            %   F: tendon force (scalar)
            %
            % OUTPUT - none (implicit class update)
            
            % Allocating space for solution vectors
            obj.X = zeros(3,2*obj.n+1);         % joint coordinates
            obj.strain = zeros(1,obj.n);        % joint strains
            obj.FOS_b = zeros(1,obj.n);           % joint buckling FOS
            
            obj.X(:,1) = obj.T_0(1:3,4);  % Initialize first joint coords
            Prev_T = obj.T_0;             % Initialize transformation matrix
            
            % Transformation of tip w.r.t base frame
            T = @(theta,s)[cos(theta) -sin(theta) 0 s*sin(theta)/theta;
                sin(theta) cos(theta) 0 s*(1-cos(theta))/theta;
                0 0 1 0;
                0 0 0 1];
            subplot(2,2,[1 3])
            
            % Draw axes
            if (obj.plotkin), obj.MakeAxes(Prev_T);hold on; end
            
            % Can introduce fudge factor so results agree with experiment
            F_adj = F;
            
            theta_last = obj.theta;
            theta_delta = 100;
            max_iter = 10;
            
            % Only run tendon model if friction is present
            if obj.mu==0
                p = max_iter-1;
            else
                p = 1;
            end
            debug_mat = [];
            debug_vector = [];
            % Fixed-point iteration
            while p<10 && theta_delta>1E-6
                
                Prev_T = obj.T_0;   % Update previous transform
                
                % Compute kinematics for flexible section
                for i = 1:2*obj.n
                    
                    % Calculate mechanics for flexible sections
                    if mod(i,2)
                        
                        % initialize modulus guess
                        idx1 = (i+1)/2;
                        obj.E_eff(idx1)=obj.E;
                        
                        % Compute distal force and moment
                        obj.F(idx1) = F_adj*exp(-obj.mu*sum(obj.theta(1:idx1)));
                        obj.M(idx1) = obj.F(idx1)*(obj.ybar(idx1)+obj.r_i);
                        if strcmp(obj.material_model,'superelastic')
 
                            pct = 100;
                            k = 1;
                            
                            % Gradient descent for non-linear modulus
                            while k<100 && pct>1E-4
                                                              
                                % Compute angular deflection of current segment
                                obj.theta(idx1) = obj.M(idx1)*obj.h(idx1)/(obj.E_eff(idx1)*obj.J(idx1));
                               
                                % Check to see if notch angle closes the notch, and
                                % if so, limit the notch angle
                                obj.CheckLimits(idx1);
                                
                                % Compute section arc length, curvature and
                                % strain
                                s1 = obj.h(idx1)-(obj.ybar(idx1))*abs(obj.theta(idx1));
                                kappa = obj.theta(idx1)/s1;
                                epsilon = (kappa.*(obj.r_o-obj.ybar(idx1)))./(1+obj.ybar(idx1)*kappa);
                                
                                % Update modulus via gradient descent
                                [stress_eff,eta] = obj.SuperElastic(abs(epsilon));
                                new_E = obj.E_eff(idx1)-eta*(obj.E_eff(idx1)-stress_eff/abs(epsilon));
                                
                                % Percent change in modulus (for convergence
                                % checking)
                                pct = abs((new_E-obj.E_eff(idx1))/obj.E_eff(idx1));
                                % NICK PRINTING STUFF %
                                debug_vec = [p, k, s1, kappa, epsilon, stress_eff,eta, new_E];
                                debug_mat = [debug_mat; debug_vec];
                                %%%%%%%
                                % Update modulus guess
                                obj.E_eff(idx1) = new_E;
                                
                                % Increment
                                k = k+1;
                            end
                            
                            % Check for convergence
                            if k==100
                                fprintf('Joint %i strain did not converge after %i cycles, Check Results:\n',idx1,k);
                                fprintf('Tube epsilon: %f pct\n',100*abs(pct));
                            end
                            
                        else
                            % Linear-elastic materials - no gradient
                            % descent - just solve kinematics with constant
                            % modulus
                            
                            % Compute angular deflection of current segment
                            obj.theta(idx1) = obj.M(idx1)*obj.h(idx1)/(obj.E*obj.J(idx1));
                            
                            % Check to see if notch angle closes the notch, and
                            % if so, limit the notch angle
                            obj.CheckLimits(idx1);
                            
                            % Compute section arc length, curvature and strain
                            s1 = obj.h(idx1)-obj.ybar(idx1)*abs(obj.theta(idx1));
                            kappa = obj.theta(idx1)/s1;
                            epsilon = (kappa.*(obj.r_o-obj.ybar(idx1)))./(1+obj.ybar(idx1)*kappa);
                            
                            stress_eff = abs(epsilon).*obj.E;
                        end
                        
                        % Compute stress
                        obj.sig(idx1) = stress_eff;
                        
                        % Compute strain factor-of-safety
                        obj.FOS(idx1) = obj.strain_max/abs(epsilon);
                        
                        % Compute buckling factor-of-safety based on Euler's
                        % critical load formula (K = 0.9)
                        obj.FOS_b(idx1) = (pi^2)*obj.E_eff(idx1)*obj.J(idx1)/(abs(F_adj)*(.9*obj.h(idx1))^2);
                        
                        % Propagate kinematics forward
                        obj.strain(idx1) = epsilon;
                        T1 = T(obj.theta(idx1),s1);
                        Xa = T1(:,4);
                        Xao = Prev_T*Xa;
                        obj.X(:,i+1) = Xao(1:3);
                        Prev_T = Prev_T*T1;
                        
                        % Update pullwire displacement
                        obj.delta_L(idx1) = obj.h(idx1)-2*(1/kappa-obj.r_i)...
                            *sin(kappa*obj.h(idx1)/(2+2*obj.ybar(idx1)*kappa));
                        
                        if (obj.plotkin), obj.MakeAxes(Prev_T); hold on; end
                        
                        % Rigid (non-bending) section
                    else
                        % Assume no bending in this section so just propagate
                        % kinematics forward
                        idx2 = (i)/2;
                        if i==2*obj.n
                            T2 = obj.HTM(0,0,0,obj.b,0,0);
                        else
                            T2 = obj.HTM(0,0,0,obj.c(idx2),0,0);
                        end
                        Xb = T2(:,4);
                        Xbo = Prev_T*Xb;
                        obj.X(:,i+1) = Xbo(1:3);
                        Prev_T =Prev_T*T2;
                        
                    end
                    
                end
                
                % compute norm difference between theta[n] and theta[n-1]
                theta_delta = sum(obj.theta)-sum(theta_last);
                p=p+1;

            end
         
            obj.T_ef = Prev_T;                 % Transformation Matrix
            obj.p_ef = obj.X(:,end);           % End effector coordinates
            obj.d_ef = obj.p_ef-obj.X(:,end-1);% End effector direction vector
            
            % Check for yielding or buckling
            if min(obj.FOS)<1, obj.yield_flag = 1; end
            if min(obj.FOS_b)<1, obj.buckle_flag = 1; end
            
            % Plot results if desired
            if obj.plotkin
                colormap jet;
                
                % Plot kinematic reconstruction
                subplot(2,2,[1 3]);
                plot3(obj.X(1,:),obj.X(2,:),obj.X(3,:),'-ko','MarkerFaceColor',[.5 .5 .5]);hold off;
                xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
                grid on;
                axis equal;
                title('Kinematic Reconstruction');
                
                % Plot stress and strain
                subplot(2,2,2)               
                yyaxis left
                cla;
                plot(1:1:obj.n,obj.strain,'k-o','MarkerFaceColor','r');hold on;
                xlabel('Joint Number'); ylabel('Strain');
                drawnow;
                yyaxis right
                cla;
                plot(1:1:obj.n,obj.sig./1000000,'k-o','MarkerFaceColor','b');hold on;
                ylabel('Stress [MPa]');
                grid on;
                legend('Strain','Stress','Location','SouthEast');
                xlim([1 obj.n]);
                ax = gca;
                ax.YColor = 'b';
                drawnow;
                
                % Plot buckling factor of safety
                subplot(2,2,4)
                cla;
                semilogy(1:1:obj.n,obj.FOS_b,'k-o','MarkerFaceColor','r');hold on;
                xlabel('Joint Number'); ylabel('Factor-Of-Safety');
                semilogy(1:1:obj.n,obj.FOS,'k-o','MarkerFaceColor','b');hold on;
                line([1 obj.n],[1 1],'Color','k','LineStyle','--'); hold off;
                grid on;
                legend('Buckling','Yielding');
                xlim([1 obj.n]);
                set(gcf,'Position',[100 100 1000 400]);
                drawnow;
                
            end
            writematrix(debug_mat,"josh_fwkin.csv")
        end
        
        function obj = CheckLimits(obj,idx)
            % DESCRIPTION: Check for joint limits at current joint index
            %
            % INPUT
            %   idx: Joint number to check (scalar)
            %
            % OUTPUT - none (implicit class update)
            
            % If deflection causes the notched section to close entirely,
            % prevent further deflection
            if obj.theta(idx)>=0 && (obj.h_c(idx)/(obj.r_o+obj.ybar(idx))<=obj.theta(idx))
                obj.theta(idx) = obj.h_c(idx)/(obj.r_o+obj.ybar(idx));
                obj.limit(idx) = 1;
            else
                obj.limit(idx) = 0;
            end
        end
   
        
        function obj = MakeAxes(obj,T)
            % DESCRIPTION: Draw local axes based on transformation matrix
            %
            % INPUT
            %   T: Homogeneous Transformation Matrix (4x4)
            %
            % OUTPUT - none (implicit class update)
            
            length = min(obj.h); % length of each axis
            
            x = T(1,4);
            y = T(2,4);
            z = T(3,4);
            
            ax1 = T*[length 0 0 1]';
            ax2 = T*[0 length 0 1]';
            ax3 = T*[0 0 length 1]';
            
            line([x ax1(1)],[y ax1(2)],[z ax1(3)],'Color','g','LineWidth',2); hold on;
            line([x ax2(1)], [y ax2(2)], [z ax2(3)],'Color','r','LineWidth',2); hold on;
            line([x ax3(1)], [y ax3(2)], [z ax3(3)],'Color','b','LineWidth',2); hold off;  
        end
        
        function [stress,eta] = SuperElastic(obj,strain)
            % DESCRIPTION: Nonlinear stress-strain characteristics of superelastic
            % nitinol
            %
            % INPUT
            %   strain: joint strain (scalar)
            % OUTPUT
            %   stress: joint stress (scalar)
            %   eta: learning rate 
            
            % Nonlinear stress-strain characteristics of Nitinol
            sigma = @(e) (e<obj.epsilon_lp).*e*obj.E+...
                (e >= obj.epsilon_lp && e < obj.epsilon_up)*((e-obj.epsilon_lp)*.08*obj.E+obj.epsilon_lp*obj.E)+...
                (e >= obj.epsilon_up)*((1.0E9)*exp(-.01/(e-obj.epsilon_up))+obj.epsilon_lp*1*obj.E+(obj.epsilon_up-obj.epsilon_lp)*.08*obj.E);
            stress = sigma(abs(strain));
            
            % Learning rate (based on where in the stress-strain curve we
            % are
            eta_fun = @(e) (e<obj.epsilon_lp).*.5+...
                 (e >= obj.epsilon_lp && e < obj.epsilon_up)*1+...
                 (e >= obj.epsilon_up).*.3;
            eta = eta_fun(strain);      
        end
        
        function [elasticModulus, eta] = GetElasticModulus(obj,strain)
           em = @(e) (e < obj.epsilon_lp).*obj.E + ...
               (e >= obj.epsilon_lp && e < obj.epsilon_up).*0.08*obj.E + ...
               (e >= obj.epsilon_up).*(1.0E9)*exp(-0.01);
           elasticModulus = em(strain);
           eta = 0.1;
        end
        
        function obj = FindMaxForce(obj,FOS_d,alpha)
            % DESCRIPTION: Iteratively find the maximum force before (1)
            % one of the notches drops below the specified
            % factor-of-safety, or (2) all notches reach their limits
            %
            % INPUT
            %   FOS_d: desired factor-of-safety [int]
            %   alpha: learning rate
            %
            % OUTPUT - none (implicit class update)
           
            F_guess = .1;        % Initial guess
            error = 100;
            max_iter = 80;
            i = 1;
            
            fprintf('=================================================\n');
            fprintf('===============Finding Max Force=================\n');
            fprintf('Desired FOS: %f\n',FOS_d);
            fprintf('Max Iteration Limit: %i\n',max_iter);
            fprintf('Verbose Output: %i\n',obj.verbose);
            fprintf('=================================================\n\n');
            % descent until FOS converges to setpoint or all of the joint
            % limits are reached
            while error>5E-3&&(any(~obj.limit))
                % Compute kinematics at guess
                
                GetKinematicsForce(obj,F_guess);
                
                % get current factors of safety
                FOS_eps = min([obj.FOS]);
                FOS_buck = min([obj.FOS_b]);
                
                % descent
                if FOS_eps<FOS_buck
                    F_new = F_guess+alpha*F_guess*(abs(obj.strain_max/FOS_d) - obj.K_sig*max(abs([obj.strain])));
                else
                    F_new = F_guess+alpha*F_guess*(min(abs([obj.FOS_b]))-abs(FOS_d));
                end
                
                % Update force and descend based on strain difference
                if i==80
                    disp('Did not converge after 80 iterations. Consider lowering learning rate.');
                    break;
                end
                error = abs((F_new-F_guess)/F_guess);
                
                % Print iterations
                if obj.verbose
                    fprintf('Iteration: %i   F: %f   Buckling FOS: %f   Strain FOS: %f   error: %f\n',i,F_guess,FOS_buck,FOS_eps,error);
                end
                i = i + 1;
                F_guess = F_new;
            end
            
            fprintf('\n=================================================\n');
            fprintf('===================Results=======================\n');
            
            % Print failure mode and max force
            if (sum(obj.limit)==obj.n)
                fprintf('Limit: Notch Closure\n');
            elseif FOS_buck<FOS_eps
                fprintf('Limit: Buckling FOS exceeded \n');
            else
                fprintf('Limit: Strain FOS exceeded\n');
            end
            fprintf('Max Force: %f N\n',F_guess);
            fprintf('=================================================\n\n');
            
            obj.F_max = F_guess;
        end
        
        function [y_bar, A_seg, Jz] = NeutralAxis(obj,ro,ri,g)
            % DESCRIPTION: Function to calculate segment cross-sectional
            % area and neutral axis location.
            %
            % INPUT
            %   r_o: outer radius [scalar]
            %   ri: inner radiu [scalar]
            %   g: notch depth [scalar]
            %
            % OUTPUT
            %   y_bar: neutral axis location (from tube center) [scalar]
            %   A_seg: segment cross-sectional area [scalar]
            %   Jz: segment second moment of area [scalar]

            t = ro-ri;  % section thickness
            
            % We're cutting down to or past the inner wall
            if g >= (2*ro-t) && g<2*ro  
                % Outer Circular Segment
                alpha_out=acos((g-ro)/ro);
                ybar_out=(2*ro/3)*(sin(alpha_out)^3)/(alpha_out-sin(alpha_out)*cos(alpha_out));
                A_seg=(ro^2)*(alpha_out-sin(alpha_out)*cos(alpha_out));
                Jz_o=(obj.r_o^4)/4*(alpha_out-sin(alpha_out)*cos(alpha_out)+2*sin(alpha_out)^3*cos(alpha_out));

                % ybar and Jz
                y_bar=ybar_out;
                Jz=Jz_o-A_seg*y_bar^2;
                
            % We're cutting down to between the radius and inner wall
            elseif g>= ro && g<(2*ro-t)  
                % Outer Circular Segment
                alpha_out=acos((g-ro)/ro);
                ybar_out=(2*ro/3)*(sin(alpha_out)^3)/(alpha_out-sin(alpha_out)*cos(alpha_out));
                A_seg_out=(ro^2)*(alpha_out-sin(alpha_out)*cos(alpha_out));
                Jz_oo=(ro^4)/4*(alpha_out-sin(alpha_out)*cos(alpha_out)+2*sin(alpha_out)^3*cos(alpha_out));

                % Inner Circular Segment
                alpha_in=acos((g-ro)/ri);
                ybar_in=(2*ri/3)*(sin(alpha_in)^3)/(alpha_in-sin(alpha_in)*cos(alpha_in));
                A_seg_in=(ri^2)*(alpha_in-sin(alpha_in)*cos(alpha_in));
                Jz_oi=(ri^4)/4*(alpha_in-sin(alpha_in)*cos(alpha_in)+2*sin(alpha_in)^3*cos(alpha_in));

                A_seg=A_seg_out-A_seg_in;
                Jz_o=Jz_oo-Jz_oi;
                
                % ybar and Jz
                y_bar=(ybar_out*A_seg_out-ybar_in*A_seg_in)/A_seg;
                Jz=Jz_o-A_seg*y_bar^2;
            
            % we're cutting out less than the radius, so we need to subtract out the segment from the tube
            elseif g>t  && g<ro  
                % Outer Circular Segment
                alpha_out=acos((ro-g)/ro);
                ybar_out=(2*ro/3)*(sin(alpha_out)^3)/(alpha_out-sin(alpha_out)*cos(alpha_out));
                A_seg_out=(ro^2)*(alpha_out-sin(alpha_out)*cos(alpha_out));
                Jz_oo=(ro^4)/4*(alpha_out-sin(alpha_out)*cos(alpha_out)+2*sin(alpha_out)^3*cos(alpha_out));
                
                % Inner Circular Segment
                alpha_in=acos((ro-g)/ri);
                ybar_in=(2*ri/3)*(sin(alpha_in)^3)/(alpha_in-sin(alpha_in)*cos(alpha_in));
                A_seg_in=(ri^2)*(alpha_in-sin(alpha_in)*cos(alpha_in));
                Jz_oi=(ri^4)/4*(alpha_in-sin(alpha_in)*cos(alpha_in)+2*sin(alpha_in)^3*cos(alpha_in));
                
                A_seg=A_seg_out-A_seg_in;
                A_tube=pi*(ro^2-ri^2);
                Jz_o=Jz_oo-Jz_oi;
                
                % ybar and Jz
                ybar_seg=(ybar_out*A_seg_out-ybar_in*A_seg_in)/A_seg;
                y_bar=ybar_seg*A_seg/(A_tube-A_seg); %b/c ybar of the tube is zero-- centered about is origin
                Jz=Jz_o-A_seg*y_bar^2;
            end

        end    

        function obj = GetInertialProperties(obj)
            % DESCRIPTION: Compute neutral axis, segment area, and second 
            % moment of area for for all cutouts
            %
            % INPUT - none (implicit class update)
            % OUTPUT - none (implicit class update)
  
            for i = 1:obj.n
                [obj.ybar(i), obj.A_seg(i), obj.J(i)] = obj.NeutralAxis(obj.r_o,obj.r_i,obj.g(i));
                obj.alpha_o(i) = acos((obj.g(i)-obj.r_o)/obj.r_o);
                
            end
        end     
        
        function obj = Angle(obj, r_o, r_i, g)
            % DESCRIPTION: Compute intermediate angle and arc length for
            % the current section
            %
            % INPUT
            %   r_o: outer radius [scalar]
            %   ri: inner radiu [scalar]
            %   g: notch depth [scalar]
            %
            % OUTPUT - none (implicit class update)
            
            % Compute angle of current section
            obj.alpha_o = acos((g-r_o)/r_o);
            obj.alpha_i = acos((g-r_i)/r_o);
            
            % Compute cross sectional arc length of current section
            obj.s_o = 2*r_o*((pi/2)-asin(r_o/(g-r_o)));
            obj.s_i = 2*r_o*((pi/2)-asin(r_o/(g-r_o)));
        end

        function T = HTM(obj,rx,ry,rz,dx,dy,dz)
            % DESCRIPTION: Computes homogeneous transformation
            %
            % INPUT
            %   rx: rotation about x [scalar]
            %   ry: rotation about y [scalar]
            %   rz: rotation about z [scalar]
            %   dx: translation along x [scalar]
            %   dy: translation along y [scalar]
            %   dz: translation along z [scalar]
            %
            % OUTPUT - none (implicit class update)
            
            Rx = @(a)[1 0 0 0;
                0 cos(rx) -sin(rx) 0;
                0 sin(rx) cos(rx) 0;
                0 0 0 1];
            
            Ry = @(b)[cos(ry) 0 sin(ry) 0;
                0 1 0 0;
                -sin(ry) 0 cos(ry) 0;
                0 0 0 1];
            
            Rz = @(c)[cos(rz) -sin(rz) 0 0;
                sin(rz) cos(rz) 0 0;
                0 0 1 0;
                0 0 0 1];
            
            % Homogeneous Transformation Matrix
            T = Rx(rx)*Ry(ry)*Rz(rz)*[1 0 0 dx;0 1 0 dy;0 0 1 dz;0 0 0 1];
        end
        
        function obj = GCode(obj,filename)
            % DESCRIPTION: Generate g-code file based on the chosen notch
            % parameters (NOTE: niti_wrist_gcode_generator() must be in the
            % same directory)
            %
            % INPUT
            %   filename: name of output file [string]
            %
            % OUTPUT - none (implicit class update)
                params1 = {obj.n,1000.*obj.r_o,1000.*fliplr(obj.h),1000.*fliplr(obj.g),1000.*fliplr(obj.c),1000.*fliplr(obj.b),2};
                niti_wrist_gcode_generator(params1,strcat(filename));
        end
    end
        
end