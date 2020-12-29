% clear
% clc
% figure(3)
% clf

%% Initializations

% E = 40e9;         %We had 75 GPa Young's Modulus in paper
% sigma_up = 400e6; %We had 600 MPa upper plateau stress in paper
% sigma_lp = 600e6;
% mus = 0.2; %coefficient of static friction. Hunter originally chose 0.4 using nonlinear least squares
h_inc=.01;

E = 40e9;         %We had 75 GPa Young's Modulus in paper
sigma_up = 40e9*0.02; %We had 600 MPa upper plateau stress in paper
sigma_lp = 40e9*0.1;
mus = 0.2; %coefficient of static friction. Hunter originally chose 0.4 using nonlinear least squares


h_vec = [0.8 0.8 0.8 0.8 0.8] .* 1e-3;
g_vec = 1.4*ones(1,5) .* 1e-3;

max_strain = [];

for i = 1:1
    %Calibrated for manufacturing tolerances - specified 0.97 mm, got 0.987 mm
    %g = 0.97e-3;     %cut depth in [m]. See Figure 4.
    %g = 0.988e-3;
    %g = g_vec(i);
    g = 1.4e-3;
    %ro = 1.16e-3/2;   %outer radius of tube in [m].
    %ri = 0.86e-3/2;   %inner radius of tube in [m].
    ro = 1.62e-3/2;   %outer radius of tube in [m].
    ri = 1.40e-3/2;   %inner radius of tube in [m].
    d = g-ro;         %intermediate variable. Depth of cut as measured from y = 0. See Figure 4.
    %h = 0.508e-3;     %height of cut in [m]
    h = 0.8e-3;
    n = 5;            %number of cuts
    
    % %Micro Wrist
    % g = 0.33e-3;     %cut depth in [m]. See Figure 4.
    % ro = 0.46e-3/2;   %outer radius of tube in [m].
    % ri = 0.28e-3/2;   %inner radius of tube in [m].
    % d = g-ro;         %intermediate variable. Depth of cut as measured from y = 0. See Figure 4.
    % h = 0.15e-3;     %height of cut in [m]
    % n = 5;            %number of cuts
    
    %Get experimental data. This is the data that Hunter took, I believe.
    %dz = [0:0.2:5]';
    % DATA = [8.13	0.003838177
    % 12.48	0.644770486
    % 17.35	1.353608926
    % 22.62	2.028785224
    % 29.15	2.570659961
    % 36.67	3.102147051
    % 47.43	3.533364544
    % 60.8	3.907519336
    % 76.08	4.220910083
    % 87.63	4.617639811
    % 100.65	4.934550852];
    
%     DATA = [5.9	0.004154105
%         11	0.357560665
%         15.4	0.776502946
%         22.9	1.20881778
%         28.8	1.622529453
%         35.3	1.947777032
%         45	2.220738834
%         55	2.409689019
%         65.6	2.541027722
%         75.1	2.66422041
%         85.1	2.816352783
%         95.1	2.949542383
%         105.6	3.118125643
%         114.9	3.247102412
%         122.9	3.405919401
%         131.6	3.632075577
%         137.7	3.880093348
%         143.7	4.202531803];
    
    theta_history = deg2rad(DATA(:,1) - DATA(1,1)); %account for initial offset.
    theta_history = DATA(:,1); %convert to degrees
    T_history = DATA(:,2);                          %get tendon force history
    
    %% Calculations
    %These formulas compute the centroid based on the integral formula, using
    % phi as acos(d/r), and simply integrating over r numerically. Note that
    % theta is being used here as a polar coordinate, and not to describe the
    % rotation of the wrist.
    fy = @(r,theta) (r.*cos(theta));
    fA = @(r,theta) (1);
    phi = @(r) acos(d./r);
    A = integral2( @(r,theta) r.*fA(r,theta), ri,ro,@(c)-phi(c),@(c)phi(c),'AbsTol',1e-12,'RelTol',1e-12 );
    ybar = 1/A*integral2( @(r,theta) r.*fy(r,theta),ri,ro,@(c)-phi(c),@(c)phi(c),'AbsTol',1e-12,'RelTol',1e-9 );
    %Instead of using the integral formula, you could use Equations (1) and
    %(2). Note that phi is theta/2, as theta is defined in:
    %http://en.wikipedia.org/wiki/Circular_segment
    
    %formula for strain as a function of distance y and curvature kappa
    %---you used this---strain = @(y,k)(k.*(y-ybar)); %Equation (13) !!!! New Equation basing curvature on the neutral axis
    strain = @(y,k)(k.*(y-ybar)/(1+ybar*k));
    
    %strain energy density as a function of strain, decomposed
    % into a linear part and a nonlinear part
    W = @(e)( (e <= sigma_up/E & e > 0).*(1/2*E*e.^2) + ...
        (e > sigma_up/E & e > 0).*(1/2*sigma_up^2/E + (e-sigma_up/E).*sigma_up) + ...
        (e <= -sigma_lp/E & e <= 0).*(1/2*sigma_lp^2/E + (-e-sigma_lp/E).*sigma_lp) + ...
        (e > -sigma_lp/E & e <= 0).*(1/2*E*e.^2));
    %Integral of Equation (15)
    
    %total strain energy stored in the wrist, as a function of kappa, computed
    %using the integral over the area, multiplied by the height, this has units
    %of [h]*[r]*[W(e)]*[dr]*[dtheta] = [m][m][ Pa ][m] = N/m^2*m^3 = Nm = J.
    %Note that theta is a dummy variable here:
    U = @(kappa) (n)*( h*integral2( ...
        @(r,theta)( r.*W(strain(r.*cos(theta),kappa)) ), ...
        ri,ro,...
        @(c)-phi(c), @(c)phi(c), ...
        'AbsTol',1e-12, ...
        'RelTol',1e-9 ...
        ) ...
        ); %Equation (16)
    
    %moment arm length
    L = (ro+ri)/2 + ybar; %Needed for Equation (17)
    
    %-------------here's the problem---------------------%
    %---you used this---kappa = @(theta)( theta/h );
    %what we had before-- kappa measured from rotation center to neutral plane. Probably wrong.
    options = optimoptions('fsolve','TolFun',1e-17,'Display','off','TolX',1e-12);
    kappa = @(theta) fsolve(@(kappa_new) theta - (h.*kappa_new)./(1+ybar*kappa_new),500,options);
    %kappa is measured from rotation center to center axis, consistent with
    %defintion in paper.
    %-----------------------------------------------------%
    
    dth = 1e-4; %set finite difference step
    gamma = @(theta)( pi - theta/2 ); %angle for calcuating friction losses. See Figure 3.
    eta = @(theta) (sin(gamma(theta)/2)-mus*cos(gamma(theta)/2))/(sin(gamma(theta)/2) + mus*cos(gamma(theta)/2));
    
    %friction model based on static balance equation at corner.
    F = @(theta)(1./L).*(1./(eta(theta).^(2*n))) ...
        .*( U(kappa(theta+dth)) - U(kappa(theta-dth)))/(n*2*dth);
    
    %Friction model based on "approximate" capstan friction,
    %in which the friction grows exponentially as the angle that the tendon
    %curves through grows
    %F = @(theta)(1./L).*(exp(mus*theta*n)) ...
    %           .*( U(kappa(theta+dth)) - U(kappa(theta-dth)))/(2*n*dth);
    
    %Force is partial of strain energy with respect to angle of single cutout, divided by the
    %moment arm length (due to Castigliano's first theorem). This is Equation (19)
    %expressed in finite differences. Note that 'n' cancels with 'n' from
    %Equation (16), above.
    
    
    %% Plot force vs. bending angle for the wrist.
    figure(3)
    hold on
    grid on
    xlabel('\theta'); ylabel('F');
    title('Tendon force vs. wrist rotation','HorizontalAlignment', 'center','FontSize',14)
    theta = linspace(0,2.5/n,50); %look at wrist bending through 2.5 radians.
    %We assume each cutout undergoes the same amount of bending.
    Ftheta = arrayfun(F,theta); %generate force required to bend each cutout by theta.
    %plot(theta*n,Ftheta,'LineWidth',3,'Color',(1/256)*[255 128 0]); %plot the model - radians
    plot(Ftheta,theta*n*180/pi,'LineWidth',3,'Color',(1/256)*[255 128 0]); %plot the model - degrees
    plot(F_vec,rad2deg(sum(theta_mat_force(:,:))),'Color','g','LineWidth',3);
    plot(T_history,theta_history,'.','MarkerSize',20,'Color',(1/256)*[96 96 96]);
    %plot the experimental data
    legend('Swaney Model','Josh Model','Experimental','Location','NorthWest')
    res = get(0,'ScreenSize');
    
    max_strain = [max_strain (ro-ybar)/(ro+ybar)];
    %set(gcf,'Position',[20 res(4)*450/1080 res(3)*655/1920 res(4)*500/1080]); %move the plot window out of the way.
    
    %%
    %matlab2tikz('fig_F_vs_Theta_NEW.tikz')
    
    
    %% fmincon to find best coefficient of static friction.
    %
    % %--Our function to minimize is "find_statics_parameters.m"---
    % a = [E,sigma_up,sigma_lp];
    % x0 = [0.01];     %initial values
    % res_initial = find_statics_parameters(x0,a)   %calculate initial residual, for reference
    % lb = [0]';       %lower bound for mus
    % ub = [1]';     %upper bound for mus
    % options = optimoptions('fmincon','Algorithm','active-set','TolCon',1e-12,'TolFun',1e-15,'Display','off');    %set fmincon options
    % mu_best = fmincon(@(x) find_statics_parameters(x,a),x0,[],[],[],[],lb,ub,[],options)               %run fmincon
    %
    % res = find_statics_parameters(x,a)    %check final residual (could probably do with an fmincon option?)
    %
    % %Note paper Hunter referenced in paper for info on Nitinol material
    % %properties:
    % %http://mrkspecialitymaterials.com/pdf/Optimisation%20of%20processing%20and%20properties%20of%20medical%20grade%20nitinol%20wire%20technical%20paper.pdf
end