%% Analyze Tip Position of image

ptom = 1.62E-3/176; %[m]/[px]
tip_pix = [250;0;743]; % [px]
tip_m = tip_pix*ptom; % [m]

theta_des = deg2rad([6.1504	7.6563	10.6232	16.4125]');
force_input = 0.5990;
tendon_disp = 1.5E-3; % [m]
% Make Wrist
wrist = MakeWrist('TipFirstTube',true);
wrist.name = 'Wrist C';

s_des = wrist.h - (wrist.ybar).*abs(theta_des);  %[m]
kappa_des = theta_des./s_des; % [1/m]
    
%% Plot Experimental Data
wrist.theta = theta_des;
wrist.s = s_des;
wrist.F = force_input;
wrist.kappa = kappa_des;
color = c(4,:);
[~,T] = wrist.robot_kin();
%     wrist.plot_stick_model('Color',color);T

error = norm(T(1:3,4) - tip_m);
fprintf("error: %f\n",error);
fprintf("x y z\n")
fprintf("%f\n",T(1:3,4) - tip_m)