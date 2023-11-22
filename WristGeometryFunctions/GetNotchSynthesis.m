function [w,h,u] = GetNotchSynthesis(maxBendingAngle,w,OD,ID,opts)
%GETNOTCHSYNTHESIS Returns cut depth, height, and uncut spacing
%   Determines the notch geometries to achieve a certain amount of bending
%   per angle based on the cut depth, outer diameter, and inner diameter of
%   the tube. Based on Alex's notchsynthesis.m script
arguments
    maxBendingAngle (1,1) double  %[deg] the maximum bending angle for the notch
    w (1,1) double  %[m] The assumed cut depth of each notch
    OD (1,1) double  %[m] outer diameter of the tube in meters
    ID (1,1) double  %[m] Inner diameter of the tube in meters
    opts.L (1,1) double = 15e-3  % [m] total length of the section
    opts.n (1,1) double = 1 % the number of notches
    opts.CutType char {mustBeMember(opts.CutType,{'on-axis','off-axis'})} = 'on-axis'  % Cutting pattern used for the notches
end

% Script to calculate the notch geometry

% Define the total notch + uncut section length:
L = opts.L; %[m]
% Define the desired maximum bending angle for the wrist
thetaMax = maxBendingAngle * pi / 180; % [radians]
%% You should not need to change anything from the code below



ro = OD/2;         % [m] tube outer radius
ri = ID/2;         % [m] tube inner radius

[ybar, ~] = GetNeutralAxis(ro, ri, w,'CutType',opts.CutType);

h = thetaMax * (ro+ybar);
u = L/n - h;
%d = h/(h+u);
disp(['Notch Width (w): ' num2str(w*1e3) ' mm']);
disp(['Notch Height (h) : ' num2str(h*1e3) ' mm']);
disp(['Notch Spacing (u): ' num2str(u*1e3) ' mm']);
end

