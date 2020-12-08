%% Comparison Script
% This script is meant to compare the outputs from force measurements of
% the nitinol notched wrist bending to the output predicted by Josh's model
clc; clear; close all;
%% Initial Conditions
% This is where we put all the tube parameters
od = 1.62E-3; % [m] - outer diameter of tube
id = 1.4E-3; % [m] - inner diameter of tube
n = 5; % number of notches
g = 0.8642*od; % [m] - Max depth which we assign to notch n.
h = 0;
u = 0;
E_lin = 40E9; % [N/m^2] - Elastic Modulus of Nitinol
E_se = 0.08*E_lin; % [N/m^2] - Slope of Super Elastic Region for Nitinol
mu = 0.2; % coefficient of friction for capstan
