function [force_vec, notch_mat] = ParseExperimentFile(file,n)
%ParseExperimentFile - parses the NxM cell passed into a force vector and notch value
%matrix
%   force vector output is Nx1 and notch_mat is Nxn where n is the number
%   of notches in the tube.
force_index = 4;
force_vec = [];
notch1_index = 5;
notch_mat = [];
[N,M] = size(file);
for i = 1:N
    if (isnumeric(file{i,force_index}) && isnumeric(file{i,notch1_index}))
        force_vec = [force_vec; file{i,force_index}];
        notch_mat = [notch_mat; file{i,notch1_index:notch1_index+n}];
    end
end

