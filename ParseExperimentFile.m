function [force_vec, notch_mat, tendon_vec] = ParseExperimentFile(file,n,options)
arguments
    file % the opened file
    n = 5 % number of notches
    options.force_index = 4 % the force index
    options.tendon_index = 3 
    options.notch1_index = 5
end
force_vec = [];
tendon_vec = [];
notch_mat = [];
[N,M] = size(file);
for v = 1:N
    if (isnumeric(file{v,options.force_index}) ...
           && isnumeric(file{v,options.notch1_index}))
        force_vec = [force_vec; file{v,options.force_index}];
        notch_mat = [notch_mat; file{v,options.notch1_index:options.notch1_index+n}];
        tendon_vec = [tendon_vec; file{v,options.tendon_index}];
    end
end
end