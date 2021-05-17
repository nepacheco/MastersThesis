function [force_vec, notch_mat, tendon_vec] = ParseExperimentFile(file,n)
force_index = 4;
force_vec = [];
tendon_index = 3;
tendon_vec = [];
notch1_index = 5;
notch_mat = [];
[N,M] = size(file);
for v = 1:N
    if (isnumeric(file{v,force_index}) && isnumeric(file{v,notch1_index}))
        force_vec = [force_vec; file{v,force_index}];
        notch_mat = [notch_mat; file{v,notch1_index:notch1_index+n}];
        tendon_vec = [tendon_vec; file{v,tendon_index}];
    end
end
end