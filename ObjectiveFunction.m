function rmse_summation = ObjectiveFunction(x,tipDeflection)
%
arguments
    x (:,4) double
    tipDeflection (1,1) logical = true
end

%%
numSets = size(x,1);
rmse_summation = ones(numSets,1);
for i = 1:size(x,1)
    load('ExperimentFiles.mat')
    E_lin = x(i,1);
    E_se = x(i,2);
    Strain_Lower = x(i,3);
    Mu = x(i,4);
    
    set = table(E_lin,E_se,Strain_Lower,Mu);
    valsTip = CompareModel('TipFirstTube',experimentFilesTip(3),true, set,'Force',3,'Plot',false);
    vals150 = CompareModel('150Tube2',experimentFiles150_2(1),true, set,'Force',3,'Plot',false);
    vals90 = CompareModel('90Tube',experimentFiles90(2),true,set,'Force',3,'Plot',false);
    
    if tipDeflection
        rmse_summation(i) = norm(valsTip(:,1,1)) + norm(vals150(:,1,1)) + norm(vals90(:,1,1));
    else
        rmse_summation(i) = norm(valsTip(1:end-1,1,1)) + norm(vals150(1:end-1,1,1)) + norm(vals90(1:end-1,1,1));
    end
end
end
