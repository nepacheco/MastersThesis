function rmse_summation = ObjectiveFunction(x)
%%
load('ExperimentFiles.mat')
E_lin = x(1);
E_se = x(2);
Strain_Lower = x(3);
Mu = x(4);

set = table(E_lin,E_se,Strain_Lower,Mu);
valsTip = CompareModel('TipFirstTube',experimentFilesTip(3),true, set,'Force',3,'Plot',false);
vals150 = CompareModel('150Tube2',experimentFiles150_2(1),true, set,'Force',3,'Plot',false);
vals90 = CompareModel('90Tube',experimentFiles90(2),true,set,'Force',3,'Plot',false);

rmse_summation = norm(valsTip(1:end-1,1,1)) + norm(vals150(1:end-1,1,1)) + norm(vals90(1:end-1,1,1));

end
