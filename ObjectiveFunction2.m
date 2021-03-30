function output = ObjectiveFunction2(x,experimentDataTip,experimentData90,experimentData150,wristTipFirst,wrist90,wrist150,tipDeflection)
%
arguments
    x (:,4) double
    experimentDataTip (1,:) cell
    experimentData90 (1,:) cell
    experimentData150 (1,:) cell
    wristTipFirst Wrist
    wrist90 Wrist
    wrist150 Wrist
    tipDeflection (1,1) logical = true
end

%%
numSets = size(x,1);
rmse_summation = ones(numSets,1);
r2_summation = ones(numSets,1);
for i = 1:size(x,1)
    E_lin = x(i,1);
    E_se = x(i,2);
    Strain_Lower = x(i,3);
    Mu = x(i,4);
    
    set = table(E_lin,E_se,Strain_Lower,Mu);
    rmseTip = CompareTipModel(wristTipFirst,experimentDataTip,set,'Plot',false);
    rmse150 = CompareTipModel(wrist150,experimentData150,set,'Plot',false);
    rmse90 = CompareTipModel(wrist90,experimentData90,set,'Plot',false);
    
    rmse_summation(i) = norm(rmseTip(:,1,1)) + norm(rmse150(:,1,1)) + norm(rmse90(:,1,1));
end
output = rmse_summation;
end
