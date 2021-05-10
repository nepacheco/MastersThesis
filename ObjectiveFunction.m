function output = ObjectiveFunction(x,experimentDataTip,experimentData90,experimentData150,wristTipFirst,wrist90,wrist150,tipDeflection)
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
%     [valsTip, r2_tip] = CompareModel(wristTipFirst,experimentDataTip,set,'Plot',false);
    valsTip = zeros(wristTipFirst.n+1, 1); r2_tip = zeros(wristTipFirst.n+1, 1);
%     [vals150, r2_150] = CompareModel(wrist150,experimentData150,set,'Plot',false);
    vals150 = zeros(wrist150.n+1,1); r2_150 = zeros(wrist150.n+1,1);
    [vals90, r2_90] = CompareModel(wrist90,experimentData90,set,'Plot',false);
    
    if tipDeflection
        rmse_summation(i) = norm(valsTip(:,1,1)) + norm(vals150(:,1,1)) + norm(vals90(:,1,1));
        r2_summation(i) = norm(r2_90(:,1,1)) + norm(r2_150(:,1,1)) + norm(r2_tip(:,1,1));
    else
        rmse_summation(i) = norm(valsTip(1:end-1,1,1)) + norm(vals150(1:end-1,1,1)) + norm(vals90(1:end-1,1,1));
        r2_summation(i) = norm(r2_90(1:end-1,1,1)) + norm(r2_150(1:end-1,1,1)) + norm(r2_tip(1:end-1,1,1));
    end
end
output = rmse_summation; %- 2*r2_summation;
end
