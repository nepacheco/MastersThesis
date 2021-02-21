%% Using MATLAB To optimize parameters
A = [];
b = [];
Aeq = [];
beq = [];
lb = [5E9,1.75E9,0.01,0.05];
ub = [30E9,4E9,0.04,0.35];
x0 = [10E9,2.25E9,0.0272,0.2389];
[x,fval,exitflag,output,lambda,grad,hessian]  = fmincon(@ObjectiveFunction,x0,A,b,Aeq,beq,lb,ub);

%% Creating a table 
x = x_withTip
set_num = PropertySets.ID(end) + 1;
E_lin = x(1);
E_se = x(2);
strain_lower = x(3);
mu = x(4);
parameter_time = string(datetime(now,'ConvertFrom','datenum'));
Precurve = true;
expFiles = "all"
Tube = "multiple"

new_row = {set_num,E_lin,E_se,strain_lower,mu,Tube,...
        Precurve,expFiles',parameter_time};
% Don't add if duplicate

PropertySets = [PropertySets; new_row];
valsTip = CompareModel('TipFirstTube',experimentFilesTip(3),true, PropertySets(end,:),'Force',3);
vals150 = CompareModel('150Tube2',experimentFiles150_2(1),true, PropertySets(end,:),'Force',3);
vals90 = CompareModel('90Tube',experimentFiles90(2),true,PropertySets(end,:),'Force',3);
