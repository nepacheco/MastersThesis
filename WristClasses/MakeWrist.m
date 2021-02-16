function wrist = MakeWrist(wristType,precurvature)
%MAKEWRIST Creates the specified wrist
%   Based on the input string creates the wrist with the desired parameters
arguments
    wristType char {mustBeMember(wristType,{'90Tube','150Tube','TipFirstTube'})} 
    precurvature logical = false
end
switch (wristType)
    case '90Tube'
        od = 1.62E-3;
        id = 1.4E-3;
        n = 5;
        h = 0.5E-3*ones(n,1);
        phi = zeros(n,1);
        c = 1.5E-3*ones(n,1);
        g = 1.4E-3*ones(n,1);
        % Precurvature for experiment on 12-19-2020
        precurve_values = deg2rad([1.521853776;1.320520452;1.255100512;1.149834263;1.336681514]); 
    case '150Tube'
        od = 1.62E-3;
        id = 1.4E-3;
        n = 5;
        h = 0.8E-3*ones(n,1);
        phi = zeros(n,1);
        c = 1.2E-3*ones(n,1);
        g = 1.4E-3*ones(n,1);        
        % this is the aversage initial reading from the 12-12 experiment 
        precurve_values = deg2rad([2.39580099;2.268315378;2.433246067;1.724263869;2.334074353]);     
    case '150Tube2'        
        od = 1.62E-3;
        id = 1.4E-3;
        n = 5;
        h = 0.8E-3*ones(n,1);
        phi = zeros(n,1);
        c = 1.2E-3*ones(n,1);
        g = 1.4E-3*ones(n,1);        
        % this is the aversage initial reading from the 12-12 experiment 
        precurve_values = deg2rad([2.536745799	1.866526998	2.31430791	2.24697836	2.114218116]);
    case 'TipFirstTube'
        od = 1.62E-3;
        id = 1.4E-3;
        n = 4;
        h = 1.0E-3*ones(n,1);
        phi = zeros(n,1);
        c = 1.0E-3*ones(n,1);
        g = [1.36,1.39,1.42,1.45].*1E-3;
        % tip first bending 01-12-20201 Experiment
        precurve_values = deg2rad([2.191;2.264;2.534;3.062]);
end
wrist = Wrist(od,id,n,h,phi,c,g,'CutType','on-axis');
if precurvature
    wrist.precurve_theta = precurve_values;
end
end

