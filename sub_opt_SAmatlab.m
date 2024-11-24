clear
close all
warning('off','all')
clc

diary('_20240523_logOptMatlabSA_30Q02_30V02_10Q01_10V01_.txt')

%% Inputs

    % Plate        
    b = 0.254;              
    tPly = 0.000150;        
    E1 =   1.81000e11;       
    E2 =   1.02700e10;      
    G12 =  7.17000e09;      
    nu12 = 0.28000;         
    r = inf;      
    
    % Stringer
    tref = 0.000;                    
    href = 0.000;            
    Eref = E1;              
    nuRef = nu12;           
    Gref = Eref/2/(1+nuRef);
    extraIyy = 0;

    % Load Factors
    fx = 1;
    fapp = 1;

    % m and n
    m = 15;
    n = 15;

    % Numerical Integration
    nIntPx = 28;
    nIntPy = 28;

    % Grouping
    inputMaterial =     {E1,E2,nu12,G12,"",Eref,Gref};
    inputGeometry =     {0,b,tPly,r,[],[],[],tref,href,extraIyy};
    inputLoads =        [fx,0,fapp];
    inputMN =           [m,n];
    inputNumIntegral =  [nIntPx,nIntPy];

    % SA
    auxFunc = @(y) y(2);
    nGroups = 4;
    nTimesPerFunc = 3;
    x = cell(nTimesPerFunc,nGroups);
    fval = zeros(nTimesPerFunc,nGroups);
    exitFlag = zeros(nTimesPerFunc,nGroups);
    nFunc = zeros(nTimesPerFunc,nGroups);
    times = zeros(nTimesPerFunc,nGroups);
    options = optimoptions('simulannealbnd','Display','diagnose');
    initGuess = [0 0 0 0];
    lwBdv = [-84 -84 -84 -84];
    upBdv = [84 84 84 84];
    group = 0;


%% 10V01 - T1.11

    matNature = "V";
    AR = 1;
    fxy = 0;
    PcrRef = 17186.2;

    a = b*AR;   
    inputMaterial{5} = matNature;
    inputGeometry{1} = a;
    inputLoads(2) = fxy;
    group = group + 1;

    for i=1:nTimesPerFunc
        rng(i)
        [x{i,group},fval(i,group),exitFlag(i,group),output] = simulannealbnd(@(dv)auxFunc(criticalStabilityLoad_forSA(dv, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral)),...
                                                                                        initGuess,lwBdv,upBdv,options);
        nFunc(i,group) = output.funccount;
        times(i,group) = output.totaltime;
    end
    for j=1:nTimesPerFunc
        fprintf("\n     10V01:  %d  ",j);
        for jj=1:length(initGuess)
            fprintf("%10.4f ",x{j,group}(jj));
        end
        fprintf("%10.2f  %5d %10.2f",fval(j,group),nFunc(j,group),times(j,group));
    end


%% 10Q01 - T1.16

    matNature = "Q";
    AR = 1;
    fxy = 0;
    PcrRef = 17186.2;

    a = b*AR;   
    inputMaterial{5} = matNature;
    inputGeometry{1} = a;
    inputLoads(2) = fxy;
    group = group + 1;

    for i=1:nTimesPerFunc
        rng(i)
        [x{i,group},fval(i,group),exitFlag(i,group),output] = simulannealbnd(@(dv)auxFunc(criticalStabilityLoad_forSA(dv, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral)),...
                                                                                        initGuess,lwBdv,upBdv,options);
        nFunc(i,group) = output.funccount;
        times(i,group) = output.totaltime;
    end
    for j=1:nTimesPerFunc
        fprintf("\n     10Q01:  %d  ",j);
        for jj=1:length(initGuess)
            fprintf("%10.4f ",x{j,group}(jj));
        end
        fprintf("%10.2f  %5d %10.2f",fval(j,group),nFunc(j,group),times(j,group));
    end


%% 30V02 - T1.13p

    matNature = "V";
    AR = 2;
    fxy = 0.5;
    PcrRef = 15866.5;

    a = b*AR;   
    inputMaterial{5} = matNature;
    inputGeometry{1} = a;
    inputLoads(2) = fxy;
    group = group + 1;

    for i=1:nTimesPerFunc
        rng(i)
        [x{i,group},fval(i,group),exitFlag(i,group),output] = simulannealbnd(@(dv)auxFunc(criticalStabilityLoad_forSA(dv, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral)),...
                                                                                        initGuess,lwBdv,upBdv,options);
        nFunc(i,group) = output.funccount;
        times(i,group) = output.totaltime;
    end
    for j=1:nTimesPerFunc
        fprintf("\n     30V02:  %d  ",j);
        for jj=1:length(initGuess)
            fprintf("%10.4f ",x{j,group}(jj));
        end
        fprintf("%10.2f  %5d %10.2f",fval(j,group),nFunc(j,group),times(j,group));
    end


%% 30Q02 - T1.14p

    matNature = "Q";
    AR = 2;
    fxy = 0.5;
    PcrRef = 15866.5;

    a = b*AR;   
    inputMaterial{5} = matNature;
    inputGeometry{1} = a;
    inputLoads(2) = fxy;
    group = group + 1;

    for i=1:nTimesPerFunc
        rng(i)
        [x{i,group},fval(i,group),exitFlag(i,group),output] = simulannealbnd(@(dv)auxFunc(criticalStabilityLoad_forSA(dv, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral)),...
                                                                                        initGuess,lwBdv,upBdv,options);
        nFunc(i,group) = output.funccount;
        times(i,group) = output.totaltime;
    end
    for j=1:nTimesPerFunc
        fprintf("\n     30Q02:  %d  ",j);
        for jj=1:length(initGuess)
            fprintf("%10.4f ",x{j,group}(jj));
        end
        fprintf("%10.2f  %5d %10.2f",fval(j,group),nFunc(j,group),times(j,group));
    end


%% end

diary off