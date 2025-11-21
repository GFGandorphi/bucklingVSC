clear
close all
warning('off','all')
clc
told = 0;
tic

% diary('_logSearch_20240926_3VR2_28.txt')

FX =  [1 0 1 1.0  1.0];
FXY = [0 1 1 0.5 -0.5];
T0vec = -84:3:84;
T1vec = -84:3:84;

for AR = 1%[1 2]
for matNature = "Q" %["V" "Q"]
for pair = 1%[1 4 5]
for T0v1 = +60%T0vec
for T1v1 = +15%T1vec
for T0v2 = T0v1%T0vec
for T1v2 = T1v1%T1vec
for iIntPoints = 28%[20 28 36 44 60 100]
for iMN = 15%5:2:17
for fapp = 1e4/0.254 %[1e-4 1e-3 1e-2 1e-1 1e0 1e1 1e2 1e3 1e4 1e5]
for hi = 0 %[0.0000,0.0030,0.0040,0.0060,0.0080,0.0100,0.0125,0.0150,0.0175,0.0200,0.0225,0.0250,0.0275,0.0300,0.0350,0.0400,0.0500,0.0600,0.0700,0.0800,0.1000,0.1600]

    % Plate Inputs
    if upper(matNature) == "Q"
        T0 = [T0v1 T0v2 T0v2 T0v1 T0v1 T0v2 T0v2 T0v1 T0v1 T0v2 T0v2 T0v1 T0v1 T0v2 T0v2 T0v1];
        T1 = [T1v1 T1v2 T1v2 T1v1 T1v1 T1v2 T1v2 T1v1 T1v1 T1v2 T1v2 T1v1 T1v1 T1v2 T1v2 T1v1];
        isMirrored =    [0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0];
    elseif upper(matNature) == "V"
        T0 = [T0v1 -T0v2 -T0v2 T0v1 T0v1 -T0v2 -T0v2 T0v1 T0v1 -T0v2 -T0v2 T0v1 T0v1 -T0v2 -T0v2 T0v1];
        T1 = [T1v1 -T1v2 -T1v2 T1v1 T1v1 -T1v2 -T1v2 T1v1 T1v1 -T1v2 -T1v2 T1v1 T1v1 -T1v2 -T1v2 T1v1];
        isMirrored = zeros(1,length(T0));
    end
    b = 0.254;                  % [m]       length in y
    a = b*AR;                   % [m]       length in x
    tPly = 0.000150;            % [m]       each ply thickness
    r = inf;                    % [m]       plate y radius     
    % FALLAHI-HAO
    E1 =   1.81000e11;          % [Pa]      modulus of elasticity in direction 1      
    E2 =   1.02700e10;          % [Pa]      modulus of elasticity in direction 2
    G12 =  7.17000e09;          % [Pa]      G modulus in direction 12
    nu12 = 0.28000;             % []        Poisson in direction 12
    % % DUTRA
    % E1 =   5.00000e10;          % [Pa]      modulus of elasticity in direction 1      
    % E2 =   2.32200e09;          % [Pa]      modulus of elasticity in direction 2
    % G12 =  6.24000e08;          % [Pa]      G modulus in direction 12
    % nu12 = 0.33300;             % []        Poisson in direction 12

    % a = 10*6*25.4e-3;                  
    % tPly = 0.01*25.4e-3/10;           
    % r = 6*25.4e-3; 
    % b = pi*r;
    % E1 =   94.8029e09;             
    % E2 =    7.1016e09;         
    % G12 =   2.8958e09;         
    % nu12 = 0.25;            

    % Stringer Inputs
    tref = 0.000;              % [m]           
    href = hi*a;               % [m]
    Eref = E1;                 % [Pa] 
    nuRef = nu12;              % []
    Gref = Eref/2/(1+nuRef);   % [Pa] 
    extraIyy = 0;%(tref*href)*(href/2)^2;
    
    % Load Factors
    fx = FX(pair);
    fxy = FXY(pair);

    % Numerical Integration Inputs
    nIntPx = iIntPoints;
    nIntPy = iIntPoints;

    % m and n
    m = iMN;
    n = iMN;

    % Evaluation
    inputMaterial =     {E1,E2,nu12,G12,matNature,Eref,Gref};
    inputGeometry =     {a,b,tPly,r,T0,T1,isMirrored,tref,href,extraIyy};
    inputLoads =        [fx,fxy,fapp];
    inputMN =           [m,n];
    inputNumIntegral =  [nIntPx,nIntPy];
    %forcedInputs = {Af,Bf,Df};

    [eigVal, eigVals, eigVecs, wbG, outputCell] = criticalStabilityLoad(inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral);
    Pcr = eigVal*fapp*b;
    
    % Additional outputs
    [ABDs, angles, C, outputGauss] = outputCell{:};

    % Prebuckling outputs - Gauss Domain
    [NxxG,NyyG,NxyG,MxxG,MyyG,MxyG,...
     exxG,eyyG,exyG,kxxG,kyyG,kxyG, ...
     umG,vmG,wmG] = outputGauss{:};
   
    % Prebuckling outputs - Uniform Domain (by interp)
    xGauss = ( (a - 0)*fun_lgwt(nIntPx,-1,1) + a + 0 ) * 0.5;
    yGauss = ( (b - 0)*fun_lgwt(nIntPy,-1,1) + b + 0 ) * 0.5;
    [Y,X] = meshgrid(yGauss,xGauss);
    xUni = (a:-a/(nIntPx-1):0)';
    yUni = (b:-b/(nIntPy-1):0)';
    [Yu,Xu] = meshgrid(yUni,xUni);
    Nxx = interp2(Y,X,NxxG,Yu,Xu);
    Nyy = interp2(Y,X,NyyG,Yu,Xu);
    Nxy = interp2(Y,X,NxyG,Yu,Xu);
    exx = interp2(Y,X,exxG,Yu,Xu);
    eyy = interp2(Y,X,eyyG,Yu,Xu);
    exy = interp2(Y,X,exyG,Yu,Xu);
    um = interp2(Y,X,umG,Yu,Xu);
    vm = interp2(Y,X,vmG,Yu,Xu);

    % w from Eigvector - Uniform Domain (by interp)
    wb = interp2(Y,X,wbG,Yu,Xu);

    % % Failure Criteria
    % Xt	= 1500e6;
    % Xc	= 1500e6;
    % Yt	= 40e6;
    % Yc	= 246e6;
    % S12	= 68e6;
    % % Xt	= 493.9e6;
    % % Xc	= 323.9e6;
    % % Yt	= 13.5e6;
    % % Yc	= 20.25e6;
    % % S12	= 35e6;
    % S23	= ((1+Yt/Yc)/(3+5*Yt/Yc) * Yt*Yc)^0.5;
    % [FYt, FYc, FXt, FXc, stress123] = fun_Hashin({E1,E2,nu12,G12}, {Xt,Xc,Yt,Yc,S12}, stressCell, tPly);
    % s11 = cellfun(@(v)v(1),stress123);
    % s22 = cellfun(@(v)v(2),stress123);
    % t12 = cellfun(@(v)v(3),stress123);

%%%% Printing Results
    fprintf("  T0_1= %6.2f", T0v1);                 fprintf("\n");
    fprintf("  T1_1= %6.2f", T1v1);                 fprintf("\n");
    fprintf("  T0_2= %6.2f", T0v2);                 fprintf("\n");
    fprintf("  T1_2= %6.2f", T1v2);                 fprintf("\n");
    fprintf("  Pcr= %9.4f kN ", Pcr/1000);          fprintf("\n");
    % fprintf("  Ncr= %9.4f N/mm ", eigVal*fapp/1000);%fprintf("\n");
    fprintf("  mat= %s ",matNature);                fprintf("\n");
    fprintf("  fx= %3.1f", fx);                     fprintf("\n");
    fprintf("  fxy= %3.1f", fxy);                   fprintf("\n");
    fprintf("  AR= %3.1f", AR);                     fprintf("\n");
    fprintf("  r= %3.1f", r);                       fprintf("\n");
    fprintf("  m= %2d", m);                         fprintf("\n");
    fprintf("  n= %2d", n);                         fprintf("\n");
    fprintf("  nx= %2d", nIntPx);                   fprintf("\n");
    fprintf("  ny= %2d", nIntPy);                   fprintf("\n");
    fprintf("  time= %6.2f s \n", toc - told);
    told = toc;
end
end
end
end
end
end
end
end
end
end
end

% diary off
