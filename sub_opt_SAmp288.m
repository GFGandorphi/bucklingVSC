clear
close all
warning('off','all')
clc
tic

% diary('_20240523_logOptMP288_30Q02_30V02_10Q01_10V01_.txt')

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
    lwBound = -84;
    upBound = 84;
    passo = 3;
    Tinit = 0.150;
    Tlow = 0.045;
    nAngles = 2;%4;
    nRuns = 3;
    nCases = 4;
    LK = 30;
    L = 15;%60;          


%% Code

    Tvec = zeros(nRuns*nCases,LK*L);
    anglesVec = cell(nRuns*nCases,LK*L);
    phiVec = 0*ones(nRuns*nCases,LK*L);
    
    anglesAstVec = cell(nRuns*nCases,LK*L);
    phiAstVec = 0*ones(nRuns*nCases,LK*L);
    
    auxVec = zeros(nRuns*nCases,nAngles+2);
    fprintf("   ")
    
    for j=1:nRuns*nCases
        if j>3*nRuns
            % 10V01 - T1.11
            jAux = j-3*nRuns;
            matNature = "V";
            AR = 1;
            fxy = 0;
            PcrRef = 17186.2;
        elseif j>2*nRuns
            % 10Q01 - T1.16
            jAux = j-2*nRuns;
            matNature = "Q";
            AR = 1;
            fxy = 0;
            PcrRef = 17186.2;
        elseif j>nRuns
            % 30V02 - T1.13p
            jAux = j-nRuns;
            matNature = "V";
            AR = 2;
            fxy = 0.5;
            PcrRef = 15866.5;
        else
            % 30Q02 - T1.14p
            jAux = j;
            matNature = "Q";
            AR = 2;
            fxy = 0.5;
            PcrRef = 15866.5;
        end
        a = b*AR;   
        inputMaterial{5} = matNature;
        inputGeometry{1} = a;
        inputLoads(2) = fxy;

        rng(jAux)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 1
        T = Tinit;
        temperatureIsHIGH = T > Tlow;
        anglesAst = generateRandomAll(nAngles, lwBound, upBound, passo);
        if nAngles > 2
            auxAngles = anglesAst;
        else
            auxAngles = [anglesAst anglesAst];
        end
        auxCr = criticalStabilityLoad_forSA(auxAngles, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral);
        phiAst = auxCr(2);
        cont = 0;
        phiBest = 0;
        nPassos0 = LK/2 + 1;
    
        figure(1)
        p2 = plot(1:LK*L, phiAstVec(j,:));
        xlim([1 LK*L]);
        ylim([-1.4 -0.2]);
    
        for K = 1:LK
            nPassos = ceil(nPassos0 - K*0.5);
            for k = 1:L
                cont = cont + 1;
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 2
                if temperatureIsHIGH
                    angles_k = generateRandomNeighbour(anglesAst, nAngles, lwBound, upBound, passo, nPassos);
                else
                    subPasso = 0;
                    if K == (LK-1); subPasso = 1;
                    elseif K == LK; subPasso = 2; end
                    angles_k = generateRandomNeighbour(anglesBest, nAngles, lwBound, upBound, passo-subPasso, nPassos);
                end
    
                if nAngles > 2
                    auxAngles = angles_k;
                else
                    auxAngles = [angles_k angles_k];
                end
                auxCr = criticalStabilityLoad_forSA(auxAngles, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral);
                phi_k = auxCr(2);
                delta = phi_k - phiAst;
            
                phiVec(j,cont) = phi_k;
                Tvec(j,cont) = T;
                anglesVec{j,cont} = angles_k;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 3
                if delta < 0
                    anglesAst = angles_k;
                    phiAst = phi_k;
                    if phi_k < phiBest
                        anglesBest = angles_k;
                        phiBest = phi_k;
                    end
                else
                    prob = exp(-delta/T)*1.0;
                    z = rand(); 
                    if z < prob
                        anglesAst = angles_k;
                        phiAst = phi_k;
                    end
                end
    
                anglesAstVec{j,cont} = anglesAst;
                phiAstVec(j,cont) = phiAst;
                if nAngles > 2
                    fprintf("\n  K=%3d || k=%3d || %6.2f ||  %10.4f %6.1f %6.1f %6.1f %6.1f %10.4f ", K, k, T, phiAst, anglesAst, -phiAst*PcrRef);
                else
                    fprintf("\n  K=%3d || k=%3d || %6.2f ||  %10.4f %6.1f %6.1f %10.4f ", K, k, T, phiAst, anglesAst, -phiAst*PcrRef);
                end
                p2.YData = phiAstVec(j,:);
                drawnow
            end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 4
            anglesAst = anglesBest;
            phiAst = phiBest;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Step 5
            rT = 0.95;
            T = rT*T;
            temperatureIsHIGH = T > Tlow;
        end
    
        if nAngles > 2
            fprintf("\n mat=%s  run=%3d  %10.4f %6.1f %6.1f %6.1f %6.1f %10.4f", matNature, j, phiBest, anglesBest, -phiBest*PcrRef);
        else
            fprintf("\n mat=%s  run=%3d  %10.4f %6.1f %6.1f %10.4f", matNature, j, phiBest, anglesBest, -phiBest*PcrRef);
        end
        fprintf("\n   ");
    end
    
    toc
    % diary off

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Aux funcitons

function angles_k = generateRandomAll(nAngles, lwBound, upBound, passo)
    angles_k = zeros(1, nAngles);
    for i = 1:nAngles
        angles_k(i) = lwBound + round(rand() * (upBound - lwBound)/passo )*passo;
    end
end

function angles_k = generateRandomNeighbour(anglesAst, nAngles, lwBound, upBound, passo, nPassos)
    angles_k = anglesAst;
    for i = 1:nAngles
        if angles_k(i)-nPassos*passo < lwBound
            lb = lwBound;
            ub = 2*nPassos*passo+lb;
        elseif angles_k(i)+nPassos*passo > upBound
            ub = upBound;
            lb = ub-2*nPassos*passo;
        else
            lb = angles_k(i)-nPassos*passo;
            ub = angles_k(i)+nPassos*passo;
        end
        angles_k(i) = lb + round(rand() * (ub - lb)/passo )*passo;
    end
end
