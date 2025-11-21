function [eigVal, eigVals, eigVecs, wb, outputCell] = criticalStabilityLoad(inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral, forcedInputs)
% -------------------------------------------------------------------------
% criticalStabilityLoad evaluates the linear critical buckling load
% -------------------------------------------------------------------------
%  INPUTS:
%      inputMaterial ----- {E1,E2,nu12,G12,matNature,Eref,Gref}
%      inputGeometry ----- {a,b,tPly,r,T0,T1,isMirrored,tref,href,extraIyy}
%      inputLoads -------- [fx,fxy,fapp]
%      inputMN ----------- [m,n]
%      inputNumIntegral -- [nIntPx,nIntPy]
%      forcedInputs ------ {Af,Bf,Df}
%      Where:
%           E1 ----------- [Pa]  plate modulus of elasticity in dir. 1
%           E2 ----------- [Pa]  plate modulus of elasticity in dir. 2
%           nu12 --------- []    Poisson in direction 12
%           G12 ---------- [Pa]  G modulus in direction 12
%           matNature ---- []    "V" for VAT, "Q" for Quasi-Parallel
%           Eref --------- [Pa]  stringer's modulus of elasticity
%           Gref --------- [Pa]  stringer's G modulus
%           a ------------ [m]   plate length in x
%           b ------------ [m]   plate length in y
%           tPly --------- [m]   thickness of single ply
%           r ------------ [m]   plate radius (in y)
%           T0 ----------- [°]   center angles (vector)
%           T1 ----------- [°]   boundary angles (vector)
%           isMirrored --- []    bool vector to mirror ply. For Q only
%           tref --------- [m]   stringer's thickness 
%           href --------- [m]   stringer's height
%           extraIyy ----- [m^4] extra inertia (due offset)
%           fx ----------- []    load factor of axial component
%           fxy ---------- []    load factor of shear component
%           fapp --------- [N/m] applied load
%           m ------------ []    number of approx. terms in x dir
%           n ------------ []    number of approx. terms in y dir
%           nIntPx ------- []    number of integration points in x
%           nIntPy ------- []    number of integration points in y
%           Af ----------- [N/m] artificial plate A matrix (3x3)
%           Bf ----------- [N]   artificial plate B matrix (3x3)
%           Df ----------- [N.m] artificial plate D matrix (3x3)
%  OUTPUTS:
%      eigVal ------------ []    minimum eigenvalue (Critical Load/Applied Load)
%      eigVals ----------- []    all eigenvalues
%      eigVecs ----------- []    all eigenvectors
%      wb ---------------- []    w buckled shape from eigvector (nIntPx x nIntPy array)
%      outputCell -------- {ABDs,angles,C,outputGauss}
%      Where:
%           ABDs --------- [+]   ABDs (nIntPx x nIntPy cell with 6x6 each)
%           angles ------- [°]   stack angles (nIntPx x nIntPy cell with nPliesx1 each)
%           C ------------ []    displacements constants vector
%           outputGauss -- {Nxx, Nyy, Nxy, Mxx, Myy, Mxy, ...
%                           exx, eyy, exy, kxx, kyy, kxy, ...
%                           um, vm, wm};
%           Nxx ---------- [N/m] x  membrane flux
%           Nyy ---------- [N/m] y  membrane flux
%           Nxy ---------- [N/m] xy membrane flux
%           Mxx ---------- [N]   x  moment flux
%           Myy ---------- [N]   y  moment flux
%           Mxy ---------- [N]   xy moment flux
%           exx ---------- [N/m] x  membrane strain
%           eyy ---------- [N/m] y  membrane strain
%           exy ---------- [N/m] xy membrane strain
%           kxx ---------- [N]   x  bending strain
%           kyy ---------- [N]   y  bending strain
%           kxy ---------- [N]   xy bending strain
%           um ----------- [m]   x prebuckling displacement
%           vm ----------- [m]   y prebuckling displacement
%           wm ----------- [m]   w prebuckling displacement
%           Note: The domain respect Gauss distribution
%           Note: All flux and displacement values are from prebuckling
%           Note: All outputGauss elements are a (nIntPx x nIntPy) array

%% 00 - INPUTS ASSIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Material - Plate
    E1 = inputMaterial{1};          % [Pa] modulus of elasticity in dir. 1
    E2 = inputMaterial{2};          % [Pa] modulus of elasticity in dir. 2
    nu12 = inputMaterial{3};        % []   Poisson in direction 12
    G12 = inputMaterial{4};         % [Pa] G modulus in direction 12
    matNature = inputMaterial{5};   % str (1ch) with plate material nature

    %%%% Material - Stringer
    Eref = inputMaterial{6};        % [Pa] stringer's modulus of elasticity
    Gref = inputMaterial{7};        % [Pa] stringer's G modulus

    %%%% Geometry - Plate
    a = inputGeometry{1};           % [m] length in x
    b = inputGeometry{2};           % [m] length in y
    tPly = inputGeometry{3};        % [m] each ply thickness
    r = inputGeometry{4};           % [m] plate radius
    T0 = inputGeometry{5};          % [°] center angle(s)
    T1 = inputGeometry{6};          % [°] boundary angle(s)
    isMirrored = inputGeometry{7};  % bool to mirror ply. For Q only

    %%%% Geometry - Stringer
    tref = inputGeometry{8};        % [m] stringer's thickness 
    href = inputGeometry{9};        % [m] stringer's height
    extraIyy = inputGeometry{10};   % [m] extra inertia (due offset)
    
    %%%% Load Factors
    fx = inputLoads(1);             % [] load factor of axial component
    fxy = inputLoads(2);            % [] load factor of shear component
    fapp = inputLoads(3);           % [N/m] applied load

    %%%% Bardell Orders and Number of Terms m and n 
    m = inputMN(1);                 % [] number of approx. terms in x dir
    n = inputMN(2);                 % [] number of approx. terms in y dir

    %%%% Numerical Integration
    nIntPx = inputNumIntegral(1);   % [] number of integration points in x
    nIntPy = inputNumIntegral(2);   % [] number of integration points in y


%% 01 - INPUTS EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Plate
    nPlies = length(T0);
    nu21 = nu12/E1 * E2;
    Q = [E1/(1-nu12*nu21)           nu12*E2/(1 - nu12*nu21)     0;
         nu12*E2/(1 - nu12*nu21)    E2/(1-nu12*nu21)            0;
         0                          0                         G12];
   
    %%%% Stringer
    yRef = b/2;
    Aref = tref * href;                 % [m^2]  
    Iyy = tref*href^3/12 + extraIyy;    % [m^4] 
    Izz = href*tref^3/12;               % [m^4] 
    GammaRef = 0;                       % [m^6]
    if href > 0
        Jcte = (1/3 - .21*tref/href*(1-tref^4/12/href^4));
    else
        Jcte = 0;
    end
    Jref = href*tref^3*Jcte;            % [m^4]
    Cref = [Eref*Aref   0           0           0           0;
            0           Eref*Izz    0           0           0;
            0           0           Eref*Iyy    0           0;
            0           0           0           Gref*Jref   0;
            0           0           0           0           Eref*GammaRef];
    
    %%%% Gauss Points
    [xAux,wx] = fun_lgwt(nIntPx,-1,1);
    [yAux,wy] = fun_lgwt(nIntPy,-1,1);
    xGauss = ( (a - 0)*xAux + a + 0 ) * 0.5;
    yGauss = ( (b - 0)*yAux + b + 0 ) * 0.5;
    

%% 02 - ABD MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if upper(matNature) == "Q"
        if size(xGauss,1)>1; xGauss=xGauss'; end
        if size(yGauss,1)>1; yGauss=yGauss'; end
        angles = zeros(nIntPx,nIntPy,nPlies);
        for i=1:nPlies
            xaux = repmat(xGauss, 1, nIntPy);
            iy = zeros(1,nIntPx*nIntPy);
            iy(1:nIntPx:end) = 1;
            iy = cumsum(iy,2);
            yaux = yGauss(iy);
            auxThetaV = fun_getQPinfo(xaux,yaux,a,b,T0(i),T1(i),isMirrored(i));
            if size(auxThetaV,2)>1; auxThetaV=auxThetaV'; end
            angles(:,:,i) = reshape(auxThetaV,nIntPx,[]);
        end
    else
        fact = (xGauss>=a/2)*-2 + 1;
        auxThetaX = fact .* 2/a .*(T0-T1) .* (xGauss - a/2) + T0;
        auxThetaX = reshape(auxThetaX,nIntPx,[],nPlies);
        angles = repmat(auxThetaX,1,nIntPy);
    end

    if exist('forcedInputs','var')
        ABDs = fun_matrixABD(angles, tPly, Q, forcedInputs);
    else
        ABDs = fun_matrixABD(angles, tPly, Q);
    end
    

%% 03 - DISPLACEMENTS DOMAINS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Initialization
    zV = zeros(1,n*m);
    etaRef = 2*yRef/b-1;

    f1pos = fun_dispBardell("f",1);
    f1neg = fun_dispBardell("f",-1);

    fref = fun_dispBardell("f",etaRef);
    dfref = fun_dispBardell("df",etaRef);
    d2fref = fun_dispBardell("d2f",etaRef);
    
    oFSm = [101 5:1:25];
    oSSm = [101 102 5:1:25];
    oSSb = [2 4 5:1:25];
    
    oFSm_x = oFSm(1:m);
    oSSm_x = oSSm(1:m);
    oSSm_y = oSSm(1:n);
    oSSb_x = oSSb(1:m);
    oSSb_y = oSSb(1:n);
    
    cX0 = 1;
    cX1 = 2/a;
    cX2 = 4/a^2;
    cY0 = 1;
    cY1 = 2/b;
    cY2 = 4/b^2;

    NuvwCell = cell(1,2);
    NuvwwCell = cell(1,2);

    % % In the case of different distributions
    % xUniform = a:-a/(nIntPx-1):0;
    % yUniform = b:-b/(nIntPy-1):0;
    % xAux = {xGauss, xUniform'};
    % yAux = {yGauss, yUniform'};
    xAux = {xGauss};
    yAux = {yGauss};

    %%% Domains Evaluation (Uniform and Gauss)
    for iDomain=1:length(xAux)
        xi = 2*xAux{iDomain}/a-1;
        eta = 2*yAux{iDomain}/b-1;
        
        fXi = fun_dispBardell("f",xi);
        fEta = fun_dispBardell("f",eta);
        
        dfXi = fun_dispBardell("df",xi);
        dfEta = fun_dispBardell("df",eta);
        
        d2fXi = fun_dispBardell("d2f",xi);
        d2fEta = fun_dispBardell("d2f",eta);
    
        %%%% Plate Displacements Fields
        Npu    = fun_getCellArray(   fXi(:,oFSm_x)*cX0,   fEta(:,oSSm_y)*cY0 );
        Npu_x  = fun_getCellArray(  dfXi(:,oFSm_x)*cX1,   fEta(:,oSSm_y)*cY0 );
        Npu_y  = fun_getCellArray(   fXi(:,oFSm_x)*cX0,  dfEta(:,oSSm_y)*cY1 );
        Npv    = fun_getCellArray(   fXi(:,oSSm_x)*cX0,   fEta(:,oSSm_y)*cY0 );
        Npv_x  = fun_getCellArray(  dfXi(:,oSSm_x)*cX1,   fEta(:,oSSm_y)*cY0 );
        Npv_y  = fun_getCellArray(   fXi(:,oSSm_x)*cX0,  dfEta(:,oSSm_y)*cY1 );
        Npw    = fun_getCellArray(   fXi(:,oSSb_x)*cX0,   fEta(:,oSSb_y)*cY0 );
        Npw_x  = fun_getCellArray(  dfXi(:,oSSb_x)*cX1,   fEta(:,oSSb_y)*cY0 );
        Npw_y  = fun_getCellArray(   fXi(:,oSSb_x)*cX0,  dfEta(:,oSSb_y)*cY1 );
        Npw_xx = fun_getCellArray( d2fXi(:,oSSb_x)*cX2,   fEta(:,oSSb_y)*cY0 );
        Npw_yy = fun_getCellArray(   fXi(:,oSSb_x)*cX0, d2fEta(:,oSSb_y)*cY2 );
        Npw_xy = fun_getCellArray(  dfXi(:,oSSb_x)*cX1,  dfEta(:,oSSb_y)*cY1 );
    
        %%%% Stringer Displacements Fields
        Nru_x   = fun_getCellArray(  dfXi(:,oFSm_x)*cX1,  fref(:,oSSm_y)*cY0 );
        Nru_xx  = fun_getCellArray( d2fXi(:,oFSm_x)*cX2,  fref(:,oSSm_y)*cY0 );
        Nrw_x   = fun_getCellArray(  dfXi(:,oSSb_x)*cX1,  fref(:,oSSb_y)*cY0 );
        Nrw_xx  = fun_getCellArray( d2fXi(:,oSSb_x)*cX2,  fref(:,oSSb_y)*cY0 );
        Nrw_xy  = fun_getCellArray(  dfXi(:,oSSb_x)*cX1, dfref(:,oSSb_y)*cY1 );
        Nrw_xyx = fun_getCellArray( d2fXi(:,oSSb_x)*cX2, dfref(:,oSSb_y)*cY1 );
        % for KG
        Nru_y   = fun_getCellArray(  fXi(:,oFSm_x)*cX0,  dfref(:,oSSm_y)*cY1 );
        Nrv_x   = fun_getCellArray( dfXi(:,oSSm_x)*cX1,   fref(:,oSSm_y)*cY0 );
        Nrv_y   = fun_getCellArray(  fXi(:,oSSm_x)*cX0,  dfref(:,oSSm_y)*cY1 );
        Nrw     = fun_getCellArray(  fXi(:,oSSb_x)*cX0,   fref(:,oSSb_y)*cY0 );
        Nrw_yy  = fun_getCellArray(  fXi(:,oSSb_x)*cX0, d2fref(:,oSSb_y)*cY2 );
        
        %%%% Edges Displacements Fields
        Npu_x0 = fun_getCellArray( f1neg(:,oFSm_x),  fEta(:,oSSm_y) );
        Npv_x0 = fun_getCellArray( f1neg(:,oSSm_x),  fEta(:,oSSm_y) );
        Npv_xa = fun_getCellArray( f1pos(:,oSSm_x),  fEta(:,oSSm_y) );
        Npu_y0 = fun_getCellArray(   fXi(:,oFSm_x), f1neg(:,oSSm_y) );
        Npu_yb = fun_getCellArray(   fXi(:,oFSm_x), f1pos(:,oSSm_y) );
    
        %%% Store values for plotting and for evaluation
        NuvwCell{iDomain} = cellfun(@(c1,c2,c3) ...
                            [c1 zV zV; ...
                             zV c2 zV; ...
                             zV zV c3], ...
                             Npu, Npv, Npw, 'UniformOutput', false);

        NuvwwCell{iDomain} = cellfun(@(c1,c2,c3,c4,c5,c6,c7,c8) ...
                             [c1            zV               zV; ...
                              zV            c2         (1/r)*c3; ...
                              c4            c5               zV; ...
                              zV            zV              -c6; ...
                              zV            c2*(1/r)        -c7; ...
                              c4*(-0.5/r)   c5*(1.5/r)    -2*c8], ...
                              Npu_x, Npv_y, Npw, Npu_y, Npv_x, Npw_xx, ...
                              Npw_yy, Npw_xy, 'UniformOutput', false);
    end
    Nuvww = NuvwwCell{1};

    %%%% For Stringer 
    Btio = cellfun(@(c1,c2,c3,c4,c5) ...
                     [c1 zV zV; ...
                      zV c2 zV; ...
                      zV zV c3; ...
                      zV zV c4; ...
                      zV zV c5], ...
                      Nru_x, Nru_xx, Nrw_xx, Nrw_xy, Nrw_xyx, ... 
                      'UniformOutput', false);

    
%% 04 - INTEGRATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Initialization
    if size(wx,2)>1; wx=wx'; end
    if size(wy,1)>1; wy=wy'; end
    intLimitsXY = (a - 0) * 0.5 * (b - 0) * 0.5;
    intLimitsX = (a - 0) * 0.5;
    intLimitsY = (b - 0) * 0.5;
    WW = wx*wy;
    WWc = num2cell(WW);
    igV = repmat(1:nIntPx,1,nIntPy);
    jgV = zeros(1,nIntPx*nIntPy);
    jgV(1:nIntPx:end) = 1;
    jgV = cumsum(jgV,2);

    %%%% Plate Stiffness
    dK = zeros(3*n*m,3*n*m);
    for I = 1:nIntPy*nIntPx
        ig = igV(I);
        jg = jgV(I);
        dK = dK + Nuvww{ig,jg}'*ABDs{ig,jg}*Nuvww{ig,jg}*WWc{ig,jg};
    end
    K = intLimitsXY * dK;

    %%%% Stringer Stiffness
    dKref = zeros(3*n*m,3*n*m);
    for ig = 1:nIntPx
        dKref = dKref + wx(ig)*Btio{ig}'*Cref*Btio{ig};
    end
    Kref = intLimitsX * dKref;

    %%%% Total Stiffness
    Ktotal = Kref + K;

    %%%% Edges Loads
    % y = 0, x = 0 to a
    dF = zeros(3*n*m,1);
    zVC = zeros(3*n*m,1);
    for ig = 1:nIntPx
        zVC(1:n*m) = wx(ig)*Npu_y0{ig}'*fxy*fapp;
        dF = dF + zVC;
    end
    FedgeY0 = intLimitsX * dF;

    % y = b, x = 0 to a
    dF = zeros(3*n*m,1);
    zVC = zeros(3*n*m,1);
    for ig = 1:nIntPx
        zVC(1:n*m) = wx(ig)*Npu_yb{ig}'*fxy*fapp;
        dF = dF + zVC;
    end
    FedgeYb = intLimitsX * dF;

    % x = 0, y = 0 to b --> X load applied
    dF = zeros(3*n*m,1);
    zVC = zeros(3*n*m,1);
    for jg = 1:nIntPy
        zVC(1:2*n*m) = wy(jg)*[Npu_x0{jg}'*fx*fapp; Npv_x0{jg}'*fxy*fapp];
        dF = dF + zVC;
    end
    FedgeX0 = intLimitsY * dF;

    % x = a, y = 0 to b --> Simply Supported on X
    dF = zeros(3*n*m,1);
    zVC = zeros(3*n*m,1);
    for jg = 1:nIntPy
        zVC(n*m+1:2*n*m) = wy(jg)*Npv_xa{jg}'*fxy*fapp;
        dF = dF + zVC;
    end
    FedgeXa = intLimitsY * dF;

    %%%% Total Edges Loads Vector
    F = (FedgeXa - FedgeX0) + (FedgeYb - FedgeY0);

    %%%% Displacements Constants Evaluation
    C = Ktotal\F;

    %%%% Prebuckling Flux Fields
    fluxes = cellfun(@(c1,c2) c2*c1*C, Nuvww, ABDs, 'UniformOutput', 0);
    fluxesA = real(cell2mat(fluxes));
    Nxx = fluxesA(1:6:end,:);
    Nyy = fluxesA(2:6:end,:);
    Nxy = fluxesA(3:6:end,:);

    %%%% Plate Geometry Stiffness
    % KGxP
    zV3 = zeros(3*n*m,3*n*m);
    dKGx = zV3;
    for I = 1:nIntPy*nIntPx
        ig = igV(I);
        jg = jgV(I);
        zV3(2*n*m+1:end, 2*n*m+1:end) = Npw_x{ig,jg}'*Nxx(ig,jg)*Npw_x{ig,jg}*WWc{ig,jg};
        dKGx = dKGx + zV3;
    end
    KGxP = intLimitsXY * dKGx;

    % KGyP
    zV3 = zeros(3*n*m,3*n*m);
    dKGy = zV3;
    for I = 1:nIntPy*nIntPx
        ig = igV(I);
        jg = jgV(I);
        zV3(2*n*m+1:end, 2*n*m+1:end) = Npw_y{ig,jg}'*Nyy(ig,jg)*Npw_y{ig,jg}*WWc{ig,jg};
        dKGy = dKGy + zV3;
    end
    KGyP = intLimitsXY * dKGy;

    % KGxyP
    zV3 = zeros(3*n*m,3*n*m);
    dKGxy = zV3;
    for I = 1:nIntPy*nIntPx
        ig = igV(I);
        jg = jgV(I);
        zV3(2*n*m+1:end, 2*n*m+1:end) = Npw_x{ig,jg}'*Nxy(ig,jg)*Npw_y{ig,jg}*WWc{ig,jg};
        dKGxy = dKGxy + zV3;
    end
    KGxyP = intLimitsXY * dKGxy;
	
    %%%% Stringer Geometry Stiffness
    % KGxR
    defsRc = cellfun(@(c1,c2,c3,c4,c5,c6,c7,c8) ...
                      [c1            zV               zV; ...
                       zV            c2         (1/r)*c3; ...
                       c4            c5               zV; ...
                       zV            zV              -c6; ...
                       zV            c2*(1/r)        -c7; ...
                       c4*(-0.5/r)   c5*(1.5/r)    -2*c8]*C, ...
                       Nru_x, Nrv_y, Nrw, Nru_y, Nrv_x, Nrw_xx, ...
                       Nrw_yy, Nrw_xy, 'UniformOutput', false);
    defsR = cell2mat(defsRc);
    exx = defsR(1:6:end)*Eref*Aref;
    zV3 = zeros(3*n*m,3*n*m);
    dKGxR = zV3;
    for ig = 1:nIntPx
        zV3(2*n*m+1:end, 2*n*m+1:end) = wx(ig)*Nrw_x{ig}'*exx(ig)*Nrw_x{ig};
        dKGxR = dKGxR + zV3;
    end
    KGxR = intLimitsX * dKGxR;

    
    %%%% Total Geometry Stiffness
    KGx = KGxP + KGxR;
    KGy = KGyP;
    KGxy = KGxyP;

    KG = KGx + KGy + 2*KGxy;% + Ksanders;

%% 05 - PREBUCKLING OUTPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    outputAux = cell(1,2);

    for iDomain=1:length(xAux)

        NuvwwOutput = NuvwwCell{iDomain};

        % Flux Fields
        fluxes = cellfun(@(c1,c2) c2*c1*C, NuvwwOutput, ABDs, 'UniformOutput', 0);
        fluxesA = real(cell2mat(fluxes));
        Nxx = fluxesA(1:6:end,:);
        Nyy = fluxesA(2:6:end,:);
        Nxy = fluxesA(3:6:end,:);
        Mxx = fluxesA(4:6:end,:);
        Myy = fluxesA(5:6:end,:);
        Mxy = fluxesA(6:6:end,:);
    
        % Strain Fields
        strains = cellfun(@(c1) c1*C, NuvwwOutput, 'UniformOutput', 0);
        strainsA = real(cell2mat(strains));
        exx = strainsA(1:6:end,:);
        eyy = strainsA(2:6:end,:);
        exy = strainsA(3:6:end,:);
        kxx = strainsA(4:6:end,:);
        kyy = strainsA(5:6:end,:);
        kxy = strainsA(6:6:end,:);

        % Displacement Fields
        NuvwOutput = NuvwCell{iDomain};
        displacements = cellfun(@(c1) c1*C, NuvwOutput, 'UniformOutput', 0);
        displacementsA = real(cell2mat(displacements));
        um = displacementsA(1:3:end,:);
        vm = displacementsA(2:3:end,:);
        wm = displacementsA(3:3:end,:);

        % Ouput Cell
        outputAux{iDomain} = {Nxx, Nyy, Nxy, Mxx, Myy, Mxy, ...
                              exx, eyy, exy, kxx, kyy, kxy, ...
                              um, vm, wm};
    end
    % % In the case of different distributions
    % outputUniform = outputAux{2};
    % outputGauss = outputAux{1};
    outputCell = {ABDs,angles,C,outputAux{1}};

%% 06 - CRITICAL LOAD EVALUATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [eigVecs, eigAux] = eig(Ktotal,KG);
    eigVals = diag(eigAux);
    [eigVal, ~] = min(eigVals(eigVals>=0));

    eigVec = eigVecs(:,eigVals==eigVal);
    cAux = eigVec(m*n*2+1:m*n*3);
    wb = real(cell2mat(cellfun(@(c1) c1*cAux, Npw, 'UniformOutput', 0)));

end