function [FYt, FYc, FXt, FXc, stress123] = fun_Hashin(inputMaterial, inputAllowable, stressCell, tPly)
% -------------------------------------------------------------------------
% fun_Hashin evaluates the failure indices of Hashin's criteria
% -------------------------------------------------------------------------
%  INPUTS:
%      inputMaterial ----- {E1,E2,nu12,G12}
%      inputAllowable ---- {Xt,Xc,Yt,Yc,S12}
%      stressCell -------- {Nxx,Nyy,Nxy,Mxx,Myy,Mxy,ABDs,angles}
%      tPly -------------- [m]    thickness of single ply
%      Where:
%           E1 ----------- [Pa]   plate modulus of elasticity in dir. 1
%           E2 ----------- [Pa]   plate modulus of elasticity in dir. 2
%           nu12 --------- []     Poisson in direction 12
%           G12 ---------- [Pa]   G modulus in direction 12
%           Xt ----------- [Pa]   tensile stength of dir. 1
%           Xc ----------- [Pa]   compression stength of dir. 1
%           Yt ----------- [Pa]   tensile stength of dir. 2
%           Yc  ---------- [Pa]   compression stength of dir. 2
%           S12  --------- [Pa]   shear stength of dir. 12
%           Nxx ---------- [N/m]  x  membrane flux (nIntPx x nIntPy array)
%           Nyy ---------- [N/m]  y  membrane flux (nIntPx x nIntPy array)
%           Nxy ---------- [N/m]  xy membrane flux (nIntPx x nIntPy array)
%           Mxx ---------- [N]    x  moment flux (nIntPx x nIntPy array)
%           Myy ---------- [N]    y  moment flux (nIntPx x nIntPy array)
%           Mxy ---------- [N]    xy moment flux (nIntPx x nIntPy array)
%           ABDs --------- [+]    ABDs (nIntPx x nIntPy cell with 6x6 each)
%           angles ------- [°]    stack angles (nIntPx x nIntPy cell with nPliesx1 each)
%  OUTPUTS:
%      FYt --------------- []     tensile failure index of dir. 2 
%      FYc --------------- []     compression failure index of dir. 2 
%      FXt --------------- []     tensile failure index of dir. 1 
%      FXc --------------- []     compression failure index of dir. 1
%      stress123 --------- [N/m²] stresses in fiber dir. (s11, s22 and t12)
%      All outputs have size of [nIntPx x nIntPy x nPlies]
%      stress123 are a cell array. Each cell is a 3x1 vector [s11;s22;t12]

%% 00 - INPUTS ASSIGNMENT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%% Material - Plate
    E1 = inputMaterial{1};
    E2 = inputMaterial{2};
    nu12 = inputMaterial{3};
    G12 = inputMaterial{4};

    %%%% Material - Allowable
    Xt	= inputAllowable{1};
    Xc	= inputAllowable{2};
    Yt	= inputAllowable{3};
    Yc	= inputAllowable{4};
    S12	= inputAllowable{5};
    S23	= ((1+Yt/Yc)/(3+5*Yt/Yc) * Yt*Yc)^0.5;

    %%%% Stresses Cells
    Nxx = num2cell(-stressCell{1});
    Nyy = num2cell(-stressCell{2});
    Nxy = num2cell(-stressCell{3});
    Mxx = num2cell(-stressCell{4});
    Myy = num2cell(-stressCell{5});
    Mxy = num2cell(-stressCell{6});

    %%%% Stacking Properties
    ABDs = stressCell{7};
    angles = stressCell{8};

    %%% Aux Values
    nu21 = nu12/E1 * E2;
    Q = [E1/(1-nu12*nu21)           nu12*E2/(1 - nu12*nu21)     0;
         nu12*E2/(1 - nu12*nu21)    E2/(1-nu12*nu21)            0;
         0                          0                         G12];
    nIntPx = size(Nxx,1);
    nIntPy = size(Nxx,2);
    nPlies = length(angles(1,1,:));
    zValues = ((1:nPlies)*2-nPlies-1)*tPly/2;
    zValues = zValues';


%% 01 - LAMINATE DEFORMATIONS AT MEAN SURFACE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    defsXYZm = cellfun(@(c1,c2,c3,c4,c5,c6,c7) ...
                         c1\[c2; c3; c4; c5; c6; c7], ...
                         ABDs, Nxx, Nyy, Nxy, Mxx, Myy, Mxy, ...
                         'UniformOutput', 0);


%% 02 - LAMINA DEFORMATIONS AND STRESSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cosksq = num2cell(cosd(angles).*cosd(angles));
    sinksq = num2cell(sind(angles).*sind(angles));
    cosksink = num2cell(cosd(angles).*sind(angles));

    defsXYZ = cell(nIntPx,nIntPy,nPlies);
    defs123 = cell(nIntPx,nIntPy,nPlies);
    for iPly=1:nPlies
        defsXYZ(:,:,iPly) = cellfun(@(c1) ...
                            [c1(1)+c1(4)*zValues(iPly); ...
                             c1(2)+c1(5)*zValues(iPly); ...
                             c1(3)+c1(6)*zValues(iPly)], defsXYZm, 'UniformOutput', 0);

        defs123(:,:,iPly) = cellfun(@(c1,c2,c3,c4) ...
                                      [c2,      c3,      c4; ...
                                       c3,      c2,     -c4; ...
                                       -2*c4, 2*c4, (c2-c3)] * c1, ...
                                      defsXYZ(:,:,iPly), cosksq(:,:,iPly), sinksq(:,:,iPly), cosksink(:,:,iPly), ...
                                      'UniformOutput', 0);
    end
    stress123 = cellfun(@(c1) Q*c1, defs123, 'UniformOutput', 0);

%% 04 - FAILURE INDICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cellFYt = cellfun(@(c1) (c1(2)/Yt)^2+(c1(3)/S12)^2, stress123, 'UniformOutput', 0);
    cellFYc = cellfun(@(c1) (c1(2)/2/S23)^2+((Yc/2/S23)^2-1)*(c1(2)/Yc)+(c1(3)/S12)^2, stress123, 'UniformOutput', 0);
    cellFXt = cellfun(@(c1) (c1(1)/Xt)^2+(c1(3)/S12)^2, stress123, 'UniformOutput', 0);
    cellFXc = cellfun(@(c1) c1(1)/-Xc, stress123, 'UniformOutput', 0);

    FYt = cell2mat(cellFYt);
    FYc = cell2mat(cellFYc);
    FXt = cell2mat(cellFXt);
    FXc = cell2mat(cellFXc);

end