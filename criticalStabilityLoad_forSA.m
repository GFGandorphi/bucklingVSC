function crVec = criticalStabilityLoad_forSA(dv, PcrRef, inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral)
  
    b = inputGeometry{2};
    fapp = inputLoads(3);
    matNature = inputMaterial{5};

    if upper(matNature) == "V"
        T0 = [dv(1) -dv(3) -dv(3) dv(1) dv(1) -dv(3) -dv(3) dv(1) dv(1) -dv(3) -dv(3) dv(1) dv(1) -dv(3) -dv(3) dv(1)];
        T1 = [dv(2) -dv(4) -dv(4) dv(2) dv(2) -dv(4) -dv(4) dv(2) dv(2) -dv(4) -dv(4) dv(2) dv(2) -dv(4) -dv(4) dv(2)];
        isMirrored = [];
    elseif upper(matNature) == "Q"
        T0 = [dv(1) dv(3) dv(3) dv(1) dv(1) dv(3) dv(3) dv(1) dv(1) dv(3) dv(3) dv(1) dv(1) dv(3) dv(3) dv(1)];
        T1 = [dv(2) dv(4) dv(4) dv(2) dv(2) dv(4) dv(4) dv(2) dv(2) dv(4) dv(4) dv(2) dv(2) dv(4) dv(4) dv(2)];
        isMirrored = [0 1 1 0 0 1 1 0 0 1 1 0 0 1 1 0];
    else
        T0 = [dv(1) -dv(2) -dv(2) dv(1) dv(1) -dv(2) -dv(2) dv(1) dv(1) -dv(2) -dv(2) dv(1) dv(1) -dv(2) -dv(2) dv(1)];
        T1 = T0;
        isMirrored = [];
    end

    inputGeometry{5} = T0;
    inputGeometry{6} = T1;
    inputGeometry{7} = isMirrored;

    eigVal = criticalStabilityLoad(inputMaterial, inputGeometry, inputLoads, inputMN, inputNumIntegral);
    eigVal = real(eigVal);

    Pcr = eigVal*fapp*b;
    phi = -Pcr/PcrRef;

    crVec = [Pcr, phi];
end