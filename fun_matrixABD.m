function ABDs = fun_matrixABD(stackingAngle, plyThickness, Q, forcedInputs)
        
    nx = size(stackingAngle,1);
    ny = size(stackingAngle,2);

    nPlies = size(stackingAngle,3);
    totalThickness = nPlies * plyThickness;
    z0 = -totalThickness / 2;

    zk_1 = z0 + plyThickness * ((1:nPlies) - 1);
    zk = zk_1 + plyThickness;

    Q11 = Q(1,1);
    Q12 = Q(1,2);
    Q22 = Q(2,2);
    Q66 = Q(3,3);
    
    cosk = cosd(stackingAngle);
    sink = sind(stackingAngle);

    Qxx = Q11*cosk.^4 + Q22*sink.^4 + 2*(Q12 + 2*Q66)*cosk.^2.*sink.^2;
    Qyy = Q11*sink.^4 + Q22*cosk.^4 + 2*(Q12 + 2*Q66)*cosk.^2.*sink.^2;
    Qxy = (Q11 + Q22 - 4*Q66)*cosk.^2.*sink.^2 + Q12*(cosk.^4 + sink.^4);
    Qss = (Q11 + Q22 - 2*Q12)*cosk.^2.*sink.^2 + Q66*(cosk.^2 - sink.^2).^2;
    Qxs = (Q11 - Q12)*cosk.^3.*sink + (Q12 - Q22)*cosk.*sink.^3 - 2*Q66*cosk.*sink.*(cosk.^2 - sink.^2);
    Qys = (Q11 - Q12)*cosk.*sink.^3 + (Q12 - Q22)*cosk.^3.*sink + 2*Q66*cosk.*sink.*(cosk.^2 - sink.^2);

    QT=zeros(nx, ny, nPlies, 3, 3);

    QT(:,:,:,1,1) = Qxx;
    QT(:,:,:,2,2) = Qyy;
    QT(:,:,:,3,3) = Qss;

    QT(:,:,:,1,2) = Qxy;
    QT(:,:,:,2,1) = Qxy;

    QT(:,:,:,1,3) = Qxs;
    QT(:,:,:,3,1) = Qxs;

    QT(:,:,:,2,3) = Qys;
    QT(:,:,:,3,2) = Qys;

    zkA = zk - zk_1;
    Aaux = QT .* reshape(zkA,1,1,[]);
    A = reshape(sum(Aaux, 3),nx, ny, 3, 3);

    zkB = 0.5*(zk.^2 - zk_1.^2);
    Baux = QT .* reshape(zkB,1,1,[]);
    B = reshape(sum(Baux, 3),nx, ny, 3, 3);

    zkD = (1/3)*(zk.^3 - zk_1.^3);
    Daux = QT .* reshape(zkD,1,1,[]);
    D = reshape(sum(Daux, 3),nx, ny, 3, 3);

    ABDzao = zeros(6*nx,6*ny);

    for j=1:3
        for i=1:3
            ABDzao(i:6:6*nx, j:6:6*ny) = A(:,:,i,j);
            ABDzao(i:6:6*nx, j+3:6:6*ny) = B(:,:,i,j);
            ABDzao(i+3:6:6*nx, j:6:6*ny) = B(:,:,i,j);
            ABDzao(i+3:6:6*nx, j+3:6:6*ny) = D(:,:,i,j);
        end
    end

    if exist('forcedInputs','var')
        ignoreFlag = -999;
        if ~isempty(forcedInputs{1})
            Af = forcedInputs{1};
            [ifA,jfA]=find(Af~=ignoreFlag);
            for iA = 1:length(ifA)
                ABDzao(ifA(iA):6:6*nx, jfA(iA):6:6*ny) = Af(ifA(iA),jfA(iA));
            end
        end
        if ~isempty(forcedInputs{2})
            Bf = forcedInputs{2};
            [ifB,jfB]=find(Bf~=ignoreFlag);
            for iB = 1:length(ifB)
                ABDzao(ifB(iB):6:6*nx, jfB(iB)+3:6:6*ny) = Bf(ifB(iB),jfB(iB));
                ABDzao(ifB(iB)+3:6:6*nx, jfB(iB):6:6*ny) = Bf(ifB(iB),jfB(iB));
            end
        end
        if ~isempty(forcedInputs{3})
            Df = forcedInputs{3};
            [ifD,jfD]=find(Df~=ignoreFlag);
            for iD = 1:length(ifD)
                ABDzao(ifD(iD)+3:6:6*nx, jfD(iD)+3:6:6*ny) = Df(ifD(iD),jfD(iD));
            end
        end
    end

    ABDs = mat2cell(ABDzao, 6*ones(1,nx), 6*ones(1,ny));

end