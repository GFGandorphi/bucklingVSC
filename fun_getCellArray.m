function cellGeneral = fun_getCellArray(xiFunc, etaFunc)

    nx = size(xiFunc,1);
    m = size(xiFunc,2);
    ny = size(etaFunc,1);
    n = size(etaFunc,2);

    auxID = zeros(1,m*n);
    auxID(1:m:m*n) = 1;
    auxID=cumsum(auxID,2);

    FXI = repmat(xiFunc,1,n*ny);
    fetaAux1 = etaFunc(:,auxID)';
    fetaAux2 = fetaAux1(:)';
    FETA = repmat(fetaAux2,nx,1);
    arrayGeneral = FXI.*FETA;

    cellGeneral = mat2cell(arrayGeneral,ones(1,nx),m*n*ones(1,ny));

end