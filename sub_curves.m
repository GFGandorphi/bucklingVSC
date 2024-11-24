clc
close all
clear

%% Inputs

    a = 0.040;
    b = 0.040;
    % T0 = randi([-85 85]);
    % T1 = randi([-85 85]);

    % Options
    QPoff = 1;
    VAToff = 2;
    QPfull = 3;
    VATfull = 4;
    QPref = 5;
    VATref = 6;
    QPpoint = 7;
    VATpoint = 8;
    VATonQPref = 9;
    joices = 2;

    % For angles
    ds = 0.0008;
    lsize = 1.1;
    nSepA = 11;
    trA = 0.85;

    % For offsets
    nSepO = 500;   % 200 para 2x1cm, 800 para 10x5cm, 1500 para 20x10cm
    nSepTH = 400;  % 300 para 2x1cm, 1000 para 10x5cm, 2500 para 20x10cm
    dNozzle = 1.2e-3;
    nOffset = 22;
    refino = 1.5;
    trO = 0.5;
    vatFac = 3.5;

    % Others
    nrPlots = 1;
    ncPlots = 1;

%% Initial Evaluation
    
    xDom = 0:a/(nSepA-1):a;
    yDom = 0:b/(nSepA-1):b;
   
    esp = a/nSepO;
    % xvec = [0:esp/refino:a/10 a/10:esp:a*(9/10)  a*(9/10):esp/refino:a];
    % yvec = [0:esp/refino:b/10 b/10:esp:b*(9/10)  b*(9/10):esp/refino:b];
    xvec = 0:esp/refino:a;
    yvec = 0:esp/refino:b;
    nx = length(xvec);
    ny = length(yvec);
    dd = a/nSepTH;

    hFigure = figure;
    t = tiledlayout(nrPlots,ncPlots);

    isMirrored = 0;
for T0 = 75%randi([-87 87]) %[-85 -75 -60:7.5:60 75 85]
for T1 = 42%randi([-87 87])%[-85 -40 0 40 85] %[-85 -75 -60:7.5:60 75 85]
    
    auxPlot = nexttile;

    xpoint = a*rand();
    ypoint = b*rand();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 01: Offsets curves of QP

if sum(joices==QPoff) > 0

    % Reference Curve
    [~, ~, yRef] = fun_getQPinfo(xvec,-1,a,b,T0,T1,isMirrored);
    for i = 1:nx
        hVec = fun_circulo(yRef(i),xvec(i),dNozzle/2,'none',[0.4 0.4 0.4 1]);
        % hVec = fun_circulo(yRef(i),xvec(i),dNozzle/2,'none',[0 0 0 1]);
    end
    hold on

    % Offsets
    xAns = repmat(xvec',1,ny);
    yAns = repmat(yvec,nx,1);
    xAns = xAns(:);
    yAns = yAns(:);
    [~, ~, ~, ~, dist] = fun_getQPinfo(xAns,yAns,a,b,T0,T1,isMirrored);
    
    change = -1;
    changeBool = 0;
    for d = dNozzle:dNozzle:dNozzle*nOffset
        change=change*-1;
        changeBool = changeBool + change;
    
        checkDist = (dist >= d-dd) & (dist <= d+dd);
        idxD = find(checkDist);
        xAnsV = xAns(idxD);
        yAnsV = yAns(idxD);
    
        if changeBool; fC=[0.3010 0.7450 0.9330 trO]; else; fC=[0 0.4470 0.7410 trO]; end
        % if changeBool; fC=[0.3010 0.7450 0.9330 trO]; else; fC=[0.4 0.4 0.4 1]; end
        % if changeBool; fC=[1 1 1 1]; else; fC=[0 0 0 1]; end
        % fC=[1 1 1 1];
        
        for i = 1:length(xAnsV)
            hVec = fun_circulo(yAnsV(i),xAnsV(i),dNozzle/2,'none',fC);
        end
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 02: Offsets curves of VAT

if sum(joices==VAToff) > 0

    [~,yref] = fun_getVATinfo(xvec,a,b,T0,T1,0);
    for ix = 1:nx
        % hVec = fun_circulo(yref(ix),xvec(ix),dNozzle/2,'none',[0.4 0.4 0.4 1]);
        % hVec = fun_circulo(yref(ix),xvec(ix),dNozzle/2,'none',[0 0 0 1]);
        hVec = fun_circulo(yref(ix),xvec(ix),dNozzle/2,'none',[0.4660 0.6740 0.1880 1]);
    end
    hold on

    xAnsV = xvec;
    change = -1;
    changeBool = 0;
    for d = vatFac*dNozzle:vatFac*dNozzle:dNozzle*vatFac*nOffset
        change=change*-1;
        changeBool = changeBool + change;
        
        if changeBool; fC=[0.1 0.6 0.2 trO]; else; fC=[0.4660 0.6740 0.1880 trO]; end
        % if changeBool; fC=[1 1 1 trO]; else; fC=[0.4660 0.6740 0.1880 trO]; end
        % if changeBool; fC=[1 1 1 1]; else; fC=[0 0 0 1]; end
        % fC=[0 0 0 1];
        
        yAnsVp = yref+d;
        yAnsVm = yref-d;
        for i = 1:length(xAnsV)
            hVec1 = fun_circulo(yAnsVp(i),xAnsV(i),dNozzle/2,'none',fC);
            hVec2 = fun_circulo(yAnsVm(i),xAnsV(i),dNozzle/2,'none',fC);
        end
        hold on
        set(gca,'xtick',[])
        set(gca,'ytick',[])
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 03: Angles of QP

if sum(joices==QPfull) > 0

    xAux = repmat(xDom,1,nSepA);
    yAux = repmat(yDom,nSepA,1);
    yAux = yAux(:)';

    auxTheta = fun_getQPinfo(xAux,yAux,a,b,T0,T1,isMirrored);
    % auxAngles = [-19.66,-25.82,-30.71,-34.64,-37.85,-40.51,-42.77,-44.70,-46.38,-47.85,-49.16,-50.34,-51.40,-52.36,-53.23,-54.04,-54.78,-55.47,-56.11,-56.71,-22.37,-28.62,-33.39,-37.12,-40.12,-42.60,-44.68,-46.46,-48.01,-49.37,-50.58,-51.66,-52.64,-53.53,-54.34,-55.09,-55.78,-56.42,-57.02,-57.58,-25.58,-31.71,-36.21,-39.65,-42.40,-44.66,-46.55,-48.18,-49.59,-50.83,-51.94,-52.93,-53.83,-54.65,-55.41,-56.10,-56.74,-57.34,-57.89,-58.41,-29.27,-34.99,-39.07,-42.16,-44.63,-46.66,-48.36,-49.83,-51.11,-52.24,-53.25,-54.16,-54.98,-55.74,-56.43,-57.07,-57.66,-58.22,-58.73,-59.22,-33.26,-38.31,-41.88,-44.60,-46.77,-48.57,-50.09,-51.41,-52.56,-53.58,-54.50,-55.33,-56.08,-56.77,-57.41,-58.00,-58.55,-59.06,-59.54,-59.99,-37.27,-41.52,-44.56,-46.91,-48.80,-50.38,-51.73,-52.90,-53.93,-54.86,-55.69,-56.44,-57.13,-57.76,-58.34,-58.89,-59.39,-59.87,-59.54,-58.82,-41.05,-44.52,-47.06,-49.06,-50.70,-52.08,-53.27,-54.31,-55.23,-56.06,-56.81,-57.49,-58.12,-58.70,-59.23,-59.73,-59.73,-59.07,-58.32,-57.48,-44.46,-47.24,-49.35,-51.05,-52.45,-53.66,-54.70,-55.63,-56.45,-57.19,-57.87,-58.49,-59.06,-59.59,-59.90,-59.29,-58.60,-57.84,-56.97,-55.97,-47.45,-49.68,-51.43,-52.86,-54.07,-55.12,-56.04,-56.85,-57.59,-58.26,-58.87,-59.43,-59.95,-59.48,-58.85,-58.15,-57.37,-56.48,-55.46,-54.26,-50.05,-51.84,-53.29,-54.51,-55.55,-56.46,-57.27,-58.00,-58.66,-59.26,-59.81,-59.65,-59.06,-58.42,-57.70,-56.91,-56.01,-54.97,-53.76,-52.30,-52.30,-53.76,-54.97,-56.01,-56.91,-57.70,-58.42,-59.06,-59.65,-59.81,-59.26,-58.66,-58.00,-57.27,-56.46,-55.55,-54.51,-53.29,-51.84,-50.05,-54.26,-55.46,-56.48,-57.37,-58.15,-58.85,-59.48,-59.95,-59.43,-58.87,-58.26,-57.59,-56.85,-56.04,-55.12,-54.07,-52.86,-51.43,-49.68,-47.45,-55.97,-56.97,-57.84,-58.60,-59.29,-59.90,-59.59,-59.06,-58.49,-57.87,-57.19,-56.45,-55.63,-54.70,-53.66,-52.45,-51.05,-49.35,-47.24,-44.46,-57.48,-58.32,-59.07,-59.73,-59.73,-59.23,-58.70,-58.12,-57.49,-56.81,-56.06,-55.23,-54.31,-53.27,-52.08,-50.70,-49.06,-47.06,-44.52,-41.05,-58.82,-59.54,-59.87,-59.39,-58.89,-58.34,-57.76,-57.13,-56.44,-55.69,-54.86,-53.93,-52.90,-51.73,-50.38,-48.80,-46.91,-44.56,-41.52,-37.27,-59.99,-59.54,-59.06,-58.55,-58.00,-57.41,-56.77,-56.08,-55.33,-54.50,-53.58,-52.56,-51.41,-50.09,-48.57,-46.77,-44.60,-41.88,-38.31,-33.26,-59.22,-58.73,-58.22,-57.66,-57.07,-56.43,-55.74,-54.98,-54.16,-53.25,-52.24,-51.11,-49.83,-48.36,-46.66,-44.63,-42.16,-39.07,-34.99,-29.27,-58.41,-57.89,-57.34,-56.74,-56.10,-55.41,-54.65,-53.83,-52.93,-51.94,-50.83,-49.59,-48.18,-46.55,-44.66,-42.40,-39.65,-36.21,-31.71,-25.58,-57.58,-57.02,-56.42,-55.78,-55.09,-54.34,-53.53,-52.64,-51.66,-50.58,-49.37,-48.01,-46.46,-44.68,-42.60,-40.12,-37.12,-33.39,-28.62,-22.37,-56.71,-56.11,-55.47,-54.78,-54.04,-53.23,-52.36,-51.40,-50.34,-49.16,-47.85,-46.38,-44.70,-42.77,-40.51,-37.85,-34.64,-30.71,-25.82,-19.6];
    % auxTheta = auxAngles;

    dy = sind(auxTheta)*ds;       
    dx = cosd(auxTheta)*ds;

    xPar(1,:) = xAux - dx;
    yPar(1,:) = yAux - dy;
    
    xPar(2,:) = xAux + dx;
    yPar(2,:) = yAux + dy;

    plot(yPar,xPar,'Color',[0.8500 0.3250 0.0980 trA],'LineWidth',lsize)
    hold on
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 04: Angles of VAT

if sum(joices==VATfull) > 0

    xAux = repmat(xDom,1,nSepA);
    yAux = repmat(yDom,nSepA,1);
    yAux = yAux(:)';

    auxTheta = fun_getVATinfo(xAux,a,b,T0,T1,0);

    dy = sind(auxTheta)*ds;       
    dx = cosd(auxTheta)*ds;

    xVat(1,:) = xAux - dx;
    yVat(1,:) = yAux - dy;
    
    xVat(2,:) = xAux + dx;
    yVat(2,:) = yAux + dy;
    
    plot(yVat,xVat,'Color',[0.9290 0.6940 0.1250 trA],'LineWidth',lsize)
    hold on
    set(gca,'xtick',[])
    set(gca,'ytick',[])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 05: Ref curve of QP

if sum(joices==QPref) > 0
    [~, ~, yRef] = fun_getQPinfo(xvec,-1,a,b,T0,T1,isMirrored);
    plot(yRef,xvec,":c",'LineWidth',2)
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 06: Ref curve of VAT

if sum(joices==VATref) > 0
    [~,yref] = fun_getVATinfo(xvec,a,b,T0,T1,0);
    plot(yref,xvec,":g",'LineWidth',2)
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 07: QP Point

if sum(joices==QPpoint) > 0
    [tt, xRef, yRef, ~, ~, xRefPair] = fun_getQPinfo(xpoint,ypoint,a,b,T0,T1,isMirrored);
    xRefPair = xRefPair{1};
    [~, ~, yRefPair] = fun_getQPinfo(xRefPair,-1,a,b,T0,T1,isMirrored);
    
    plot([ypoint yRefPair(1)], [xpoint xRefPair(1)], 'color',0.5*[1 1 1], 'LineWidth', 1.5)
    plot([ypoint yRefPair(2)], [xpoint xRefPair(2)], 'color',0.5*[1 1 1], 'LineWidth', 1.5)
    plot([ypoint yRef], [xpoint xRef], "b", 'LineWidth', 1.5)
    hold on
    plot(ypoint, xpoint, "ok", 'LineWidth', 1.5, 'MarkerSize', 6)
    plot(ypoint, xpoint, ".k", 'LineWidth', 1.5, 'MarkerSize', 6)
    text(ypoint+b/10,xpoint-a/40,num2str(tt,3),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',14,'FontWeight','bold')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 08: VAT Point

if sum(joices==VATpoint) > 0
    tt = fun_getVATinfo(xpoint,a,b,T0,T1,0);
    plot(ypoint, xpoint, "ok", 'LineWidth', 1.5, 'MarkerSize', 6)
    plot(ypoint, xpoint, ".k", 'LineWidth', 1.5, 'MarkerSize', 6)
    text(ypoint+b/15,xpoint-a/35,num2str(tt,3),'VerticalAlignment','bottom','HorizontalAlignment','right','FontSize',14,'FontWeight','bold')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 09: Angles of VAT on QP ref curve

if sum(joices==VATonQPref) > 0

    % QP Reference Curve
    [~, ~, yRef] = fun_getQPinfo(xvec,-1,a,b,T0,T1,isMirrored);
    plot(yRef,xvec,":c",'LineWidth',2)
    hold on
    [~, ~, yRef] = fun_getQPinfo(xDom,-1,a,b,T0,T1,isMirrored);

    % VAT angles
    auxTheta = fun_getVATinfo(xDom,a,b,T0,T1,0);

    dy = sind(auxTheta)*ds;       
    dx = cosd(auxTheta)*ds;

    xVat = zeros(2,nSepA);
    yVat = zeros(2,nSepA);

    xVat(1,:) = xDom - dx;
    yVat(1,:) = yRef - dy;
    
    xVat(2,:) = xDom + dx;
    yVat(2,:) = yRef + dy;
    
    plot(yVat,xVat,'Color',[0.9290 0.6940 0.1250 trA],'LineWidth',lsize)
    hold on
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot

axis equal
axis([0 b 0 a])
title(auxPlot, "T0 = "+T0+" | T1 = "+T1+" | mirrored = "+isMirrored) 
% set(gca,'XTick',[],'YTick',[])
hold off

end
end

set(hFigure, 'MenuBar', 'none');
set(hFigure, 'ToolBar', 'none');
set(hFigure, 'color','w');
% set(hFigure, 'Position', [0 50 1000*27/21 1000]);
t.Padding = 'compact';
t.TileSpacing = 'compact';
axis off

