function [thetaRef, xRef, yRef, cfs, dist, xRefPair, imin] = fun_getQPinfo(xp,yp,a,b,T0,T1,isMirrored)
% -------------------------------------------------------------------------
% fun_getQPinfo evaluates several parameters for Quasi-Parallel approach
% -------------------------------------------------------------------------
%   INPUTS:
%       xp ------------ x point to have infos evaluated. If yp=-1, xp=xRefs
%       yp ------------ y point to have infos evaluated. If yp=-1, xp=xRefs
%       a ------------- plate length (x-direction)
%       b ------------- plate width  (y-direction)
%       T0 ------------ angle at plate's center
%       T1 ------------ angle at plate's x boundaries
%       isMirrored ---- bool if the ref. curve is mirrored regarding x-axis

%% 00 - Initializations
    xp = xp(:)';
    yp = yp(:)';
    nVals = length(xp);

    thetaRef = zeros(1, nVals);
    xRef = zeros(1, nVals);
    yRef = zeros(1, nVals);
    dist = zeros(1, nVals);
    xRefPair = cell(1, nVals);
    imin = zeros(1, nVals);

%% 01 - Reference curve definition. Two 3rd order polynomial interpolation
    
    % syms As Bs Cs Ds Es Fs Gs Hs as bs tgT0s tgT1s x
    % 
    % eqn1 = As*(as/2)^3 + Bs*(as/2)^2 + Cs*(as/2) + Ds == bs/2;
    % eqn2 = As*(as)^3 + Bs*(as)^2 + Cs*(as) + Ds == bs;
    % eqn3 = 3*As*(as/2)^2 + 2*Bs*(as/2) + Cs == tgT0s;
    % eqn4 = 3*As*(as)^2 + 2*Bs*(as) + Cs == tgT1s;
    % sol1 = solve([eqn1, eqn2, eqn3, eqn4], [As, Bs, Cs, Ds]);
    % 
    % eqn1 = Es*(as/2)^3 + Fs*(as/2)^2 + Gs*(as/2) + Hs == bs/2;
    % eqn2 = 3*Es*(as/2)^2 + 2*Fs*(as/2) + Gs == tgT0s;
    % eqn3 = Hs == 0;
    % eqn4 = Gs == tgT1s;
    % sol2 = solve([eqn1, eqn2, eqn3, eqn4], [Es, Fs, Gs, Hs]);
    
    Av = 4*(a*tand(T0) - 2*b + a*tand(T1))/a^3;
    Bv = -2*(5*a*tand(T0) - 9*b + 4*a*tand(T1))/a^2;
    Cv = (8*a*tand(T0) - 12*b + 5*a*tand(T1))/a;
    Dv = 3*b - 2*a*tand(T0) - a*tand(T1);
    Ev = 4*(a*tand(T0) - 2*b + a*tand(T1))/a^3;
    Fv = -2*(a*tand(T0) - 3*b + 2*a*tand(T1))/a^2;
    Gv = tand(T1);
    Hv = 0;

    cs = [Av, Bv, Cv, Dv, a/2, inf;
          Ev, Fv, Gv, Hv, -inf, a/2];

    if isMirrored
        cfs = [-Av, -Bv, -Cv, -Dv+b, -Ev, -Fv, -Gv, -Hv+b];
    else
        cfs = [Av, Bv, Cv, Dv, Ev, Fv, Gv, Hv];
    end


%% 02 - Evaluation of yRef and/or min distance 
    % If yp == -1: xp=xRef and yRef will be evaluated with xp directly
    % If yp != -1: Min. distance evaluation between (xp,yp) and ref. curve
    
    if length(yp) == 1 && yp(1) == -1
        for i=1:nVals
            % Curve definition
            if xp(i)>a/2; ic=1; else; ic=2; end
            [A, B, C, D] = deal(cs(ic,1),cs(ic,2),cs(ic,3),cs(ic,4));
            yRef(i) = polyval([A, B, C, D],xp(i));
        end
    
        % Outputs
        xRef = xp;
        if isMirrored
            yRef = -yRef + b;
        end
    else
        for i=1:nVals
    
            xRefPairVec = zeros(2,1);
            minDist2Vec = zeros(2,1);
    
            for ic=1:2
                
                % Curve definition
                [A, B, C, D, xlim1, xlim2] = deal(cs(ic,1),cs(ic,2),cs(ic,3),cs(ic,4),cs(ic,5),cs(ic,6));
    
                % System to retrieve minimum distance per curve
                % syms x A B C D xp yp
                % aux1 = (A*(x)^3 + B*(x)^2 + C*(x) + D)*(3*A*(x)^2 + 2*B*(x) + C);
                % aux2 = yp*(3*A*(x)^2 + 2*B*(x) + C);
                % aux3 = x - xp + aux1 - aux2;
                % aux3 = expand(aux3);
                % aux3 = collect(aux3);
                if isMirrored
                    auxPolyMinDist2 = [3*A^2, 5*A*B, (4*A*C+2*B^2), (3*B*C+3*A*D+3*A*yp(i)-3*A*b), (1+C^2+2*B*D+2*B*yp(i)-2*B*b), (C*D+C*yp(i)-C*b-xp(i))];
                else
                    auxPolyMinDist2 = [3*A^2, 5*A*B, (4*A*C+2*B^2), (3*B*C+3*A*D-3*A*yp(i)), (1+C^2+2*B*D-2*B*yp(i)), (C*D-C*yp(i)-xp(i))];
                end
    
                auxResultsMinDist2 = roots(auxPolyMinDist2); 
                auxResultsMinDist2 = auxResultsMinDist2(imag(auxResultsMinDist2)==0);
                auxResultsMinDist2 = auxResultsMinDist2(auxResultsMinDist2>=xlim1);
                auxResultsMinDist2 = auxResultsMinDist2(auxResultsMinDist2<=xlim2);
                xRefsVec = [auxResultsMinDist2; xlim1; xlim2];
    
                yRefsVec = polyval([A, B, C, D],xRefsVec);
                if isMirrored
                    yRefsVec = -yRefsVec + b;
                end
    
                d2Vec = (xRefsVec - xp(i)).^2 + (yRefsVec - yp(i)).^2;
                [d2,id2] = min(d2Vec);
                minDist2Vec(ic) = d2;
                xRefPairVec(ic) = xRefsVec(id2);
            end
            [dist2,imind] = min(minDist2Vec);
            [A, B, C, D] = deal(cs(imind,1),cs(imind,2),cs(imind,3),cs(imind,4));
    
            % Outputs
            imin(i) = imind;
            xRefPair{i} = xRefPairVec;
            xRef(i) = xRefPairVec(imind);
            yRef(i) = polyval([A, B, C, D],xRef(i));
            if isMirrored
                yRef(i) = -yRef(i) + b;
            end
            dist(i) = dist2^(0.5);
        end
    end

%% 03 - Theta evaluation

    for i=1:nVals
        
        if xRef(i)>a/2; ic=1; else; ic=2; end
        [A, B, C] = deal(cs(ic,1),cs(ic,2),cs(ic,3));
        
        thetaRef(i) = atand(polyval([3*A, 2*B, C],xRef(i)));
        if isMirrored
            thetaRef(i) = -thetaRef(i);
        end 
    end

end

