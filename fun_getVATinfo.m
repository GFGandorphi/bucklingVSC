function [thetaRef, yRef, xRef] = fun_getVATinfo(xp,a,b,T0,T1,offsetY)
% -------------------------------------------------------------------------
% fun_getVATinfo evaluates yRef, xRef for VAT approach
% -------------------------------------------------------------------------
%   INPUTS:
%       xp ------------ x point(s) to have yRef, xRef evaluated
%       a ------------- plate length (x-direction)
%       b ------------- plate width  (y-direction)
%       T0 ------------ angle at plate's center
%       T1 ------------ angle at plate's x boundaries
%       offsetY ------- offset from center curve
    
    T0r = deg2rad(T0);
    if T0==T1
        T1r = deg2rad(T1+1e-5);
    else
        T1r = deg2rad(T1);
    end

    yRef = zeros(1,length(xp));
    thetaRef = zeros(1,length(xp));

    for i=1:length(xp)
        if xp(i) >= a/2
            delta = a-xp(i);
            fact = -1;
        else
            delta = xp(i);
            fact = 1;
        end
    
        xpAux = a - delta;
    
        xRefAux = xpAux - a/2;
        yRefAux = a/2/(T1r-T0r)*( -log(cos(T0r + 2/a*xRefAux*(T1r-T0r))) + log(cos(T0r)) );
        yRefAux = yRefAux + b/2;
    
        if xp(i) < a/2
            yRefAux = -yRefAux + b;
        end
    
        yRef(i) = yRefAux + offsetY;

        thetaRef(i) = fact * 2/a *(T0-T1) * (xp(i) - a/2) + T0;
    end
    xRef = xp;

end

