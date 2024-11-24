function h = fun_circulo(x,y,r,eC,fC)
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position',[px py d d],'Curvature',[1,1],'EdgeColor',eC,'FaceColor',fC);
    %daspect([1,1,1])
end