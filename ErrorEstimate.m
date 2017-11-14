function [error,H1norm,L2ele] = ErrorEstimate(u,xarray,nInt,Grid_Size,a,kappa,l1,l2,f,g2)
global h

nElement = Grid_Size-1;
error = 0;
[xi,weight] = GaussQuad(nInt, 1);
L2ele = zeros(nElement,1);
H1norm = 0;
for ele = 1:nElement
    uE = [u(ele),u(ele+1)];
    len = length(xi);
    ShapeFunc = zeros(2,len);
    %uA = zeros(1,len);
    for n = 1:len
        ShapeFunc(:,n) = [1-xi(n);xi(n)];
    end
    uH = uE*ShapeFunc;
    uH1 = uE*[-1;1]*(1/h(ele));
    
    x =  xarray(ele)+xi*h(ele);
    [uA,uAx] = Exact(x,a,kappa,l1,l2,f,g2);
    L2ele(ele) = (((uH-uA').^2*weight*h(ele)));
    error = ((uH-uA').^2*weight*h(ele))+error;
    H1norm = ((uH1-uAx').^2*weight*h(ele))+H1norm;
    
end
error = error^0.5;
H1norm = H1norm^0.5;
end

