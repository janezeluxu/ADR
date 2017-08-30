function [c1] = solvec1HO( u,x,nPoints,k0,a0,f0,c2h,c2H,BCval)

c1 = zeros(nPoints-1,1);
nEle = nPoints-1;
for ele = 2:nEle-1
    %[a,kappa,~]=ProblemDefinition(a0,k0,f0,x(node),bcValue);
    J = x(ele+1)-x(ele);
    InvJ = 1/J;
    [uH] = quadraticFit(u,x,ele,nPoints,BCval);
    uh = [u(ele);u(ele+1)];
    c1(ele) = Elementc1(a0,k0,f0,J,InvJ,x,ele,c2h,c2H,uH,uh,BCval);
end    
    
for ele = 1:nEle
     if c1(ele)<0
         %c1(ele) = 0;
         c1(ele) = -c1(ele);
     end
end
end

function [uH] = quadraticFit(u,x,ele,nPoints,bcValue)
global h
[SFH0,~] = ShapeFuncDiv(2,0);
[SFH1,~] = ShapeFuncDiv(2,0.5);
[SFH2,~] = ShapeFuncDiv(2,1);
M = [SFH0';SFH1';SFH2']^-1;
if ele ==1
    ufit = [u(ele),u(ele+1),u(ele+2)];
    xfit = [x(ele),x(ele+1),x(ele+2)];
    
    p = polyfit(xfit,ufit,2);
    x1 = [x(ele),x(ele)+h(ele)/2,x(ele+1)];
    uList = polyval(p,x1);
    uH = M*uList';
    uH(1) = bcValue(1);
elseif ele == nPoints-1
    ufit = [u(ele-1),u(ele),u(ele+1)];
    xfit = [x(ele-1),x(ele),x(ele+1)];
    
    p = polyfit(xfit,ufit,2);
    x1 = [x(ele),x(ele)+h(ele)/2,x(ele+1)];
    uHx = p(1)*2*x1(2)+p(2);
    
    uList = polyval(p,x1);
    uH = M*uList';
    uH(3) = bcValue(2);
else
    ufit = [u(ele-1),u(ele),u(ele+1),u(ele+2)];
    xfit = [x(ele-1),x(ele),x(ele+1),x(ele+2)];
    
    p = polyfit(xfit,ufit,2);
    x1 = [x(ele),x(ele)+h(ele)/2,x(ele+1)];
    uList = polyval(p,x1);
    uH = M*uList';
    
end


end

function [c1ele] = Elementc1(a0,k0,f0,J,InvJ,xarray,ele,c2h,c2H,uH,uh,BCval)
global h
n = nIntergerPoints(2,1);
[xi,w] = GaussQuad(n,1);
[SFh,DivSFh] = ShapeFuncDiv(1,xi);
[SFH,DivSFH] = ShapeFuncDiv(2,xi);
L1 = 0; L2 = 0; M1 = 0; M2 = 0;
for i = 1:length(w)
    NxH = DivSFH(:,i)*InvJ;
    NH = SFH(:,i);
    
    Nxh = DivSFh(:,i)*InvJ;
    Nh = SFh(:,i);

    x =  xarray(ele)+xi(i)*h(ele);
    [a,kappa,l, da]=ProblemDefinition(a0,k0,f0,x,BCval);
    
    wH = NH'*uH;
    wHx = NxH'*uH;
    whx = wHx;
    wh = wH;
    H = h(ele);
    L1 = L1+whx*a*a*h(ele)^c2h*(Nxh'*uh)*w(i)*J-whx*a*l*h(ele)^c2h;
    L2 = L2+wHx*a*a*H^c2H*(NxH'*uH)*w(i)*J-wHx*a*l*H^c2H;
    M1 = M1+(whx*(-a)*(NH'*uH)*w(i)*J)+wh*(-da)*(NH'*uH)*w(i)*J+(kappa*(NxH'*uH)*w(i)*J);
    M2 = M2+(wHx*(-a)*(Nh'*uh)*w(i)*J)+wH*(-da)*(Nh'*uh)*w(i)*J+(kappa*(Nxh'*uh)*w(i)*J);
end
 L = L1-L2;
 M = M1-M2;
 c1ele = M/L;
end