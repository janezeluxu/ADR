function [Ga,GaL2] = GaCalculate(nPoints,x,BCval,c1,itau,kappa,a,B,l,c2,u)

%source
nElement = nPoints-1;
%l = zeros(nElement,1)*0.2;

IEN(:,1) = 1:nPoints-1;
IEN(:,2) = 2:nPoints;
Ga = zeros(nPoints,1);

for element = 1:nElement
    J = x(element+1)-x(element);
    InvJ = 1/J;
    %element
    u1h = u(element);
    u2h = u(element+1);
    [ga] = ElementGa(a,kappa,B,itau,J,InvJ,l,x,BCval,element,c1,c2,u1h,u2h);
    
    %assemble        
    Ga(IEN(element,1)) = ga(1)+Ga(IEN(element,1));
    Ga(IEN(element,2)) = ga(2)+ Ga(IEN(element,2));
end
GaL2 = (sum(Ga(2:end-1).^2))^0.5/(nPoints)^0.5;
end

function [ga] = ElementGa(a0,k0,B,itau,J,InvJ,f0,xarray,BCval,ele,c1,c2,u1,u2)
global h;
global stabilizer;
p=1;
n = nIntergerPoints(p,60);
[xi,w] = GaussQuad(n,1);
[ShapeFunc,DivSF] = ShapeFuncDiv(p,xi);
ga = zeros(p+1,1);
c1ele = max(c1(ele),c1(ele+1));
for i = 1:length(w)
    Nax = DivSF(:,i)*InvJ;
    Nbx = Nax';
    Na = ShapeFunc(:,i);
    Nb = Na';
    x =  xarray(ele)+xi(i)*h(ele);
    %c1ele = [c1(ele),c1(ele+1)]*Na;
    [a,kappa,l]=ProblemDefinition(a0,k0,f0,x,BCval);
    [tau] = tauFunc(ele,c1ele,itau,a,kappa,InvJ,c2);
    if strcmp(stabilizer, 'GLS') == 1
    ga = ga+((Na*a*(Nax(1)*u1+Nax(2)*u2))+...
        (Nax*kappa*(Nax(1)*u1+Nax(2)*u2))+...
        (Na*B*(Na(1)*u1+Na(2)*u2))+...
        (a*Nax+B*Na)*tau*(a*(Nax(1)*u1+Nax(2)*u2)+B*(Na(1)*u1+Na(2)*u2))...
        -(Na*l)-((a*Nax+B*Na)*tau*l))*J*w(i);
    elseif strcmp(stabilizer, 'SUPG') == 1
        ga = ga+((Na*a*(Nax(1)*u1+Nax(2)*u2))+...
        (Nax*kappa*(Nax(1)*u1+Nax(2)*u2))+...
        (Na*B*(Na(1)*u1+Na(2)*u2))+...
        (a*Nax)*tau*(a*(Nax(1)*u1+Nax(2)*u2)+B*(Na(1)*u1+Na(2)*u2))...
        -(Na*l)-((a*Nax)*tau*l))*J*w(i);
    elseif strcmp(stabilizer, 'VMS') == 1
    ga = ga+((Na*a*(Nax(1)*u1+Nax(2)*u2))+...
        (Nax*kappa*(Nax(1)*u1+Nax(2)*u2))+...
        (Na*B*(Na(1)*u1+Na(2)*u2))+...
        (a*Nax-B*Na)*tau*(a*(Nax(1)*u1+Nax(2)*u2)+B*(Na(1)*u1+Na(2)*u2))...
        -(Na*l)-((a*Nax-B*Na)*tau*l))*J*w(i);
    else
        fprintf('Should specify stabilizer');
    end
end
end

function [tau] = tauFunc(element,c1element,itau,a,kappa,InvJ,c2)
global h;
if(itau == 1)
    alpha = norm(a)*h(element)/(2*norm(kappa));
    zi = coth(alpha)-1/alpha;
    tau = h(element)*zi/(2*norm(a));
elseif(itau == 2)  
    tau1 = h(element)/(2*a);
    tau2 = h(element)^2/(4*kappa);
    tau = 1/sqrt(tau1^-2 + 9*tau2^-2);
elseif(itau == 3)                    
    tau = c1element*(h(element)^c2); 
elseif(itau == 4)
    tau = (2*InvJ*norm(a)+2*InvJ^2*kappa)^-1;
elseif(itau == 0)
    tau = 0;
end

end

function [ShapeFunc,DivSF] = ShapeFuncDiv(p,x)
L = length(x);
ShapeFunc = zeros(p+1,L);
DivSF = zeros(p+1,L);
for v = 0:p
    coeff = nchoosek(p,v);
    ShapeFunc(v+1,:) = coeff*x.^v.*(1-x).^(p-v);
    DivSF(v+1,:) = coeff*v*x.^(v-1).*(1-x).^(p-v)-coeff*x.^(v).*(p-v).*(1-x).^(p-v-1);
end
end