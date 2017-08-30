function [k,f] = ElementStiff(a0,k0,B,itau,J,InvJ,f0,xarray,BCval,ele,c1,c2)
global h
global stabilizer
p=1;
n = nIntergerPoints(p,20);
[xi,w] = GaussQuad(n,1);
[ShapeFunc,DivSF] = ShapeFuncDiv(p,xi);
k = zeros(p+1,p+1);
f = zeros(p+1,1);
%fileID = fopen('x.txt','a+');
%c1ele = 0.5*(c1(ele)+c1(ele+1));
c1ele = max(c1(ele),c1(ele+1));
for i = 1:length(w)
    Nax = DivSF(:,i)*InvJ;
    Nbx = Nax';
    Na = ShapeFunc(:,i);
    Nb = Na';
    x =  xarray(ele)+xi(i)*h(ele);
    [a,kappa,l]=ProblemDefinition(a0,k0,f0,x,BCval);
    %c1ele = [c1(ele),c1(ele+1)]*Na;
    [tau] = tauFunc(ele,c1ele,itau,kappa,a,B,InvJ,c2);
    k = k+Nax*kappa*Nbx*J*w(i)+Na*a*Nbx*J*w(i)+Na*B*Nb*J*w(i);
    f = f+Na*l*J*w(i);
    if strcmp(stabilizer, 'GLS') == 1
        k = k+(a*Nax+B*Na)*tau*(a*Nbx+B*Nb)*J*w(i);
        f = f+(a*Nax+B*Na)*tau*l*J*w(i);
    elseif strcmp(stabilizer, 'SUPG') == 1
        k = k+(a*Nax)*tau*(a*Nbx+B*Nb)*J*w(i);
        f = f+(a*Nax)*tau*l*J*w(i);
    elseif strcmp(stabilizer, 'VMS') == 1
        k = k+(a*Nax-B*Na)*tau*(a*Nbx+B*Nb)*J*w(i);
        f = f+(a*Nax-B*Na)*tau*l*J*w(i);
    else
        fprintf('Should specify stabilizer');
    end

end
%fclose(fileID);
end

function [tau] = tauFunc(element,c1element,itau,kappa,a,B,InvJ,c2)
global h;
if(itau == 1)
    alpha = norm(a)*h(element)/(2*norm(kappa));
    zi = coth(alpha)-1/alpha;
    tau = h(element)*zi/(2*norm(a));
elseif(itau == 2) 
    tau1 = h(element)/(2*a);
    tau2 = h(element)^2/(4*kappa);
    tau = 1/sqrt(tau1^(-2) + 9*tau2^(-2)+B^2);
    %tau = 10e-5;
elseif(itau == 3)                    
    tau = c1element*(h(element)^c2); 
elseif(itau == 4)
    tau = (2*InvJ*norm(a)+2*InvJ^2*kappa)^-1;
elseif(itau == 0)
    tau = 0;
end

end