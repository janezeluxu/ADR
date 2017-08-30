function [u] = testRun()
clc
clear all
nPoints = 41;
c1 = ones(nPoints,1);
global stabilizer
stabilizer = 'VMS';
global h
h = ones(nPoints,1)*(1/(nPoints-1));

x = 0:(1/40):1;
% (global) Peclet number
PeG = 20;
% (global) Damkohler numner
DaG = PeG*500;
% sign of reaction term
bsign = 1;
% num. of points

L=1;
a=1;
kappa = abs(a)*L/PeG; % diffusion
B = bsign*abs(a)*DaG/L; % reaction
l=0;
[u] = solveUtest(nPoints,x,[0,1],c1,2,kappa,a,B,l,0);
plot(x,u,'ro-')
end

function [u] = solveUtest(nPoints,x,BCval,c1,itau,kappa,a,B,l,c2)
global h;
%source
nElement = nPoints-1;
IEN(:,1) = 1:nPoints-1;
IEN(:,2) = 2:nPoints;
K = zeros(nPoints,nPoints);
F = zeros(nPoints,1);
for element = 1:nElement
    J = h(element);
    InvJ = 1/J; 
    [k,f] = ElementStiff(a,kappa,B,itau,J,InvJ,l,x,BCval,element,c1,c2);
    
    %assemble    
    K(IEN(element,1),IEN(element,1)) = k(1,1)+K(IEN(element,1),IEN(element,1));
    K(IEN(element,1),IEN(element,2)) = k(1,2)+K(IEN(element,1),IEN(element,2));
    K(IEN(element,2),IEN(element,1)) = k(2,1)+K(IEN(element,2),IEN(element,1));
    K(IEN(element,2),IEN(element,2)) = k(2,2)+K(IEN(element,2),IEN(element,2));
    
    F(IEN(element,1)) = f(1)+F(IEN(element,1));
    F(IEN(element,2)) = f(2)+ F(IEN(element,2));
end
%% apply Strong boundary condition
for i = 1:nPoints
    F(i) = F(i)-K(i,1)*BCval(1)-K(i,nPoints)*BCval(2);
end
F(1) = BCval(1);
F(nPoints) = BCval(2);
K(1,:) = 0;
K(:,1) = 0;
K(1,1) = 1;
K(nPoints,:) = 0;
K(:,nPoints) = 0;
K(nPoints,nPoints) = 1;
%cond(K)
u = K\F;
end

function [k,f] = ElementStiff(a0,k0,B,itau,J,InvJ,f0,xarray,BCval,ele,c1,c2)
global h
global stabilizer
p=1;
n = 20;
[xi,w] = GaussQuad(n,1);
[ShapeFunc,DivSF] = ShapeFuncDiv(p,xi);
k = zeros(p+1,p+1);
f = zeros(p+1,1);
c1ele = max(c1(ele),c1(ele+1));
for i = 1:length(w)
    Nax = DivSF(:,i)*InvJ;
    Nbx = Nax';
    Na = ShapeFunc(:,i);
    Nb = Na';
    x =  xarray(ele)+xi(i)*h(ele);
    a =a0; kappa=k0;l=f0;
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
end

function [tau] = tauFunc(element,c1element,itau,kappa,a,B,InvJ,c2)
global h;
if(itau == 2) 
    tau1 = h(element)/(2*a);
    tau2 = h(element)^2/(4*kappa);
    tau = 1/sqrt(tau1^(-2) + 9*tau2^(-2)+B^2);
    tau = 10e-5;
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

function [xx, ww] = GaussQuad(n, a, b)

   % Check number of input arguments.
   narginchk(1, 3);

   % Assign default values to missing arguments.
   switch nargin
      case 1                    % GAUSSQUAD(N)
         b = 1;
         a = -1;
      case 2                    % GAUSSQUAD(N, C)
         b = a;
         a = 0;
   end

   u = 1 : n-1;
   u = u ./ sqrt(4*u.^2 - 1);

   % Same as A = diag(u, -1) + diag(u, 1), but faster (no addition).
   A = zeros(n, n);
   A( 2 : n+1 : n*(n-1) ) = u;
   A( n+1 : n+1 : n^2-1 ) = u;

   % Find the base points X and weight factors W for the interval [-1,1].
   [v, x] = eig(A);
   [x, k] = sort(diag(x));
   w = 2 * v(1,k)'.^2;

   % Linearly transform from [-1,1] to [a,b].
   x = (b - a) / 2 * x + (a + b) / 2;   % transform base points X
   w = (b - a) / 2 * w;                 % adjust weigths

   % If output arguments are given, return output and exit.
   if nargout
      xx = x;
      ww = w;
      return
   end
end