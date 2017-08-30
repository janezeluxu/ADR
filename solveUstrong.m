function [u] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,B,l,c2)
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
    %element
    %tau = tauFunc(element,c1ele,itau,kappa,a,l,x,BCval,InvJ,c2);    
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
% for i = 1:nPoints
%     tempK = K(i,i);
%        K(i,:) = K(i,:)/tempK; 
%        
%        F(i) = F(i)/tempK; 
% end
K(nPoints,:) = 0;
K(:,nPoints) = 0;
K(nPoints,nPoints) = 1;
%cond(K)
u = K\F;
dlmwrite('testsKF/Pe1000g0b0K.txt',K,'delimiter','\t','precision',16)
dlmwrite('testsKF/Pe1000g0b0F.txt',F,'delimiter','\t','precision',16)
end

