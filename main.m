function [u,c1,c1new] = main(itau,nPoints,x,uTol,maxIter,dt,c1,kappa,a,B,f,BCval,FigNum,lineType,c2)
global h;
[u] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,B,f,c2);
%uGa = GaCalculate(nPoints,x,BCval,c1,itau,kappa,a,B,f,c2,u)
c1new = 0;
G = 0;
Q=0;
if itau == 3
    for iteration = 1:maxIter
        iteration
        [c1new] = solvec1Project( u,x,nPoints,kappa,a,f,c2,BCval);
        c1 = (1/(1+dt))*c1+(dt/(1+dt))*c1new;
        [unew] = solveUstrong(nPoints,x,BCval,c1,itau,kappa,a,B,f,c2);
        
        normu = norm(u-unew)/norm(unew);
        [Ga] = GaCalculate(nPoints,x,BCval,c1,itau,kappa,a,B,f,c2,u)/(nPoints)^0.5;
        [Qa] = QaCalculate(nPoints,x,BCval,c1,kappa,a,f,c2,u)/(nPoints)^0.5;
        G(iteration) = Ga;
        Q(iteration) = Qa;
        u = unew;
        if  (Ga <uTol) 
            break
        end
        
    end
        
end
figure (FigNum)
plot(x,u,lineType)

Pe = norm(a)*h(1)/(2*norm(kappa))
title(horzcat('Pe = 10000 Da = Pe/500'),'FontSize',20);
hold on
end

