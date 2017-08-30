function [Qa] = QaCalculate(nPoints,x,BCval,c1,k0,a0,f0,c2,u)
[~,l1,m1] = solvec1Project( u,x,nPoints,k0,a0,f0,c2,BCval);
[a,kappa,l]=ProblemDefinition(a0,k0,f0,x,BCval);
Qa = l1.*c1-m1;
Qa = (sum(Qa.^2))^0.5;
end