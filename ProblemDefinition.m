function [a,kappa,f,da]=ProblemDefinition(a0,k0,f0,x,BCval)
global option
if (option==1)
    a = a0;
    da = 0;
    kappa = k0;
    f = f0;
elseif (option==2)
    a = a0;
    da = 0;
    kappa = k0;
    f = f0;
end
end