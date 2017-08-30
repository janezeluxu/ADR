function [ue] = Exact(x,a,kappa,l1,l2,f,BCval)
global option
if (option ==1)
ue = (exp(-l1*(1-x))-exp(-l2*(1-x)))./((l1/l2)^(-l1/(l1-l2))-(l1/l2)^(-l2/(l1-l2)));
elseif (option==2)
    g1=BCval(1);
    g2=BCval(2);
ue = g1*(exp(-l1*(1-x))-exp(-l2*(1-x)))./(exp(-l1)-exp(-l2))+...
    g2*(-exp(-(l1+l2)+l1*x)+exp(-(l1+l2)+l2*x))./(exp(-l1)-exp(-l2));
elseif (option ==3)
    g2 = BCval(2);
    r = -a/kappa;
    [ue] = -(exp(r/2*(1-x))+exp(r/2))./(exp(r/2)+1).*(exp(r/2*(1-x))-exp(r/2))./(exp(r/2)-1);
    ue = (ue*g2+(f/a)*(x-ue)).*exp(x);
elseif (option==4)
    ue =(erf(x./(2*kappa)^0.5))./(erf(1/(2*kappa)^0.5));
end