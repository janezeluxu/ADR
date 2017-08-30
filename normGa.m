function [normGa,normQa] = normGa(c,u,nPoints,f, a, k, c2)
G = zeros(nPoints,1);
global h;
[GA, QA]=middleNodeCalculate();
%% loop through each node, direct put into global matrix
for i = 2:nPoints-1
%for i = nPoints-1
    hh = h(i);
    H = 2*h(i);
    u1h=u(i-1);
    u2h=u(i);
    u3h=u(i+1);
    
    c1h=c(i-1);
    c2h=c(i);
    c3h=c(i+1);
    
    %each middle node, pre-calculated
    %GA = eval(GA)
    G(i) = eval(GA);  
    Q(i) = eval(QA);
end
%lastG = G(nPoints-1)
normGa = sum(G.^2);
normGa = normGa^0.5;

normQa = sum(Q.^2);
normQa = normQa^0.5;
end