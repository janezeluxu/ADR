function [x]=meshDefinition(nPoints,domain,q)
global h
nElements = nPoints-1;
h = zeros(1,nElements);
x = zeros(1,nPoints);

if q == 1
    h1 = (domain(2)-domain(1))/(nPoints-1);
else
    h1 = (1-q)/(1-q^(nElements));
end

for i = 1:nElements
    h(i) = h1*q^(i-1);
end

for i = 1:nPoints
    x(i) = domain(1)+ sum(h(1:i-1));
end
end