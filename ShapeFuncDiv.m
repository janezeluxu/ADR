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