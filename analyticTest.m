Pe = 1000;
Da = Pe/5; % less than Pe/4

kappa = 1/Pe;
b = -Da;

% assume a/advection to be positive and b/reaction to be negative
l1 = (Pe+sqrt(Pe)*sqrt(Pe-4*Da))/2;
l2 = (Pe-sqrt(Pe)*sqrt(Pe-4*Da))/2;
g1 = (exp(-l1)-exp(-l2))/((l1/l2)^(-l1/(l1-l2))-(l1/l2)^(-l2/(l1-l2)));

npoints=1000;

for i=1:npoints
  x(i)=1.0*(i-1)/(npoints-1);
  u(i) = (exp(-l1*(1-x(i)))-exp(-l2*(1-x(i))))/((l1/l2)^(-l1/(l1-l2))-(l1/l2)^(-l2/(l1-l2)));
end

last_npoints=300;

plot(x(end-last_npoints:end),u(end-last_npoints:end))
%plot(x,u)