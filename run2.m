%% clear all
global option
global stabilizer
stabilizer = 'VMS';
%nPoints = 41;
r=1;
q = 1/r;
domain = [0,1];
option =1;
global h;
%% problem setup
% (global) Peclet number
PeG = 100;
% (global) Damkohler numner
DaG = PeG/500;
% sign of reaction term
bsign = -1;
% num. of points
nPoints = 41;

L=1;
a=1;
kappa = abs(a)*L/PeG; % diffusion
B = bsign*abs(a)*DaG/L; % reaction
%B=0
% assume a/advection to be positive and b/reaction to be negative
l1 = (PeG+sqrt(PeG)*sqrt(PeG+sign(B)*4*DaG))/2;
l2 = (PeG-sqrt(PeG)*sqrt(PeG+sign(B)*4*DaG))/2;
g1 = (exp(-l1)-exp(-l2))/((l1/l2)^(-l1/(l1-l2))-(l1/l2)^(-l2/(l1-l2)));

f = 0;
bcValue = [g1,0];
c2 = 2;

%% iteration setup
uTol = 1e-3;
BurnIter = 10;
Burndt  = 3/7;

[xarray]=meshDefinition(nPoints,domain,q);
c0 = ones(nPoints,1)*0;

for i = 1:nPoints
    x = xarray(i);
    [a0,k0,l0]=ProblemDefinition(a,kappa,f,x,bcValue);
    c0(i) = abs(h(1)/kappa);
    %c0(i) = abs(1/a0);
    
end
nInt = 60;
[u2,c1,c1new] = main(2,nPoints,xarray,uTol,BurnIter,Burndt,c0, kappa,a,B,f,bcValue,1,'bs-',c2);
%[error] = ErrorEstimate(u,xarray,nInt,nPoints,a,kappa,f,bcValue)
[u3,c1,c1new] = main(3,nPoints,xarray,uTol,BurnIter,Burndt,c0, kappa,a,B,f,bcValue,1,'rs-',c2);
%ylim([0,1.1])
%xlim([0.8,1])
ue = Exact(xarray,a,kappa,l1,l2,f,bcValue);
plot(xarray,ue,'c*-')
%figure(2)
%semilogy(xarray,abs(ue'-u2),'bs-')
%hold on
%semilogy(xarray,abs(ue'-u3),'rs-')

%% error calculation