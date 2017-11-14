clear all
global option
global stabilizer
stabilizer = 'VMS';
r=1;
q = 1/r;
domain = [0,1];
option =1;
caseNumber = 12;
%% problem setup
% (global) Peclet number
PeG = 10;
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

% assume a/advection to be positive and b/reaction to be negative
l1 = (PeG+sqrt(PeG)*sqrt(PeG+sign(B)*4*DaG))/2;
l2 = (PeG-sqrt(PeG)*sqrt(PeG+sign(B)*4*DaG))/2;
g1 = (exp(-l1)-exp(-l2))/((l1/l2)^(-l1/(l1-l2))-(l1/l2)^(-l2/(l1-l2)));

f = 0;
bcValue = [g1,0];
c2 = 2;

%% iteration setup
uTol = 1e-3;
BurnIter = 20;
Burndt  = 3/7;

[xarray]=meshDefinition(nPoints,domain,q);
c0 = ones(nPoints,1)*0;

for i = 1:nPoints
    x = xarray(i);
    [a0,k0,l0]=ProblemDefinition(a,kappa,f,x,bcValue);
    c0(i) = 0;
    %c0(i) = abs(1/a0);
end
nInt = 60;
[u1] = main(1,nPoints,xarray,uTol,BurnIter,Burndt,c0, kappa,a,B,f,bcValue,1,'y*-',c2);
[error1,H1norm,L2ele1] = ErrorEstimate(u1,xarray,nInt,nPoints,a,kappa,l1,l2,f,bcValue);
error1
H1norm
[u2] = main(2,nPoints,xarray,uTol,BurnIter,Burndt,c0, kappa,a,B,f,bcValue,1,'bs-',c2);
[error2,H2norm,L2ele2] = ErrorEstimate(u2,xarray,nInt,nPoints,a,kappa,l1,l2,f,bcValue);
error2
H2norm
[u3,c3,Ga3,G3] = main(3,nPoints,xarray,uTol,BurnIter,Burndt,c0, kappa,a,B,f,bcValue,1,'rs-',c2);
[error3,H3norm,L2ele3] = ErrorEstimate(u3,xarray,nInt,nPoints,a,kappa,l1,l2,f,bcValue);
error3
H3norm
xx=0:0.01:1;
ue = Exact(xx,a,kappa,l1,l2,f,bcValue);
plot(xx,ue,'c-')

[u0] = main(0,nPoints,xarray,uTol,BurnIter,Burndt,c0, kappa,a,B,f,bcValue,1,'k-',c2);
[error0,H0norm,L2ele0] = ErrorEstimate(u0,xarray,nInt,nPoints,a,kappa,l1,l2,f,bcValue);
error0
H0norm
%xlim([0.8 1])
%% write solution to file
writeSolution(caseNumber,ue,u1,u2,u3,u0,c3,Ga3,G3,error0,error2,error3, ...
    L2ele0, L2ele2,L2ele3,H0norm,H2norm,H3norm)