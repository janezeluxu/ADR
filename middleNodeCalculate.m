function [GA, QA]=middleNodeCalculate()
syms x A u1h u2h u3h c1h c2h c3h f  hh a k c2
x1 = (A-2)*hh;
x2 = (A-1)*hh;
x3 = A*hh;

N1 = (x2-x)/hh;
N21 = (x-x1)/hh;
N22 = (x3-x)/hh;
N3 = (x-x2)/hh;

N1x = -(1/hh);
N21x = 1/hh;
N22x = -(1/hh);
N3x = (1/hh);

GA1 = int(N21*a*(N1x*u1h+N21x*u2h), x,x1, x2)+int(N22*a*(N22x*u2h+N3x*u3h), x,x2, x3);
GA2 = int(N21x*k*(N1x*u1h+N21x*u2h), x, x1, x2)+int(N22x*k*(N22x*u2h+N3x*u3h), x,x2, x3);
GA3 = int(a*N21x*0.5*(c1h+c2h)*a*hh^c2*(N1x*u1h+N21x*u2h), x,x1, x2)+int(a*N22x*0.5*(c2h+c3h)*a*hh^c2*(N22x*u2h+N3x*u3h), x,x2, x3);
GA4 = int(N21x*0.5*(c1h+c2h)*a*hh^c2*f, x,x1, x2)+int(N22x*0.5*(c2h+c3h)*a*hh^c2*f, x,x2, x3);
GA5 = int(N21*f, x,x1, x2)+int(N22*f, x,x2, x3);
GA = GA1 +GA2+GA3-GA4-GA5;

syms H 
M = [1.25 0.25;0.25,1.25]^-1*[1 0.5 0;0 0.5 1];

u1H = M(1,1)*u1h+M(1,2)*u2h+M(1,3)*u3h;
u2H = M(2,1)*u1h+M(2,2)*u2h+M(2,3)*u3h;

N1H=(x3-x)/(2*hh);
N2H = (x-x1)/(2*hh);
N1xH = -1/(2*hh);
N2xH = 1/(2*hh);
WX = 1;
wx1 = 1;
wx2 = 1;

QA1 = a^2*int(wx1*c2h*hh^c2*(N1x*u1h+N21x*u2h), x,x1, x2)+...
    a^2*int(wx2*c2h*hh^c2*(N22x*u2h+N3x*u3h), x,x2, x3);
QA2 = a^2*int(WX*c2h*H^c2*(N1xH*u1H+N2xH*u2H), x,x1, x3);

QA3 = int(-a*WX*(N1H*u1H+N2H*u2H), x,x1, x3);
QA4 = int(k*WX*(N1xH*u1H+N2xH*u2H), x,x1, x3);

QA5 = int(-a*wx1*(N1*u1h+N21*u2h), x,x1, x2)+int(-a*wx2*(N22*u2h+N3*u3h), x,x2, x3);
QA6 = int(k*wx1*(N1x*u1h+N21x*u2h), x,x1, x2)+int(k*wx2*(N22x*u2h+N3x*u3h), x,x2, x3);

QA = QA1-QA2-(QA3+QA4-(QA5+QA6));

GAu1h = diff(GA,u1h);
GAu2h = diff(GA,u2h);
GAu3h = diff(GA,u3h);
GAc1h = diff(GA,c1h);
GAc2h = diff(GA,c2h);
GAc3h = diff(GA,c3h);

QAu1h = diff(QA,u1h);
QAu2h = diff(QA,u2h);
QAu3h = diff(QA,u3h);
QAc1h = diff(QA,c1h);
QAc2h = diff(QA,c2h);
QAc3h = diff(QA,c3h);
end