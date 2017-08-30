function [c1,l1,m1] = solvec1Project( u,x,nPoints,k0,a0,f0,c2,bcValue)

global h 
global stabilizer;
c1 = zeros(nPoints,1);
l1 = zeros(nPoints,1);
m1 = zeros(nPoints,1);
for node = 2:nPoints-1
    H = h(node-1)+h(node);
    
    uh = [u(node-1);u(node);u(node+1)];
    CordMatrix = [1,x(node-1);1,x(node);1,x(node+1)];
    
    [dUx, uHa,uHb, wHa, wHb] = LeastSquare( x,CordMatrix, uh,node,nPoints,bcValue);
    
    %% add quadrature rules, avoid x=0 point    
     n = nIntergerPoints(1,10);
    [xi,w] = GaussQuad(n,1);
    [ShapeFunc,DivSF] = ShapeFuncDiv(1,xi);
    A=0; B=0; C=0; D=0; E=0; F=0;G=0;I = 0;
    
    for i = 1:length(w)
        J1 = h(node-1);
        InvJ1 = 1/J1;
        Nax1 = DivSF(:,i)*InvJ1;        
        x1 =  x(node-1)+xi(i)*h(node);
        [a1,kappa1,f1,da1]=ProblemDefinition(a0,k0,f0,x1,bcValue);
        duha = [u(node-1),u(node)]*Nax1;
        
        J2 = h(node);
        InvJ2 = 1/J2;
        Nax2 = DivSF(:,i)*InvJ2;
        x2 =  x(node)+xi(i)*h(node);        
        [a2,kappa2,f2,da2]=ProblemDefinition(a0,k0,f0,x2,bcValue);
        duhb = [u(node),u(node+1)]*Nax2;
        
        xH =  x(node-1)+xi(i)*H;
        [a,kappa,f,da]=ProblemDefinition(a0,k0,f0,xH,bcValue);

        Na = ShapeFunc(:,i);
        uha = [u(node-1),u(node)]*Na;
        uhb = [u(node),u(node+1)]*Na;
        uH  =[uHa,uHb]*Na;
        
        wH = [wHa,wHb]*Na;
        wHx = (wHb-wHa)/H;

        whx1 = wHx;whx2 = wHx;
        wh1 = wH; wh2 = wH;
        
        if strcmp(stabilizer, 'GLS') == 1
        A = A+(a1*whx1+B*wh1)*h(node-1)^(c2)*(a1*duha+B*uha)*J1*w(i)...
            +(whx2*a2+B*wh2)*h(node)^(c2)*(a2*duhb+B*uhb)*J2*w(i);
        B = B+(a*wHx+B*wH)*H^(c2)*(a*dUx+B*uH)*H*w(i);
        
        C = C-wHx*a*(uH*H*w(i))-wH*da*(uH*H*w(i));
        D = D+wHx*kappa*dUx*H*w(i)+wH*B*uH*H*w(i);
        
        E = E+whx1*a1*(-uha)*J1*w(i)+whx2*a2*(-uhb)*J2*w(i)+...
            wh1*da1*(-uha)*J1*w(i)+wh2*da2*(-uhb)*J2*w(i);
        F = F+whx1*kappa1*(duha)*J1*w(i)+whx2*kappa2*(duhb)*J2*w(i)+...
            wh1*B*uha*J1*w(i)+wh2*B*uhb*J2*w(i);
        
        G = G+(-(a*wHx+B*wH)*H^(c2)*f*H*w(i))+...
            -(-(a1*whx1+B*wh1)*h(node-1)^(c2)*f1*J1*w(i))...
            -(-(whx2*a2+B*wh2)*h(node)^(c2)*f2*J2*w(i));
        I = I+wH*f*H*w(i)-wh1*f1*J1*w(i)-wh2*f2*J2*w(i);
        elseif strcmp(stabilizer, 'SUPG') == 1
            
            A = A+(a1*whx1)*h(node)^(c2)*(a1*duha+B*uha)*J1*w(i)...
                +(whx2*a2)*h(node-1)^(c2)*(a2*duhb+B*uhb)*J2*w(i);
            B = B+(a*wHx)*H^(c2)*(a*dUx+B*uH)*H*w(i);
            
            C = C-wHx*a*(uH*H*w(i))-wH*da*(uH*H*w(i));
            D = D+wHx*kappa*dUx*H*w(i)+wH*B*uH*H*w(i);
            
            E = E+whx1*a1*(-uha)*J1*w(i)+whx2*a2*(-uhb)*J2*w(i)+...
                wh1*da1*(-uha)*J1*w(i)+wh2*da2*(-uhb)*J2*w(i);
            F = F+whx1*kappa1*(duha)*J1*w(i)+whx2*kappa2*(duhb)*J2*w(i)+...
                wh1*B*uha*J1*w(i)+wh2*B*uhb*J2*w(i);
            
            G = G+(-(a*wHx)*H^(c2)*f*H*w(i))+...
                -(-(a1*whx1)*h(node-1)^(c2)*f1*J1*w(i))...
                -(-(a2*whx2)*h(node)^(c2)*f2*J2*w(i));
            I = I+wH*f*H*w(i)-wh1*f1*J1*w(i)-wh2*f2*J2*w(i);
        
        elseif strcmp(stabilizer, 'VMS') == 1
        A = A+(a1*whx1-B*wh1)*h(node)^(c2)*(a1*duha+B*uha)*J1*w(i)...
            +(whx2*a2-B*wh2)*h(node-1)^(c2)*(a2*duhb+B*uhb)*J2*w(i);
        B = B+(a*wHx-B*wH)*H^(c2)*(a*dUx+B*uH)*H*w(i);
        
        C = C-wHx*a*(uH*H*w(i))-wH*da*(uH*H*w(i));
        D = D+wHx*kappa*dUx*H*w(i)+wH*B*uH*H*w(i);
        
        E = E+whx1*a1*(-uha)*J1*w(i)+whx2*a2*(-uhb)*J2*w(i)+...
            wh1*da1*(-uha)*J1*w(i)+wh2*da2*(-uhb)*J2*w(i);
        F = F+whx1*kappa1*(duha)*J1*w(i)+whx2*kappa2*(duhb)*J2*w(i)+...
            wh1*B*uha*J1*w(i)+wh2*B*uhb*J2*w(i);
        
        G = G+(-(a*wHx-B*wH)*H^(c2)*f*H*w(i))+...
            -(-(a1*whx1-B*wh1)*h(node-1)^(c2)*f1*J1*w(i))...
            -(-(whx2*a2-B*wh2)*h(node)^(c2)*f2*J2*w(i));
        I = I+wH*f*H*w(i)-wh1*f1*J1*w(i)-wh2*f2*J2*w(i);   
        else
            fprintf('Should specify stabilizer');
        end
    end
    %%compute c1
    l1(node) = (A-B)+G;
    m1(node) = (C+D-(E+F));%+I;
    
end  
  for node = 2:nPoints-1
     c1(node) = (m1(node)/l1(node));
     if abs(l1(node))<1e-12
         c1(node) = 0;
     elseif c1(node)<0 
         c1(node) = 0;
         %c1(node) = -c1(node);
     end
  end
 
 %c1(nPoints) = c1(nPoints-1);
 
end
function [dUx, uHa,uHb, wHa,wHb] = LeastSquare( x,CordMatrix, uh,node,nPoints,bcValue)
global h
if node ==2
    value = bcValue(1);
    A = CordMatrix\uh;
    a2 = A(2);
    a1 = value-a2*CordMatrix(1,2);
    dUx = a2;
    uHa = a1+a2*x(node-1);
    uHb = a1+a2*x(node+1);
    wHa = 0;
    wHb = uHb;
elseif node == nPoints-1
    value = bcValue(2);
    A = CordMatrix\uh;
    a2 = A(2);
    a1 = value-a2*CordMatrix(3,2);
    
    dUx = a2;
    uHa = a1+a2*x(node-1);
    uHb = a1+a2*x(node+1);
    
    H = 2*h(node);
    uHa = uh(1);
    uHb = uh(3);    
    dUx = (uHb-uHa)/H;
    
    wHa = uHa;
    wHb = 0;
else
    A = CordMatrix\uh;
    a1 = A(1);
    a2 = A(2);
    
    dUx = a2;
    uHa = a1+a2*x(node-1);
    uHb = a1+a2*x(node+1);
    
%     H = 2*h(node);
%     uHa = uh(1);
%     uHb = uh(3);
%     dUx = (uHb-uHa)/H;
    
    wHa = uHa;
    wHb = uHb;
end
end