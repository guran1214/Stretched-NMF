function [x] = lsqnonneg(C,d,x,nZeros) 
% tol = 10*eps*norm(C,1)*length(C); 
tol = 1e-14;
n = size(C,2);
if nargin<4
nZeros = zeros(n,1); 
end
count=0;
wz = nZeros; 
if nargin<3
    x=nZeros;
end
Z=x<=tol;
P=~Z;
z = nZeros;

outeriter = 0; 
iter = 0; 
itmax = 3*n; 


 z(P) = C(:,P)\d; count=count+1;
while 1
    while any(z(P) <= 0) 
        iter = iter + 1; 
        Q = (z <= 0) & P; 
        alpha = min(x(Q)./(x(Q) - z(Q)));
        x = x + alpha*(z - x); 
        Z = ((abs(x) < tol) & P) | Z; 
        P = ~Z;
        z = nZeros;
        z(P) = C(:,P)\d; count=count+1;
    end
    x = z;
    resid = d - C*x; 
    w = C'*resid;
    
    if ~(any(Z) && any(w(Z) > tol))
        break;
    end
    outeriter = outeriter + 1; 
    z = nZeros;
    wz(P) = -Inf; 
    wz(Z) = w(Z); 
    [~,t] = max(wz); 
    P(t) = true; 
    Z(t) = false; 
    
    z(P) = C(:,P)\d; count=count+1;
end
% lambda = w;          
% resnorm = resid'*resid;
    
end





