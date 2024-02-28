function [Y] = matchweight(X,T)

% min_{u,v} ||X.*(u*v')-T||^2
% ||U*X*V-T||^2
maxiter=1e4;
[n,k] = size(T);
%% method 1
u = ones(n,1);
v = ones(k,1);
for iter=1:maxiter
    u=((X.*T)*v)./(X.^2*v.^2);
    if norm(v.*(u'.^2*X.^2)'-(u'*(X.*T))')<1e-6
        break;
    end
    v=(u'*(X.*T))'./(u'.^2*X.^2)';
end
Y=X.*(u*v');
end