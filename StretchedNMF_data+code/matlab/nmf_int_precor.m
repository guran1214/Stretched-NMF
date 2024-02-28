function [ X,Y,histtime,histnorm,histf] = nmf_int_precor( M,X0,Y0,time0 )
% NMF solver: interior point method
% min ||M-XY||_F, X>=0, Y>=0
% Citation: Gu R, Billinge S J L, Du Q. 
% A fast two-stage algorithm for non-negative matrix factorization in smoothly varying data[J]. 
% Acta Crystallographica Section A: Foundations and Advances, 2023, 79(2).
[n,m]=size(M);
k=size(X0,2);
maxiter=200;
tol=1e-6;
rho=tol;
eta=0;
histtime=[];
histnorm=[];
histf=[];
Y=Y0;
X=X0;
resM=(X*Y-M);
    graX=resM*Y';
    graY=X'*resM;

    r=1*ones(k*n,1)*max(abs(graX(:)));
    s=1*ones(k*m,1)*max(abs(graY(:)));

mu=inf;sigma=inf;
for iter=1:maxiter
    resM=(X*Y-M);
    %% balance    
    graX=resM*Y';
    graY=X'*resM;
    grax=reshape(graX',[],1);
    gray=graY(:);
    t=sqrt((norm(gray-s)/sqrt(m*k))/(norm(grax-r)/sqrt(n*k)));
    X=X/t;Y=Y*t;r=r*t;s=s/t;
    %%
    graX=resM*Y';
    graY=X'*resM;
    grax=reshape(graX',[],1);
    gray=graY(:);
    x=reshape(X',[],1);
    y=Y(:);
    fprintf('%i\t%f\t%f\t%f\t%f\t%f\n',iter,norm(resM,'fro')^2/2,norm([grax-r;gray-s]),...
        norm([x.*r;y.*s]),mu,sigma);
    histtime=[histtime,cputime-time0];
    histnorm=[histnorm,max(norm([grax-r;gray-s]),norm([x.*r;y.*s]))];
    histf=[histf,norm(resM,'fro')^2/2];
    if max(norm([grax-r;gray-s]),norm([x.*r;y.*s]))<tol||(iter>2&&abs(histf(end)-histf(end-1))/abs(histf(end))<tol*1e-4)
        break;
    end
    %% Xrow Ycolumn
    mu=0;
    nu=0;
    eta=sigma<1e-2;
%     eta=0;
    
    Z=Y*Y';
    W=X'*X;

    if iter==1
    Q2=zeros(k*m);
    end
       
    for i=1:m
        Q2((i-1)*k+1:i*k,(i-1)*k+1:i*k)=W+diag(s((i-1)*k+1:i*k)./y((i-1)*k+1:i*k)+rho);
    end
    
    
    if iter==1
        CP1i=zeros(m*k,n*k);Ri=zeros(n*k,k);Rit=Ri;
    end
    Rx=r./x;
    for j=1:n
        T=inv(chol(Z+diag(Rx((j-1)*k+1:j*k)+rho)));
        Ri((j-1)*k+1:j*k,:)=T;Rit((j-1)*k+1:j*k,:)=T';
        xx=X(j,:);
        F=repmat(xx',m,k).* reshape(repmat(T'*Y,k,1),k,m*k)'+...
            reshape(reshape(T',[],1)*(eta*resM(j,:)),k,[])';
        CP1i(:,(j-1)*k+1:j*k)=F;
    end
    
    %% predict
    beq1=-grax;
    beq2=-gray;
%     P1ib1=P1i*beq1;
    P1ib1=sum(Rit.*reshape(repmat(reshape(beq1,k,n),k,1),k,k*n)',2);
    [L,U]=lu(Q2-CP1i*CP1i');
    solu2=U\(L\(beq2-CP1i*P1ib1));
    solu1=sum(Ri.*reshape(repmat(reshape(P1ib1-CP1i'*solu2,k,n),k,1),k,k*n)',2);
    dr=mu./x-r-r.*solu1./x;
    ds=nu./y-s-s.*solu2./y;
    step1=min( max(-[solu1;solu2]./[x;y])^-1,1);
    step2=min( max(-[dr;ds]./[r;s])^-1,1);
    muaff=([x;y]+step1*[solu1;solu2])'*([r;s]+step2*[dr;ds])/(n*k+m*k);
    mu=([x;y]'*[r;s])/(n*k+m*k);
    sigma=(muaff/mu)^3;
    sigma=min(sigma,0.99);
    %% correct
    aff1=solu1.*dr;aff2=solu2.*ds;
    beq1=(mu*sigma-aff1)./x-grax;
    beq2=(mu*sigma-aff2)./y-gray;
    
    P1ib1=sum(Rit.*reshape(repmat(reshape(beq1,k,n),k,1),k,k*n)',2);
    solu2=U\(L\(beq2-CP1i*P1ib1));
    solu1=sum(Ri.*reshape(repmat(reshape(P1ib1-CP1i'*solu2,k,n),k,1),k,k*n)',2);
    dr=(mu*sigma-aff1)./x-r-r.*solu1./x;
    ds=(mu*sigma-aff2)./y-s-s.*solu2./y;
    
    step1=min( 0.9*max(-[solu1;solu2]./[x;y])^-1,1);
    step2=min( 0.9*max(-[dr;ds]./[r;s])^-1,1);
    
    dX=reshape(solu1,k,n)';
    dY=reshape(solu2,k,m);
    %%
    X1=X+dX*step1;
    Y1=Y+dY*step1;
    r1=r+dr*step2;%
    s1=s+ds*step2;%
    
    %%
    f0=norm(X*Y-M,'fro')^2/2-mu*sigma*sum(log(X(:)))-mu*sigma*sum(log(Y(:)));
    gra0=[graX(:)-mu*sigma*X(:).^-1;graY(:)-mu*sigma*Y(:).^-1];
    d0=[dX(:);dY(:)];
    if gra0'*d0>=0
        continue;
    end

    %%
    X=X1;Y=Y1;
    r=r1;s=s1;
end


end
