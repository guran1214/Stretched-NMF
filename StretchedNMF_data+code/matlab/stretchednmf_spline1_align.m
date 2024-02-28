function [ output,rate,Y,A ] = stretchednmf_spline1_align( MM,Y0,X0,A )

% align Y0 to MM
% min ||M(r)-X(ar)Y||_F, s.t. Y>=0
% variables: Y, a
% X(0) is in X
% written by Ran Gu (rgu@nankai.edu.cn)

[N,M]=size(MM);
K=size(Y0,1);
maxiter=300;
index=1:N;
X = X0; Y=Y0;

R=getR(X,Y,A);
fun=objfun(R,A);
histF=[inf*ones(1,0),fun];
numupdate=0;
for outiter=1:maxiter

    updateY2();    
    updateA2();
    numupdate=numupdate+1;
    R=getR(X,Y,A);
    fun=objfun(R,A);histF=[histF,fun];
    diffun=histF(end-1)-fun;
    if diffun<fun*1e-6&&outiter>=20
        break;
    end

end
output=[MM,MM+R,-R];% in reality, it is [X,targ, X-targ]
figure;plot([output(:,1:2),output(:,3)+min(output(:,1:2),[],'all')-max(output(:,3))])
rate=norm(R)/norm(MM);
fprintf('%f\n',rate)

%%

    function updateY2()
        y=zeros(K,1);
    for m=1:M
        T=zeros(N,K);
        for k=1:K
            T(:,k)=getAfun(A(k,m),X(:,k));
        end
        y=lsqnonneg(T(index,:),MM(index,m),y);
        Y(:,m)=y;
    end        
    end
    function updateA2()        
        options=optimoptions(@fmincon,'Algorithm','trust-region-reflective','Display','off','SpecifyObjectiveGradient',true,'CheckGradients',false,'HessianFcn','objective');
        A=fmincon(@funregu,A,[],[],[],[],0.1*ones(K,M),[],[],options);
    end
%%
    function [fun,gra,hess]=funregu(A)
        if nargout==1
            [AX]=getAfun_matrix(A,X,Y);
            RA=reshape(sum(reshape(AX,N*M,K),2),N,M)-MM;
            fun=objfun(RA,A);
        elseif nargout==2
            [AX,TX]=getAfun_matrix(A,X,Y);
            RA=reshape(sum(reshape(AX,N*M,K),2),N,M)-MM;
            fun=objfun(RA,A);
            gra=reshape(sum(TX(index,:).*repmat(RA(index,:),1,K)),M,K)';
        else
            [AX,TX,HX]=getAfun_matrix(A,X,Y);
            RA=reshape(sum(reshape(AX,N*M,K),2),N,M)-MM;
            fun=objfun(RA,A);
            gra=reshape(sum(TX(index,:).*repmat(RA(index,:),1,K)),M,K)';
            hess = zeros(M*K);
            for m=1:M
                Tx=TX(index,m+M*(0:K-1));hess((m-1)*K+1:m*K,(m-1)*K+1:m*K)=Tx'*Tx;
            end
            hess= hess +spdiags(reshape(reshape(sum(HX(index,:).*repmat(RA(index,:),1,K)),M,K)',M*K,1),0,M*K,M*K);
        end
    end
    function [Ax,Tx,Hx]=getAfun(a,x)
        ii=(0:N-1)'/a;
        II=floor(ii);
        I=II(II<N);
        i=ii(II<N);
%         Ax2=sparse([(1:length(I))';(1:length(I))'],[min(I+2,N);I+1],[(i-I).*(I+2<=N);1-i+I],length(I),N)*x;
        Ax2=x(I+1).*(1-i+I)+x(min(I+2,N)).*(i-I).*(I+2<=N);
        Ax=[Ax2;zeros(N-length(I),1)];
        %-------------------------------------------
        if nargout>1
            di=(-i/a);
        Tx2=x(I+1).*(-di)+x(min(I+2,N)).*di.*(I+2<=N);
        Tx=[Tx2;zeros(N-length(I),1)];
        if nargout>2
            ddi=-di/a+i*a^-2;
            Hx2=x(I+1).*(-ddi)+x(min(I+2,N)).*ddi.*(I+2<=N);
            Hx=[Hx2;zeros(N-length(I),1)];
        end
        end
    end
    function [Ax,Tx,Hx]=getAfun_matrix(A,X,Y)
        
        AA=repmat(reshape(A',1,M*K).^-1,N,1);
        ii=repmat((0:N-1)',1,K*M).*AA;
        YY=repmat(reshape(Y',1,M*K),N,1);
%         bias=repmat((0:K-1)*(N+1),N,M);
        bias=kron((0:K-1)*(N+1),ones(N,M));
        % X1,X2,XK
        X1=[X;zeros(1,K)];
        II=floor(ii);
        II1=min(II+1,N+1);
        II2=min(II1+1,N+1);
        iI=ii-II;
        II1_=II1+bias;
        II2_=II2+bias;
        XI1=reshape( X1(II1_), N,K*M);
        XI2=reshape( X1(II2_), N,K*M);
        Ax2=XI1.*(1-iI)+XI2.*(iI);
        Ax=Ax2;
        Ax=Ax.*YY;
        %-------------------------------------------
        if nargout>1
            di=(-ii.*AA);
        Tx2=XI1.*(-di)+XI2.*di;
        Tx=Tx2;
        Tx=Tx.*YY;
        if nargout>2
            ddi=-di.*AA*2;
            Hx2=XI1.*(-ddi)+XI2.*ddi;
            Hx=Hx2;
            Hx=Hx.*YY;
        end
        end
    end
    function R=getR(X,Y,A)
        R=-MM;
        for m=1:M
            r=R(:,m);
            for k=1:K
                r=r+Y(k,m)*getAfun(A(k,m),X(:,k));
            end
            R(:,m)=r;
        end
    end
    function fun=objfun(R,A)
        fun=0.5*norm(R(index,:),'fro')^2;
    end

end
