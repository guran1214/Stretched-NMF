function [ X,Y,A,fun ] = stretchednmf_sparselh_spline1_lbcd( MM,Y0,X0,A,rho,index,eta,maxiter )

%
% min ||M(r)-X(ar)Y||_F, X>=0, Y>=0
% variables: X, Y, a
% X(0) is in X
% if fillzero, fill zeros, else, fill X(end,:)
% may use mkr_box(Hungerlander and Rendl, SiOpt 2015) for updating Y
% written by Ran Gu (rgu@nankai.edu.cn)
% citation: Gu R, Rakita Y, Lan L, et al. 
% Stretched Non-negative Matrix Factorization[J]. 
% arXiv preprint arXiv:2311.15173, 2023.

isplot=1;
isdisp=1;
fillzero=0;
[N,M]=size(MM);
K=size(Y0,1);
if nargin<8
    maxiter=300;
end
if nargin<7
    eta=1;
end
if nargin<6 || isempty(index)
    index=1:N;
end
Nindex=length(index);
if nargin<5
rho=1e18;
end
if isplot
    figure;
end
fprintf('\n rho = %g , K = %i \n',rho,K);

P=0.25*sparse([1:M-2,1:M-2,1:M-2],[1:M-2,2:M-1,3:M],[ones(1,M-2),-2*ones(1,M-2),ones(1,M-2)],M-2,M); % 2nd order spline
% P=0.5*sparse([1:M-1,1:M-1],[1:M-1,2:M],[ones(1,M-1),-ones(1,M-1)],M-1,M); % 1st order spline
PP=P'*P;PPPP=sparse(M*K);
for k=1:K
    PPPP(M*(k-1)+1:M*k,M*(k-1)+1:M*k)=PP;
end
seq=reshape(reshape(1:M*K,M,K)',M*K,1);
PPPP=PPPP(seq,seq);
if nargin<4
    A=ones(K,M)+randn(K,M)*1e-3;
end
if nargin<3
    X0=rand(N,K);
end
Y=max(0,Y0);X=max(0,X0);


R=getR(X,Y,A);
fun=objfun(R,A);
histF=[inf*ones(1,0),fun];
if isdisp
fprintf('%i\t%g\t%g\t%i\n',0,fun,fun-0.5*rho*norm(P*A','fro')^2-eta*sum(X(:).^0.5),0);
end
if isplot
clf;drawfig;
end
numupdate=0;
preX=[];curX=X;pregraX=[];curgraX=[];
for outiter=1:maxiter
    for iter=1:4
    updateX();
    numupdate=numupdate+1;
    R=getR(X,Y,A);
    fun=objfun(R,A);histF=[histF,fun];
    if outiter==1&&iter==1
        diffun=histF(end-1)-fun;
    end
    preX=curX;curX=X;pregraX=curgraX;
    
    updateY2();
    numupdate=numupdate+1;
    R=getR(X,Y,A);
    fun=objfun(R,A);histF=[histF,fun];
    if histF(end-2)-fun<diffun*1e-3
        break;
    end
    end
    
    updateA2();
    numupdate=numupdate+1;
    R=getR(X,Y,A);
    fun=objfun(R,A);histF=[histF,fun];
    diffun=histF(end-1)-fun;
    if isdisp
    fprintf('%i\t%g\t%g\t%i\n',numupdate,fun,fun-0.5*rho*norm(P*A','fro')^2-eta*sum(X(:).^0.5),outiter);
    end
    if isplot
    clf;drawfig;
    end
    if diffun<fun*1e-6&&outiter>=20
        break;
    end

end

%%
    function drawfig()
        maxy=max(Y,[],2);X_=X.*repmat(maxy',N,1);
        subplot('position',[0.55 0.05 0.4 0.9]);plot(repmat([0,cumsum(max(X_(:,1:end-1)))],N,1)+X_);xlim([0,N]);
        subplot('position',[0.05 0.55 0.4 0.4]);plot(A');
        legend(arrayfun(@(mode) sprintf('stretch %d', mode), 1:K, 'UniformOutput', false))
        subplot('position',[0.05 0.05 0.4 0.4]);plot((Y./repmat(maxy,1,M))');
        legend(arrayfun(@(mode) sprintf('weight %d', mode), 1:K, 'UniformOutput', false))
        pause(1e-4)
    end
    function updateX()
        [AX]=getAfun_matrix(A,X,Y);
        RA=reshape(sum(reshape(AX,N*M,K),2),N,M);
        RR=RA-MM;
        GraX=full(getATfun_matrix(A,RR,Y));
        curgraX=GraX;
        L0=max(eig(Y'*Y))*max([A(:);A(:).^-1]);
        if outiter==1&&iter==1
            L=L0;
        else
            L=sum((GraX-pregraX).*(curX-preX),'all')/norm(curX-preX,'fro')^2;
            if L<=0,L=L0; end
        end
        while 1
            X_=curX-GraX*L^-1;
            X=rooth(-X_,eta/(2*L)).^2;
            X=(X.^2*L/2-L*X.*X_+eta*X.^0.5<0).*X;
        if histF(end)-objfun(getR(max(0,X),Y,A),A)>0
            break;
        end
        L=L*2;
        if isinf(L)
            break;
        end
        end
    end
    function updateY2()
        y=zeros(K,1);
    for m=1:M
        T=zeros(N,K);
        for k=1:K
            T(:,k)=getAfun(A(k,m),X(:,k));
        end
        y = mkr_box( T(index,:)'*T(index,:), -T(index,:)'*MM(index,m), zeros(K,1), 1*ones(K,1), 1, [], [], [], [], 0, 10^-10, 0, 0, 0, 1);
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
            gra=reshape(sum(TX(index,:).*repmat(RA(index,:),1,K)),M,K)'+rho*A*P'*P;
        else
            [AX,TX,HX]=getAfun_matrix(A,X,Y);
            RA=reshape(sum(reshape(AX,N*M,K),2),N,M)-MM;
            fun=objfun(RA,A);
            gra=reshape(sum(TX(index,:).*repmat(RA(index,:),1,K)),M,K)'+rho*A*P'*P;
            hess = zeros(M*K);
            for m=1:M
                Tx=TX(index,m+M*(0:K-1));hess((m-1)*K+1:m*K,(m-1)*K+1:m*K)=Tx'*Tx;
            end
            hess= hess +spdiags(reshape(reshape(sum(HX(index,:).*repmat(RA(index,:),1,K)),M,K)',M*K,1),0,M*K,M*K)+rho*PPPP;
        end
    end
    function AT=getATfun_matrix(A,R,Y)
        AA=repmat(reshape(A,1,M*K).^-1,Nindex,1);
        ii=repmat((index-1)',1,K*M).*AA;
        YY=repmat(reshape(Y,1,M*K),Nindex,1);
        II=floor(ii);
        II1=min(II+1,N+1);
        II1_=II1;
        II2=min(II1+1,N+1);
        II2_=II2;
        iI=ii-II;
        repm=repmat(1:K,Nindex,M);
        kro=kron(R(index,:),ones(1,K));
        kroiI=kro.*(iI);
        iIYY=(iI-1).*YY;
        x2=sparse(II1_,repm,kro.*-iIYY,(N+1),K);
        x3=sparse(II2_,repm,kroiI.*YY,(N+1),K);
        if fillzero==0
            x2(N,:)=x2(N,:)+x2((N+1),:);
            x3(N,:)=x3(N,:)+x3((N+1),:);
        end
        x2((N+1),:)=[];
        x3((N+1),:)=[];
        AT=x2+x3;
    end
    function [Ax,Tx,Hx]=getAfun(a,x)
        ii=(0:N-1)'/a;
        II=floor(ii);
        I=II(II<N);
        i=ii(II<N);
        if fillzero
            Ax2=x(I+1).*(1-i+I)+x(min(I+2,N)).*(i-I).*(I+2<=N);
            Ax=[Ax2;zeros(N-length(I),1)]; % fill zeros for the tail
        else
            Ax2=x(I+1).*(1-i+I)+x(min(I+2,N)).*(i-I);
            Ax=[Ax2;Ax2(end)*ones(N-length(I),1)]; % fill constant for the tail
        end
        %-------------------------------------------
        if nargout>1
            di=(-i/a);
            if fillzero
                Tx2=x(I+1).*(-di)+x(min(I+2,N)).*di.*(I+2<=N);
            else
                Tx2=x(I+1).*(-di)+x(min(I+2,N)).*di;
            end
            Tx=[Tx2;zeros(N-length(I),1)];
        if nargout>2
            ddi=-di/a+i*a^-2;
            if fillzero
                Hx2=x(I+1).*(-ddi)+x(min(I+2,N)).*ddi.*(I+2<=N);
            else
                Hx2=x(I+1).*(-ddi)+x(min(I+2,N)).*ddi;
            end
            Hx=[Hx2;zeros(N-length(I),1)];
        end
        end
    end
    function [Ax,Tx,Hx]=getAfun_matrix(A,X,Y)
        
        AA=repmat(reshape(A',1,M*K).^-1,N,1);
        ii=repmat((0:N-1)',1,K*M).*AA;
        YY=repmat(reshape(Y',1,M*K),N,1);
        bias=kron((0:K-1)*(N+1),ones(N,M));
        % X1,X2,XK
        if fillzero
        X1=[X;zeros(1,K)];
        else
            X1=[X;X(end,:)];
        end
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
        fun=0.5*norm(R(index,:),'fro')^2+0.5*rho*norm(P*A','fro')^2+eta*sum(X(:).^0.5);
    end

end

function y =  rooth(p,q)
% x^3+px+q=0 the maximal root
if q==0
    y=max(0,-p).^0.5;
else
delta=(q/2).^2+(p/3).^3;
d=(delta).^0.5;
a1=(-q/2+d).^(1/3);
a2=(-q/2-d).^(1/3);
w=(3^0.5*1i-1)/2;
y1=a1+a2;
y2=w*a1+w^2*a2;
y3=w^2*a1+w*a2;
y=max(real(y1),max(real(y2),real(y3))).*(delta<0);
end
end