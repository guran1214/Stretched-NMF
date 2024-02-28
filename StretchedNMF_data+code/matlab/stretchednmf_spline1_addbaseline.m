function [ X,Y,A ] = stretchednmf_spline1_addbaseline( MM,Y0,X0,A,rho,index )

%
% min ||M(r)-X(ar)Y||_F, X>=0
% variable: X
% X(0) is in X
% written by Ran Gu (rgu@nankai.edu.cn)
isplot = 1;
[N,M]=size(MM);
K=size(Y0,1);
maxiter=200;
if nargin<6
    index=1:N;
end
Nindex=length(index);
if nargin<5
rho=1e18;
end
if isplot
figure;
fprintf('\n rho = %g , K = %i \n',rho,K);
end
P=0.25*sparse([1:M-2,1:M-2,1:M-2],[1:M-2,2:M-1,3:M],[ones(1,M-2),-2*ones(1,M-2),ones(1,M-2)],M-2,M); % 2nd order spline
% P=0.5*sparse([1:M-1,1:M-1],[1:M-1,2:M],[ones(1,M-1),-ones(1,M-1)],M-1,M); % 1st order spline

PP=P'*P;PPPP=sparse(M*K);
for k=1:K
    PPPP(M*(k-1)+1:M*k,M*(k-1)+1:M*k)=PP;
end
seq=reshape(reshape(1:M*K,M,K)',M*K,1);
PPPP=PPPP(seq,seq);

X = X0; Y = Y0;
R=getR(X,Y,A);
fun=objfun(R,A);
histF=[inf*ones(1,0),fun];
if isplot
fprintf('%i\t%g\t%g\t%i\n',0,fun,fun-0.5*rho*norm(P*A','fro')^2,0);
clf;drawfig;
end
numupdate=0;
preX=[];curX=X;pregraX=[];curgraX=[];
    for iter=1:maxiter
    updateX();
    numupdate=numupdate+1;
    R=getR(X,Y,A);
    fun=objfun(R,A);histF=[histF,fun];
    
    preX=curX;curX=X;pregraX=curgraX;
    if isplot
        fprintf('%i\t%g\t%g\t%i\n',numupdate,fun,fun-0.5*rho*norm(P*A','fro')^2,1);
        clf;drawfig;
    end
    if iter>1&&(histF(end-1)-fun)/fun<1e-6
        break;
    end
    end


%%
    function drawfig()
        maxy=max(Y,[],2);X_=X.*repmat(maxy',N,1);
        subplot('position',[0.55 0.05 0.4 0.9]);plot(repmat([0,cumsum(max(X_(:,1:end-1))-min(X_(:,2:end)))],N,1)+X_);xlim([0,N]);
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
        if iter==1
            L=L0;
        else
            L=sum((GraX-pregraX).*(curX-preX),'all')/norm(curX-preX,'fro')^2;
            if L<=0,L=L0; end
        end
        while 1
            X=curX-GraX*L^-1;
        if histF(end)-objfun(getR(X,Y,A),A)>0
            break;
        end
        L=L*2;
        if isinf(L)
            break;
        end
        end
    end

%%

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
        x2((N+1),:)=[];
        x3=sparse(II2_,repm,kroiI.*YY,(N+1),K);
        x3((N+1),:)=[];
        AT=x2+x3;
    end
    function [Ax,Tx,Hx]=getAfun(a,x)
        ii=(0:N-1)'/a;
        II=floor(ii);
        I=II(II<N);
        i=ii(II<N);
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
        fun=0.5*norm(R(index,:),'fro')^2+0.5*rho*norm(P*A','fro')^2;
    end

end