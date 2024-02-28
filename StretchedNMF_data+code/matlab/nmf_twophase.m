function [X,Y,histtime,histnorm,histf]=nmf_twophase( M,X0,Y0)
% NMF solver
% min ||M-XY||_F^2, X>=0, Y>=0
% Citation: Gu R, Billinge S J L, Du Q. 
% A fast two-stage algorithm for non-negative matrix factorization in smoothly varying data[J]. 
% Acta Crystallographica Section A: Foundations and Advances, 2023, 79(2).
M=max(M,0);
[n,m]=size(M);
k=size(Y0,1);
maxiter=200;
Y=Y0;X=X0;

elapse=cputime;
% 1st phase
anlsas();

rho=1e-3;
X=max(X,rho);
Y=max(Y,rho);
% 2nd phase
[X,Y,histt2,histg2,histf2]=nmf_int_precor( M,X,Y,elapse);
histtime=[histtime,histt2];
histnorm=[histnorm,histg2];
histf=[histf,histf2];

    function anlsas()
        
        resM=(X*Y-M);
        graX=resM*Y';
        graY=X'*resM;
        gnorm=max(norm([graX(graX<0);graY(graY<0)]),norm([max(graX(:),0).*X(:);max(graY(:),0).*Y(:)]));
        histnorm=gnorm;
        histf=norm(resM,'fro')^2/2;
        histtime=0;
        
        flag=0;
        iter=1;
        while iter<=maxiter
            if iter>1
                preX=X;
                preY=Y;
            end
            %update X
            [X]=getW(M',Y');X=X';
            %updata Y
            [Y]=getW(M,X);
            
            % balance
            scal=(sum(Y,2)./sum(X)').^0.5;
            X=repmat(scal',n,1).*X;
            Y=repmat(scal.^-1,1,m).*Y;
            
            resM=(X*Y-M);
            graX=resM*Y';
            graY=X'*resM;
            
            gnorm=max(norm([graX(graX<0);graY(graY<0)]),norm([max(graX(:),0).*X(:);max(graY(:),0).*Y(:)]));
            histtime=[histtime,cputime-elapse];
            histnorm=[histnorm,gnorm ];
            histf=[histf,norm(resM,'fro')^2/2];
            
            fprintf('%i\t%f\t%f\n',iter,norm(resM,'fro')^2/2,gnorm);
            if (iter>1)&&(all((preX(:)>0)==(X(:)>0))&&all((preY(:)>0)==(Y(:)>0))) || iter==50
                flag=flag+1;
                if flag==2
                    break;
                end
            else
                flag=0;
            end
            iter=iter+1;
        end
    end

end

function [W]=getW(Z,X)
[n,m]=size(Z);
R=chol(X'*X);
k=size(X,2);
D=R'\(X'*Z);
W=zeros(k,m);
w=zeros(k,1);
seq=1:m;
nZeros=zeros(k,1);

[w]=lsqnonneg(R,D(:,seq(1)),w,nZeros);
W(:,seq(1))=w;
for i=2:m
    [w]=lsqnonneg(R,D(:,seq(i)),w,nZeros);
    W(:,seq(i))=w;
end

end
