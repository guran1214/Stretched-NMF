% Section 5.2 xrd
% Tab. 3, Fig. 4
% Gu R, Rakita Y, Lan L, et al. 
% Stretched Non-negative Matrix Factorization[J]. 
% arXiv preprint arXiv:2311.15173, 2023.
% XRD_mix_ZnSe-w_0.2_a_0.1__BaTiO3_c-w_0.8_a_0.2.csv
% BaTiO3_c: stretch 0 - 0.2, weight 0.8
% ZnSe: stretch 0 - 0.1 , weight 0.2  
% Data form https://github.com/yevgenyr/diffpysim
load('test_section5_2_xrd.mat')

%% stretch NMF
rho=1e13;
eta=0;
K=2;
% W0=rand(K,size(Data_input,2));
% X0=rand(size(Data_input,1),K);
% A0=1+1e-2*randn(K,size(Data_input,2));
%%
[ X,W,A ] = stretchednmf_sparselh_spline1_lbcd( Data_input,W0,X0,A0,rho,[],0 );

%% align
j=2;targ=xrdsolu1;list=zeros(1,size(targ,2));
for i=1:20
    list(i)=norm(targ(:,i)'*X(:,j)/norm(targ(:,i))^2*targ(:,i)-X(:,j))/norm(X(:,j));
end
[fmin,ind]=min(list);
i=ind;
% align1=stretchednmf_spline1_align( X0(:,j),targ(:,i)'*X0(:,j)/norm(targ(:,i))^2,targ(:,i),1.01);
align1=stretchednmf_spline1_align( targ(:,i),targ(:,i)'*X(:,j)/norm(X(:,j))^2,X(:,j),0.99);
pearson(align1(:,1:2))
%% align
j=1;targ=xrdsolu2;list=zeros(1,size(targ,2));
for i=1:20
    list(i)=norm(targ(:,i)'*X(:,j)/norm(targ(:,i))^2*targ(:,i)-X(:,j))/norm(X(:,j));
end
[fmin,ind]=min(list);
i=ind;
% align1=stretchednmf_spline1_align( X0(:,j),targ(:,i)'*X0(:,j)/norm(targ(:,i))^2,targ(:,i),1.01);
align2=stretchednmf_spline1_align( targ(:,i),targ(:,i)'*X(:,j)/norm(X(:,j))^2,X(:,j),0.99);
pearson(align2(:,1:2))




