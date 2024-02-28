% Section 5.1 pdf
% Tab. 2, Fig. 3
% Gu R, Rakita Y, Lan L, et al. 
% Stretched Non-negative Matrix Factorization[J]. 
% arXiv preprint arXiv:2311.15173, 2023.
% PDF_mix_ZnSe-w_0.5_a_0.1__BaTiO3_c-w_0.5_a_0.2.csv
% BaTiO3_c: stretch 0 - 0.2, weight 0.8
% ZnSe: stretch 0 - 0.1 , weight 0.2 
% Data form https://github.com/yevgenyr/diffpysim
load('testsimupdflift20230320');

%% remove baseline through lifting
data_input = Data - min(Data(:));

%% stretch nmf
K=2; rho=1e11;
% W0=rand(K,size(data_input,2));
% X0=rand(size(data_input,1),K);
% A0=1+1e-2*randn(K,size(data_input,2));
[ Xstr,Wstr,Astr ] = stretchednmf_sparselh_spline1_lbcd( data_input,W0,X0,A0,rho,[],0 );
Wstr = diag(max(Wstr,[],2))\Wstr; Astr = diag(max(Astr,[],2))\Astr;
[ Xstr_add,Wstr_add,Astr_add ] =stretchednmf_spline1_addbaseline( Data,Wstr,Xstr,Astr,rho);

%% stretch NMF plot1
j=1;targ=pdfsolu1;list=zeros(1,size(targ,2));
for i=1:20
    if targ(:,i)'*Xstr_add(:,j)<0
        list(i) = 1;
    else
        list(i)=norm(targ(:,i)'*Xstr_add(:,j)/norm(targ(:,i))^2*targ(:,i)-Xstr_add(:,j))/norm(Xstr_add(:,j));
    end
end
[fmin,ind]=min(list);
i=ind;
% align1=stretchednmf_spline1_align( X0(:,j),targ(:,i)'*X0(:,j)/norm(targ(:,i))^2,targ(:,i),1.01);
align3=stretchednmf_spline1_align( targ(:,i),targ(:,i)'*Xstr_add(:,j)/norm(Xstr_add(:,j))^2,Xstr_add(:,j),1.01);
pearson(align3(:,1:2))
%% stretch NMF plot2
j=2;targ=pdfsolu2;list=zeros(1,size(targ,2));
for i=1:20
    if targ(:,i)'*Xstr_add(:,j)<0
        list(i)=1;
    else
        list(i)=norm(targ(:,i)'*Xstr_add(:,j)/norm(targ(:,i))^2*targ(:,i)-Xstr_add(:,j))/norm(Xstr_add(:,j));
    end
end
[fmin,ind]=min(list);
i=ind;
align4=stretchednmf_spline1_align( targ(:,i),targ(:,i)'*Xstr_add(:,j)/norm(Xstr_add(:,j))^2,Xstr_add(:,j),1.01);
pearson(align4(:,1:2))


