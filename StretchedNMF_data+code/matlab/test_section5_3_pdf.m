% Section 5.3 pdf
% Tab. 4, Fig. 5
% Gu R, Rakita Y, Lan L, et al. 
% Stretched Non-negative Matrix Factorization[J]. 
% arXiv preprint arXiv:2311.15173, 2023.
% PDF_mix_ZnSe-w_0.5_a_0.02__BaTiO3_c-w_0.5_a_0.04.csv
% BaTiO3_c: stretch 0 - 0.04, weight 0.5
% ZnSe: stretch 0 - 0.02 , weight 0.5
% Data form https://github.com/yevgenyr/diffpysim
load('test_section5_3_pdf.mat')

%% plot raw data
figure;plot(rr, Data,'k');
figure;plot(rr, Data+repmat((0:size(Data,2)-1)*10,size(Data,1),1),'k');
%% remove baseline (lift)
data_input = Data - min(Data(:));
figure;plot(rr,data_input+repmat(10*(0:size(data_input,2)-1),size(data_input,1),1),'k');
figure;plot(rr,data_input)

%% stretch nmf
% K=2; rho=1e11;
% W0=rand(K,size(data_input,2));
% X0=rand(size(data_input,1),K);
% A0=1+1e-2*randn(K,size(data_input,2));
[ Xstr,Wstr,Astr ] = stretchednmf_sparselh_spline1_lbcd( data_input,W0,X0,A0,rho,[],0 );
Wstr = diag(max(Wstr,[],2))\Wstr; Astr = diag(max(Astr,[],2))\Astr;
[ Xstr_add,Wstr_add,Astr_add ] =stretchednmf_spline1_addbaseline( Data,Wstr,Xstr,Astr,rho);

%% stretch NMF plot1
j=2;targ=pdfsolu1;list=zeros(1,size(targ,2));
for i=1:20
    if targ(:,i)'*Xstr_add(:,j)<0
        list(i) = 1;
    else
        list(i)=norm(targ(:,i)'*Xstr_add(:,j)/norm(targ(:,i))^2*targ(:,i)-Xstr_add(:,j))/norm(Xstr_add(:,j));
    end
end
[fmin,ind]=min(list);
i=ind;
align3=stretchednmf_spline1_align( targ(:,i),targ(:,i)'*Xstr_add(:,j)/norm(Xstr_add(:,j))^2,Xstr_add(:,j),1.01);
pearson(align3(:,1:2))
%% stretch NMF plot2
j=1;targ=pdfsolu2;list=zeros(1,size(targ,2));
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


