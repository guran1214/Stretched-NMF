% Section 5.4 xrd
% Tab. 7, Fig. 8 & 9
% Gu R, Rakita Y, Lan L, et al. 
% Stretched Non-negative Matrix Factorization[J]. 
% arXiv preprint arXiv:2311.15173, 2023.
% XRD_MgMnO_YCl real data
% Data from O’Nolan D, Huang G, Kamm GE, et al (2020) 
% A thermal-gradient approach to variable-temperature measurements resolved in space. 
% Journal of Applied Crystallography 53(3):662–670.
load('test_section5_4_xrd.mat'); 

%% data 668K K=2  sparse stretch nmf
K=2;rho=1e12;eta=610;
% W0=rand(K,size(Data668_input,2));
% X0=rand(size(Data668_input,1),K);
% A0=1+1e-2*randn(K,size(Data668_input,2));
% J=[1,4];b=perms(1:2);
[ X,W,A ] = stretchednmf_sparselh_spline1_lbcd( Data668_input,W0,X0,A0,rho,i1+2:size(Data668_input,1),eta );
% normalization
Wspstr=W;Aspstr=A;
Wspstr = diag(max(Wspstr,[],2))\Wspstr;
Aspstr = diag(max(Aspstr,[],2))\Aspstr;
[ Xspstr,Wspstr,Aspstr ] =stretchednmf_sparselh_spline1_normx( Data668_input,Wspstr,X,Aspstr,rho,i1+2:size(Data668_input,1),eta  );
X=Xspstr;
%% align
j=2;
i=1;[align1,r1]=stretchednmf_spline1_align( Xsim_q(:,J(i)),Xsim_q(:,J(i))'*X(:,b(j,i))/norm(X(:,b(j,i)))^2,X(:,b(j,i)),0.999);
pearson(align1(:,1:2))
i=2;[align2,r2]=stretchednmf_spline1_align( Xsim_q(:,J(i)),Xsim_q(:,J(i))'*X(:,b(j,i))/norm(X(:,b(j,i)))^2,X(:,b(j,i)),0.999);
pearson(align2(:,1:2))

%% conventional nmf
[ Xnmf,Wnmf ] = nmf_twophase( Data668_input,X0,W0);
% normalization
Wnmf = diag(max(Wnmf,[],2))\Wnmf;
[ Xnmf_add,Wnmf_add,Anmf_add ] = stretchednmf_spline1_addbaseline( Data668_input,Wnmf,Xnmf,ones(size(Wnmf)), 0);
%% align
j=2;i=1;[align3,r3]=stretchednmf_spline1_align( Xsim_q(:,J(i)),Xsim_q(:,J(i))'*Xnmf_add(:,j)/norm(Xnmf_add(:,j))^2,Xnmf_add(:,j),0.999);
pearson(align3(:,1:2))
j=1;i=2;[align4,r4]=stretchednmf_spline1_align( Xsim_q(:,J(i)),Xsim_q(:,J(i))'*Xnmf_add(:,j)/norm(Xnmf_add(:,j))^2,Xnmf_add(:,j),0.999);
pearson(align4(:,1:2))

%% stretch nmf  without sparse
[ Xstr,Wstr,Astr ] = stretchednmf_sparselh_spline1_lbcd( Data668_input,W0,X0,A0,rho,i1+2:size(Data668_input,1),0 );
% normalization
Wstr = diag(max(Wstr,[],2))\Wstr;
Astr = diag(max(Astr,[],2))\Astr;
[ Xstr_add,Wstr_add,Astr_add ] =stretchednmf_sparselh_spline1_normx( Data668_input,Wstr,Xstr,Astr,rho,i1+2:size(Data668_input,1),0  );
X=Xstr_add;
%% align
j=1;i=1;[align5,r5]=stretchednmf_spline1_align( Xsim_q(:,J(i)),Xsim_q(:,J(i))'*X(:,j)/norm(X(:,j))^2,X(:,j),0.999);
pearson(align5(:,1:2))
j=2;i=2;[align6,r6]=stretchednmf_spline1_align( Xsim_q(:,J(i)),Xsim_q(:,J(i))'*X(:,j)/norm(X(:,j))^2,X(:,j),0.999);
pearson(align6(:,1:2))

%% compare weights
Wmod_spstr=matchweight(Wspstr',Wtrue);
Wmod_str=matchweight(Wstr',Wtrue);
Wmod_nmf=matchweight(Wnmf([2,1],:)',Wtrue);
figure;plot(tgrid(5:end),[Wtrue(:,1),Wmod_nmf(:,1),Wmod_str(:,1),Wmod_spstr(:,1)]);ylim([0,100]);xlim([368,668]);
legend('Rietveld fit','Conventional NMF','Stretched NMF','Sparse Stretched NMF');title('percentage of weight of tYOCl')

figure;plot(tgrid(5:end),[Wtrue(:,2),Wmod_nmf(:,2),Wmod_str(:,2),Wmod_spstr(:,2)]);ylim([0,100]);xlim([368,668]);
legend('Rietveld fit','Conventional NMF','Stretched NMF','Sparse Stretched NMF');title('percentage of weight of MgMn2O4')

    