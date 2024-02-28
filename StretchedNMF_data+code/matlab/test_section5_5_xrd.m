% Section 5.5 xrd
% Fig. 10
% Gu R, Rakita Y, Lan L, et al. 
% Stretched Non-negative Matrix Factorization[J]. 
% arXiv preprint arXiv:2311.15173, 2023.
% XRD_CuCl2_Na2Se2 real data
% Data form Martinolich AJ, Kurzman JA, Neilson JR (2015)
% Polymorph Selectivity of Superconducting CuSe2 Through Kinetic Control of Solid-State Metathesis. 
% Journal of the American Chemical Society 137(11):3827–3833.
load('test_section5_5_xrd.mat'); 

%% sparse streched nmf
K=6;
% X0=rand(size(Data,1),K);
% W0=rand(K,size(Data,2));
% A0= 1+randn(K,size(Data,2))*1e-3;
eta=1.2e6; 
rho=1e21; 
[ X3,W3,A3 ] = stretchednmf_sparselh_spline1_lbcd( Data,W0,X0,A0,rho,[],eta ); 
% normalization
X3=X3/diag(W3'\ones(length(id),1));
W3 = diag(W3'\ones(length(id),1))*W3;  
seq = [1 2 5 6 3 4]; X3=X3(:,seq); W3=W3(seq,:); A3=A3(seq,:);
% show solution
figure;subplot('position',[0.05 0.05+0.10*6 0.3 0.13]);plot(temp(id),W3');xlim([min(temp(id)),max(temp(id))]);ylim([0,1]);title('weight')
subplot('position',[0.05 0.05+0.10*6+0.04+0.13 0.3 0.13]);plot(temp(id),A3');xlim([min(temp(id)),max(temp(id))]);title('stretch')
subplot('position',[0.4 0.05 0.55 0.9]);plot(X3+repmat(cumsum([0,max(X3(:,1:end-1),[],1)]),length(rr),1));xlim([1000,length(rr)]);
title('Sparse Stretched NMF with background removal')
for i =1:K
    subplot('position',[0.05 0.05+0.1*(i-1) 0.3 0.1]);plot(temp(id),W3(i,:)');xlim([min(temp(id)),max(temp(id))]);ylim([0,1]);
end
%% align xrd
move=19;
j = 2; targ = tr1x; targ = [targ(move+1:end);zeros(move,1)]; % align 1 NaCl
align1=stretchednmf_spline1_align( targ,targ'*X3(:,j)/norm(X3(:,j))^2,X3(:,j),1.01);
pearson(align1(:,1:2))
move=18;
j = 4; targ = tr2x; targ = [targ(move+1:end);zeros(move,1)];% align 2 CuSe
align2=stretchednmf_spline1_align( targ,targ'*X3(:,j)/norm(X3(:,j))^2,X3(:,j),0.99);
pearson(align2(:,1:2))
move=14;
j = 5; targ = tr3x; targ = [targ(move+1:end);zeros(move,1)];% align 3 Cu2Se
align3=stretchednmf_spline1_align( targ,targ'*X3(:,j)/norm(X3(:,j))^2,X3(:,j),0.99);
pearson(align3(:,1:2))
move=0;
j = 1; targ = tr4x; targ = [targ(move+1:end);zeros(move,1)];% align 4 Se
align4=stretchednmf_spline1_align( targ,targ'*X3(:,j)/norm(X3(:,j))^2,X3(:,j),0.95);
pearson(align4(:,1:2))
move=22;
j = 3; targ = tr5x; targ = [targ(move+1:end);zeros(move,1)];% align 5 pyrite
align5=stretchednmf_spline1_align( targ,targ'*X3(:,j)/norm(X3(:,j))^2,X3(:,j),1.01);
pearson(align5(:,1:2))
move=0;
j = 6; targ = tr6x; targ = [targ(move+1:end);zeros(move,1)];% align 6 marcasite
align6=stretchednmf_spline1_align( targ,targ'*X3(:,j)/norm(X3(:,j))^2,X3(:,j),1.01);
pearson(align6(:,1:2))

%% align weight 
id1=3:215;
j = 2; targ = true1w;% align 1 NaCl
align1w=[targ,targ'*W3(j,id1)'/norm(W3(j,id1)')^2*W3(j,id1)'];
figure;plot(truetemp,align1w);
[norm(align1w(:,2)-align1w(:,1))/norm(align1w(:,1)),pearson(align1w(:,1:2))]
j = 4; targ = true2w;% align 2 CuSe
align2w=[targ,targ'*W3(j,id1)'/norm(W3(j,id1)')^2*W3(j,id1)'];
figure;plot(truetemp,align2w);
[norm(align2w(:,2)-align2w(:,1))/norm(align2w(:,1)),pearson(align2w(:,1:2))]
j = 5; targ = true3w;% align 3 Cu2Se
align3w=[targ,targ'*W3(j,id1)'/norm(W3(j,id1)')^2*W3(j,id1)'];
figure;plot(truetemp,align3w);
[norm(align3w(:,2)-align3w(:,1))/norm(align3w(:,1)),pearson(align3w(:,1:2))]
j = 1; targ = true4w;% align 4 Se
align4w=[targ,targ'*W3(j,id1)'/norm(W3(j,id1)')^2*W3(j,id1)'];
figure;plot(truetemp,align4w);
[norm(align4w(:,2)-align4w(:,1))/norm(align4w(:,1)),pearson(align4w(:,1:2))]
j = 3; targ = true5w;% align 5 pyrite
align5w=[targ,targ'*W3(j,id1)'/norm(W3(j,id1)')^2*W3(j,id1)'];
figure;plot(truetemp,align5w);
[norm(align5w(:,2)-align5w(:,1))/norm(align5w(:,1)),pearson(align5w(:,1:2))]
j = 6; targ = true6w;% align 6 marcasite
align6w=[targ,targ'*W3(j,id1)'/norm(W3(j,id1)')^2*W3(j,id1)'];
figure;plot(truetemp,align6w);
[norm(align6w(:,2)-align6w(:,1))/norm(align6w(:,1)),pearson(align6w(:,1:2))]






