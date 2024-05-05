%The function comp_CI is the main function to compute the confidence region, simutaneous
%confidence intervals for z0, and individual confidence intervals for z0
%and x0. 

%The method of computing the confidence region is based on the paper 
%"Symmetric confidence regions and confidence intervals for normal map
%formualtions of stochastic variational inequalities, S Lu, SIOPT vol.
%24(3), pp 1458-1484."

%The method of computing individual confidence intervals for z0 and x0 is
%based on the method in the paper "Statistical Inference for Piecewise Normal Distributions and
%Stochastic Variational Inequalities, S Lu."

%Input: 
%n: the dimension of the problem
%N: sample size
%zN: the SAA solution for the normal map formulation
%d_fN_xN: the derivative of fN at xN. For an LCP, it is equal to A_saa
%Sigma_N: the sample covariance matrix for F(xN,\xi)
%alpha: 1-alpha is the confidence level for the individual confidence intervals
%alpha1: 1-alpha1 is the confidence level for the confidence region and
%simultaneous confidence intervals for z0
%I: we use an array I to denote an index set that is a subset of {1,2,...n} 
%that defines the set S=\{x\in R^n | x_i >=0, i\in I \}
%For example, I=[3,4] means that I={3,4} is the index set

%Output:
%inv_Lambda_N: the matrix defining the confidence region of z0, the inverse
%of Lambda_N
%r_conf_reg: the radius for the confidence region of z0
%The confidence region for z0 is 
%{z| (z-zN)^T * inv_Lambda_N*(z-zN)) <= r}
%inv_Lambda_N is the inverse of the matrix Lambda_N
%Lambda_N: the matrix \Lambda_N in the paper
%simCI_z0: (n,1) vector consisting of half-lengths of simultenous confidence
%intervals for z0 at confidence level 1-alpha1. In other words, the box
%[zN - simCI_z0,zN + simCI_z0] is the minimal bounding box of the confdience
%region of z0
%C_iN: a tuple with values 0, 1 and -1 that indicates the cell $C_{iN}$
%which is the estimator of $C_{i0}$, the cell containing z0 in its relative
%interior
%indCI_z0:   half lengths of  individual confidence intervals for z0 at 
%confidence level 1-alpha. In other words, the confidence interval for (z0)_j
%is [(tilde_zN)_j - indCI_zN_j,(tilde_zN)_j + indCI_zN_j].
%indCI_x0: half lengths of  individual confidence intervals for x0
%tilde_zN: center of individual confidence intervals for z0
%tilde_xN: center of individual confidence intervals for x0

function [inv_Lambda_N, r_conf_reg, Lambda_N, simCI_z0,C_iN, indCI_z0,indCI_x0,tilde_zN,tilde_xN]=comp_CI(n,N,zN,d_fN_xN,Sigma_N,alpha,alpha1,I)

%p: the number of indices in the index set I
p=size(I,2);


%Compute M_N = d f^{nor}_{N,S}(z_N)= dfN(xN) * d \Pi_S(z_N) + I - d \Pi_S(z_N)
%This is the derivative of f^{nor}_{N,S} at zN

%first, compute d_PiS_zN: the (n,n) matrix representing d \Pi_S(z_N)
% S=\{x\in R^n | x_i >=0, i\in I \}

d_PiS_zN = eye(n);
for i=1:1:p
    j=I(i);
    if zN(j)<0; d_PiS_zN(j,j)=0;end;
end


%Then, compute M_N 
M_N=d_fN_xN * d_PiS_zN + eye(n)-d_PiS_zN;

%Compute Lambda_N = MN^{-1} * Sigma_N * MN^{-T}. The matrix that represents
%the covariance matrix of zN locally
invMN=inv(M_N);
Lambda_N = invMN*Sigma_N*invMN';


%Identify C_{iN} and its parallel space as in Lemma 9 of the paper by using
%simultaneous confidence intervals of z0

%The confidence region QN at level 1-alpha1 is defined as
%{z| N* (M_N(z-zN))^T * inv(Sigma_N)*(M_N(z-zN)) <= chi_n_alpha1}
%We use the formula below Theorem 7. 

chi_n_alpha1=chi2inv(1-alpha1,n);%chi_n_alpha1 is a number such that the probability 
%for a chi2 r.V. with n degrees of freedom to be <= this number equals 1-alpha1

%omega: (n,1) vector consisting of half-lengths of simultenous confidence
%intervals for z0 at confidence level 1-alpha1. In other words, the box
%[zN(j)-omega(j),zN(j)+omega(j)]_j is the minimal bounding box of QN


omega = zeros(n,1);
for j=1:1:n
    omega(j)=sqrt(chi_n_alpha1/N* invMN(j,:) *Sigma_N* invMN(j,:)');
end

%inv_Lambda_N=M_N^T * inv(Sigma_N) * M_N: the inverse of Lambda_N


inv_Lambda_N=M_N'* inv(Sigma_N)*M_N;

r_conf_reg=chi_n_alpha1/N;



%Identify C_iN: C_{iN} in Lemma 9 using the minimal bounding box instead of QN


%With S=\{x\in R^n | x_i >=0, i\in I \}, each cell in its normal manifold is represented by a triple
%(I0, I+, I-), with C(I0,I+,I-)={x\in R^n| x=0, i\in I0; x>=0, i\in I+;
%x<=0, i\in I-}, and I0 \cup I+ \cup I- = I
%Here, C_iN is a p-tuple with values 0, 1 and -1; indices of 0 components in C_iN
%correspond to indices in I0, and so on.
%It is easy to see that the C_iN from codes below meets the box
%[zN(j)-omega(j),zN(j)+omega(j)]_j
%It is also of the least dimension of all such cells because
%it has as many zeros as possible

C_iN=sign(zN(I));
for i=1:1:p
    j=I(i);
    if C_iN(i)<0 &&  zN(j)+omega(j)>=0;  C_iN(i)=0; end;
    if C_iN(i)>0 &&  zN(j)-omega(j)<=0;  C_iN(i)=0; end;
end 

%Because C_{iN} always contains 0, par(C_{iN}) is the same as aff(C_{iN})
%W1: a matrix whose columns form a basis of tilde_E=par(C_{iN})
%Here, it is not necessary to try to represent tilde_E itself in the codes
%It is enough to represent it using W1
%W2: a matrix whose columns form a basis for a complementary space of
%tilde_E
tmp=logical(ones(n,1));

for i=1:1:p
    if C_iN(i)==0; tmp(I(i))=0; end;
end
    
I_n=eye(n);
W1=I_n(:,tmp);
tmp=~tmp;
W2=I_n(:,tmp);

%Construct individual confidence intervals for z0 using Theorem 8

%For the set S=\{x\in R^n | x_i >=0, i\in I \}, all cells in the normal manifold of S consists of 0, so it is most
%convenient to let a0=0.
a0=0;

%Compute tPi_tE_LambdaN= \tilde{\Pi}_{\tilde_E}(Lambda_N)

tPi_tE_LambdaN= tPi(W1,W2,Lambda_N);

%tilde_zN: \tilde{z}_N defined in Theorem 8
tilde_zN=tPi_tE_LambdaN * zN; 

%delta_alpha2_z: \delta_\alpha2(Lambda_N, \tilde{\Pi}_{\tilde{E}})
%tmp: (n,n) matrix used in Definition 2
tmp=tPi_tE_LambdaN * Lambda_N * tPi_tE_LambdaN';
delta_alpha2_z=zeros(n,1);

alpha2=alpha-alpha1;

chi_1_alpha2=chi2inv(1-alpha2,1);
for j=1:1:n
    delta_alpha2_z(j)=sqrt(chi_1_alpha2*tmp(j,j));
end

%Construct individual confidence intervals for x0 using Theorem 9

%tilde_F is used to represent the set \tilde{F}=C_{iN} \cap S


%S = \{x\in R^n | x_i >=0, i\in I \}
%C_{iN}={x\in R^n| x=0, i\in I0; x>=0, i\in I+; x<=0, i\in I-}, then
%So the set \tilde{F}=C_{iN}\cap S ={x\in R^n| x=0, i\in I0 \cup I-; x>=0, i\in I+}

%tilde_F is a p-tuple with values 0 and 1; indices of 0 components
%correspond to indices in I0 and I-, and indices of 1 components 
%correspond to indices in I+


tilde_F=C_iN;
tmp=C_iN <0;
tilde_F(tmp)=0;

%When S = \{x\in R^n | x_i >=0, i\in I \}, \tilde{F} always contains the origin, so 
%\tilde{H} = par(\tilde{F}) is the same as aff(\tilde{F}).
%\tilde{H}=par(\tilde{F})={x\in R^n| x=0, i\in I0\cup I-}
%So the Euclidean projector onto \tilde{H} is just to reduce (I0 \cup I-)-elements of
%a vector to 0

%tilde_xN=Pi_{aff(\tilde{F})}(tilde_zN)
%Pi_tH: an (n,n) matrix representing the Euclidean projector Pi_{\tilde{H}

tilde_xN=tilde_zN;
Pi_tH=eye(n);

for i=1:1:p
    if tilde_F(i)==0
        j=I(i);
        tilde_xN(j)=0;
        Pi_tH(j,j)=0;
    end
end


%delta_alpha2_x: \delta_\alpha2(Lambda_N, \Pi_{\tilde{H}}\circ \tilde{\Pi}_{\tilde{E}})
%tmp: Pi_{\tilde{H}} *  \tilde{\Pi}_{\tilde_E}(Lambda_N)
%tmp2: (n,n) matrix as used in Definition 3
tmp=Pi_tH * tPi_tE_LambdaN;
tmp2=tmp * Lambda_N * tmp';
delta_alpha2_x=zeros(n,1);

for j=1:1:n
    delta_alpha2_x(j)=sqrt(chi_1_alpha2*tmp2(j,j));
end

%Summary of all simultaneous and individual confidence intervals

simCI_z0=omega;

indCI_z0=delta_alpha2_z/sqrt(N);

indCI_x0=delta_alpha2_x/sqrt(N);

end






