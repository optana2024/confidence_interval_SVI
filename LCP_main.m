%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Confidence region for linear complementarity problems
%Shu Lu, July-November 2020, University of North Carolina at Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%F(x,\xi)= A(\xi)x+b(\xi)
%A(\xi) is an n*n matrix
%Entries of A(\xi) and b(\xi) follow certain distributions

%The SAA problem:
%0 <=  A_saa x+b_saa \perp x >= 0
%The true problem:
%0 <=  A_true x+b_true \perp x >= 0

%Generate SAA and true problems

%A_allsamples is a matrix of dimension (n,n,N,numRep),
%b_allsamples is a matrix of dimension (n,N,numRep)
%A_true is a matrix of dimension (n,n),
%b_true is a matrix of dimension (n,1)
%A_saa is a matrix of dimension (n,n,numRep),
%b_saa is a matrix of dimension (n,numRep)

%LCP: n=2
%[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('2d');

%LCP: n=10, example fixedA1
[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('10d_fixedA1');


%LCP: n=10, example randomA1
%[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('10d_randomA1');

%LCP: n=10, example fixedA2
%[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('10d_fixedA2');

%LCP: n=10, example randomA2
%[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('10d_randomA2');


%LCP: n=30, example 1
%[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('30d_1');

%LCP: n=30, example 2
%[A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim('30d_2');


%numRep: the number of replications
%N: the sample size in each replication
%n: the dimension of each LCP

[n,N,numRep]=size(b_allsamples);

%Compute the true solution.

%The LCP: 0 \in A_true*x+b_true + N_{R^2_+}(x)
%i.e., 0 \le A_true*x+b_true \perp x \ge 0

[y0,x0,retcode] = LCPSolve(A_true,b_true,1e-8,1e4);
if retcode(1)==2 || retcode(2)==1e4
    error('The true LCP is not solved successfully'); 
    return
end
z0=x0-y0; 

%Compute the SAA solution for each replication.

%The LCP: 0 \in A_saa*x+b_saa + N_{R^2_+}(x)
%i.e., 0 \le A_saa*x+b_saa \perp x \ge 0

xN=zeros(n,numRep);
yN=zeros(n,numRep);
zN=zeros(n,numRep);

for k=1:1:numRep
    [yN_tmp,xN_tmp,retcode] = LCPSolve(A_saa(:,:,k),b_saa(:,k),1e-8,1e4);
    if retcode(1)==2 || retcode(2)==1e4
        error('The %dth SAA problem is not solved successfully',k); 
        return
    end
    xN(:,k)=xN_tmp;
    yN(:,k)=yN_tmp;
    zN(:,k)=xN(:,k) - yN(:,k); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%For each replication, compute individual confidence intervals for z0 and x0 of confidence leval
%1-alpha, using (1-alpha1)-confidence regions (simultaneous confidence intervals)
%for z0 to identify the piecewise structure.
%The main method of computing individual confidence intervals for true solutions 
%is in the function comp_CI.m, which is based on the paper 
%"Statistical Inference for Piecewise Normal Distributions and
%Stochastic Variational Inequalities" by Shu Lu.


%For the purpose of checking the validity of this method, also compute
%individual confidence intervals for z0 and x0 based on true spaces E and H

%For the purpose of comparison, also use the naive method (that ignores the
%piecewise structure) to compute individual confidence intervals for z0 and
%x0. This method is known to be asymptotic exact when the true asymptotic distribution has
%at most two pieces, so it works for the example in which n=2. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


alpha=0.1;
alpha1=0.05;
alpha2=alpha-alpha1;

simCI_z0=zeros(n,numRep);
C_iN=int8(zeros(n,numRep));
indCI_z0=zeros(n,numRep);
indCI_x0=zeros(n,numRep);
tilde_zN=zeros(n,numRep);
tilde_xN=zeros(n,numRep);
indCI_z0_check=zeros(n,numRep);
indCI_x0_check=zeros(n,numRep);
bar_zN=zeros(n,numRep);
bar_xN=zeros(n,numRep);
indCI_z0_naive=zeros(n,numRep);
indCI_x0_naive=zeros(n,numRep);


%z0_inConfGeg(k)=1 if z0 is included in the confidence region of the kth
%replication
z0_inConfReg=logical(zeros(numRep,1));

%z0_inIndCI(j,k)=1 if z0(j) is included in the individual confidence intervals for z0(j) 
%of the kth replication
%z0_inIndCI_check(j,k)=1 if z0(j) is included in the individual confidence
%intervals for z0(j) using the true space E of the kth replication
%z0_inIndCI_naive(j,k)=1 if z0(j) is included in the individual confidence intervals for z0(j) 
%of the kth replication computed using the naive method


%similarly for x0_inIndCI, x0_inIndCI_check, x0_inIndCI_naive
%tol: a parameter used in checking if components of z0 or x0 are included in
%their corresponding confidence intervals by taking account for numerical inacurracies 

z0_inIndCI=logical(zeros(n,numRep));
x0_inIndCI=logical(zeros(n,numRep));
z0_inIndCI_check=logical(zeros(n,numRep));
x0_inIndCI_check=logical(zeros(n,numRep));
z0_inIndCI_naive=logical(zeros(n,numRep));
x0_inIndCI_naive=logical(zeros(n,numRep));

tol=1e-6;


for k=1:1:numRep
    %Compute the covariance matrix Sigma_N, the sample covariance matrix of 
    %F(xN,\xi) for all N samples \xi
    Fxn=zeros(n,N);
    for i=1:1:N
        Fxn(:,i)=A_allsamples(:,:,i,k)*xN(:,k) +b_allsamples(:,i,k);
    end
    Sigma_N = cov(Fxn');
    %All LCP examples considered here have S=R^n_+, 
    %The function comp_CI considers the case in which S = {x\in R^n | x_i >=0, i\in I}
    %Hence, we let I={1,2,...n} which is done by letting I=[1,2,...n]
    I=[1:1:n];
    %call the function comp_CI to compute confidence regions, confidence
    %intervals for this replication, using the method in the paper
    %"Statistical Inference for Piecewise Normal Distributions and
    %Stochastic Variational Inequalities" by Shu Lu.
    
    [inv_Lambda_N, r_conf_reg, Lambda_N, simCI_z0(:,k),C_iN(:,k),indCI_z0(:,k),indCI_x0(:,k),tilde_zN(:,k),tilde_xN(:,k)]=comp_CI(n,N,zN(:,k),A_saa(:,:,k),Sigma_N,alpha,alpha1,I);
    tmp=(z0-zN(:,k))' * inv_Lambda_N *(z0-zN(:,k));
    if tmp <= r_conf_reg
        z0_inConfReg(k)=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %begin of temporary codes for graphing specific replications
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if n==2 && k==2
        graph2d_oneRep(zN(:,k),inv_Lambda_N,r_conf_reg,simCI_z0(:,k),tilde_zN(:,k),indCI_z0(:,k),z0);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %end of temporary codes for graphing specific replications
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %call the function comp_CI_check to compute individual confidence
    %intervals using the true spaces E and H for this replication
    
    
    [indCI_z0_check(:,k),indCI_x0_check(:,k),bar_zN(:,k),bar_xN(:,k)]=comp_CI_check(n,N,zN(:,k),Lambda_N, alpha2,z0,I); %n,N,zN,Lambda_N, alpha2,z0
    
    
    %call the function comp_CI_naive to compute individual confidence
    %intervals using the naive method
    
    [indCI_z0_naive(:,k),indCI_x0_naive(:,k)]=comp_CI_naive(n,N,zN(:,k),A_saa(:,:,k),Sigma_N,alpha,I)
    
    
    for j=1:1:n
        z0_inIndCI(j,k)=isIncluded(z0(j),tilde_zN(j,k),indCI_z0(j,k),tol);
        x0_inIndCI(j,k)=isIncluded(x0(j),tilde_xN(j,k),indCI_x0(j,k),tol);
        z0_inIndCI_check(j,k)=isIncluded(z0(j),bar_zN(j,k),indCI_z0_check(j,k),tol);
        x0_inIndCI_check(j,k)=isIncluded(x0(j),bar_xN(j,k),indCI_x0_check(j,k),tol);
        z0_inIndCI_naive(j,k)=isIncluded(z0(j),zN(j,k),indCI_z0_naive(j,k),tol);
        x0_inIndCI_naive(j,k)=isIncluded(x0(j),xN(j,k),indCI_x0_naive(j,k),tol);

    end
    
    
        
end

%rate of coverages

%rate_cov_region: the ratio of replications that include z0 in the
%confidence region

%rate_cov_z0_ind(j): the ratio of replications that include z0(j) in the
%individual confidence intervals for z0(j)

%rate_cov_x0_ind(j): the ratio of replications that include x0(j) in the
%individual confidence intervals for x0(j)

%rate_cov_z0_ind_check(j): the ratio of replications that include z0(j) in the
%individual confidence intervals for z0(j) obtained using the true space E

%rate_cov_x0_ind_check(j): the ratio of replications that include x0(j) in the
%individual confidence intervals for x0(j) obtained using the trus spaces E
%and H

%rate_cov_z0_ind_naive(j): the ratio of replications that include z0(j) in the
%individual confidence intervals for z0(j) obtained using the naive method

%rate_cov_x0_ind_naive(j): the ratio of replications that include x0(j) in the
%individual confidence intervals for x0(j) obtained using the naive method


rate_cov_region=sum(z0_inConfReg)/numRep;

rate_cov_z0_ind=sum(z0_inIndCI,2)/numRep;
rate_cov_x0_ind=sum(x0_inIndCI,2)/numRep;

rate_cov_z0_ind_check=sum(z0_inIndCI_check,2)/numRep;
rate_cov_x0_ind_check=sum(x0_inIndCI_check,2)/numRep;

rate_cov_z0_ind_naive=sum(z0_inIndCI_naive,2)/numRep;
rate_cov_x0_ind_naive=sum(x0_inIndCI_naive,2)/numRep;


%average lengths of indivivual confidence intervals

%ave_indCI_x0(j): the average length of indivividual confidence intervals
%for x0(j) across all replications

%ave_indCI_z0(j): the average length of indivividual confidence intervals
%for z0(j) across all replications

%ave_indCI_x0_check(j): the average length of indivividual confidence intervals
%for x0(j) across all replications obtained using the true space E

%ave_indCI_z0_check(j): the average length of indivividual confidence intervals
%for z0(j) across all replications obtained using the true space E

%ave_indCI_x0_naive(j): the average length of indivividual confidence intervals
%for x0(j) across all replications obtained using the naive method

%ave_indCI_z0_naive(j): the average length of indivividual confidence intervals
%for z0(j) across all replications obtained using the naive method


ave_indCI_x0=mean(indCI_x0,2);
ave_indCI_z0=mean(indCI_z0,2);
ave_indCI_x0_check=mean(indCI_x0_check,2);
ave_indCI_z0_check=mean(indCI_z0_check,2);
ave_indCI_x0_naive=mean(indCI_x0_naive,2);
ave_indCI_z0_naive=mean(indCI_z0_naive,2);






















