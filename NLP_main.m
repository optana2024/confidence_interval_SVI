%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Confidence region for a stochastic nonlinear program
%Shu Lu, July-November 2020, University of North Carolina at Chapel Hill
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%the true problem

%min true_obj(x1,x2)=x1^3+x2^3+x1^2+x2^2-x1x2
%s.t. - x1 <=0, -x2<=0

%solution: x0=[0,0], lambda0=[0,0]

%options = optimoptions('fmincon','Algorithm','interior-point','SpecifyObjectiveGradient',true);

%g0=@(x)x(1)^3+x(2)^3+x(1)^2+x(2)^2 - x(1)*x(2);

%xstart=[0; 0];

%[x0,fval,exitflag,output] = fmincon(g0,xstart,[], [], [],[], [0;0], [inf; inf]);
%Instead of using the function fmincon to find the true solutions, I just
%input them here. 
x0=[0; 0];
lambda0=[0; 0];
z0=[0;0;0;0];

%specify dimension n, sample size N and number of replications numRep
n=4;

N=100;
numRep=50000;

%generate data for SAA problems
[xi1_allsamples,  xi2_allsamples,xi3_allsamples,xi4_allsamples,xi5_allsamples, xi1_SAA, xi2_SAA, xi3_SAA, xi4_SAA, xi5_SAA]=NLP_SAA_sim(N,numRep);

%solve all SAA problems

%min x_1^3 + x_2^3 +  x_1^2 + x_2^2 + N^{-1} \sum_{i=1}^N \xi^i_1 x_1x_2 + 
%N^{-1} \sum_{i=1}^N \xi^i_2 x_1 + N^{-1} \sum_{i=1}^N \xi^i_3 x_2 \\
%s.t. x_1 >= N^{-1} \sum_{i=1}^N \xi^i_4
% x_2 \ge N^{-1} \sum_{i=1}^N \xi^i_5

xN=zeros(2,numRep);
lambdaN=zeros(2,numRep);
zN=zeros(4,numRep);
exitflag=zeros(1,numRep);

for k=1:1:numRep
    gN=@(x)x(1)^3+x(2)^3+x(1)^2+x(2)^2+xi1_SAA(k)*x(1)*x(2)+xi2_SAA(k)*x(1)+xi3_SAA(k)*x(2);
    xstart=[0;0];
    lb=[xi4_SAA(k); xi5_SAA(k)];
    ub=[inf;inf];
    options = optimoptions('fmincon','OptimalityTolerance',1e-8);
    [xN_tmp,fval,exitflag_tmp,output] = fmincon(gN,xstart,[], [], [],[],lb,ub,[],options);
    exitflag(k)=exitflag_tmp;
    lambdaN_tmp=zeros(2,1);
    lambdaN_tmp(1)=3*xN_tmp(1)^2 + 2*xN_tmp(1) + xi1_SAA(k)*xN_tmp(2)+xi2_SAA(k);
    lambdaN_tmp(2)=3*xN_tmp(2)^2 + 2*xN_tmp(2) + xi1_SAA(k)*xN_tmp(1)+xi3_SAA(k);
    %check if KKT conditions hold at (xN_tmp, lambdaN_tmp)
    yN_tmp=zeros(2,1);
    yN_tmp(1)=xN_tmp(1)-xi4_SAA(k);
    yN_tmp(2)=xN_tmp(2)-xi5_SAA(k);
    
    zN_tmp=zeros(4,1);
    zN_tmp(1)=xN_tmp(1);
    zN_tmp(2)=xN_tmp(2);
    zN_tmp(3)=lambdaN_tmp(1) - yN_tmp(1);
    zN_tmp(4)=lambdaN_tmp(2) - yN_tmp(2);
    
    xN(:,k)=xN_tmp;
    lambdaN(:,k)=lambdaN_tmp;
    
    zN(:,k)=zN_tmp;

end

%For each replication, compute individual confidence intervals for z0 and x0 of confidence leval
%1-alpha, using (1-alpha1)-confidence regions (simultaneous confidence intervals)
%for z0 to identify the piecewise structure

%For the purpose of checking the validity of this method, also compute
%individual confidence intervals for z0 and x0 based on true spaces E and H

%We choose not to compare the intervals from comp_CI with the naive method in comp_CI_naive, because this example
%has 4 pieces, so the naive method does not have theoretical asymptotic exactness guarantee 

alpha=0.1;
alpha1=0.05;
alpha2=alpha-alpha1;

simCI_z0=zeros(n,numRep);

%We will use the function comp_CI to compute individual confidence intervals.
%This fucntion deals with a general set S = {x \in R^n | x_i
%>=0, i \in I}
%For the example here,  S=R^2 \times R^2_+, so I={3, 4}.
%We use the array I=[3,4] to denote the index set {3,4}
%p is the cardinality of the index set I    

p=2;
I=[3,4];

C_iN=int8(zeros(p,numRep));
indCI_z0=zeros(n,numRep);
indCI_x0_lambda0=zeros(n,numRep);
tilde_zN=zeros(n,numRep);
tilde_xN_lambdaN=zeros(n,numRep);
indCI_z0_check=zeros(n,numRep);
indCI_x0_lambda0_check=zeros(n,numRep);
bar_zN=zeros(n,numRep);
bar_xN_lambdaN=zeros(n,numRep);

%z0_inConfGeg(k)=1 if z0 is included in the confidence region of the kth
%replication
z0_inConfReg=logical(zeros(numRep,1));

%z0_inIndCI(j,k)=1 if z0(j) is included in the individual confidence intervals for z0(j) 
%of the kth replication
%z0_inIndCI_check(j,k)=1 if z0(j) is included in the individual confidence
%intervals for z0(j) using the true space E of the kth replication

%similarly for x0_lambda0_inIndCI and x0_lambda0_inIndCI_check
%tol: a parameter used in checking if components of z0 or x0 are included in
%their corresponding confidence intervals by taking account for numerical inacurracies 

z0_inIndCI=logical(zeros(n,numRep));
x0_lambda0_inIndCI=logical(zeros(n,numRep));
z0_inIndCI_check=logical(zeros(n,numRep));
x0_lambda0_inIndCI_check=logical(zeros(n,numRep));
tol=1e-6;


for k=1:1:numRep
    %Compute the covariance matrix Sigma_N, the sample covariance matrix of 
    %F(xN, lambdaN, \xi) for all N samples \xi
    
    %F(x,lambda,\xi)=
    %[3x1^2+2x1+xi1 x2 + xi2 - lambda1, 
    % 3x2^2+2x2+xi1 x1 + xi3 - lambda2,
    % x1 - xi4,
    % x2 - xi5]
    
    %To simplify computation, we don't need to use terms in
    %F(xN,lambdaN,\xi) that do not depend on xi
    %Hence, the terms we use are
    %[xi1 x2 + xi2 , 
    % xi1 x1 + xi3 ,
    % - xi4,
    % - xi5]
    
    Fxn_reduced=zeros(n,N);
    for i=1:1:N
        %Fxn(1,i)=3*xN(1,k)^2 + 2*xN(1,k) + xi1_allsamples(i,k)*xN(2,k)+xi2_allsamples(i,k)-lambdaN(1,k);
        Fxn_reduced(1,i)= xi1_allsamples(i,k)*xN(2,k)+xi2_allsamples(i,k);
        Fxn_reduced(2,i)= xi1_allsamples(i,k)*xN(1,k)+xi3_allsamples(i,k);
        Fxn_reduced(3,i)= -xi4_allsamples(i,k);
        Fxn_reduced(4,i)= -xi5_allsamples(i,k);
    end
    Sigma_N = cov(Fxn_reduced');
    
    d_fN_xN_lambdaN=zeros(4,4);
    d_fN_xN_lambdaN(1,1)=6*xN(1,k)+2;
    d_fN_xN_lambdaN(1,2)=xi1_SAA(k);
    d_fN_xN_lambdaN(1,3)=-1;
    d_fN_xN_lambdaN(2,1)=xi1_SAA(k);
    d_fN_xN_lambdaN(2,2)=6*xN(2,k)+2;
    d_fN_xN_lambdaN(2,4)=-1;
    d_fN_xN_lambdaN(3,1)=1;
    d_fN_xN_lambdaN(4,2)=1;
    
    
    %call the function comp_CI to compute confidence regions, confidence
    %intervals for this replication
    
    %The function comp_CI deals with a general set S = {x \in R^n | x_i
    %>=0, i \in I}
    %Here, S=R^2 \times R^2_+, so I={3, 4}.
    %We use the array I=[3,4] to denote the index set {3,4}
    I=[3,4];
    
    [inv_Lambda_N, r_conf_reg, Lambda_N, simCI_z0(:,k),C_iN(:,k),indCI_z0(:,k),indCI_x0_lambda0(:,k),tilde_zN(:,k),tilde_xN_lambdaN(:,k)]=comp_CI(n,N,zN(:,k),d_fN_xN_lambdaN,Sigma_N,alpha,alpha1,I);
    tmp=(z0-zN(:,k))' * inv_Lambda_N *(z0-zN(:,k));
    if tmp <= r_conf_reg
        z0_inConfReg(k)=1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %begin of temporary codes for graphing specific replications
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %end of temporary codes for graphing specific replications
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %call the function comp_CI_check to compute individual confidence
    %intervals using the true spaces E and H for this replication
    
    
    [indCI_z0_check(:,k),indCI_x0_lambda0_check(:,k),bar_zN(:,k),bar_xN_lambdaN(:,k)]=comp_CI_check(n,N,zN(:,k),Lambda_N, alpha2,z0,I); %n,N,zN,Lambda_N, alpha2,z0
    
    
    for j=1:1:n
        z0_inIndCI(j,k)=isIncluded(z0(j),tilde_zN(j,k),indCI_z0(j,k),tol);
        z0_inIndCI_check(j,k)=isIncluded(z0(j),bar_zN(j,k),indCI_z0_check(j,k),tol);
        if j<=2 %when j=1,2, the jth elements in tilde_xN_lambdaN(:,k) and bar_xN_lambdaN are estimators for x0(j)
            x0_lambda0_inIndCI(j,k)=isIncluded(x0(j),tilde_xN_lambdaN(j,k),indCI_x0_lambda0(j,k),tol);
            x0_lambda0_inIndCI_check(j,k)=isIncluded(x0(j),bar_xN_lambdaN(j,k),indCI_x0_lambda0_check(j,k),tol);
        else %when j=3,4, the jth elements in tilde_xN_lambdaN(:,k) is an estimator for lambda0(j-2)
            x0_lambda0_inIndCI(j,k)=isIncluded(lambda0(j-2),tilde_xN_lambdaN(j,k),indCI_x0_lambda0(j,k),tol);
            x0_lambda0_inIndCI_check(j,k)=isIncluded(lambda0(j-2),bar_xN_lambdaN(j,k),indCI_x0_lambda0_check(j,k),tol);
        end
    end
end

%rate of coverages

%rate_cov_region: the ratio of replications that include z0 in the
%confidence region

%rate_cov_z0_ind(j): the ratio of replications that include z0(j) in the
%individual confidence intervals for z0(j)

%rate_cov_x0_ind(j): the ratio of replications that include x0(j) in the
%individual confidence intervals for x0(j)

%rate_cov_lambda0_ind(j): the ratio of replications that include lambda0(j) in the
%individual confidence intervals for lambda0(j)


%rate_cov_z0_ind_check(j): the ratio of replications that include z0(j) in the
%individual confidence intervals for z0(j) obtained using the true space E

%rate_cov_x0_ind_check(j): the ratio of replications that include x0(j) in the
%individual confidence intervals for x0(j) obtained using the trus spaces E
%and H

%rate_cov_lambda0_ind_check(j): the ratio of replications that include lambda0(j) in the
%individual confidence intervals for lambda0(j) obtained using the trus spaces E
%and H


rate_cov_region=sum(z0_inConfReg)/numRep;

rate_cov_z0_ind=sum(z0_inIndCI,2)/numRep;
rate_cov_x0_ind=sum(x0_lambda0_inIndCI(1:2,:),2)/numRep;
rate_cov_lambda0_ind=sum(x0_lambda0_inIndCI(3:4,:),2)/numRep;

rate_cov_z0_ind_check=sum(z0_inIndCI_check,2)/numRep;
rate_cov_x0_ind_check=sum(x0_lambda0_inIndCI_check(1:2,:),2)/numRep;
rate_cov_lambda0_ind_check=sum(x0_lambda0_inIndCI_check(3:4,:),2)/numRep;


%average lengths of indivivual confidence intervals

%ave_indCI_x0(j): the average length of indivividual confidence intervals
%for x0(j) across all replications

%ave_indCI_lambda0(j): the average length of indivividual confidence intervals
%for lambda0(j) across all replications


%ave_indCI_z0(j): the average length of indivividual confidence intervals
%for z0(j) across all replications

%ave_indCI_x0_check(j): the average length of indivividual confidence intervals
%for x0(j) across all replications obtained using the true space E

%ave_indCI_lambda0_check(j): the average length of indivividual confidence intervals
%for lambda0(j) across all replications obtained using the true space E


%ave_indCI_z0_check(j): the average length of indivividual confidence intervals
%for z0(j) across all replications obtained using the true space E

ave_indCI_x0=mean(indCI_x0_lambda0(1:2,:),2);
ave_indCI_lambda0=mean(indCI_x0_lambda0(3:4,:),2);
ave_indCI_z0=mean(indCI_z0,2);

ave_indCI_x0_check=mean(indCI_x0_lambda0_check(1:2,:),2);
ave_indCI_lambda0_check=mean(indCI_x0_lambda0_check(3:4,:),2);
ave_indCI_z0_check=mean(indCI_z0_check,2);


function [xi1_allsamples,  xi2_allsamples,xi3_allsamples,xi4_allsamples,xi5_allsamples, xi1_SAA, xi2_SAA, xi3_SAA, xi4_SAA, xi5_SAA]=NLP_SAA_sim(N,numRep); 

rng(123);

%\xi_1: uniform distribution on [-2, 0]
xi1_allsamples=rand(N,numRep)*2-2;

%\xi_2, \xi_3, \xi_4, \xi_5: uniform distribution on [-1, 1]
xi2_allsamples=rand(N,numRep)*2-1;
xi3_allsamples=rand(N,numRep)*2-1;
xi4_allsamples=rand(N,numRep)*2-1;
xi5_allsamples=rand(N,numRep)*2-1;

xi1_SAA=mean(xi1_allsamples,1);
xi2_SAA=mean(xi2_allsamples,1);
xi3_SAA=mean(xi3_allsamples,1);
xi4_SAA=mean(xi4_allsamples,1);
xi5_SAA=mean(xi5_allsamples,1);
end
