
%Generate SAA and true problems for an LCP

%F(x,\xi)= A(\xi)x+b(\xi)
%A(\xi) is an n*n matrix
%Entries of A(\xi) and b(\xi) follow certain distributions

%The SAA problem:
%0 <=  A_saa x+b_saa \perp x >= 0
%The true problem:
%0 <=  A_true x+b_true \perp x >= 0


%Here, A(\xi) and b(\xi) follow uniform distributions represented by 
%rand(n,n)*scale_A + constant_A
%and 
%rand(n,1).*scale_b+constant_b
%where rand(n,n) is an (n,n) matrix whose entries are iid [0,1] uniform RVs
%and rand(n,1) is an (n,1) matrix whose entries are iid [0,1] uniform RVs

%scale_A: (n,n) matrix used to scale up random numbers in A_allsamples
%constant_A: (n,n) matrix used to shift the random numbers in A_allsamples
%scale_b: (n,1) matrix used to scale up random numbers in b_allsamples
%constant_b: (n,1) matrix used to shift random numbers in b_allsamples

%numRep: the number of replications
%N: the sample size in each replication
%n: the dimension of each LCP

%inputs:
%exm_code: a text that indicates which example 

%outputs:
%A_allsamples is a matrix of dimension (n,n,N,numRep),
%b_allsamples is a matrix of dimension (n,N,numRep)
%A_true is a matrix of dimension (n,n),
%b_true is a matrix of dimension (n,1)
%A_saa is a matrix of dimension (n,n,numRep),
%b_saa is a matrix of dimension (n,numRep)



function [A_allsamples,b_allsamples, A_true, b_true, A_saa, b_saa] = LCP_SAA_sim(exm_code)

[rng_code,n,N,numRep,scale_A,constant_A,scale_b,constant_b]=parse_code(exm_code);

rng(rng_code);

A_allsamples=zeros(n,n,N,numRep);
b_allsamples=zeros(n,N,numRep);
A_saa=zeros(n,n,numRep);
b_saa=zeros(n,numRep);

for k=1:1:numRep
    for i=1:1:N
        A_allsamples(:,:,i,k)=rand(n,n).*scale_A+constant_A;
        b_allsamples(:,i,k)=rand(n,1).*scale_b+constant_b;
    end
    A_saa(:,:,k)=mean(A_allsamples(:,:,:,k),3);
    b_saa(:,k)=mean(b_allsamples(:,:,k),2);
end 



A_true=scale_A/2+constant_A;
b_true=scale_b/2+constant_b;

%%begin of temparary testing codes for testing if A_allsamples and b_allsamples are uniform
%%distributed by scattle plots

%scatter(squeeze(A_allsamples(1,1,:,1)),squeeze(A_allsamples(1,2,:,1)));
%cov(squeeze(A_allsamples(1,1,:,1)),squeeze(A_allsamples(1,2,:,1)));

%End of temparary testing codes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This function accepts the exm_code, and then generate the right
%data that will be used to simulate A_allsamples, b_allsamples

function [rng_code,n,N,numRep,scale_A,constant_A,scale_b,constant_b]=parse_code(exm_code)
if strcmp(exm_code,'2d')
    rng_code=123;
    %rng_code=84396;
    n=2;
    N=100;
    numRep=600;
    scale_A=[2, 1; 2, 4]; constant_A=[0, 0; 0, 0];
    scale_b=[2; 1]; constant_b=[-1; 0];
elseif strcmp(exm_code,'10d_fixedA1')
    rng_code=123;
    %rng_code=84396;
    n=10;
    N=500;
    numRep=500;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    constant_b=-ones(n,1);
    for i=1:1:n
        constant_A(i,i)=2;
        for j=1:1:n
            if i<j
                constant_A(i,j)=1.5;
            elseif i>j
                constant_A(i,j)=1;
            end
        end
    end
elseif strcmp(exm_code,'10d_randomA1')
    rng_code=123;
    %rng_code=84396;
    n=10;
    N=2000;
    numRep=500;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    constant_b=-ones(n,1);
    for i=1:1:n
        scale_A(i,i)=4;
        for j=1:1:n
            if i<j
                scale_A(i,j)=3;
            elseif i>j
                scale_A(i,j)=2;
            end
        end
    end
elseif strcmp(exm_code,'10d_fixedA2')
    rng_code=123;
    %rng_code=84396;
    n=10;
    N=2000;
    numRep=500;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    scale_b(5:10)=0.8*ones(6,1);
    constant_b=-ones(n,1);
    constant_b(1)=-2;constant_b(2)=-2;
    for i=1:1:n
        constant_A(i,i)=2;
        for j=1:1:n
            if i<j
                constant_A(i,j)=1.5;
            elseif i>j
                constant_A(i,j)=1;
            end
        end
    end
elseif strcmp(exm_code,'10d_randomA2')
    rng_code=123;
    %rng_code=84396;
    n=10;
    N=2000;
    numRep=500;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    scale_b(5:10)=0.8*ones(6,1);
    constant_b=-ones(n,1);
    constant_b(1)=-2;constant_b(2)=-2;
    for i=1:1:n
        scale_A(i,i)=4;
        for j=1:1:n
            if i<j
                scale_A(i,j)=3;
            elseif i>j
                scale_A(i,j)=2;
            end
        end
    end
elseif strcmp(exm_code,'10d_2randomA_tmp')
    %rng_code=123;
    rng_code=84396;
    n=10;
    N=2000;
    numRep=2000;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    scale_b(5:10)=0.8*ones(6,1);
    constant_b=-ones(n,1);
    constant_b(1)=-2;constant_b(2)=-2;
    for i=1:1:n
        scale_A(i,i)=4;
        for j=1:1:n
            if i<j
                scale_A(i,j)=3;
            elseif i>j
                scale_A(i,j)=2;
            end
        end
    end    
elseif strcmp(exm_code,'30d_1')
    rng_code=123;
    %rng_code=84396;
    n=30;
    N=500;
    numRep=500;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    constant_b=-ones(n,1);
    for i=1:1:n
        scale_A(i,i)=4;
        for j=1:1:n
            if i<j
                scale_A(i,j)=3;
            elseif i>j
                scale_A(i,j)=2;
            end
        end
    end
elseif strcmp(exm_code,'30d_2')
    %With the following data, most solutions zN have zero in them,
    %I suspect that this is becasue the sameple size N=500 is not large
    %enough to represent the true distribution of zN. 
    %When I change N to N=5000, however, Matlab has a problem with 
    %storage. So I decide to switch to d=10.
    
    rng_code=123;
    %rng_code=84396;
    n=30;
    N=5000;
    numRep=5000;
    %specify values of scale_A, constant_A, scale_b and constant_b
    scale_A=zeros(n,n); 
    constant_A=zeros(n,n); 
    scale_b=2*ones(n,1);
    scale_b(11:30,1)=0.8*ones(20,1);
    constant_b=-ones(n,1);
    constant_b(1)=-2;constant_b(2)=-2;
    for i=1:1:n
        scale_A(i,i)=4;
        for j=1:1:n
            if i<j
                scale_A(i,j)=3;
            elseif i>j
                scale_A(i,j)=2;
            end
        end
    end

end
