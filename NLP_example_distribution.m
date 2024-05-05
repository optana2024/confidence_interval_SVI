%NLP_example_distribution 

%This script is written specifically for the NLP example.
%To run this code, first run NLP_main.m, and then run this script

%%%%Below is the limiting distribution for N^1/2(zN-x0)


Pi1=eye(4);%for points z in R^2\times R^2_+
Pi2=Pi1;%for points z in R^2\times R_+ \times R_-
Pi2(4,4)=0;
Pi3=Pi1;
Pi3(3,3)=0;%for points z in R^2\times R_- \times R_+
Pi4=Pi1; %for points z in R^2\times R^2_-
Pi4(3,3)=0;
Pi4(4,4)=0;
L=[2 -1 -1 0; -1 2 0 -1; 1  0  0  0;  0  1  0  0];
M1=L*Pi1+eye(4)-Pi1;
M2=L*Pi2+eye(4)-Pi2;
M3=L*Pi3+eye(4)-Pi3;
M4=L*Pi4+eye(4)-Pi4;
Sigma0=1/3*eye(4);
invM1=inv(M1);
invM2=inv(M2);
invM3=inv(M3);
invM4=inv(M4);
Q1=invM1*Sigma0*invM1';
Q2=invM2*Sigma0*invM2';
Q3=invM3*Sigma0*invM3';
Q4=invM4*Sigma0*invM4';

%%%Below is the limiting distribution for N^1/2(bar{z}_N-z0)


W1=[1 0; 0 1; 0 0; 0 0];
W2=[0 0; 0 0; 1 0; 0 1];

tPi_E_Q1=tPi(W1,W2,Q1);
tPi_E_Q2=tPi(W1,W2,Q2);
tPi_E_Q3=tPi(W1,W2,Q3);
tPi_E_Q4=tPi(W1,W2,Q4);

tmp1=tPi_E_Q1*Q1*tPi_E_Q1';
tmp2=tPi_E_Q2*Q2*tPi_E_Q2';
tmp3=tPi_E_Q3*Q3*tPi_E_Q3';
tmp4=tPi_E_Q4*Q4*tPi_E_Q4';


%the limiting distribution for N^1/2(bar{x}_N-x0)
PiH=[1 0 0 0; 0 1 0 0; 0 0 0 0; 0 0 0 0];
PiH*tmp1*PiH'

%%%Below is the limiting distribution for N^1/2(xN-x0)
N1=Pi1*Q1*Pi1';
N2=Pi2*Q2*Pi2';
N3=Pi3*Q3*Pi3';
N4=Pi4*Q4*Pi4';

%%%Below is to compute the sample covariance matrix of
%%%N^1/2(\bar{x}_N-x0)
tmp=bar_xN_lambdaN(1:2,:);
data=sqrt(N)*(tmp-x0);
cov(data');

%%%Below is to compute the sample covariance matrices of
%%%N^1/2(x_N-x0)
data=sqrt(N)*(xN-x0);
cov(data');

%%Below is to compute the density of the limiting distribution of N^1/2(xN-x0). 
%Note that this is just the margional density function of the first two components
%of the limiting distribution of N^{1/2}(zN-z0). 

%On K1= R^2 \times R^2_+, this limiting distribution coincides with N(0,Q1)
%On K2= R^2 \times R_+ \times R_-, this limiting distribution coincides with N(0,Q2)
%On K3= R^2 \times R_- \times R_+, this limiting distribution coincides with N(0,Q3)
%On K4= R^2 \times R^2_-, this limiting distribution coincides with N(0,Q4)

%x_points: we will evaluate the density of this distribution at
%x_points(:,i) for each i=1,2, ..., num_x_points

%num_x_points: the number of x_points 

%the density function of a multinormal distribution N(mu, Sigma) is
%exp(-1/2 (x-mu)^T inv(Simga)(x-mu))/(sqrt((2pi)^k det(Sigma))


x_points=[-0.05 0 0 0.025 0.025 0.05 0.05 0.075 0.075 0.1 ; -0.05 -0.04 0 -0.02 0.02 0 0.05 0 0.05 0.05    ]

num_x_points=size(x_points,2);
density_at_points=zeros(1,num_x_points);

for i=1:1:num_x_points
    %compute the density of this distribution at (x1,x2)
    x1=x_points(1,i);
    x2=x_points(2,i);
    
    %integration over R^2_+
    det_Q1=det(Q1);
    inv_Q1=inv(Q1);
    c7=sqrt((2*pi)^4*det_Q1);
    
    
    fun1=@(u,v) exp(-1/2* [x1 x2 u v]*inv_Q1*[x1;x2;u;v])/c7;
 %   q1=integral2(fun,0,inf,0,inf);
    
    
    %integration over K1
    %the integrand is 
    %exp(-1/2 (x1, x2, u, v)^T * inv(Q1) * (x1,x2,u,v))/sqrt((2pi)^k det(Q1))
    %as a function of (u,v), it is of the form
    %exp(c1*u^2 + c2*u v + c3*v^2 + c4*u+ c5*v + c6)/c7
    %[c1,c2,c3,c4,c5,c6,c7]=compute_coef(Q1,x1,x2);   
    %fun=@(u,v) exp(c1*u.^2+c2*u.*v+c3*v.^2+c4*u+c5*v+c6)/c7;
    %q1=integral2(fun,0,inf,0,inf);
    
     %integration over R_+ \times R_-     
    det_Q2=det(Q2);
    inv_Q2=inv(Q2);
    c7=sqrt((2*pi)^4*det_Q2);
    fun2=@(u,v) exp(-1/2* [x1 x2 u v]*inv_Q2*[x1;x2;u;v])/c7;
 %   q2=integral2(fun,0,inf,-inf,0);
    
    
     %integration over R_- \times R_+     
    det_Q3=det(Q3);
    inv_Q3=inv(Q3);
    c7=sqrt((2*pi)^4*det_Q3);
    fun3=@(u,v) exp(-1/2* [x1 x2 u v]*inv_Q3*[x1;x2;u;v])/c7;
    %q3=integral2(fun,-inf,0,0,inf);
    
    %integration over R^2_-     
    det_Q4=det(Q4);
    inv_Q4=inv(Q4);
    c7=sqrt((2*pi)^4*det_Q4);
    fun4=@(u,v) exp(-1/2* [x1 x2 u v]*inv_Q4*[x1;x2;u;v])/c7;
    %q4=integral2(fun,-inf,0,-inf,0);
    
    %density_at_points(1,i)=q1+q2+q3+q4;
    
   
end





% function [c1,c2,c3,c4,c5,c6,c7]=compute_coef(Qi,x1,x2)
%     det_Qi=det(Qi);
%     inv_Qi=inv(Qi);
%     %compute c6
%     tmp=inv_Qi(1:2,1:2);
%     c6= -[x1 x2]*tmp*[x1; x2]/2;
%     
%     %c7:sqrt((2pi)^k det_Qi)
%     c7=sqrt((2*pi)^4*det_Qi);
%    
%     %c1: the coeffient of u^2 in (x1,x2,u,v)^T * inv_Q1*(x1,x2,u,v) is
%     %inv_Qi(3,3). This needs to be multiplied with -1/2.
%     
%     c1= -inv_Qi(3,3)/2;
%     %c3: the coeffient of v^2 in (x1,x2,u,v)^T * inv_Q1*(x1,x2,u,v) is
%     %inv_Qi(4,4). This needs to be multiplied with -1/2.
%     c3=-inv_Qi(4,4)/2;
%     %c2: the coeffient of u v in (x1,x2,u,v)^T * inv_Q1*(x1,x2,u,v) is
%     %2*inv_Qi(3,4). This needs to be multiplied with -1/2.
%     c2=-inv_Qi(3,4);
%     %c4: the terms for x1*u and x2*u in (x1,x2,u,v)^T * inv_Q1*(x1,x2,u,v)
%     %are 2*inv_Qi(1,3)x1 u + 2* inv_Qi(2,3)x2 u. Mutliply this by -1/2.   
%     c4=-(inv_Qi(1,3)*x1+ inv_Qi(2,3)*x2);
%     %c5: the terms for x1*v and x2*v in (x1,x2,u,v)^T * inv_Q1*(x1,x2,u,v)
%     %are 2*inv_Qi(1,4)x1 v + 2* inv_Qi(2,4)x2 v. Mutliply this by -1/2.   
%     c5=-(inv_Qi(1,4)*x1+ inv_Qi(2,4)*x2);
% end
