%The function tPi is to construct the matrix 
%\tilde{\Pi}_E(\Sigma) using Definition 1 of the paper
%Input: W1: a matrix whose columns form a basis of E
%W2: a matrix whose columns form a basis for a complementary subspace of E
%Sigma: a symmetirc, positive definite matrix

%Output: tPi_E_Sigma=\tilde{\Pi}_E(\Sigma)

function tPi_E_Sigma=tPi(W1,W2,Sigma)
k=size(W1,2);
n=size(W1,1);
if k>0
    W=[W1, W2];
    invW=inv(W);

    tilde_Sigma=invW*Sigma* invW';
    tilde_Sigma_12=tilde_Sigma(1:k,k+1:n);
    tilde_Sigma_22=tilde_Sigma(k+1:n, k+1:n);
    inv_tSigma_22=inv(tilde_Sigma_22);

    %tilde_W_2=W1* tilde_Sigma_12*inv_tSigma_22+W2;
    %tilde_W=[W1, tilde_W_2];


    %note that the definition in equation (8) can be rewritten as
    %\tilde{\Pi}_E(\Sigma) = [W_1, -W_1 \tilde{\Sigma}_{12}*
    %\tilde{\Sigma}_{22}^{-1}] W^{-1} in view of the formula for \tilde{W}^{-1}
    %above equation (11)

    tmp=W1*tilde_Sigma_12*inv_tSigma_22;

    tPi_E_Sigma=[W1, -tmp]*invW; 
else
    tPi_E_Sigma = zeros(n,n);
end
end