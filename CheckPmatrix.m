%Check if a matrix is a P matrix by enumerating all principal minors
%Copyright (c) 2020  Shu Lu. All rights reserved.
%Written on June 21, 2020

%If the matirx is a P matrix, the output IsP is 1. Otherwise it is 2. 

function [IsP] = CheckPmatrix(M)

%n: the number of rows in M
n=size(M,1);

%checking if M is a square matrix
if size(M)~=[n n]; error('Matrices are not compatible'); end;

%Num_matrices: the number of sub_matrices to enumerate
Num_matrices= 2^n-1;

IsP=1;

for count=1:1:Num_matrices
    %convert count into a binary string
    binStr = dec2bin(count,n);
    %convert binStr into a logical array
    binArray = binStr == '1'; 
    %take a square, principal, submatrix of M by using indices with element 1 in binArray
    sub_M=M(binArray,binArray);
    if det(sub_M)<0
        IsP=2;
        break;
    end
end


