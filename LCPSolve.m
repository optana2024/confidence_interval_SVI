% % LCPSolve(M,q) solves the linear complementarity problem:
% 
%       w = M*z + q,   w and z >= 0,   w'*z = 0
% 
%   The function takes the matrix M and the vector q as arguments.  The
%   function has three return variables. The first the vectors w and the
%   second is the vector z, found by complementary pivoting.  The third 
%   return is a 1 by 2 vector.  The first component is a 1 if the algorithm
%   was successful, and a 2 if a ray termination resulted.  The second 
%   component is the number of iterations performed in the outer loop.
%
% % Example problem
% >> M=[2,-1;-1,1]; q=[-3;1];
% 
% >> [w,z,retcode] = LCPSolve(M,q);
%
% % Printed solution
% w =
%      0
%      0
%
% z =
%     2.0000
%     1.0000
% 
% retcode =
%      1     2
%
% % Copyright (c) 2013 Andreas Almqvist, Andrew Spencer and Peter Wall
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are 
% met: 
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer. 
% 2. Redistributions in binary form must reproduce the above copyright 
%    notice, this list of conditions and the following disclaimer in the
%    documentation and/or other materials provided with the distribution. 
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS 
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR 
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR 
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% We acknowledge the work LCPSolve.py in the OpenOpt python package by
% Rob Dittmar, Enzo Michelangeli and IT Vision Ltd

%comments by Shu Lu on June 18, 2020
%I don't understand why it is necessary to define dimen=size(M,1)
%Isn't dimen always equal to n since M is a square matrix?

%comments by Shu Lu on June 19, 2020
%Also, it seems to me that this method is an implementation of Lemke's
%method as outlined in [Ferris et al. Linear Programming with Matlab,
%Chapter 7]. I did not check carefully whether this method is exactly
%Lemke's method; the method in [Ferris et al.] is explained using Jordan
%exchanges, while the method here is implemented by elementary row operations 
%in the "full" tableau. (Of course, these are just differences on the
%surface.) It is worthy noting that the identification of a termination ray 
%(i.e retcode(1)=2) does not completely rule out the possibility that the 
%LCP has a solution, because the theory in p. 182 of the book is only for
%special matrices. 
%I am a little afraid that this code would not produce a solution when the
%matrix M is a P matrix. Should be OK. See Adler_Lemke_method_solvability,
%page 3. 

function [w,z,retcode] = LCPSolve(M,q,pivtol,maxits)

if nargin<3, pivtol = 1e-8; maxits = 1e4; end;
if nargin<4, maxits = 1e3; end;
n = length(q);
if size(M)~=[n n]; error('Matrices are not compatible'); end;

rayTerm = false;
loopcount = 0;
if min(q)>=0 % If all elements are positive a trivial solution exists
    %  As w - Mz = q, if q >= 0 then w = q a, and z = 0
    w = q;
    z = zeros(size(q));
else 
    dimen = size(M,1); % Number of rows
    % Create initial tableau
    tableau = [eye(dimen), -M, -ones(dimen, 1), q];
    % Let artificial variable enter the basis
    basis = 1:dimen; % A set of row indices in the tableau
    %comments by Shu Lu on June 18, 2020
    %I don't understand the comment. tableau has 2*dimen+2 columns instead
    %of 2*dimen+1.
    [~,locat] = min(tableau(:,end)); % Row of minimum element in column 
                                     % 2*dimen+1 (last of tableau)
    basis(locat) = 2*dimen+1; % Replace that index with the column
    cand = locat + dimen;
    %comments by Shu Lu on June 18, 2020
    %It is a little strange, because tableau(locat,2*dimen+1)=-1 here.
    %So pivot is just -tableau(locat,:)
    
    pivot = tableau(locat,:)/tableau(locat,2*dimen+1);
    % From each column subtract the column 2*dimen+1, multiplied by pivot
    
     %comments by Shu Lu on June 18, 2020
    %This is just to subtract the locat'th row from each row of tableau
    
    tableau = tableau - tableau(:,2*dimen+1)*pivot; 
    
    %comments by Shu Lu on June 18, 2020
    %This is just to negate the signs of all elements of the locat'th row
    %So, after this step, the 2*dimen+1'th column becomes a unit vector
    %and the original locat'th column becomes all -1
    
    tableau(locat,:) = pivot; % set all elements of row locat to pivot
    % Perform complementary pivoting
    while max(basis) == 2*dimen+1 && loopcount < maxits
        loopcount = loopcount + 1;
        eMs = tableau(:,cand); % This is used to check convergence (pivtol)
        missmask = eMs <= 0;  % Check if elements of eMs are less than zero
        quots = tableau(:,2*dimen+2)./eMs;
        quots(missmask) = Inf;
        [~,locat] = min(quots);
        % Check if at least one element is not missing
        
        %Comments by Shu Lu on June 19
        %It seems that the condition "abs(eMs(locat)) > pivtol"
        %is to detect whether the pivot element is
        %very close to 0. I guess, if the pivot element is very close to 0,
        %then it is similar to the case when all elements of the pivot
        %column are nonpositive.
        
        if  sum(missmask)~=dimen && abs(eMs(locat)) > pivtol 
            % Reduce tableau
            pivot = tableau(locat,:)/tableau(locat,cand);
            tableau = tableau - tableau(:,cand)*pivot;
            tableau(locat,:) = pivot;
            oldVar = basis(locat);
            % New variable enters the basis
            basis(locat) = cand;
            % Select next candidate for entering the basis
            
            
            
            if oldVar > dimen
                cand = oldVar - dimen;
            else
                cand = oldVar + dimen;
            end
        else
            %Comments by Shu Lu on June 19
            %rayTerm = true means that a ray is found. This implies that
            %the LCP does not have a solution if M is positive
            %semidefinite. See page 182 of [Ferris et al].       
            
            rayTerm = true; % Problem was solved
            break % Break out of the while loop
        end
    end
    % Return the solution to the LCP
    vars = zeros(2*dimen+1,1);
    vars(basis) = tableau(:,2*dimen+2).';
    w = vars(1:dimen,1);
    z = vars(dimen+1:2*dimen,1);
end

if rayTerm
    retcode = [2, loopcount];  % Ray termination
else
    retcode = [1, loopcount];  % Success
end