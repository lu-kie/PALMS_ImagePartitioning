function [C_lin,S_lin,C_const,S_const] = calcGivensAngles(n,eta)
% [C_lin,S_lin,C_const,S_const] = calcGivensAngles(n,eta)
% computes the recurrence coefficients necessary for the solver of the
% univariate subproblems in the ADMM scheme. It performs effectively a QR
% decomposition of the system of a least squares problem and saves the
% necessary rotation angles to obtain the upper triangluar matrix R
%  
% Inputs 
%	n: (maximum) size of 1D problems
%   eta: (positive) weight of the offset/pixel values in the subproblem 
% Output: recurrence coefficients/Givens rotation angles
%
% Copyright (c) 2020 Lukas Kiefer <lukas.kiefer2@gmail.com>
% 
% This program is free software: you can redistribute it and/or modify it 
% under the terms of the GNU Affero General Public License as published by 
% the Free Software Foundation, either version 3 of the License, 
% or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but 
% WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero General Public 
% License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%
% System matrices:
% The matrix A_q of eq. (27) can be treated separately w.r.t. the 2nd column
% removed and rows 3i+1 and 3i+3 (i=0...q-1) and the 2nd column with rows
% 3i+2 (i=0,...,q-1).
% In the following lines A corresponds to the first submatrix
% and C_lin,S_lin to the according Givens coefficients and B corresponds
% to the second submatrix and C_const,S_const to the according Givens coefficients
A = zeros(2*n,2);
A(1:2:end,1) = eta*[1:n]';
A(2:2:end,1) = 1;
A(1:2:end,2) = eta;

B = ones(n,1);

C_lin = zeros(size(A));
S_lin = zeros(size(A));
C_const = zeros(size(B));
S_const = zeros(size(B));
% Givens Rotations Coefficients
% linear part
for i = 2:size(A,1)
    for j=1:min(2,i-1)
        if A(i,j)==0
            continue
        end
        % Update Givens rotation coefficients
        rho = sign(A(j,j))*sqrt(A(j,j)^2+A(i,j)^2);
        c = A(j,j)/rho;
        s = A(i,j)/rho;
        C_lin(i,j) = c;
        S_lin(i,j) = s;
        % Update A
        Aj_old = A(j,:);
        Ar_old = A(i,:);
	    A(j,:) = c*Aj_old+s*Ar_old;
        A(i,:) = -s*Aj_old+c*Ar_old;
    end
end
% pure constant part
for i = 2:size(B,1)
    % Update Givens rotation coefficients
    rho = sign(B(1,1))*sqrt(B(1,1)^2+B(i,1)^2);
    c = B(1,1)/rho;
    s = B(i,1)/rho;
    C_const(i,1) = c;
    S_const(i,1) = s;
    % Update B
    B(1,:) = c*B(1,:)+s*B(i,:);
end
end
