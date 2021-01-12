function L = getPartitioningFromJetField(a,b,c)
% L = getPartitioningFromJetField(a,b,c) 
% generates a partition L image from the piecewise constant field of
% first order polynomials corresponding to its slopes a,b and offsets c
% computed by affineLinearPotts_ADMM.
% Each segment is identified with an integer stored in L
% Input:
%   piecewise constant field of first order polynomials a,b,c
% Output: 
%   L: integer image encoding the partition
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
[m,n,~] = size(c);
tol = 10^-2; % allowed relative difference between neighboring first coefficients

% Get adjacency matrix of generic 8-connected graph
I_1 = (find ( mod([1:m*n],m) ~= 0))';
J_1 = I_1 + 1;
I_2 = (find ([1:m*n] + m <= m*n))';
J_2 = I_2 + m;

I_3 = (find ( mod([1:m*n],m) ~= 0 & [1:m*n] + m <= m*n))';
J_3 = I_3 + m + 1;

I_4 = find ( mod([1:m*n],m) ~= 1 & [1:m*n] + m <= m*n)';
J_4 = I_4 + m - 1;

I_8 = [I_1;I_2;I_3;I_4];
J_8 = [J_1;J_2;J_3;J_4];
% Find segment edges
dirs = [1 0 1 -1;
        0 1 1  1];
M = repmat([1:m]',[1,n]);
N = repmat([1:n],[m,1]);
I = [];
J = [];
for s = 1:size(dirs,2)
    d_s = dirs(:,s);
    % Find edges for direction d_s by checking the respective jet values
    a_shift = circshift(a,-d_s);
    b_shift = circshift(b,-d_s);
    c_shift = circshift(c,-d_s);
    
    B = max(abs(a-a_shift) ./ (abs(a)+abs(a_shift)),[],3) > tol |...
        max(abs(b-b_shift) ./ (abs(b)+abs(b_shift)),[],3) > tol |...
        max(abs(c-c_shift) ./ (abs(c)+abs(c_shift)),[],3) > tol &...
        M+d_s(1) >= 1 & M+d_s(1) <= m & N+d_s(2) <= n;% Exclude image boundary
    
    
    % Store edge locations
    I_curr = find(B);
    I = [I;I_curr];
    J = [J;I_curr + d_s(1) + m*d_s(2)];
    
end
% Exclude segment boundary edges from 8-neighborhood edges
[~,idx] = intersect([I_8 J_8],[I J],'rows');
I_8(idx)=[];
J_8(idx)=[];

% Create sparse adjacency matrix
S = sparse(I_8,J_8,ones(size(I_8)),m*n, m*n);
% Get connected components from the associated undirected graph
G = graph(S,'upper');
L = reshape(conncomp(G),[m n]);
end
