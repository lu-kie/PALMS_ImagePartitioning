function [dirs, omegas] = getDirsAndWeights(nr_dirs)
% This function returns the directions of the discrete gradient and the
% corresponding weights
% Input:
%   nr_dirs: 2 directions amount to the anisotropic 4-neighborhood and 4
%   directions amount to the near-isotropic 8-neighborhood
% Output: 
%   dirs: 2 x nr_dirs matrix whose columns denote the directions
%   weights: corresponding weights
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
switch nr_dirs
    case 2
        dirs = [0 1;
                1 0];
        omegas = ones(2,1);
    case 4
        dirs = [0 1 1 1;
                1 0 1 -1];
        omegas = [sqrt(2)-1 sqrt(2)-1 1-sqrt(2)/2  1-sqrt(2)/2]';
end
