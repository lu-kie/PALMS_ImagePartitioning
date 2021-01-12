function [u_up,a_up,b_up,c_up,partition_up]=scaleUpResult(u,a,b,c,partition,scale_factor)
% [u_up,a_up,b_up,c_up,partition_up]=scaleUpResult(u,a,b,c,partition,scale_factor)
% produces an up-scaled version partition_up of the input partition,
% the corresponding jet (a_up,b_up,c_up)
% and the piecewise affine-linear image u_up
% Input:
%	'a,b,c': piecewise constant field of first order polynomials computed for
%   a down-scaled input image with factor 1/scale_factor
%   'partition': label image of integers corresponding to a,b,c
%   'scale_factor': upscaling factor (=2 or =4)
% Output:
%   'a_up,b_up,c_up': up-scaled version of the field  a,b,c
%   'u_up': corresponding up-scaled  piecewise affine-linear image
%   'partition_up': corresponding scaled up partition
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
% Get image dimensions
[m,n,nr_channels] = size(u);
m_up = scale_factor*m;
n_up = scale_factor*n;
% Sort the segments of the input partition
L = unique(partition);
label_counts = hist(partition(:),L);
[~,J] = sort(label_counts,'descend');
% Get basic up-scaled partition
partition_up = imresize(partition,scale_factor,'nearest');
% Assign integer labels in terms of the segment size
partition_up = changem(partition_up,L,J);
% Smooth the boundaries of the up-scaled partition and favor large segments
se = strel('disk',scale_factor);
partition_up = imopen(partition_up,se);
% Get the corresponding up-scaled jet
a_up = zeros(m_up,n_up,nr_channels);
b_up = zeros(m_up,n_up,nr_channels);
c_up = zeros(m_up,n_up,nr_channels);

for j = 1:numel(J) % go through all segments
    % If label j has no segment in the up-scaled partition anymore, continue
    if isempty(find(partition_up == j , 1))
        continue;
    end
    % Get the jet element of the current segment
    jj = find(partition == J(j),1) + [0:nr_channels-1]*m*n;
    a_j = a(jj);
    b_j = b(jj);
    c_j = c(jj);
    % Channel-wise assign the jet element to the up-scaled jet 
    for ch = 1:nr_channels
        idx_chan = find((partition_up == j)) + (ch-1)*m_up*n_up;
        
        a_up(idx_chan) = a_j(ch);
        b_up(idx_chan) = b_j(ch);
        c_up(idx_chan) = c_j(ch);
    end
end
% Obtain the up-scaled piecewise affine-linear image
II = 1/scale_factor:1/scale_factor:n;
JJ = [1/scale_factor:1/scale_factor:m]';
u_up = c_up+...
    repmat(II,m_up,1,nr_channels).*a_up+...
    repmat(JJ,1,n_up,nr_channels).*b_up;
end