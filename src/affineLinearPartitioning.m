function [u,partition,a,b,c] = affineLinearPartitioning(f,varargin)
% [u,partition,a,b,c] = affineLinearPartitioning(f,parameters)
% computes the result of the ADMM approach to the piecewise affine-linear
% Mumford-Shah model based on a Taylor jet formulation [1],[2]
%
% Inputs:
%   f: Gray, Color or vector-valued double image to be partitioned
%   varargin:
%   'downScale': determines if the input image is scaled down before the computations
%   and the results are scaled up afterwards to accelerate computation time, the factor is 2 for downScale=1
%   and 4 for downScale=2 (default=0, no scaling)
%   'gamma': positive user-driven parameter that penalizes the total length
%   of the segment boundaries. Larger choices yield fewer segments
%   'maxIter' is the maximum number of allowed ADMM iterations
%   (default: 500)
%   'muNuStep' is the progression parameter which determines how fast the
%   splitting variables have to become equal
%   (corresponds to the factor phi in the paper) (default: 1.3)
%   'isotropic' determines if the near-isotropic (1) or anisotropic
%   discretization (0) is used (default: 1)
%   'splitTol' parameter for the relative difference stopping criterion (default:10^-2)
%   'u_0,a_0,b_0': initializations (default: f,0,0)
%   'nr_threads': number of threads for OpenMP multicore support (default: 32)
%   'verbose': toggles wether iterations and total iteration number are displayed (default: true)
%
%
% Outputs:
%   u: piecewise affine-linear approximation of the input image f
%   partition: integer-valued label image; each integer represents a
%   segment of the partitioning corresponding to u
%   a: horizontal slopes of the affine-linear polynomials in each pixel
%   corresponding to u
%   b: vertical slopes of the affine-linear polynomials in each pixel
%   corresponding to u
%   c: offsets (matrix origin) of the affine-linear polynomials in each pixel
%   corresponding to u
%
%
% Ref:
% [1] Kiefer,Storath,Weinmann - PALMS Image Partitioning -
%	A new Parallel Algorithm for the
%	Piecewise Affine-Linear Mumford-Shah Model, IPOL
%
% [2] Kiefer,Storath,Weinmann - An Efficient Algorithm for the Piecewise
% Affine-Linear Mumford-Shah Model Based on a Taylor Jet Splitting, IEEE
% Transactions on Image Processing
%
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
% Parse options
ip = inputParser;
addParameter(ip,'gamma', 1.0);
addParameter(ip,'downScale', 0);
addParameter(ip,'maxIter', 500);
addParameter(ip,'muNuStep', 1.3);
addParameter(ip,'isotropic', '1');
addParameter(ip,'splitTol', 10^-2);
addParameter(ip,'u_0', f);
addParameter(ip,'a_0', zeros(size(f)));
addParameter(ip,'b_0', zeros(size(f)));
addParameter(ip,'nr_threads', 32);
addParameter(ip,'verbose', true);

parse(ip, varargin{:});
par = ip.Results;

if ~isa(par.gamma,'numeric')
    par.gamma=str2num(par.gamma);
end
if ~isa(par.muNuStep,'numeric')
    par.muNuStep=str2num(par.muNuStep);
end
if ~isa(par.isotropic,'integer') && ~isa(par.isotropic,'numeric')
    par.isotropic=str2num(par.isotropic);
end
if ~isa(par.downScale,'integer') && ~isa(par.downScale,'numeric')
    par.zoomDown=str2num(par.downScale);
end


if par.isotropic == 0
    nr_dirs = 2;
elseif par.isotropic == 1
    nr_dirs = 4;
else
    disp(['Invalid isotropic parameter, assert 8 neighborhood'])
    nr_dirs = 4;
end

% Check input arguments
assert(par.gamma > 0, 'Edge penalty gamma must be > 0.');
assert(par.muNuStep > 1, 'Progression muNuStep must be > 1.');
assert(par.maxIter > 1, 'Number of max. Iterations maxIter must be > 1.');
assert(par.splitTol > 0, 'Stopping parameter etaIter must be > 0.');
assert(par.nr_threads > 0, 'Number of threads must be > 0.');
assert(isequal(size(par.u_0) , size(f)) && isequal(size(par.a_0) , size(f)) && isequal(size(par.b_0) , size(f)),...
    'Initializations must have the same dimensions as input image.');

% Scale down the input image if requested
if par.downScale > 0
    if par.downScale == 2
        scale_factor = 4;
    else
        scale_factor = 2;
    end
    % adjust the model penalty to the smaller image size
    par.gamma = par.gamma / scale_factor;
    % scale down the input image (cut off surplus rows and columns)
    m_new = size(f,1)-mod(size(f,1),scale_factor);
    n_new = size(f,2)-mod(size(f,2),scale_factor);
    f = f(1:m_new,1:n_new,:);
    f = imresize(f,1/scale_factor);
    % scale down the initial values accordingly
    par.u_0 = par.u_0(1:m_new,1:n_new,:);
    par.a_0 = par.a_0(1:m_new,1:n_new,:);
    par.b_0 = par.b_0(1:m_new,1:n_new,:);
    par.u_0 = imresize(par.u_0,1/scale_factor);
    par.a_0 = imresize(par.a_0,1/scale_factor);
    par.b_0 = imresize(par.b_0,1/scale_factor);
    if par.verbose
        disp(['The result is computed for a scaled down version of the input image by factor ' num2str(scale_factor) ...
            ' and gamma is adjusted accordingly.'])
    end
end
% Peform ADMM strategy
[u,a,b,c] = affineLinearMS_ADMM(f,par.gamma,nr_dirs,par.maxIter,par.splitTol,par.muNuStep,par.u_0,par.a_0,par.b_0,par.nr_threads,par.verbose);
% Get the associated partition from the computed (pcw constant) jet field
partition = getPartitioningFromJetField(a,b,c);
% Scale up the result if the input was scaled down before
if par.downScale > 0
    [u,a,b,c,partition]=scaleUpResult(u,a,b,c,partition,scale_factor);
end

end
