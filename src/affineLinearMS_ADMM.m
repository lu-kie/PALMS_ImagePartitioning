function [u,a,b,c] = affineLinearMS_ADMM(f,gamma,nr_dirs,MAX_IT,splitTol,prog,u_0,a_0,b_0,nr_threads,verbose)
% [u,a,b,c] = affineLinearMS_ADMM(f,gamma,nr_dirs,MAX_IT,split_tol,prog,u_0,a_0,b_0)
% computes the result of the ADMM approach to the piecewise affine-linear
% Mumford-Shah model based on a Taylor jet splitting. This function should
% not be called directly, but via the function affineLinearPartitioning
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
% Initialization
[~,omegas]  = getDirsAndWeights(nr_dirs); % weights of directions
[m,n,nr_channels] = size(f);
max_stripe_length = max(m,n); % max size of 1D subproblems
% Allocate splitting variables and multipliers
us(1:nr_dirs,:) = {u_0};
as(1:nr_dirs,:) = {a_0};
bs(1:nr_dirs,:) = {b_0};

lambdas(1:nr_dirs,1:nr_dirs) = {zeros(m,n,nr_channels)};
taus(1:nr_dirs,1:nr_dirs)    = {zeros(m,n,nr_channels)};
rhos(1:nr_dirs,1:nr_dirs)    = {zeros(m,n,nr_channels)};
% Initial coupling penalties
mu = 10^-3;
nu = min(450*gamma*mu,1);
% ADMM Iterations
for i = 1:MAX_IT
    % Data weight of subproblems
    eta = sqrt((2+mu*nr_dirs*(nr_dirs-1))/(nu*nr_dirs*(nr_dirs-1)));
    % Calculate recurrence coefficients for the subproblems, i.e., the Givens rotation angles
    [C_lin,S_lin,C_const,S_const] = calcGivensAngles(max_stripe_length,eta);
    % Solve linewise jet problems for each direction (corresponding to the first line of equation (16))
    for s = 1:nr_dirs
        % jump penalty of univariate subproblems (corresponds to gamma' after eq. (20) in section 2.2)
        gamma_s = (2*omegas(s)*gamma) / ((nr_dirs-1)*nu);
        % calc input data for subproblem solver
        [us_data,as_data,bs_data] = computeLinewiseData(f,lambdas,taus,rhos,us,as,bs,mu,nu,nr_dirs,s);
        % Call the solver of the subproblems
        [us{s},as{s},bs{s}] = LinewiseSolver(us_data,as_data,bs_data,s,gamma_s,eta,C_lin,S_lin,C_const,S_const,nr_threads);
    end
    % Update Lagrange multipliers
    [lambdas,taus,rhos]= updateMultipliers(us,as,bs,lambdas,taus,rhos,mu,nu,nr_dirs);
    % Check stopping criterion
    stop_bool = partitioningStopCrit(us,as,bs,nr_dirs,splitTol);
    if stop_bool
        if verbose
            fprintf(['\nTotal number iterations: ' num2str(i) '\n']);
        end
        break;
    end
    % Update coupling penalties
    mu = mu*prog;
    nu = nu*prog;
    if verbose
        fprintf('*');
    end
end
if i == MAX_IT
    fprintf(['\nWarning: Max number of iterations (' num2str(MAX_IT) ') reached\n'])
end
% Output u,a,b,c
[u,a,b] = splittingMeans(us,as,bs,nr_dirs,nr_channels);
c = u-repmat(1:n,m,1,nr_channels).*a-repmat([1:m]',1,n,nr_channels).*b;

end


%% Auxiliary functions
function [u_mean,a_mean,b_mean] = splittingMeans(us,as,bs,nr_dirs,nr_channels)
% This function returns the means of all splitting variables
m = size(us{1},1);
n = size(us{2},2);
u_mean = zeros(m,n,nr_channels);
a_mean = zeros(m,n,nr_channels);
b_mean = zeros(m,n,nr_channels);

for s = 1:nr_dirs
    u_mean = u_mean + us{s};
    a_mean = a_mean + as{s};
    b_mean = b_mean + bs{s};
end
u_mean = u_mean / nr_dirs;
a_mean = a_mean / nr_dirs;
b_mean = b_mean / nr_dirs;
end

function [us_data,as_data,bs_data] = computeLinewiseData(f,lambdas,taus,rhos,us,as,bs,mu,nu,nr_dirs,s)
% This function computes the data for the subproblems as in the lines 6-8 of Algorithm 1
[m,n,nr_channels] = size(f);

% Allocation
w_s = zeros(m,n,nr_channels);
y_s = zeros(m,n,nr_channels);
z_s = zeros(m,n,nr_channels);
% Already updated splitting variables
for r = 1:s-1
    w_s = w_s + us{r}...
        + lambdas{r,s}/mu;
    
    y_s = y_s + as{r}...
        + taus{r,s}/nu;
    
    z_s = z_s + bs{r}...
        + rhos{r,s}/nu;
end
% Not yet updated splitting variables
for t = s+1:nr_dirs
    w_s = w_s + us{t}...
        - lambdas{s,t} /mu;
    
    y_s = y_s + as{t}...
        - taus{s,t} /nu;
    
    z_s = z_s + bs{t}...
        - rhos{s,t} /nu;
end
% Gather the above
% Offset data
us_data = (2*f+mu*nr_dirs*w_s) / (2+mu*nr_dirs*(nr_dirs-1));
% Vertical slope data
as_data = y_s/(nr_dirs-1);
% Horizontal slope data
bs_data = z_s/(nr_dirs-1);
end

function [lambdas,taus,rhos]=...
    updateMultipliers(us,as,bs,lambdas,taus,rhos,mu,nu,nr_dirs)
% This function performs the gradient ascents of the Lagrange multipliers in the ADMM scheme
% (corresponding to eq. (16) and lines 11-15 of Algorithm 1, respectively.)

for s = 1:nr_dirs
    for t = s+1:nr_dirs
        lambdas{s,t} = lambdas{s,t} + mu*(us{s}-us{t});
        taus{s,t} = taus{s,t} + nu*(as{s}-as{t});
        rhos{s,t} = rhos{s,t} + nu*(bs{s}-bs{t});
        
    end
end
end

function stop_bool = partitioningStopCrit(us,as,bs,nr_dirs,splitTol)
% This function checks if the relative deviation between consecutive splitting
% variables (offsets and slopes) is smaller than split_tol
% (corresponding to line 17 of Algorithm 1)
u_1 = us{1};
u_2 = us{2};
a_1 = as{1};
a_2 = as{2};
b_1 = bs{1};
b_2 = bs{2};

stop_bool = false;
if max(abs(u_1(:)-u_2(:))./(abs(u_1(:))+abs(u_2(:)))) > splitTol
    return;
end

if max(abs(a_1(:)-a_2(:))./(abs(a_1(:))+abs(a_2(:)))) > splitTol
    return;
end

if max(abs(b_1(:)-b_2(:))./(abs(b_1(:))+abs(b_2(:)))) > splitTol
    return;
end

if nr_dirs > 2
    u_3 = us{3};
    u_4 = us{4};
    a_3 = as{3};
    a_4 = as{4};
    b_3 = bs{3};
    b_4 = bs{4};
    
    if max(abs(u_3(:)-u_4(:))./(abs(u_3(:))+abs(u_4(:)))) > splitTol
        return;
    end
    
    if max(abs(a_3(:)-a_4(:))./(abs(a_4(:))+abs(a_4(:)))) > splitTol
        return;
    end
    
    if max(abs(b_3(:)-b_4(:))./(abs(b_3(:))+abs(b_4(:)))) > splitTol
        return;
    end
    
end
stop_bool = true;
end
