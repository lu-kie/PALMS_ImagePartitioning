function [u_s,a_s,b_s] = LinewiseSolver(us_data,as_data,bs_data,s,gamma_s,eta,C_lin,S_lin,C_const,S_const,nr_threads)
% [u_s,a_s,b_s] = LinewiseSolver(us_data,as_data,bs_data,s,gamma_s,eta_s,C_lin,S_lin,C_const,S_const)
% is a wrapper to call the ADMM subproblem solver in C++
% Inputs:
%   us_data: Offset data
%   as_data: Horizontal slope data
%   bs_data: Vertical slope data
%   s: indicator of current direction
%   gamma_s: jump penalty of 1D problems
%   eta: data weight of 1D problems
%   C_lin,S_lin,C_const,S_const: Recurrence coefficients/Givens rotation
%   angles for the solver
%   nr_threads: number of threads for OpenMP multicore support (default: 32)
% Outputs:
%   u_s: linewise regularized functional values
%   a_s: linewise regularized horizontal slopes
%   b_s: linewise regularized vertical slopes
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
% Transform slope data according to the s-th direction
switch s
    case 1
        x_slope_data = bs_data;
        y_slope_data = as_data;
    case 2
        x_slope_data = as_data;
        y_slope_data = bs_data;
    case 3
        x_slope_data = as_data + bs_data;
        y_slope_data =  as_data - bs_data;
    case 4
        x_slope_data = as_data - bs_data;
        y_slope_data =  as_data + bs_data;
end
% Get current direction
D = getDirsAndWeights(4);
dir_s = D(:,s);
% Call the linewise Solver in C++
[u_s,x,y] = LinewiseSolver_mexWrapper(us_data,x_slope_data,y_slope_data,dir_s,gamma_s,eta,C_lin,S_lin,C_const,S_const,nr_threads);

% Back transform the slopes x and y to a_s and b_s
switch s
    case 1
        a_s = y;
        b_s = x;
    case 2
        a_s = x;
        b_s = y;
    case 3
        a_s = (x+y)/2;
        b_s = (x-y)/2;
    case 4
        a_s = (x+y)/2;
        b_s = (y-x)/2;
end
