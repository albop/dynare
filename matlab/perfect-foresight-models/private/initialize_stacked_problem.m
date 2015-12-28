function [options, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, dynamicmodel] = initialize_stack_solve_algo_7(endogenousvariables, options, M)

% Copyright (C) 2015 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.
    
if ~isfield(options.lmmcp,'lb')
    [lb, ub, pfm.eq_index] = get_complementarity_conditions(M, options.ramsey_policy);
    options.lmmcp.lb = repmat(lb, options.periods, 1);
    options.lmmcp.ub = repmat(ub, options.periods, 1);
end

y0 = endogenousvariables(:,1);
yT = endogenousvariables(:,options.periods+2);
z = endogenousvariables(:,2:options.periods+1);
illi = M.lead_lag_incidence';
[i_cols, junk,i_cols_j] = find(illi(:));
illi = illi(:,2:3);
[i_cols_J1, junk,i_cols_1] = find(illi(:));
i_cols_T = nonzeros(M.lead_lag_incidence(1:2,:)');
dynamicmodel = str2func([M.fname,'_dynamic']);