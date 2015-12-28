function [endogenousvariables, info] = solve_stacked_problem(endogenousvariables, exogenousvariables, steadystate, M, options);

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
    
[options, y0, yT, z, i_cols, i_cols_J1, i_cols_T, i_cols_j, i_cols_1, dynamicmodel] = ...
    initialize_stacked_problem(endogenousvariables, options, M);

[y, check] = dynare_solve(@perfect_foresight_problem,z(:),options, ...
                         dynamicmodel, y0, yT, ...
                         exogenousvariables, M.params, steadystate, ...
                         M.maximum_lag, options.periods, M.endo_nbr, i_cols, ...
                         i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                         M.NNZDerivatives(1));

if all(imag(y)<.1*options.dynatol.x)
    if ~isreal(y)
        y = real(y);
    end
else
    check = 1;
end

endogenousvariables = [y0 reshape(y, M.endo_nbr, options.periods) yT];

if check
    info.status = false;
else
    info.status = true;
end