function a = isestimation()

% Returns 1 if we are currently estimating a model, 0 otherwise.

% Copyright (C) 2014 Dynare Team
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

a = 0;

tmp = struct2cell(dbstack);
list_of_previously_called_routines = transpose(tmp(2,:));

if ~isempty(strmatch('dsge_likelihood', list_of_previously_called_routines, 'exact')) || ...
    ~isempty(strmatch('dsge_var_likelihood', list_of_previously_called_routines, 'exact')) || ...
    ~isempty(strmatch('non_linear_dsge_likelihood', list_of_previously_called_routines, 'exact')) || ...
    ~isempty(strmatch('simulated_moments_estimation', list_of_previously_called_routines, 'exact'))
    a = 1;
end