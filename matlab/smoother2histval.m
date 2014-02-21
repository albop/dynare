function smoother2histval(file, period, invars, outvars)
% This function takes values from oo_.SmoothedVariables and copies them into
% M_.histval.
% 
% INPUTS
%    file:        An optional *_results MAT file created by Dynare.
%                 If present, oo_.SmoothedVariables is read from there.
%                 Otherwise, it is read from the global workspace.
%    period:      An optional period number to use as the starting point
%                 for subsequent simulations. It should be between 1 and
%                 the number of observations that were used to produce the
%                 smoothed values. If absent, the last observation is used.
%    invars:      An optional cell array listing variables to read in
%                 oo_.SmoothedVariables. If absent, all the endogenous
%                 variables of the current model are used (from M_.endo_names)
%    outvars:     An optional cell array listing variables to be written in
%                 M_.histval. This cell must be of same length than invars,
%                 and there is a mapping between the input variable at the
%                 i-th position in invars, and the output variable at the
%                 i-th position in outvars. If absent, then taken as equal
%                 to invars.

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

global M_ options_ oo_

if nargin == 0
    if ~isfield(oo_, 'SmoothedVariables')
        error('Could not find smoothed variables; did you set the "smoother" option?')
    end
    smoothedvals = oo_.SmoothedVariables;
else
    S = load(file);
    if ~isfield(S, 'oo_') || ~isfield(S.oo_, 'SmoothedVariables')
        error('Could not find smoothed variables in file; is this a Dynare results file, and did you set the "smoother" option when producing it?')
    end
    smoothedvals = S.oo_.SmoothedVariables;
end

n = size(getfield(smoothedvals, fieldnames(smoothedvals){1}));

if n < M_.maximum_endo_lag
    error('Not enough observations to create initial conditions')
end

if nargin < 2
    period = n;
else
    if period > n
        error('The period that you indicated is beyond the data sample')
    end
    if period < M_.maximum_endo_lag
        error('The period that you indicated is too small to construct initial conditions')
    end
end

if nargin < 3
    invars = cellstr(M_.endo_names);
end

if nargin < 4
    outvars = invars;
else
    if length(invars) ~= length(outvars)
        error('The number of input and output variables is not the same')
    end
end

M_.endo_histval = repmat(oo_.steady_state, 1, M_.maximum_endo_lag);

for i = 1:length(invars)
    s = getfield(smoothedvals, invars{i});
    j = strmatch(outvars{i}, M_.endo_names, 'exact');
    if isempty(j)
        error(['Ouput variable ' outvars{i} ' does not exist'])
    end
    M_.endo_histval(j, :) = s((period-M_.maximum_endo_lag+1):period);
end

end
