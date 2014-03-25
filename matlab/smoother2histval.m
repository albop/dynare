function smoother2histval(infile, invars, period, outfile, outvars)
% This function takes values from oo_.SmoothedVariables and copies them into
% M_.histval.
% 
% INPUTS
%    infile:      An optional *_results MAT file created by Dynare.
%                 If present, oo_.SmoothedVariables is read from there.
%                 Otherwise, it is read from the global workspace.
%    invars:      An optional cell array listing variables to read in
%                 oo_.SmoothedVariables. If absent, all the endogenous
%                 variables present in oo_.SmoothedVariables are used.
%    period:      An optional period number to use as the starting point
%                 for subsequent simulations. It should be between 1 and
%                 the number of observations that were used to produce the
%                 smoothed values. If absent, the last observation is used.
%    outfile:     An optional MAT file in which to save the histval structure.
%                 If absent, the output will be written in M_.endo_histval
%    outvars:     An optional cell array listing variables to be written in
%                 outfile or M_.endo_histval. This cell must be of same
%                 length than invars, and there is a mapping between the input
%                 variable at the i-th position in invars, and the output
%                 variable at the i-th position in outvars. If absent, then
%                 taken as equal to invars.
%
% The function also uses the value of option_.parameter_set

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

if nargin == 0 || isempty(infile)
    if ~isfield(oo_, 'SmoothedVariables')
        error('Could not find smoothed variables; did you set the "smoother" option?')
    end
    smoothedvals = oo_.SmoothedVariables;
else
    S = load(infile);
    if ~isfield(S, 'oo_') || ~isfield(S.oo_, 'SmoothedVariables')
        error('Could not find smoothed variables in file; is this a Dynare results file, and did you set the "smoother" option when producing it?')
    end
    smoothedvals = S.oo_.SmoothedVariables;
end

% Hack to determine if oo_.SmoothedVariables was computed after a Metropolis
if isstruct(getfield(smoothedvals, fieldnames(smoothedvals){1}))
    post_metropolis = 1;
else
    post_metropolis = 0;
end

% If post-Metropolis, select the parameter set
if isempty(options_.parameter_set)
    if post_metropolis
        smoothedvals = smoothedvals.Mean;
    end
else
    switch options_.parameter_set
      case 'calibration'
        if post_metropolis
            error('Option parameter_set=calibration is not consistent with computed smoothed values.')
        end
      case 'posterior_mode'
        if post_metropolis
            error('Option parameter_set=posterior_mode is not consistent with computed smoothed values.')
        end
      case 'posterior_mean'
        if ~post_metropolis
            error('Option parameter_set=posterior_mean is not consistent with computed smoothed values.')
        end
        smoothedvals = smoothedvals.Mean;
      case 'posterior_median'
        if ~post_metropolis
            error('Option parameter_set=posterior_median is not consistent with computed smoothed values.')
        end
        smoothedvals = smoothedvals.Median;
      otherwise
        error([ 'Option parameter_set=' options_.parameter_set ' unsupported.' ])
    end
end

% Determine number of periods
n = size(getfield(smoothedvals, fieldnames(smoothedvals){1}));

if n < M_.maximum_endo_lag
    error('Not enough observations to create initial conditions')
end

if nargin < 2 || isempty(invars)
    invars = fieldnames(smoothedvals);
end

if nargin < 3 || isempty(period)
    period = n;
else
    if period > n
        error('The period that you indicated is beyond the data sample')
    end
    if period < M_.maximum_endo_lag
        error('The period that you indicated is too small to construct initial conditions')
    end
end

if nargin < 5 || isempty(outvars)
    outvars = invars;
else
    if length(invars) ~= length(outvars)
        error('The number of input and output variables is not the same')
    end
end

endo_histval = repmat(oo_.steady_state, 1, M_.maximum_endo_lag);

for i = 1:length(invars)
    s = getfield(smoothedvals, invars{i});
    j = strmatch(outvars{i}, M_.endo_names, 'exact');
    if isempty(j)
        if strncmp('AUX_', outvars{i}, 4)
            warning(['smoother2histval: output auxiliary variable ' outvars{i} ' does not exist, ignoring.'])
        else
            error(['smoother2histval: output variable ' outvars{i} ' does not exist.'])
        end
    else
        endo_histval(j, :) = s((period-M_.maximum_endo_lag+1):period);
    end
end

if nargin < 4 || isempty(outfile)
    M_.endo_histval = endo_histval;
else
    save(outfile, 'endo_histval')
end

end
