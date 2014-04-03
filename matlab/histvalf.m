function histvalf(fname)

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

global M_ oo_

if ~exist(fname)
    error(['Can''t find datafile: ' fname ]);
end

M_.endo_histval = repmat(oo_.steady_state, 1, M_.maximum_endo_lag);

S = load(fname);

outvars = fieldnames(S);

for i = 1:length(outvars)
    ov_ = outvars{i};
    if ov_(end) == '_'
        ov = ov_(1:end-1);
        j = strmatch(ov, M_.endo_names, 'exact');
        if isempty(j)
            error(['smoother2histval: output variable ' ov ' does not exist.'])
        end
    else
        % Lagged endogenous, search through aux vars
        z = strsplit(ov_, '_');
        ov = z{1};
        lead_lag = str2num(z{2});
        j = [];
        for i = 1:length(M_.aux_vars)
            if M_.aux_vars(i).type ~= 1 
                continue
            end
            orig_var = deblank(M_.endo_names(M_.aux_vars(i).orig_index, :));
            if strcmp(orig_var, ov) && M_.aux_vars(i).orig_lead_lag == lead_lag
                j = M_.aux_vars(i).endo_index;
            end
        end
        if isempty(j)
            continue
        end
    end
    M_.endo_histval(j, :) = getfield(S, ov_);
end
