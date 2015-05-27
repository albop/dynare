function oo_ = convert_oo_(oo_, ver)
%function oo_ = convert_oo_(oo_, ver)
% Converts oo_ from oo_.dynare_version to ver
%
% INPUTS
%    oo_    [struct]    dynare output struct
%    ver    [string]    desired oo_ output version
%
% OUTPUTS
%    oo_    [struct]    dynare output struct
%
% SPECIAL REQUIREMENTS
%    none

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

check_valid_ver(ver);

if isfield(oo_, 'dynare_version')
    ver_orig = oo_.dynare_version;
else
    ver_orig = '4.4.3';
end

if strcmp(ver_orig, ver)
    return;
end

if ver_less_than(ver_orig, '4.5.0') && ver_greater_than_equal(ver, '4.5.0')
    oo_.exo_simul = oo_.exo_simul';
end
end
