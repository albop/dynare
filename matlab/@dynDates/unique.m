function dd = unique(dd)
% unique method for dynDates class.

%@info:
%! @deftypefn {Function File} {@var{a} =} unique (@var{a})
%! @anchor{dynDates/unique}
%! @sp 1
%! Unique method for the Dynare dates class (removes repetitions).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Object instantiated by @ref{dynDates}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item a
%! Object instantiated by @ref{dynDates}.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2012-2013 Dynare Team
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

if ~isa(dd,'dynDates')
    error(['dynDates::unique: Input argument ' inputname(dd) ' has to be a dynDates object.'])
end

if dd.ndat==1
    return
end

[tmp,id,jd] = unique(dd.time,'rows');
dd.time = dd.time(sort(id),:);
dd.ndat = size(dd.time,1);

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950q1';
%$ B4 = '1945Q3';
%$ B5 = '1950Q2'; 
%$
%$ % Define expected results.
%$ e.time = [1953 4; 1950 1; 1945 3; 1950 2];
%$ e.freq = 4;
%$ e.ndat = 4;
%$
%$ % Call the tested routine.
%$ d = dynDates(B1,B2,B3,B4,B5);
%$ d = d.unique;
%$ 
%$ % Check the results.
%$ t(1) = dyn_assert(d.time,e.time);
%$ t(2) = dyn_assert(d.freq,e.freq);
%$ t(3) = dyn_assert(d.ndat,e.ndat);
%$ T = all(t);
%@eof:1