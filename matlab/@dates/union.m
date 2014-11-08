function D = union(varargin) % --*-- Unitary tests --*--

% Overloads union function for dates objects (removes repetitions if any).
%
% INPUTS 
%  o A    dates object.
%  o B    dates object.
%  o C    dates object.
%  o ...
%
% OUPUTS 
%  o D    dates object (elements are sorted by increasing order).

% Copyright (C) 2013 Dynare Team
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

if isequal(nargin,1)
    D = sort(unique(varargin{1}));
    return;
end

D = sort(unique(horzcat(varargin{:})));

%@test:1
%$ % Define some dates objects
%$ d1 = dates('1950Q1'):dates('1959Q4') ;
%$ d2 = dates('1960Q1'):dates('1969Q4') ;
%$ d3 = dates('1970Q1'):dates('1979Q4') ;
%$
%$ % Call the tested routine.
%$ e1 = union(d1);
%$ e2 = union(d1,d2);
%$ e3 = union(d1,d2,d3);
%$ e4 = union(d1,d2,d3,d2+d3);
%$ e5 = union(d1,d2,d3,d2);
%$
%$ % Check the results.
%$ t(1) = dassert(isequal(e1,d1),1);
%$ t(2) = dassert(isequal(e2,d1+d2),1);
%$ t(3) = dassert(isequal(e3,d1+d2+d3),1);
%$ t(4) = dassert(isequal(e4,d1+d2+d3),1);
%$ t(5) = dassert(isequal(e5,d1+d2+d3),1);
%$ T = all(t);
%@eof:1