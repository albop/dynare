function C = eq(A,B) % --*-- Unitary tests --*--

% Overloads == operator for dynDates objects.
%
% INPUTS 
%  o A    dynDates object with n or 1 elements.
%  o B    dynDates object with n or 1 elements.
%
% OUTPUTS 
%  o C    column vector of max(n,1) elements (zeros or ones).

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

if ~isequal(nargin,2)
    error('dynDates::eq: I need exactly two input arguments!')
end

if ~isa(A,'dynDates') || ~isa(B,'dynDates')
    error(['dynDates::eq: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' have to be a dynDates objects!'])
end

if ~isequal(A.freq,B.freq)
    C = 0;
    return
end

if isequal(A.ndat, B.ndat)
    C = isequal(A.time, B.time);
else
    if isequal(A.ndat,1) || isequal(B.ndat,1)
        C = transpose(all(transpose(bsxfun(@eq,A.time,B.time))));
    else
        C = 0;
    end
end

%@test:1
%$ % Define some dynDates objects
%$ d1 = dynDate('1950Q1'):dynDate('1959Q4') ;
%$ d2 = dynDate('1960Q1'):dynDate('1979Q4') ;
%$ d3 = dynDate('1970M1'):dynDate('1979M12') ;
%$
%$ % Call the tested routine.
%$ t1 = d1==d1;
%$ t2 = d1==d2;
%$ t3 = d1==d3;
%$
%$ % Check the results.
%$ t(1) = dyn_assert(t1,1);
%$ t(2) = dyn_assert(t2,0);
%$ t(2) = dyn_assert(t3,0);
%$ T = all(t);
%@eof:1