function B = mtimes(A,n)

% Overloads the times operator (*). Returns dynDates object A replicated n times.
%
% INPUTS 
%  o A    dynDates object with m elements.
%
% OUTPUTS 
%  o B    dynDates object with m*n elements.
%
% EXAMPLE 1
%  If A = dynDates('2000Q1'), then B=A*3 is a dynDates object equal to dynDates('2000Q1','2000Q1','2000Q1')  
%
% EXAMPLE 2
%  If A = dynDates('2003Q1','2009Q2'), then B=A*2 is a dynDates object equal to dynDates('2003Q1','2009Q2','2003Q1','2009Q2')

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

if ~(isscalar(n) && isint(n))
    error('dynDates::m: First and second input arguments have to be a dynDates object and a scalar integer!')
end

A.time = repmat(A.time,n,1);
A.ndat = A.ndat*n;