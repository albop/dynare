function B = uplus(A)

% Overloads the unary plus operator for dynDates objects. Shifts all the elements by one period.
%
% INPUTS 
%  o A    dynDates object with n elements.
%
% OUTPUTS 
%  o B    dynDates object with n elements.

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

B = dynDates(A);
B.time(:,2) = B.time(:,2)-1;
idx = find(B.time(:,2)==0);
B.time(idx,1) = B.time(idx,1)-1;
B.time(idx,2) = B.freq;

%@test:1
%$ % Define some dates
%$ date_1 = '1950Y';
%$ date_2 = '1950Q2';
%$ date_3 = '1950Q1';
%$ date_4 = '1950M2';
%$ date_5 = '1950M1';
%$
%$ % Call the tested routine.
%$ d1 = dynDates(date_1); d1 = -d1;
%$ d2 = dynDates(date_2); d2 = -d2;
%$ d3 = dynDates(date_3); d3 = -d3;
%$ d4 = dynDates(date_4); d4 = -d4;
%$ d5 = dynDates(date_5); d5 = -d5;
%$ i1 = (d1==dynDates('1949Y'));
%$ i2 = (d2==dynDates('1950Q1'));
%$ i3 = (d3==dynDates('1949Q4'));
%$ i4 = (d4==dynDates('1950M1'));
%$ i5 = (d5==dynDates('1949M12'));
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i1,1);
%$ t(2) = dyn_assert(i2,1);
%$ t(3) = dyn_assert(i3,1);
%$ t(4) = dyn_assert(i4,1);
%$ t(5) = dyn_assert(i5,1);
%$ T = all(t);
%@eof:1