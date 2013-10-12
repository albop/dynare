function C = min(varargin)
    
% Overloads the min function for dynDates objects.
    
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

switch nargin
  case 1
    A = varargin{1};
    switch length(A)
      case 0
        C= dynDates();
      case 1
        C = A;
      otherwise
        tmp = sortrows(A.time);
        C = dynDates();
        C.freq = A.freq;
        C.ndat = 1;
        C.time = tmp(1,:);
    end
  case 2
    A = varargin{1};
    if ~isa(A,'dynDates')
        error('dynDates::min: All inputs must be dynDates objects!')
    end
    B = varargin{2};
    if ~isa(B,'dynDates')
        error('dynDates::min: All inputs must be dynDates objects!')
    end
    C = min(A+B);
  otherwise
    if ~isa(varargin{1},'dynDates')
        error('dynDates::min: All inputs must be dynDates objects!')
    end
    A = varargin{1};
    for i=2:nargin
        if ~isa(varargin{i},'dynDates')
            error('dynDates::min: All inputs must be dynDates objects!')
        end
        A = A + varargin{i};
    end
    C = min(A);
end

%@test:1
%$ % Define some dates
%$ d3 = dynDates('1950q2');
%$ d4 = dynDates('1950Q3');
%$ d5 = dynDates('1950m1');
%$ d6 = dynDates('1948M6');
%$ m2 = min(d3,d4);
%$ i2 = (m2==d3);
%$ m3 = min(d5,d6);
%$ i3 = (m3==d6);
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i2,1);
%$ t(2) = dyn_assert(i3,1);
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ d = dynDates('1950Q2','1951Q3','1949Q1','1950Q4');
%$ m = min(d);
%$ i = (m==dynDates('1949Q1'));
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i,1);
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dates
%$ m = min(dynDates('1950Q2','1951Q3'),dynDates('1949Q1'),dynDates('1950Q4'));
%$ i = (m==dynDates('1949Q1'));
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i,1);
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define some dates
%$ m = min(dynDates('1950Q2'),dynDates('1951Q3'),dynDates('1949Q1'),dynDates('1950Q4'));
%$ i = (m==dynDates('1949Q1'));
%$
%$ % Check the results.
%$ t(1) = dyn_assert(i,1);
%$ T = all(t);
%@eof:4