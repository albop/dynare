function C = plus(A,B) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function File} {@var{C} =} plus (@var{A},@var{B})
%! @anchor{@dynDates/plus}
%! @sp 1
%! Overloads the plus (addition) operator for the @ref{dynDates} class. Combines two dynDates objects, A and B, without removing repetitions
%! if A and B are not disjoints.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! @ref{dynDates} object.
%! @item B
%! @ref{dynDates} object.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item C
%! @ref{dynDates} object.
%! @end table
%! @end deftypefn
%@eod:

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

if isa(A,'dynDates') && isa(B,'dynDates')
    % Concatenate dynDates objects without removing repetitions if A and B are not disjoint sets of dates.
    if ~isequal(A.freq,B.freq)
        error(['dynDates::plus: Input arguments ''' inputname(1) ''' and ''' inputname(2) ''' must have common frequencies!'])
    end
    if isempty(B)
        C = A;
        return
    end
    if isempty(A)
        C = B;
        return
    end
    C = dynDates();
    C.freq = A.freq;
    C.time = [A.time; B.time];
    C.ndat = A.ndat+B.ndat;
elseif isa(A,'dynDates') && ( (isvector(B) && isequal(length(B),A.ndat) && all(isint(B))) || isscalar(B) && isint(B) || isequal(length(A),1) && isvector(B) && all(isint(B)))
    C.time = add_periods_to_array_of_dates(A.time, A.freq, B);
elseif isa(B,'dynDates') && ( (isvector(A) && isequal(length(A),B.ndat) && all(isint(A))) || isscalar(A) && isint(A) )
    C.time = add_periods_to_array_of_dates(B.time, B.freq, A);
else
    error('dynDates::plus: I don''t understand what you want to do! Check the manual.')
end

%@test:1
%$ % Define some dynDates objects
%$ d1 = dynDates('1950Q1','1950Q2') ;
%$ d2 = dynDates('1950Q3','1950Q4') ;
%$ d3 = dynDates('1950Q1','1950Q2','1950Q3','1950Q4') ;
%$
%$ % Call the tested routine.
%$ try
%$   e1 = d1+d2;
%$   e2 = d1+d2+d3;
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1) 
%$   t(2) = dyn_assert(e1==d3,1);
%$   t(3) = dyn_assert(e2==dynDates('1950Q1','1950Q2','1950Q3','1950Q4','1950Q1','1950Q2','1950Q3','1950Q4'),1);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dynDates objects
%$ d1 = dynDates('1950Q1');
%$ e1 = dynDates('1950Q2');
%$ e2 = dynDates('1950Q3');
%$ e3 = dynDates('1950Q4');
%$ e4 = dynDates('1951Q1');
%$ e5 = dynDates('1950Q2','1950Q3','1950Q4','1951Q1');
%$
%$ % Call the tested routine.
%$ try
%$   f1 = d1+1
%$   f2 = d1+2
%$   f3 = d1+3
%$   f4 = d1+4
%$   f5 = d1+transpose(1:4)
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dyn_assert(e1==f1,1);
%$   t(3) = dyn_assert(e2==f2,1);
%$   t(4) = dyn_assert(e3==f3,1);
%$   t(5) = dyn_assert(e4==f4,1);
%$   t(6) = dyn_assert(e5==f5,1);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define some dynDates objects
%$ d1 = dynDates('1950Q1');
%$ e1 = dynDates('1949Q4');
%$ e2 = dynDates('1949Q3');
%$ e3 = dynDates('1949Q2');
%$ e4 = dynDates('1949Q1');
%$ e5 = dynDates('1948Q4');
%$
%$ % Call the tested routine.
%$ try
%$   f1 = d1+(-1)
%$   f2 = d1+(-2)
%$   f3 = d1+(-3)
%$   f4 = d1+(-4)
%$   f5 = d1+(-5)
%$   t(1) = 1;
%$ catch
%$   t(1) = 0;
%$ end
%$
%$ if t(1)
%$   t(2) = dyn_assert(e1==f1,1);
%$   t(3) = dyn_assert(e2==f2,1);
%$   t(4) = dyn_assert(e3==f3,1);
%$   t(5) = dyn_assert(e4==f4,1);
%$   t(6) = dyn_assert(e5==f5,1);
%$ end
%$ T = all(t);
%@eof:3