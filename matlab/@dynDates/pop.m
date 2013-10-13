function B = pop(A,a) % --*-- Unitary tests --*--

% pop method for dynDates class (removes a date).
%
% INPUTS 
%  o A     dynDates object.
%  o a     dynDates object with one element, string which can be interpreted as a date or integer scalar.
%
% OUTPUTS 
%  o B     dynDates object (with B.ndat==A.ndat-1).
%
% REMARKS 
%  If a is a date appearing more than once in A, then only the last occurence is removed. If one wants to
%  remove all the occurences of a in A, the setdiff function should be used instead.

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

if nargin<2
    % Remove last date
    B = dynDates();
    B.ndat = A.ndat-1;
    B.freq = A.freq;
    B.time = NaN(B.ndat,2);
    B.time = A.time(1:end-1,:);
    return
end

if ~(isa(a,'dynDates') || (ischar(a) && isdate(a)) || (isscalar(a) && isint(a)))
    error(['dynDates::pop: Input argument ' inputname(2) ' has to be a dynDates object with a single element, a string (which can be interpreted as a date) or an integer.'])
end

if ischar(a)
    a = dynDates(a);
end

B = dynDates();
B.ndat = A.ndat-1;
B.freq = A.freq;
B.time = NaN(B.ndat,2);

if isscalar(a) && isint(a)
    idx = find(transpose(1:A.ndat)~=a);
    B.time = A.time(idx,:);
else
    if ~isequal(A.freq,a.freq)
        error('dynDates::pop: Inputs must have common frequency!')
    end
    idx = find(A==a);
    jdx = find(transpose(1:A.ndat)~=idx(end));
    B.time = A.time(jdx,:);
end

%@test:1
%$ % Define some dates
%$ B1 = '1953Q4';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009Q2';
%$
%$ % Define expected results
%$ e.time = [1945 3; 1950 1; 1950 2; 1953 4; 2009 2];
%$ e.freq = 4;
%$ e.ndat = 5;
%$
%$ % Call the tested routine
%$ d = dynDates(B4,B3,B2,B1);
%$ d = d.append(dynDates(B5));
%$ f = d.pop();
%$ t(1) = dyn_assert(f.time,e.time(1:end-1,:));
%$ t(2) = dyn_assert(f.freq,e.freq);
%$ t(3) = dyn_assert(f.ndat,e.ndat-1);
%$ f = d.pop(B1);
%$ t(4) = dyn_assert(f.time,[1945 3; 1950 1; 1950 2; 2009 2]);
%$ t(5) = dyn_assert(f.freq,e.freq);
%$ t(6) = dyn_assert(f.ndat,e.ndat-1);
%$ f = d.pop(dynDates(B1));
%$ t(7) = dyn_assert(f.time,[1945 3; 1950 1; 1950 2; 2009 2]);
%$ t(8) = dyn_assert(f.freq,e.freq);
%$ t(9) = dyn_assert(f.ndat,e.ndat-1);
%$
%$ % Check the results.
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define some dates
%$ B1 = '1950Q1';
%$ B2 = '1950Q2';
%$ B3 = '1950Q1';
%$ B4 = '1945Q3';
%$ B5 = '2009Q2';
%$
%$ % Call the tested routine
%$ d = dynDates(B1,B2,B3,B4);
%$ d = d.append(dynDates(B5));
%$ f = d.pop();
%$ t(1) = dyn_assert(isequal(f,dynDates(B1,B2,B3,B4)),1);
%$ f = d.pop(B1);
%$ t(2) = dyn_assert(isequal(f,dynDates(B1,B2,B4,B5)),1);
%$ g = f.pop(1);
%$ t(3) = dyn_assert(isequal(g,dynDates(B2,B4,B5)),1);
%$ T = all(t);
%@eof:2