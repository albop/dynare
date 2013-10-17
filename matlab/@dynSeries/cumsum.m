function B = cumsum(varargin) % --*-- Unitary tests --*--

% Overloads matlab's cumsum function for dynSeries objects.
%
% INPUTS 
%  o A     dynSeries object [mandatory].
%  o d     dates object [optional]
%  o v     dynSeries object with one observation [optional]
%
% OUTPUTS 
%  o B     dynSeries object.

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

% Get indices of the columns without NaNs
idx = find(~any(isnan(varargin{1}.data)));

if isempty(idx)
    error('dynSeries::cumsum: All the variables have NaNs. The cumulated sum cannot be computed!')
end

if ~isequal(idx(:),transpose(1:varargin{1}.vobs))
    warning('dynSeries::cumsum: The cumulated sum is not computed for some variables because they have NaNs!')
end

switch nargin
    case 1
      % Initialize the output.
      B = varargin{1};
      % Perform the cumulated sum
      B.data(:,idx) = cumsum(B.data(:,idx));
      % Change the name of the variables
      for i=1:B.vobs
          B.name(i) = {['cumsum(' B.name{i} ')']};
          B.tex(i) = {['\sum_t ' B.tex{i}]};
      end
  case 2
    if isa(varargin{2},'dynSeries')
        if ~isequal(varargin{1}.vobs, varargin{2}.vobs)
            error('dynSeries::cumsum: First and second input arguments must be dynSeries objects with the same number of variables!')
        end
        if ~isequal(varargin{1}.name, varargin{2}.name)
            warning('dynSeries::cumsum: First and second input arguments must be dynSeries objects do not have the same variables!')
        end
        if ~isequal(varargin{2}.nobs,1)
            error('dynSeries::cumsum: Second input argument must be a dynSeries object with only one observation!')
        end
        B = cumsum(varargin{1});
        B.data = bsxfun(@plus,B.data,varargin{2}.data);
    elseif isdates(varargin{2})
        B = cumsum(varargin{1});
        t = find(B.time==varargin{2});
        if isempty(t)
            if varargin{2}==(B.init-1)
                return
            else
                error(['dynSeries::cumsum: date ' date2string(varargin{2}) ' is not in the sample!'])
            end
        end
        B.data = bsxfun(@minus,B.data,B.data(t,:));
    else
        error('dynSeries::cumsum: Second input argument must be a dynSeries object or a dates object!')
    end
  case 3
    if ~isdates(varargin{2})
        error('dynSeries::cumsum: Second input argument must be a dates object!')
    end
    if ~isa(varargin{3},'dynSeries')
        error('dynSeries::cumsum: Third input argument must be a dynSeries object!')
    end
    if ~isequal(varargin{1}.vobs, varargin{3}.vobs)
        error('dynSeries::cumsum: First and third input arguments must be dynSeries objects with the same number of variables!')
    end
    if ~isequal(varargin{1}.name, varargin{3}.name)
        warning('dynSeries::cumsum: First and third input arguments must be dynSeries objects do not have the same variables!')
    end
    if ~isequal(varargin{3}.nobs,1)
        error('dynSeries::cumsum: Third input argument must be a dynSeries object with only one observation!')
    end
    B = cumsum(varargin{1});
    t = find(B.time==varargin{2});
    if isempty(t)
        if varargin{2}==(B.init-1)
            B.data = bsxfun(@plus,B.data,varargin{3}.data);
            return
        else
            error(['dynSeries::cumsum: date ' date2string(varargin{2}) ' is not in the sample!'])
        end
    end
    B.data = bsxfun(@plus,B.data,varargin{3}.data-B.data(t,:));
  otherwise
    error('dynSeries::cumsum: Wrong number of input arguments!')
end

%@test:1
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ ts1 = cumsum(ts1);
%$
%$ % Expected results.
%$ ts2 = dynSeries(transpose(1:10), [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dyn_assert(isequal(ts1,ts2),1);
%$ warning on
%$ T = all(t);
%@eof:1

%@test:2
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$
%$ % Call the tested method.
%$ ts1 = ts1.cumsum();
%$
%$ % Expected results.
%$ ts2 = dynSeries(transpose(1:10), [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dyn_assert(isequal(ts1,ts2),1);
%$ warning on
%$ T = all(t);
%@eof:2

%@test:3
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ 
%$ % Call the tested method.
%$ ts1 = cumsum(ts1,dates('3Y'));
%$
%$ % Expected results.
%$ ts2 = dynSeries([-2; -1; 0; 1; 2; 3; 4; 5; 6; 7], [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dyn_assert(isequal(ts1,ts2),1);
%$ warning on
%$ T = all(t);
%@eof:3

%@test:4
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(pi, [], A_name, []);
%$
%$ % Call the tested method.
%$ ts3 = cumsum(ts1,dates('3Y'),ts2);
%$
%$ % Expected results.
%$ ts4 = dynSeries([-2; -1; 0; 1; 2; 3; 4; 5; 6; 7]+pi, [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dyn_assert(isequal(ts3,ts4),1);
%$ warning on
%$ T = all(t);
%@eof:4

%@test:4
%$ % Define a data set.
%$ A = ones(10,1);
%$
%$ % Define names
%$ A_name = {'A1'};
%$
%$ % Instantiate a time series object.
%$ ts1 = dynSeries(A,[],A_name,[]);
%$ ts2 = dynSeries(pi, [], A_name, []);
%$
%$ % Call the tested method.
%$ ts3 = ts1.cumsum(dates('3Y'),ts2);
%$
%$ % Expected results.
%$ ts4 = dynSeries([-2; -1; 0; 1; 2; 3; 4; 5; 6; 7]+pi, [], A_name, []);
%$
%$ % Check the results.
%$ warning off, % Because the names of the variables are not the same...
%$ t(1) = dyn_assert(isequal(ts3,ts4),1);
%$ warning on
%$ T = all(t);
%@eof:4