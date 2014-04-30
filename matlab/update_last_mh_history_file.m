function update_last_mh_history_file(MetropolisFolder, ModelName, record)

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

BaseName = [MetropolisFolder filesep ModelName];

% Get the list of all the mh_history files.
mh_history_files = dir([BaseName '_mh_history_*.mat']);

% Check the existence of mh-files (assuming version 2, ie dynare version greater than 4.3.x).
if isequal(length(mh_history_files),0)
    error(['Estimation::update_mh_file: I cannot find any mh-history file in ' MetropolisFolder '!'])
end

BaseName = [BaseName '_mh_history_'];

save([BaseName num2str(length(mh_history_files)-1) '.mat'],'record','-append');