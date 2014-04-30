function i = write_mh_history_file(MetropolisFolder, ModelName, record)

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

% Set base name for the mh-history files.
BaseName = [MetropolisFolder filesep ModelName '_mh_history_'];

% Get the list of all the mh-history files.
mh_history_files = dir([BaseName '*.mat']);

i = length(mh_history_files);

% Save record structure.
save([BaseName num2str(i) '.mat'],'record');

% Add version number.
version = 2;
save([BaseName num2str(i) '.mat'],'version','-append');