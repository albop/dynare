function info = load_last_mh_history_file(MetropolisFolder, ModelName)
    
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

% Consistency with older versions of Dynare.
if isequal(length(mh_history_files),0)
    if exist([BaseName '_mh_history.mat'])
        copyfile([BaseName '_mh_history.mat'],[BaseName '_mh_history_0.mat'])
        version = 1;
        save([BaseName '_mh_history_0.mat'],'version','-append');
        mh_history_files = dir([BaseName '_mh_history_*.mat']);
    else
        error(['Estimation::load_mh_file: I cannot find any mh-history file in ' MetropolisFolder '!'])
    end
end

BaseName = [BaseName '_mh_history_'];

% ... And read the last one.
mh_history_info = load([BaseName num2str(length(mh_history_files)-1) '.mat']);

% Solve version issues.
if mh_history_info.version<2
    mh_history_info.record.LastLogPost = mh_history_info.record.LastLogLiK;
    mh_history_info.record.LastSeeds = mh_history_info.record.Seeds;
    mh_history_info.record.InitialSeeds = NaN; % This information is forever lost...
    mh_history_info.record = rmfield(mh_history_info.record,'LastLogLiK');
    mh_history_info.record = rmfield(mh_history_info.record,'Seeds');
    mh_history_info.version = 2;
    record = mh_history_info.record;
    version = mh_history_info.record;
    save([BaseName '_mh_history_0.mat'],'record'); clear record;
    save([BaseName '_mh_history_0.mat'],'version','-append'); clear version;
end

if isequal(nargout,0)
    assignin('caller', 'record', mh_history_info.record);
else
    info = mh_history_info.record;
end