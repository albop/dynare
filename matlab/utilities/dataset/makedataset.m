function [DynareDataset, DatasetInfo] = makedataset(DynareOptions)

% Initialize a dataset as a dseries object.
%
%
% INPUTS
% ======
%
%     DynareOptions [struct] Structure of options built by Dynare's preprocessor.
%
%
% OUTPUTS
% =======
%
%     DynareDataset [dseries]  The dataset.
%     DatasetInfo   [struct]   Various informations about the dataset (descriptive statistics and missing observations).
%
% EXAMPLE
% =======
%
%     [dataset_, dataset_info] = makedataset(options_) ;
%
%
% See also dynare_estimation_init

if isempty(DynareOptions.datafile) && isempty(DynareOptions.dataset.file)
    if gsa_flag
        DynareDataset = dseries();
        return
    else
        error('datafile option is missing')
    end
end

if isempty(DynareOptions.datafile) && ~isempty(DynareOptions.dataset.file)
    datafile = DynareOptions.dataset.file;
    newdatainterface = 1;
elseif ~isempty(DynareOptions.datafile) && isempty(DynareOptions.dataset.file)
    datafile = DynareOptions.datafile;
    newdatainterface = 0;
elseif isempty(DynareOptions.datafile) && ~isempty(DynareOptions.dataset.file)
    error('You cannot use simultaneously the data command and the datafile option (in the estimation command)!')
else
    error('You have to specify the datafile!')
end

% Check extension.
allowed_extensions = {'m','mat','csv','xls','xlsx'};
datafile_extension = get_file_extension(datafile);
if isempty(datafile_extension)
    available_extensions = {}; j = 1;
    for i=1:length(allowed_extensions)
        if exist([datafile '.' allowed_extensions{i}])
            available_extensions(j) = {allowed_extensions{i}};
            j = j+1;
        end
    end
    if isempty(available_extensions)
        error(['I can''t find a datafile (with allowed extension)!'])
    end
    if length(available_extensions)>1
        error(sprintf(['You did not specify an extension for the datafile, but more than one candidate ' ...
                       'are available in the designed folder!\nPlease, add an extension to the datafile ' ...
                       '(m, mat, csv, xls or xlsx are legal extensions).']));
    end
    datafile = [datafile '.' available_extensions{1}];
end

% Load the data in a dseries object.
DynareDataset = dseries(datafile);

% Select a subset of the variables.
DynareDataset = DynareDataset{DynareOptions.varobs{:}};

% Apply log function if needed.
if DynareOptions.loglinear && ~DynareOptions.logdata
    DynareDataset = DynareDataset.log();
end

% Test if an initial period (different from its default value) is explicitely defined in the datafile.
if isequal(DynareDataset.init, dates(1,1))
    dataset_default_initial_period = 1;
else
    dataset_default_initial_period = 0;
end

%  Test if an initial period (different from its default value) is explicitely defined in the mod file with the set_time command.
if isequal(DynareOptions.initial_period, dates(1,1))
    set_time_default_initial_period = 1;
else
    set_time_default_initial_period = 0;
end

if ~set_time_default_initial_period && dataset_default_initial_period
    % Overwrite the initial period in dataset (it was set to default).
    % Note that the updates of freq and time members are auto-magically
    % done by dseries::subsasgn overloaded method.
    DynareDataset.init = DynareOptions.initial_period;
end

if set_time_default_initial_period && ~dataset_default_initial_period
    % Overwrite the global initial period defined by set_time (it was set to default).
    DynareOptions.initial_period = DynareDataset.init;
end

if ~set_time_default_initial_period && ~dataset_default_initial_period
    % Check if dataset.init and options_.initial_period are identical.
    if  ~isequal(DynareOptions.initial_period, DynareDataset.init)
        error('The date as defined by the set_time command is not consistent with the initial period in the database!')
    end
end

% Set firstobs, lastobs and nobs
if newdatainterface
    if isempty(DynareOptions.dataset.firstobs)
        % first_obs option was not used in the data command.
        firstobs = DynareDataset.init;
    else
        firstobs = DynareOptions.dataset.firstobs;
    end
    if isnan(DynareOptions.dataset.nobs)
        % nobs option was not used in the data command.
        if isempty(DynareOptions.dataset.lastobs)
            % last_obs option was not used in the data command.
            nobs = DynareDataset.nobs;
            lastobs = DynareDataset.dates(end);
        else
            lastobs = DynareOptions.dataset.lastobs;
            nobs = lastobs-firstobs+1;
        end
    else
        nobs = DynareOptions.dataset.nobs;
        if isempty(DynareOptions.dataset.lastobs)
            % last_obs option was not used in the data command.
            lastobs = firstobs+(nobs-1);
        else
            % last_obs and nobs were used in the data command. Check that they are consistent (with firstobs).
            if ~isequal(lastobs,firstobs+(nobs-1))
                error(sprintf('Options last_obs (%s), first_obs (%s) and nobs (%s) are not consistent!',char(lastobs),char(firstobs),num2str(nobs)));
            end
        end
    end
else
    if isnan(DynareOptions.first_obs)
        firstobs = DynareDataset.init;
    else
        firstobs = DynareDataset.dates(DynareOptions.first_obs);
    end
    if isnan(DynareOptions.nobs)
        lastobs = DynareDataset.dates(end);
        nobs = lastobs-firstobs+1;
    else
        nobs = DynareOptions.nobs;
        lastobs = firstobs+(nobs-1);
    end
end

% Check that firstobs belongs to DynareDataset.dates
if firstobs<DynareDataset.init
    error(sprintf('first_obs (%s) cannot be less than the first date in the dataset (%s)!',char(firstobs),char(DynareDataset.init)))
end

% Check that lastobs belongs to DynareDataset.dates...
if newdatainterface
    if lastobs>DynareDataset.dates(end)
        error(sprintf('last_obs (%s) cannot be greater than the last date in the dataset (%s)!',char(lastobs),char(DynareDataset.dates(end))))
    end
else
    % ...  or check that nobs is smaller than the number of observations in DynareDataset.
    if nobs>DynareDataset.nobs
        error(sprintf('nobs (%s) cannot be greater than the last date in the dataset (%s)!', num2str(nobs), num2str(nobs)))
    end
end

% Select a subsample.
DynareDataset = DynareDataset(firstobs:lastobs);

% Initialize DatasetInfo structure.
DatasetInfo = struct('missing', struct('state', NaN, 'aindex', [], 'vindex', [], 'number_of_observations', NaN, 'no_more_missing_observations', NaN), ...
                     'descriptive', struct('mean', [], 'covariance', [], 'correlation', [], 'autocovariance', []));

% Fill DatasetInfo.missing if some observations are missing
DatasetInfo.missing.state = isanynan(DynareDataset.data);
if DatasetInfo.missing.state
    [DatasetInfo.missing.aindex, DatasetInfo.missing.number_of_observations, DatasetInfo.missing.no_more_missing_observations, DatasetInfo.missing.vindex] = ...
        describe_missing_data(DynareDataset.data);
else
    DatasetInfo.missing.aindex = num2cell(transpose(repmat(1:DynareDataset.vobs,DynareDataset.nobs,1)),1);
    DatasetInfo.missing.no_more_missing_observations = 1;
end

% Compute the empirical mean of the observed variables.
DatasetInfo.descriptive.mean = nanmean(DynareDataset.data);

% Compute the empirical covariance matrix of the observed variables.
DatasetInfo.descriptive.covariance = nancovariance(DynareDataset.data);

% Compute the empirical correlation matrix of the observed variables.
normalization_matrix = diag(1./sqrt(diag(DatasetInfo.descriptive.covariance)));
DatasetInfo.descriptive.correlation = normalization_matrix*DatasetInfo.descriptive.covariance*normalization_matrix;

% Compute autocorrelation function.
DatasetInfo.descriptive.autocovariance = nanautocovariance(DynareDataset.data, DynareOptions.ar);