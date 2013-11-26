function [freq, init, data, varlist] = load_csv_file_data(file, withtime, withnames, noemptycell)

% Loads data in a csv file.
%
% INPUTS 
%  o file        string, name of the csv file (with path).
%  o withtime    integer scalar, nonzero iff the first column is for the dates of the observations.
%  o withnames   integer scalar, nonzero iff the first row is for the names of the variables.
%  o noemptycell integer scalar, nonzero the csv file does not have empty cells.
%
% OUTPUTS 
%  o freq        integer scalar equal to 1, 4, 12 or 52 (for annual, quaterly, monthly or weekly frequencies).
%  o init        dates object, initial date in the dataset.
%  o data        matrix of doubles, the data.
%  o varlist     cell of strings, names of the variables.

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

% Set defaults.
if nargin<4
    noemptycell = 1;
    if nargin<3
        withnames = 1;
        if nargin <2
            withtime = 1;
            if nargin<1
                error('load_csv_file_data:: I need at least one input (name of the csv file)!')
            end
        end
    end
end

if ~withtime && ~withnames && noemptycell
    % Use matlab builtin routine!
    data = csvread(file);
end

if ~( isequal(withtime,0) || isequal(withtime,1) )
    error('load_csv_file_data:: Second input argument has to be equal to 1 or 0!')
end

if ~( isequal(withnames,0) || isequal(withnames,1) )
    error('load_csv_file_data:: Third input argument has to be equal to 1 or 0!')
end

% Output initialization
freq = 1;                  % Default frequency is annual.
init = dates(1,1);         % Default initial date is year one.
varlist = [];
if ~isoctave
    % Under Matlab, save time by using importdata
    assert(exist(file, 'file') == 2, ['load_csv_file_data: I can''t find file ' file '!']);
    A = importdata(file, ',', withnames);
    if withnames && withtime
        if size(A.textdata, 1) == 1
            % year dates confused for data
            varlist = A.textdata(1, 2:end);
            init = dates([num2str(A.data(1, 1)) 'Y']);
            data = A.data(:, 2:end);
        else
            varlist = A.textdata(1, 2:end);
            init = dates(A.textdata{2, 1});
            data = A.data;
        end
    elseif withnames && ~withtime
        varlist = A.textdata;
        data = A.data;
    elseif ~withnames && withtime
        if ~isstruct(A)
            % year dates confused for data
            init = dates([num2str(A.data(1, 1)) 'Y']);
            data = A(:, 2:end);
        else
            init = dates(A.textdata{1, 1});
            data = A.data;
        end
    else
        error('load_csv_file_data:: Shouldn''t arrive here');
    end
    freq = init.freq;
    varlist = transpose(varlist);
    return
end

% Check if file exists.
if check_file_extension(file,'csv')
    try
        fid = fopen(file,'r');
    catch
        error(['load_csv_file_data: I can''t find file ' file '!'])
    end
else
    error('load_csv_file_data: Wrong file extension!')
end

% bfile contains a vector of ascii codes.
bfile = fread(fid);

% Close (csv) file. 
fclose(fid);

% Set newline code (ok for *nix, check for mac and windows)
if isunix
    newline_code = 10;
elseif ispc
    newline_code = 13;
elseif ismac
    newline_code = 10;
else
    error('load_csv_file_data:: Not implemented for your OS!')
end

% Get the positions of the end-of-line code;
end_of_line_locations = find(bfile==newline_code);
if ispc && isempty(end_of_line_locations)
    newline_code=10;
    end_of_line_locations = find(bfile==newline_code);
end;
tmp = find(bfile==newline_code);

% Get the number of lines in the file.
ndx = length(tmp);                       

% Create a cell of indices for each line.
b = [1; end_of_line_locations+1];
c = [end_of_line_locations-1; length(bfile)+1];
b = b(1:end-1);
c = c(1:end-1);

linea = 1;

if withnames
    % Get the first line of the csv file (names of the variables).
    linee = char(transpose(bfile(b(linea):c(linea))));
    % Get the content of the first line and determine the number of variables and their names.
    [B,C] = get_cells_id(linee,',');
    if withtime
        B = B(2:end);
        C = C(2:end);
    end
    varlist = cell(length(B),1);
    number_of_variables = length(varlist);
    for i=1:number_of_variables
        varlist(i) = {linee(B(i):C(i))};
    end
    varlist = strtrim(varlist);
    linea = linea+1;
end

% Get following line (number 1 or 2 depending on withnames flag)
linee = char(transpose(bfile(b(linea):c(linea))));
comma_locations = transpose(strfind(linee,','));
B = 1;
C = comma_locations(1)-1;

if withtime
    tmp = linee(B:C);
    % Check the dates formatting
    if isnumeric(tmp) && isint(tmp)
        tmp = [num2str(tmp) 'Y'];
    end
    if ~isdate(tmp)
        error('load_csv_file_data:: Formatting error. I can''t read the dates!')
    end
    init = dates(tmp);
    freq = init.freq;
    first = 2;
else
    first = 1;
end

if ~withnames
    number_of_variables = length(tmp)-withtime;
end

% Initialization of matrix data.
data = zeros(ndx,number_of_variables);

% Populate data.
for linea = 1+withnames:ndx
    linee = char(transpose(bfile(b(linea):c(linea))));
    [B,C] = get_cells_id(linee,',');
    for i=first:length(B)
        data(linea,i-withtime) = str2double(linee(B(i):C(i)));
    end
end

% Remove first line if withnames
data = data(1+withnames:ndx,:);