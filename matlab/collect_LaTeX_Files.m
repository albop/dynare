function collect_LaTeX_Files(M_)
% function collect_LaTeX_Files(M_);
% Creates TeX-File embedding all eps-loaders created for current mod-file
% 
% Inputs:
%   o M_                    model structure
% 
% Notes: 
%   - The packages loaded enable pdflatex to run
%   - The _dynamic and _static TeX-model files are not included as they are standalone TeX-files

% Copyright (C) 2015 Dynare Team
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

%% Write header
f_name_binder=[M_.fname,'_TeX_binder.TeX'];
fid=fopen(f_name_binder,'w+');
fprintf(fid,'%s \n','\documentclass[12pt]{article}');
fprintf(fid,'%s \n','\usepackage{psfrag}');
fprintf(fid,'%s \n','\usepackage{graphicx}');
fprintf(fid,'%s \n','\usepackage{epstopdf}');
fprintf(fid,'%s \n','\usepackage{longtable}');
fprintf(fid,'%s \n','\begin{document}');

%% Root directory
TeX_Files=dir([M_.fname,'*.TeX']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder) && ...
        ~strcmp(TeX_Files(ii).name,[M_.fname,'_dynamic.tex']) && ...
        ~strcmp(TeX_Files(ii).name,[M_.fname,'_static.tex']) && ...
        ~strcmp(TeX_Files(ii).name,[M_.fname,'_original.tex'])
        fprintf(fid,'%s \n',['\include{',f_name,'}']);    
    end
end

%% Output directory
TeX_Files=dir([M_.dname filesep 'Output' filesep  M_.fname '*.TeX']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/Output' '/',f_name,'}']);    
    end
end

%5 graphs directory
TeX_Files=dir([M_.dname filesep 'graphs' filesep  M_.fname '*.TeX']);
for ii=1:length(TeX_Files)
    [pathstr,f_name,ext] = fileparts(TeX_Files(ii).name);
    if ~strcmp(TeX_Files(ii).name,f_name_binder)
        fprintf(fid,'%s \n',['\include{', M_.dname '/graphs' '/',f_name,'}']);    
    end
end

%% Write footer
fprintf(fid,'%s \n','\end{document}');

fclose(fid);