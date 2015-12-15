function dyn_latex_table(M_,title,LaTeXtitle,headers,labels,values,label_width,val_width,val_precis,optional_header)

OutputDirectoryName = CheckPath('Output',M_.dname);

%% get width of label column
if ~isempty(label_width)
    label_width = max(size(deblank(char(headers(1,:),labels)),2)+2, ...
        label_width);
else %use default length
    label_width = max(size(deblank(char(headers(1,:),labels)),2))+2;
end
label_format_leftbound  = sprintf('$%%-%ds$',label_width);

%% get width of label column
if all(~isfinite(values))
    values_length = 4;
else
    values_length = max(ceil(max(max(log10(abs(values(isfinite(values))))))),1)+val_precis+1;
end
if any(values) < 0 %add one character for minus sign
    values_length = values_length+1;
end

%% get width of header strings
headers_length = max(size(deblank(headers(2:end,:)),2));
if ~isempty(val_width)
    val_width = max(max(headers_length,values_length)+4,val_width);
else
    val_width = max(headers_length,values_length)+4;
end
value_format  = sprintf('%%%d.%df',val_width,val_precis);
header_string_format  = sprintf('$%%%ds$',val_width);

%Create and print header string
if length(headers) > 0
    header_string = sprintf(label_format_leftbound ,deblank(headers(1,:)));
    header_code_string='l';
    for i=2:size(headers,1)
        header_string  = [header_string '\t & \t ' sprintf(header_string_format,strrep(deblank(headers(i,:)),'\','\\'))];
        header_code_string= [header_code_string 'c'];
    end
end
header_string=[header_string '\\\\\n'];

filename = [OutputDirectoryName '/' M_.fname '_' LaTeXtitle '.TeX'];
fidTeX = fopen(filename,'w');
fprintf(fidTeX,['%% ' datestr(now,0)]);
fprintf(fidTeX,' \n');
fprintf(fidTeX,' \n');
fprintf(fidTeX,'\\begin{center}\n');
fprintf(fidTeX,['\\begin{longtable}{%s} \n'],header_code_string);
fprintf(fidTeX,['\\caption{',title,'}\\\\\n ']);

fprintf(fidTeX,['\\label{Table:',LaTeXtitle,'}\\\\\n']);
fprintf(fidTeX,'\\toprule \n');
if nargin==10
    for ii=1:size(optional_header,1)
        fprintf(fidTeX,'%s\n',optional_header{ii});
    end
end
fprintf(fidTeX,header_string);
fprintf(fidTeX,'\\midrule \\endfirsthead \n');
fprintf(fidTeX,'\\caption{(continued)}\\\\\n ');
fprintf(fidTeX,'\\toprule \\\\ \n');
if nargin==10
    for ii=1:size(optional_header,1)
        fprintf(fidTeX,'%s\n',optional_header{ii});
    end
end
fprintf(fidTeX,header_string);
fprintf(fidTeX,'\\midrule \\endhead \n');
fprintf(fidTeX,['\\bottomrule \\multicolumn{',num2str(size(headers,1)),'}{r}{(Continued on next page)} \\\\ \\hline \\endfoot \n']);
fprintf(fidTeX,'\\bottomrule \\endlastfoot \n');
for i=1:size(values,1)
    fprintf(fidTeX,label_format_leftbound,deblank(labels(i,:)));
    fprintf(fidTeX,['\t & \t' value_format],values(i,:));
    fprintf(fidTeX,' \\\\ \n');
end

fprintf(fidTeX,'\\end{longtable}\n ');
fprintf(fidTeX,'\\end{center}\n');
fprintf(fidTeX,'%% End of TeX file.\n');
fclose(fidTeX);