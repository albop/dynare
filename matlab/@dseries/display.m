function display(A)
%@info:
%! @deftypefn {Function File} display (@var{A})
%! @anchor{@dseries/display}
%! @sp 1
%! Overloads the disp method for the Dynare time series class (@ref{dseries}).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Dynare time series object instantiated by @ref{dseries}.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! None
%! @end deftypefn
%@eod:

vspace = ' ';
TABLE = ' ';

if A.vobs<=10
    if A.nobs<=40
        separator = repmat(' | ',A.nobs+1,1);
        for t=1:A.nobs
            TABLE = char(TABLE, date2string(A.dates(t)));
        end
        for i = 1:A.vobs
            TABLE = horzcat(TABLE,separator);
            tmp = A.name{i};
            for t=1:A.nobs
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            TABLE = horzcat(TABLE, tmp);
        end
    else
        n = 10;
        separator = repmat(' | ',2*n+3,1);
        for t=1:n
            TABLE = char(TABLE, date2string(A.dates(t)));
        end
        TABLE = char(TABLE,vspace);
        for t = A.nobs-n:A.nobs
            TABLE = char(TABLE, date2string(A.dates(t)));
        end
        for i=1:A.vobs
            TABLE = horzcat(TABLE,separator);
            tmp = A.name{i};
            for t=1:10
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            tmp = char(tmp,vspace);
            for t=A.nobs-10:A.nobs
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            TABLE = horzcat(TABLE, tmp);
        end
    end
else
    m = 4;
    if A.nobs<=40
        separator = repmat(' | ',A.nobs+1,1);
        for t=1:A.nobs
            TABLE = char(TABLE, date2string(A.dates(t)));
        end
        for i = 1:m
            TABLE = horzcat(TABLE,separator);
            tmp = A.name{i};
            for t=1:A.nobs
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            TABLE = horzcat(TABLE, tmp);
        end
        TABLE = horzcat(TABLE, separartor, repmat(' ... ', A.nobs+1,1));
        for i = A.vobs-m+1:A.vobs
            TABLE = horzcat(TABLE,separator);
            tmp = A.name{i};
            for t=1:A.nobs
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            TABLE = horzcat(TABLE, tmp);
        end
    else
        n = 10;
        separator = repmat(' | ',2*n+3,1);
        for t=1:n
            TABLE = char(TABLE, date2string(A.dates(t)));
        end
        TABLE = char(TABLE,vspace);
        for t = A.nobs-n:A.nobs
            TABLE = char(TABLE, date2string(A.dates(t)));
        end
        for i=1:m
            TABLE = horzcat(TABLE,separator);
            tmp = A.name{i};
            for t=1:10
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            tmp = char(tmp,vspace);
            for t=A.nobs-10:A.nobs
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            TABLE = horzcat(TABLE, tmp);
        end
        TABLE = horzcat(TABLE, separator, repmat(' ... ', 2*n+3,1));
        for i=A.vobs-m+1:A.vobs
            TABLE = horzcat(TABLE,separator);
            tmp = A.name{i};
            for t=1:10
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            tmp = char(tmp,vspace);
            for t=A.nobs-10:A.nobs
                tmp = char(tmp,num2str(A.data(t,i)));
            end
            TABLE = horzcat(TABLE, tmp);
        end
    end
end
disp(vspace)
disp([inputname(1) ' is a dseries object:'])
disp(vspace);
if ~isempty(strtrim(TABLE))
    disp(TABLE);
    disp(vspace);
end
