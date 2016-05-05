function [res,fjac,domer] = mcp_func(x,jacflag)
global mcp_data

domer = 0;
if jacflag
    [res,fjac] = mcp_data.func(x,mcp_data.args{:});
else
    res = mcp_data.func(x,mcp_data.args{:});
    fjac = [];
end
disp(norm(res))