function [res,fjac,domer] = mcp_func(x,jacflag)
global mcp_data

if jacflag
    [res,fjac] = mcp_data.func(x,mcp_data.args{:});
    fjac = sparse(fjac);
else
    res = mcp_data.func(x,mcp_data.args{:});
    fjac = [];
end
if isreal(res)
    domer = 0;
else
    domer = 1;
end
