function [] = portfolios(varargin)
    
    % options : set options to proceed with portfolios solving

    global options_    
    
    if nargin > 1

        s = options_.portfolios;
        t = reshape(varargin,2,nargin/2);
        t = t';
        for i=1:size(t,1);
            s = setfield(s,cell2mat(t(i,1)),cell2mat(t(i,2)));
            %s = cell2struct(t(:,2),t(:,1),1)
        end
        options_.portfolios = s;
        %error('Options for command ''portfolios'' not implemented yet')
    end

    if ~(isfield(options_.portfolios,'method'))
        options_.portfolios.method = 'devereux-sutherland';
    end

    if isfield(options_.portfolios,'dynamic') && strcmp(options_.portfolios.dynamic,'true')
        [alphas,dalphas] = solve_for_portfolios;
    else
        alphas = solve_for_portfolios;
    end
end
