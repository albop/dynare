function [] = portfolios_options(varargin)

    global M_ options_ oo_
    
    portfolios_init;

    % default options
    defaults = struct;
    defaults.solver = 'csolve';    
    defaults.method = 'iterative';
    defaults.dynamic = 'false';
    defaults.only_covariances = 1; 
    defaults.keep_deltax = 0;
    options_.check_solution = true;

    solver_options = optimset();
    solver_options.TolFun = 1e-10;
    solver_options.Display = 'iter';
    solver_options.MaxFunEvals = 10000;
    solver_options.Jacobian = 'off';

    options_.solver_options = solver_options;


    % user-supplied options
    options_.portfolios = defaults;

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

end

