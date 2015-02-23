function [opt_par_values,fval,exitflag,hessian_mat,options_,Scale]=dynare_minimize_objective(objective_function,start_par_value,minimizer_algorithm,options_,bounds,parameter_names,prior_information,Initial_Hessian,varargin)

% Calls a minimizer
%
% INPUTS
%   objective_function  [function handle]                   handle to the objective function
%   start_par_value     [n_params by 1] vector of doubles   starting values for the parameters
%   minimizer_algorithm [scalar double]                     code of the optimizer algorithm
%   options_            [matlab structure]                  Dynare options structure
%   bounds              [n_params by 2] vector of doubles   2 row vectors containing lower and upper bound for parameters
%   parameter_names     [n_params by 1] cell array          strings containing the parameters names   
%   prior_information   [matlab structure]                  Dynare prior information structure (bayestopt_) provided for algorithm 6
%   Initial_Hessian     [n_params by n_params] matrix       initial hessian matrix provided for algorithm 6
%   varargin            [cell array]                        Input arguments for objective function
%    
% OUTPUTS
%   opt_par_values      [n_params by 1] vector of doubles   optimal parameter values minimizing the objective
%   fval                [scalar double]                     value of the objective function at the minimum
%   exitflag            [scalar double]                     return code of the respective optimizer
%   hessian_mat         [n_params by n_params] matrix       hessian matrix at the mode returned by optimizer
%   options_            [matlab structure]                  Dynare options structure (to return options set by algorithms 5)
%   Scale               [scalar double]                     scaling parameter returned by algorith 6
%    
% SPECIAL REQUIREMENTS
%   none.
%  
% 
% Copyright (C) 2014-2015 Dynare Team
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


%% set bounds and parameter names if not already set
n_params=size(start_par_value,1);
if isempty(bounds)
    if minimizer_algorithm==10 
        error('Algorithm 10 (simpsa) requires upper and lower bounds')
    else
        bounds=[-Inf(n_params,1) Inf(n_params,1)];
    end
end

if minimizer_algorithm==10 && any(any(isinf(bounds)))
    error('Algorithm 10 (simpsa) requires finite upper and lower bounds')
end

if isempty(parameter_names)
    parameter_names=[repmat('parameter ',n_params,1),num2str((1:n_params)')];
end

%% initialize function outputs
hessian_mat=[];
Scale=[];
exitflag=1;
fval=NaN;
opt_par_values=NaN(size(start_par_value));


switch minimizer_algorithm
  case 1
    if isoctave
        error('Optimization algorithm 1 is not available under Octave')
    elseif ~user_has_matlab_license('optimization_toolbox')
        error('Optimization algorithm 1 requires the Optimization Toolbox')
    end
    % Set default optimization options for fmincon.
    optim_options = optimset('display','iter', 'LargeScale','off', 'MaxFunEvals',100000, 'TolFun',1e-8, 'TolX',1e-6);
    if isfield(options_,'optim_opt')
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if options_.analytic_derivation,
        optim_options = optimset(optim_options,'GradObj','on','TolX',1e-7);
    end
    [opt_par_values,fval,exitflag,output,lamdba,grad,hessian_mat] = ...
        fmincon(objective_function,start_par_value,[],[],[],[],bounds(:,1),bounds(:,2),[],optim_options,varargin{:});
  case 2
    error('Optimization algorithm 1 option (Lester Ingber''s Adaptive Simulated Annealing) is no longer available')
  case 3
    if isoctave && ~user_has_octave_forge_package('optim')
        error('Optimization algorithm 3 requires the optim package')
    elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
        error('Optimization algorithm 3 requires the Optimization Toolbox')
    end
    % Set default optimization options for fminunc.
    optim_options = optimset('display','iter','MaxFunEvals',100000,'TolFun',1e-8,'TolX',1e-6);
    if isfield(options_,'optim_opt')
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    if options_.analytic_derivation,
        optim_options = optimset(optim_options,'GradObj','on');
    end
    if ~isoctave
        [opt_par_values,fval,exitflag] = fminunc(objective_function,start_par_value,optim_options,varargin{:});
    else
        % Under Octave, use a wrapper, since fminunc() does not have a 4th arg
        func = @(x) objective_function(x,varargin{:});
        [opt_par_values,fval,exitflag] = fminunc(func,start_par_value,optim_options);
    end
  case 4
    % Set default options.
    H0 = 1e-4*eye(n_params);
    crit = options_.csminwel.tolerance.f;
    nit = options_.csminwel.maxiter;
    numgrad = options_.gradient_method;
    epsilon = options_.gradient_epsilon;
    % Change some options.
    if isfield(options_,'optim_opt')
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                nit = options_list{i,2};
              case 'InitialInverseHessian'
                H0 = eval(options_list{i,2});
              case 'TolFun'
                crit = options_list{i,2};
              case 'NumgradAlgorithm'
                numgrad = options_list{i,2};
              case 'NumgradEpsilon'
                epsilon = options_list{i,2};
              otherwise
                warning(['csminwel: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    % Set flag for analytical gradient.
    if options_.analytic_derivation
        analytic_grad=1;
    else
        analytic_grad=[];
    end
    % Call csminwell.
    [fval,opt_par_values,grad,hessian_mat,itct,fcount,exitflag] = ...
        csminwel1(objective_function, start_par_value, H0, analytic_grad, crit, nit, numgrad, epsilon, varargin{:});
  case 5
    if options_.analytic_derivation==-1 %set outside as code for use of analytic derivation
        analytic_grad=1;
        crit = options_.newrat.tolerance.f_analytic;
        newratflag = 0; %analytical Hessian
    else
        analytic_grad=0;
        crit=options_.newrat.tolerance.f;
        newratflag = options_.newrat.hess; %default
    end
    nit=options_.newrat.maxiter;
    
    if isfield(options_,'optim_opt')
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                nit = options_list{i,2};
              case 'Hessian'
                flag=options_list{i,2};
                if options_.analytic_derivation && flag~=0
                    error('newrat: analytic_derivation is incompatible with numerical Hessian. Using analytic Hessian')
                else
                    newratflag=flag; 
                end
              case 'TolFun'
                crit = options_list{i,2};
              otherwise
                warning(['newrat: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    [opt_par_values,hessian_mat,gg,fval,invhess] = newrat(objective_function,start_par_value,analytic_grad,crit,nit,0,varargin{:});
    if options_.analytic_derivation %Hessian is already analytic one, reset option
        options_.analytic_derivation = ana_deriv;
    elseif ~options_.analytic_derivation && newratflag ==0 %Analytic Hessian wanted, but not computed yet
            hessian_mat = reshape(mr_hessian(0,opt_par_values,objective_function,1,crit,varargin{:}), n_params, n_params);
    end
  case 6
    % Set default options
    gmhmaxlikOptions = options_.gmhmaxlik;
    if ~isempty(Initial_Hessian);
        gmhmaxlikOptions.varinit = 'previous';
    else
        gmhmaxlikOptions.varinit = 'prior';
    end
    if isfield(options_,'optim_opt')
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'NumberOfMh'
                gmhmaxlikOptions.iterations = options_list{i,2};
              case 'ncov-mh'
                gmhmaxlikOptions.number = options_list{i,2};
              case 'nscale'
                gmhmaxlikOptions.nscale = options_list{i,2};
              case 'nclimb'
                gmhmaxlikOptions.nclimb = options_list{i,2};
              case 'InitialCovarianceMatrix'
                switch options_list{i,2}
                  case 'previous'
                    if isempty(Initial_Hessian)
                        error('gmhmaxlik: No previous estimate of the Hessian matrix available! You cannot use the InitialCovarianceMatrix option!')
                    else
                        gmhmaxlikOptions.varinit = 'previous';
                    end
                  case {'prior', 'identity'}
                    gmhmaxlikOptions.varinit = options_list{i,2};
                  otherwise
                    error('gmhmaxlik: Unknown value for option ''InitialCovarianceMatrix''!')
                end
              case 'AcceptanceRateTarget'
                gmhmaxlikOptions.target = options_list{i,2};
                if gmhmaxlikOptions.target>1 || gmhmaxlikOptions.target<eps
                    error('gmhmaxlik: The value of option AcceptanceRateTarget should be a double between 0 and 1!')
                end
              otherwise
                warning(['gmhmaxlik: Unknown option (' options_list{i,1}  ')!'])
            end
        end
    end
    % Evaluate the objective function.
    fval = feval(objective_function,start_par_value,varargin{:});
    OldMode = fval;
    if ~exist('MeanPar','var')
        MeanPar = start_par_value;
    end
    switch gmhmaxlikOptions.varinit
      case 'previous'
        CovJump = inv(Initial_Hessian);
      case 'prior'
        % The covariance matrix is initialized with the prior
        % covariance (a diagonal matrix) %%Except for infinite variances ;-)
        stdev = prior_information.p2;
        indx = find(isinf(stdev));
        stdev(indx) = ones(length(indx),1)*sqrt(10);
        vars = stdev.^2;
        CovJump = diag(vars);
      case 'identity'
        vars = ones(length(prior_information.p2),1)*0.1;
        CovJump = diag(vars);
      otherwise
        error('gmhmaxlik: This is a bug! Please contact the developers.')
    end
    OldPostVar = CovJump;
    Scale = options_.mh_jscale;
    for i=1:gmhmaxlikOptions.iterations
        if i == 1
            if gmhmaxlikOptions.iterations>1
                flag = '';
            else
                flag = 'LastCall';
            end
            [opt_par_values,PostVar,Scale,PostMean] = ...
                gmhmaxlik(objective_function,start_par_value,bounds,gmhmaxlikOptions,Scale,flag,MeanPar,CovJump,varargin{:});
            fval = feval(objective_function,opt_par_values,varargin{:});
            mouvement = max(max(abs(PostVar-OldPostVar)));
            skipline()
            disp('========================================================== ')
            disp(['   Change in the covariance matrix = ' num2str(mouvement) '.'])
            disp(['   Mode improvement = ' num2str(abs(OldMode-fval))])
            disp(['   New value of jscale = ' num2str(Scale)])
            disp('========================================================== ')
            OldMode = fval;
        else
            OldPostVar = PostVar;
            if i<gmhmaxlikOptions.iterations
                flag = '';
            else
                flag = 'LastCall';
            end
            [opt_par_values,PostVar,Scale,PostMean] = ...
                gmhmaxlik(objective_function,opt_par_values,[bounds bounds],...
                          gmhmaxlikOptions,Scale,flag,PostMean,PostVar,varargin{:});
            mouvement = max(max(abs(PostVar-OldPostVar)));
            fval = feval(objective_function,opt_par_values,varargin{:});
            skipline()
            disp('========================================================== ')
            disp(['   Change in the covariance matrix = ' num2str(mouvement) '.'])
            disp(['   Mode improvement = ' num2str(abs(OldMode-fval))])
            disp(['   New value of jscale = ' num2str(Scale)])
            disp('========================================================== ')
            OldMode = fval;
        end
        hessian_mat = inv(PostVar);
    end
    skipline()
    disp(['Optimal value of the scale parameter = ' num2str(Scale)])
    skipline()
  case 7
    % Matlab's simplex (Optimization toolbox needed).
    if isoctave && ~user_has_octave_forge_package('optim')
        error('Option mode_compute=7 requires the optim package')
    elseif ~isoctave && ~user_has_matlab_license('optimization_toolbox')
        error('Option mode_compute=7 requires the Optimization Toolbox')
    end
    optim_options = optimset('display','iter','MaxFunEvals',1000000,'MaxIter',6000,'TolFun',1e-8,'TolX',1e-6);
    if isfield(options_,'optim_opt')
        eval(['optim_options = optimset(optim_options,' options_.optim_opt ');']);
    end
    [opt_par_values,fval,exitflag] = fminsearch(objective_function,start_par_value,optim_options,varargin{:});
  case 8
    % Dynare implementation of the simplex algorithm.
    simplexOptions = options_.simplex;
    if isfield(options_,'optim_opt')
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                simplexOptions.maxiter = options_list{i,2};
              case 'TolFun'
                simplexOptions.tolerance.f = options_list{i,2};
              case 'TolX'
                simplexOptions.tolerance.x = options_list{i,2};
              case 'MaxFunEvals'
                simplexOptions.maxfcall = options_list{i,2};
              case 'MaxFunEvalFactor'
                simplexOptions.maxfcallfactor = options_list{i,2};
              case 'InitialSimplexSize'
                simplexOptions.delta_factor = options_list{i,2};
              otherwise
                warning(['simplex: Unknown option (' options_list{i,1} ')!'])
            end
        end
    end
    [opt_par_values,fval,exitflag] = simplex_optimization_routine(objective_function,start_par_value,simplexOptions,parameter_names,varargin{:});
  case 9
    % Set defaults
    H0 = 1e-4*ones(n_params,1);
    cmaesOptions = options_.cmaes;
    % Modify defaults
    if isfield(options_,'optim_opt')
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                cmaesOptions.MaxIter = options_list{i,2};
              case 'TolFun'
                cmaesOptions.TolFun = options_list{i,2};
              case 'TolX'
                cmaesOptions.TolX = options_list{i,2};
              case 'MaxFunEvals'
                cmaesOptions.MaxFunEvals = options_list{i,2};
              otherwise
                warning(['cmaes: Unknown option (' options_list{i,1}  ')!'])
            end
        end
    end
    warning('off','CMAES:NonfinitenessRange');
    warning('off','CMAES:InitialSigma');
    [x, fval, COUNTEVAL, STOPFLAG, OUT, BESTEVER] = cmaes(func2str(objective_function),start_par_value,H0,cmaesOptions,varargin{:});
    opt_par_values=BESTEVER.x;
    fval=BESTEVER.f;
  case 10
    simpsaOptions = options_.simpsa;
    if isfield(options_,'optim_opt')
        options_list = read_key_value_string(options_.optim_opt);
        for i=1:rows(options_list)
            switch options_list{i,1}
              case 'MaxIter'
                simpsaOptions.MAX_ITER_TOTAL = options_list{i,2};
              case 'TolFun'
                simpsaOptions.TOLFUN = options_list{i,2};
              case 'TolX'
                tolx = options_list{i,2};
                if tolx<0
                    simpsaOptions = rmfield(simpsaOptions,'TOLX'); % Let simpsa choose the default.
                else
                    simpsaOptions.TOLX = tolx;
                end
              case 'EndTemparature'
                simpsaOptions.TEMP_END = options_list{i,2};
              case 'MaxFunEvals'
                simpsaOptions.MAX_FUN_EVALS = options_list{i,2};
              otherwise
                warning(['simpsa: Unknown option (' options_list{i,1}  ')!'])
            end
        end
    end
    simpsaOptionsList = options2cell(simpsaOptions);
    simpsaOptions = simpsaset(simpsaOptionsList{:});
    [opt_par_values, fval, exitflag] = simpsa(func2str(objective_function),start_par_value,bounds(:,1),bounds(:,2),simpsaOptions,varargin{:});
  case 11
     options_.cova_compute = 0 ;
     [opt_par_values,stdh,lb_95,ub_95,med_param] = online_auxiliary_filter(start_par_value,varargin{:}) ;
  case 101
    myoptions=soptions;
    myoptions(2)=1e-6; %accuracy of argument
    myoptions(3)=1e-6; %accuracy of function (see Solvopt p.29)
    myoptions(5)= 1.0;
    [opt_par_values,fval]=solvopt(start_par_value,objective_function,[],myoptions,varargin{:});
  case 102
    %simulating annealing
    LB = start_par_value - 1;
    UB = start_par_value + 1;
    neps=10;
    %  Set input parameters.
    maxy=0;
    epsilon=1.0e-9;
    rt_=.10;
    t=15.0;
    ns=10;
    nt=10;
    maxevl=100000000;
    idisp =1;
    npar=length(start_par_value);

    disp(['size of param',num2str(length(start_par_value))])
    c=.1*ones(npar,1);
    %*  Set input values of the input/output parameters.*

    vm=1*ones(npar,1);
    disp(['number of parameters= ' num2str(npar) 'max= '  num2str(maxy) 't=  ' num2str(t)]);
    disp(['rt_=  '  num2str(rt_) 'eps=  '  num2str(epsilon) 'ns=  '  num2str(ns)]);
    disp(['nt=  '  num2str(nt) 'neps= '   num2str(neps) 'maxevl=  '  num2str(maxevl)]);
    disp '  ';
    disp '  ';
    disp(['starting values(x)  ' num2str(start_par_value')]);
    disp(['initial step length(vm)  '  num2str(vm')]);
    disp(['lower bound(lb)', 'initial conditions', 'upper bound(ub)' ]);
    disp([LB start_par_value UB]);
    disp(['c vector   ' num2str(c')]);

    [opt_par_values, fval, nacc, nfcnev, nobds, ier, t, vm] = sa(objective_function,xparstart_par_valueam1,maxy,rt_,epsilon,ns,nt ...
                                                          ,neps,maxevl,LB,UB,c,idisp ,t,vm,varargin{:});
  otherwise
    if ischar(minimizer_algorithm)
        if exist(options_.mode_compute)
            [opt_par_values, fval, exitflag] = feval(options_.mode_compute,objective_function,start_par_value,varargin{:});
        else
            error('No minimizer with the provided name detected.')
        end
    else
        error(['Optimization algorithm ' int2str(minimizer_algorithm) ' is unknown!'])
    end
end
