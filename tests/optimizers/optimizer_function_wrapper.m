function [opt_par_values,fval,exitflag]=optimizer_function_wrapper(objective_function_handle,start_par_value,varargin)
% function [opt_par_values,fval,exitflag]=optimizer_function_wrapper(objective_function_handle,start_par_value,varargin)
% Demonstrates how to invoke external optimizer for mode_computation

%set options of optimizer
H0 = 1e-4*eye(length(start_par_value),length(start_par_value));
nit=1000;
crit = 1e-7;
numgrad = 2; 
epsilon = 1e-6; 
analytic_grad=[];

%call optimizer
[fval,opt_par_values,grad,hessian_mat,itct,fcount,exitflag] = ...
    csminwel1(objective_function_handle, start_par_value, H0, analytic_grad, crit, nit, numgrad, epsilon, varargin{:});
end