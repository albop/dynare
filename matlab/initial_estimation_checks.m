function initial_estimation_checks(xparam1,gend,data)
global bayestopt_ estim_params_ exo_nbr
  
nv = size(data,1);
  
if nv > exo_nbr + estim_params_.nvn
  error(['Estimation can''t take place because there are less shocks than' ...
	 'observed variables'])
end
  
r = rank(data);
if r < nv
  error(['Estimation can''t take place because the data are perfectly' ...
	 ' correlated']);
end

[fval,cost_flag,ys,trend_coeff,info] = DsgeLikelihood(xparam1,gend,data);

if info(1) > 0
  disp('Error in computing likelihood for initial parameter values')
  print_info(info)
end
