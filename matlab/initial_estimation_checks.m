function initial_estimation_checks(xparam1,gend,data)
global dr1_test_ bayestopt_ estim_params_ exo_nbr
  
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

fval = DsgeLikelihood(xparam1,gend,data);

if exist('dr1_test_') & dr1_test_(1)>0
  % disp(dr1_test_)
  switch(dr1_test_(1))
   case 1
    error('The steady state can''t be found');
   case 2
    error(['Estimation can''t take place because there are an infinity of' ...
	   ' stable solutions']);
   case 3
    error(['Estimation can''t take place because there is no stable' ...
	   ' solution']);
   case 4
    error(['Estimation can''t take place because of singularity in Kalman' ...
	   ' filter']);
   otherwise
  end
end