function PosteriorOddsTable = model_comparison(ModelNames,ModelPriors)
% 05-30-2005
%
% type is a string  = Laplace
%                   = ModifiedHarmonicMean
% ModelPriors is a m*1 column vector
% ModelNames is m*1 cell array

global oo_ options_ fname_

tmp_oo = oo_;
if isfield(options_,'model_comparison_approximation')
  type = options_.model_comparison_approximation;
else
  type = 'LaplaceApproximation';
end

if strcmp(type,'Laplace')
  type = 'LaplaceApproximation';
end

NumberOfModels = size(ModelNames,1);
MarginalLogDensity = zeros(NumberOfModels,1);
empty_prior = 1;
if isempty(ModelPriors)
  empty_prior = 0;
  ModelPriors = ones(NumberOfModels,1)/NumberOfModels;
end

if abs(sum(ModelPriors)-1) > 1e-6
  disp('MODEL_COMPARISON: model priors renormalized so as to sum to 1 !');
  ModelPriors = ModelPriors/sum(ModelPriors);
end

% Get the estimates of the log marginal densities

for i = 1:NumberOfModels
  if strcmp(ModelNames{i},fname_)
    oo_ = tmp_oo;
  else
    load([ModelNames{i} '_results.mat' ],'oo_');
  end
  try
    eval(['MarginalLogDensity(i) = oo_.MarginalDensity.' type ';']) 
  catch
    if strcmpi(type,'LaplaceApproximation')
      disp(['MODEL_COMPARISON: I cant''t find the Laplace approximation associated to model ' ModelNames{i}])
      oo_ = tmp_oo;
      return
    elseif strcmpi(type,'ModifiedHarmonicMean')
      disp(['MODEL_COMPARISON: I cant''t find the modified harmonic mean' ...
	    ' estimate associated to model ' ModelNames{i}])
      oo_ = tmp_oo;
      return
    end
  end
end

%in order to avoid overflow, we divide the numerator and the denominator
%of th Posterior Odds Ratio by the largest Marginal Posterior Density
lmpd = log(ModelPriors)+MarginalLogDensity;
[maxval,k] = max(lmpd);
elmpd = exp(lmpd-maxval);

% Now I display the posterior probabilities
title = 'Model Comparison'; 
headers = strvcat('Model',ModelNames{:});
if empty_prior
  labels = strvcat('Priors','Log Marginal Density','Bayes Ratio', ...
		   'Posterior Model Probability');
  values = [ModelPriors';MarginalLogDensity';exp(lmpd-lmpd(1))'; ...
	    elmpd'/sum(elmpd)];
else
  labels = strvcat('Priors','Log Marginal Density','Bayes Ratio','Posterior Odds Ratio', ...
		   'Posterior Model Probability');
  values = [ModelPriors';MarginalLogDensity'; exp(MarginalLogDensity-MarginalLogDensity(1))'; ...
	    exp(lmpd-lmpd(1))'; elmpd'/sum(elmpd)];
end
  
table(title,headers,labels,values, 0, 15, 6)

oo_ = tmp_oo;