function [forecast,info] = dyn_forecast(var_list,M,options,oo,task)
% function dyn_forecast(var_list,task)
%   computes mean forecast for a given value of the parameters
%   computes also confidence band for the forecast    
%
% INPUTS
%   var_list:    list of variables (character matrix)
%   task:        indicates how to initialize the forecast
%                either 'simul' or 'smoother'
% OUTPUTS
%   nothing is returned but the procedure saves output
%   in oo_.forecast.Mean
%      oo_.forecast.HPDinf
%      oo_.forecast.HPDsup
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2015 Dynare Team
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

info = 0;

make_ex_;

maximum_lag = M.maximum_lag;

endo_names = M.endo_names;
if isempty(var_list)
    var_list = endo_names(1:M.orig_endo_nbr, :);
end
i_var = [];
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t and endogenous variable'])
    end
    i_var = [i_var; tmp];
end

n_var = length(i_var);

trend = 0;
switch task
  case 'simul'
    horizon = options.periods;
    if horizon == 0
        horizon = 5;
    end
    if isempty(M.endo_histval)
        y0 = repmat(oo.dr.ys,1,maximum_lag);
    else
        y0 = M.endo_histval;
    end
  case 'smoother'
    horizon = options.forecast;
    y_smoothed = oo.SmoothedVariables;
    y0 = zeros(M.endo_nbr,maximum_lag);
    for i = 1:M.endo_nbr
        v_name = deblank(M.endo_names(i,:));
        y0(i,:) = y_smoothed.(v_name)(end-maximum_lag+1:end)+oo.dr.ys(i);
    end
    gend = options.nobs;
    if isfield(oo.Smoother,'TrendCoeffs')
        var_obs = options.varobs;
        endo_names = M.endo_names;
        order_var = oo.dr.order_var;
        i_var_obs = [];
        trend_coeffs = [];
        for i=1:length(var_obs)
            tmp = strmatch(var_obs{i},endo_names(i_var,:),'exact');
            if ~isempty(tmp)
                i_var_obs = [ i_var_obs; tmp];
                trend_coeffs = [trend_coeffs; oo.Smoother.TrendCoeffs(i)];
            end
        end
        if ~isempty(trend_coeffs) 
          trend = trend_coeffs*(gend+(1-M.maximum_lag:horizon)); 
        end
    end
    global bayestopt_
    if isfield(bayestopt_,'mean_varobs')
        trend = trend + repmat(bayestopt_.mean_varobs,1,horizon+M.maximum_lag);
    end
  otherwise
    error('Wrong flag value')
end 

if M.exo_det_nbr == 0
    [yf,int_width] = forcst(oo.dr,y0,horizon,var_list);
else
    exo_det_length = size(oo.exo_det_simul,1)-M.maximum_lag;
    if horizon > exo_det_length
        ex = zeros(horizon,M.exo_nbr);
        oo.exo_det_simul = [ oo.exo_det_simul;...
                            repmat(oo.exo_det_steady_state',...
                                   horizon- ... 
                                   exo_det_length,1)];
    elseif horizon < exo_det_length 
        ex = zeros(exo_det_length,M.exo_nbr); 
    end
    [yf,int_width] = simultxdet(y0,ex,oo.exo_det_simul,...
                                options.order,var_list,M,oo,options);
end

if ~isscalar(trend)
    yf(i_var_obs,:) = yf(i_var_obs,:) + trend;
end

for i=1:n_var
    vname = deblank(var_list(i,:));
    forecast.Mean.(vname) = yf(i,maximum_lag+(1:horizon))';
    forecast.HPDinf.(vname)= yf(i,maximum_lag+(1:horizon))' - int_width(1:horizon,i);
    forecast.HPDsup.(vname) = yf(i,maximum_lag+(1:horizon))' + int_width(1:horizon,i);
end

for i=1:M.exo_det_nbr
    forecast.Exogenous.(deblank(M.exo_det_names(i,:))) = oo.exo_det_simul(maximum_lag+(1:horizon),i);
end

if options.nograph == 0
    oo.forecast = forecast;
    forecast_graphs(var_list,M, oo,options)
end
