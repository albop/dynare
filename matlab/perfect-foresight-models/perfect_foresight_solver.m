function perfect_foresight_solver()
% Computes deterministic simulations
%  
% INPUTS
%   None
%  
% OUTPUTS
%   none
%    
% ALGORITHM
%   
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 1996-2015 Dynare Team
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

global M_ options_ oo_

check_input_arguments(options_, M_, oo_);

if isempty(options_.scalv) || options_.scalv == 0
    options_.scalv = oo_.steady_state;
end

options_.scalv= 1;

if options_.debug
    model_static = str2func([M_.fname,'_static']);
    for ii=1:size(oo_.exo_simul,1)
        [residual(:,ii)] = model_static(oo_.steady_state, oo_.exo_simul(ii,:),M_.params);
    end
    problematic_periods=find(any(isinf(residual)) | any(isnan(residual)))-M_.maximum_endo_lag;
    if ~isempty(problematic_periods) 
        period_string=num2str(problematic_periods(1));
        for ii=2:length(problematic_periods)
            period_string=[period_string, ', ', num2str(problematic_periods(ii))];
        end
        fprintf('\n\nWARNING: Value for the exogenous variable(s) in period(s) %s inconsistent with the static model.\n',period_string);   
        fprintf('WARNING: Check for division by 0.\n')
    end
end

initperiods = 1:M_.maximum_lag;
lastperiods = (M_.maximum_endo_lag+options_.periods+1):(M_.maximum_endo_lag+options_.periods+M_.maximum_endo_lead);


[oo_.endo_simul,oo_.deterministic_simulation.status] = perfect_foresight_solver_core(M_,options_,oo_);

% If simulation failed try homotopy.
if ~oo_.deterministic_simulation.status && ~options_.no_homotopy
    skipline()
    disp('Simulation of the perfect foresight model failed!')
    skipline()
    
    % Disable warnings if homotopy
    warning off all
    % Do not print anything
    oldverbositylevel = options_.verbosity;
    options_.verbosity = 0;
    
    exosim = oo_.exo_simul;
    exoinit = repmat(oo_.exo_steady_state',M_.maximum_lag+options_.periods+M_.maximum_lead,1);
    
    endosim = oo_.endo_simul;
    endoinit = repmat(oo_.steady_state, 1,M_.maximum_lag+options_.periods+M_.maximum_lead);

    current_weight = 0; % Current weight of target point in convex combination
    step = .5;
    success_counter = 0;
    iteration = 0;

    fprintf('Iter. \t | Lambda \t | status \t | Max. residual\n')
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    
    while (step > options_.dynatol.x)

        if ~isequal(step,1)
            options_.verbosity = 0;
        end

        iteration = iteration+1;
        new_weight = current_weight + step; % Try this weight, and see if it succeeds

        if new_weight >= 1
            new_weight = 1; % Don't go beyond target point
            step = new_weight - current_weight;
        end

        % Compute convex combination for exo path and initial/terminal endo conditions
        % But take care of not overwriting the computed part of oo_.endo_simul
        oo_.exo_simul = exosim*new_weight + exoinit*(1-new_weight);
        
        oo_.endo_simul(:,[initperiods, lastperiods]) = ...
            new_weight*oo_.endo_simul(:,[initperiods, lastperiods])+(1-new_weight)*endoinit(:,[initperiods, lastperiods]);
        
        path_with_nans = any(any(isnan(oo_.endo_simul)));
        path_with_cplx = any(any(~isreal(oo_.endo_simul)));
        
        if isequal(iteration,1)
            oo_.endo_simul(:,M_.maximum_lag+1:end-M_.maximum_lead) = endoinit(:,1:options_.periods);
        elseif path_with_nans || path_with_cplx
            oo_.endo_simul(:,M_.maximum_lag+1:end-M_.maximum_lead) = saved_endo_simul(:,1+M_.maximum_endo_lag:end-M_.maximum_endo_lead);
        end
        
        saved_endo_simul = oo_.endo_simul;

        [oo_.endo_simul,oo_.deterministic_simulation.status,me] = perfect_foresight_solver_core(M_,options_,oo_);

        if oo_.deterministic_simulation.status == 1
            current_weight = new_weight;
            if current_weight >= 1
                fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_weight, 'succeeded', me)
                break
            end
            success_counter = success_counter + 1;
            if success_counter >= 3
                success_counter = 0;
                step = step * 2;
            end
            fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_weight, 'succeeded', me)
        else
            oo_.endo_simul = saved_endo_simul;
            success_counter = 0;
            step = step / 2;
            if isreal(me)
                fprintf('%i \t | %1.5f \t | %s \t | %e\n', iteration, new_weight, 'failed', me)
            else
                fprintf('%i \t | %1.5f \t | %s \t | %s\n', iteration, new_weight, 'failed', 'Complex')
            end
        end
    end
    fprintf('++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n')
    skipline()
    options_.verbosity = oldverbositylevel;
    warning on all
end

if oo_.deterministic_simulation.status == 1
    disp('Perfect foresight solution found.')
    skipline()
else
    warning('Failed to solve perfect foresight model')
end

dyn2vec;

if ~isdates(options_.initial_period) && isnan(options_.initial_period)
    initial_period = dates(1,1);
else
    initial_period = options_.initial_period;
end

ts = dseries(transpose(oo_.endo_simul),initial_period,cellstr(M_.endo_names));
assignin('base', 'Simulated_time_series', ts);
