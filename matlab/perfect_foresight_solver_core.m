function [endo_simul, status] = perfect_foresight_solver_core(M,oo,options)
% Core function to compute deterministic simulations
%  
% INPUTS
%   M: model structure
%   oo: output structure
%   options: options structure
%  
% OUTPUTS
%   endo_simul: matrix endogenous variables
%   deterministic_simulation: simulation status
%    
% ALGORITHM
%   
%   various
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

endo_simul = oo.endo_simul;
status = 0;
deterministic_simulation = struct();

if(options.block)
    if(options.bytecode)
        [info, endo_simul] = bytecode('dynamic');
        if info == 1
            status = 0;
        else
            status = 1;
        end
        mexErrCheck('bytecode', info);
    else
        eval([M.fname '_dynamic']);
    end
else
    if(options.bytecode)
        [info, endo_simul]=bytecode('dynamic');
        if info == 1
            status = 0;
        else
            status = 1;
        end;
        mexErrCheck('bytecode', info);
    else
        if M.maximum_endo_lead == 0 
            % Purely backward model
            global oo_ options_
            oo_ = oo;
            options_ = options;
            sim1_purely_backward;
            endo_simul = oo_.endo_simul;
            if oo_.deterministic_simulation.status == 1
                status = 0;
            end
        elseif M.maximum_endo_lag == 0 
            % Purely forward model
            global oo_ options_
            oo_ = oo;
            options_ = options;
            sim1_purely_forward;
            endo_simul = oo_.endo_simul;
            if oo_.deterministic_simulation.status == 1
                status = 0;
            end
        else 
            % General case
            if options.stack_solve_algo == 0
                oo = sim1(M,options,oo);
                endo_simul = oo.endo_simul;
                if oo.deterministic_simulation.status == 1
                    status = 0;
                end
            elseif options.stack_solve_algo == 6
                global oo_ options_
                oo_ = oo;
                options_ = options;
                sim1_lbj;
                endo_simul = oo_.endo_simul;
                if oo_.deterministic_simulation.status == 1
                    status = 0;
                end
            elseif options.stack_solve_algo == 7
                periods = options.periods;
                if ~isfield(options.lmmcp,'lb')
                    [lb,ub,pfm.eq_index] = get_complementarity_conditions(M);
                    options.lmmcp.lb = repmat(lb,periods,1);
                    options.lmmcp.ub = repmat(ub,periods,1);
                end

                y = endo_simul;
                y0 = y(:,1);
                yT = y(:,periods+2);
                z = y(:,2:periods+1);
                illi = M.lead_lag_incidence';
                [i_cols,~,i_cols_j] = find(illi(:));
                illi = illi(:,2:3);
                [i_cols_J1,~,i_cols_1] = find(illi(:));
                i_cols_T = nonzeros(M.lead_lag_incidence(1:2,:)');
                [y,info] = dynare_solve(@perfect_foresight_problem,z(:),1, ...
                                        str2func([M.fname '_dynamic']),y0,yT, ...
                                        oo.exo_simul,M.params,oo.steady_state, ...
                                        options.periods,M.endo_nbr,i_cols, ...
                                        i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                                        M.NNZDerivatives(1));
                endo_simul = [y0 reshape(y,M.endo_nbr,periods) yT];
                if info == 1
                    status = 1;
                else
                    status = 0;
                end;
            end
        end
    end
end


