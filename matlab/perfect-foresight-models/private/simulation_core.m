function [oo_, maxerror] = simulation_core(options_, M_, oo_)
%function [oo_, maxerror] = simulation_core(options_, M_, oo_)

% Copyright (C) 2015 Dynare Team
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

if options_.block
    if options_.bytecode
        try
            [info, tmp] = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state,1,options_.periods+2), options_.periods);
        catch
            info = 0;
        end
        if info
            oo_.deterministic_simulation.status = false;
        else
            oo_.endo_simul = tmp;
            oo_.deterministic_simulation.status = true;
        end
        if options_.no_homotopy
            mexErrCheck('bytecode', info);
        end
    else
        oo_ = feval([M_.fname '_dynamic'], options_, M_, oo_);
    end
else
    if options_.bytecode
        try
            [info, tmp] = bytecode('dynamic', oo_.endo_simul, oo_.exo_simul, M_.params, repmat(oo_.steady_state,1,options_.periods+2), options_.periods);
        catch
            info = 0;
        end
        if info
            oo_.deterministic_simulation.status = false;
        else
            oo_.endo_simul = tmp;
            oo_.deterministic_simulation.status = true;
        end
        if options_.no_homotopy
            mexErrCheck('bytecode', info);
        end
    else
        if M_.maximum_endo_lead == 0 % Purely backward model
            oo_ = sim1_purely_backward(options_, M_, oo_);
        elseif M_.maximum_endo_lag == 0 % Purely forward model
            oo_ = sim1_purely_forward(options_, M_, oo_);
        else % General case
            if options_.stack_solve_algo == 0
                oo_ = sim1(options_, M_, oo_);
            elseif options_.stack_solve_algo == 6
                oo_ = sim1_lbj(options_, M_, oo_);
            elseif options_.stack_solve_algo == 7
                periods = options_.periods;
                if ~isfield(options_.lmmcp,'lb')
                    [lb,ub,pfm.eq_index] = get_complementarity_conditions(M_);
                    options_.lmmcp.lb = repmat(lb,periods,1);
                    options_.lmmcp.ub = repmat(ub,periods,1);
                end
                y = oo_.endo_simul;
                y0 = y(:,1);
                yT = y(:,periods+2);
                z = y(:,2:periods+1);
                illi = M_.lead_lag_incidence';
                [i_cols,junk,i_cols_j] = find(illi(:));
                illi = illi(:,2:3);
                [i_cols_J1,junk,i_cols_1] = find(illi(:));
                i_cols_T = nonzeros(M_.lead_lag_incidence(1:2,:)');
                [y,info] = dynare_solve(@perfect_foresight_problem,z(:),options_, ...
                                 str2func([M_.fname '_dynamic']),y0,yT, ...
                                 oo_.exo_simul,M_.params,oo_.steady_state, ...
                                 options_.periods,M_.endo_nbr,i_cols, ...
                                 i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                                 M_.NNZDerivatives(1));
                if all(imag(y)<.1*options_.dynatol.f)
                    if ~isreal(y)
                        y = real(y);
                    end
                else
                    info = 1;
                end
                oo_.endo_simul = [y0 reshape(y,M_.endo_nbr,periods) yT];
                if info == 1
                    oo_.deterministic_simulation.status = false;
                else
                    oo_.deterministic_simulation.status = true;
                end
            end
        end
    end
end

if nargout>1
    y0 = oo_.endo_simul(:,1);
    yT = oo_.endo_simul(:,options_.periods+2);
    yy  = oo_.endo_simul(:,2:options_.periods+1);
    if ~exist('illi')
        illi = M_.lead_lag_incidence';
        [i_cols,junk,i_cols_j] = find(illi(:));
        illi = illi(:,2:3);
        [i_cols_J1,junk,i_cols_1] = find(illi(:));
        i_cols_T = nonzeros(M_.lead_lag_incidence(1:2,:)');
    end
    if options_.block && ~options_.bytecode
        maxerror = oo_.deterministic_simulation.error;
    else
        if options_.bytecode
            [chck, residuals, junk]= bytecode('dynamic','evaluate', oo_.endo_simul, oo_.exo_simul, M_.params, oo_.steady_state, 1);
        else
            residuals = perfect_foresight_problem(yy(:),str2func([M_.fname '_dynamic']), y0, yT, ...
                                                  oo_.exo_simul,M_.params,oo_.steady_state, ...
                                                  options_.periods,M_.endo_nbr,i_cols, ...
                                                  i_cols_J1, i_cols_1, i_cols_T, i_cols_j, ...
                                                  M_.NNZDerivatives(1));
        end
        maxerror = max(max(abs(residuals)));
    end
end
