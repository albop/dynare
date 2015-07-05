function oo_ = sim1(options_, M_, oo_)
% function sim1
% Performs deterministic simulations with lead or lag on one period.
% Uses sparse matrices.
%
% INPUTS
%   ...
% OUTPUTS
%   ...
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   None.

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

verbose = options_.verbosity;

lead_lag_incidence = M_.lead_lag_incidence;

ny = M_.endo_nbr;
nx = M_.exo_nbr;

maximum_lag = M_.maximum_lag;
max_lag = M_.maximum_endo_lag;

nyp = nnz(lead_lag_incidence(1,:)) ;
ny0 = nnz(lead_lag_incidence(2,:)) ;
nyf = nnz(lead_lag_incidence(3,:)) ;

nd = nyp+ny0+nyf; % size of y (first argument passed to the dynamic file).

periods = options_.periods;
y_steady_state = oo_.steady_state;
x_steady_state = oo_.exo_steady_state;

params = M_.params;
endo_simul = oo_.endo_simul;
exo_simul = oo_.exo_simul;

% Indices in A.
ip   = find(lead_lag_incidence(1,:)');
ic   = find(lead_lag_incidence(2,:)');
in   = find(lead_lag_incidence(3,:)');
icn  = find(lead_lag_incidence(2:3,:)');
ipcn = find(lead_lag_incidence');

% Indices in y.
jp  = nonzeros(lead_lag_incidence(1,:)');
jc  = nonzeros(lead_lag_incidence(2,:)');
jn  = nonzeros(lead_lag_incidence(3,:)');
jpc = [jp; jc];
jcn = [jc; jn];

jexog = transpose(nd+(1:nx));
jendo = transpose(1:nd);

i_upd = maximum_lag*ny+(1:periods*ny);

% Center the endogenous and exogenous variables around the deterministic steady state.
endo_simul = bsxfun(@minus,endo_simul,y_steady_state);
exo_simul =bsxfun(@minus,exo_simul,transpose(x_steady_state));

Y = endo_simul(:);

if verbose
    skipline()
    printline(80)
    disp('MODEL SIMULATION:')
    skipline()
end

model_dynamic = str2func([M_.fname,'_dynamic']);
z = y_steady_state([ip; ic; in]);

% Evaluate the Jacobian of the dynamic model at the deterministic steady state.
[d1,jacobian] = model_dynamic(z,transpose(x_steady_state), params, y_steady_state,1);

% Check that the dynamic model was evaluated at the steady state.
if max(abs(d1))>1e-12
    error('Jacobian is not evaluated at the steady state!')
end

[rT,cT,vT] = find(jacobian(:,jpc));
[r1,c1,v1] = find(jacobian(:,jcn));
[rr,cc,vv] = find(jacobian(:,jendo));

ivT = 1:length(vT);
iv1 = 1:length(v1);
iv  = 1:length(vv);

% Initialize the vector of residuals.
res = zeros(periods*ny,1);

% Initialize the sparse Jacobian.
iA = zeros(periods*M_.NNZDerivatives(1),3);

h2 = clock ;
i_rows = (1:ny)';
i_cols_A = ipcn;
i_cols = ipcn+(maximum_lag-1)*ny;
m = 0;
for it = (maximum_lag+1):(maximum_lag+periods)
    if it == maximum_lag+periods
        nv = length(vT);
        iA(ivT+m,:) = [i_rows(rT),i_cols_A(jpc(cT)),vT];
    elseif it == maximum_lag+1
        nv = length(v1);
        iA(iv1+m,:) = [i_rows(r1),icn(c1),v1];
    else
        nv = length(vv);
        iA(iv+m,:) = [i_rows(rr),i_cols_A(cc),vv];
    end
    z(jendo) = Y(i_cols);
    z(jexog) = transpose(exo_simul(it,:));
    res(i_rows) = jacobian*z;
    m = m + nv;
    i_rows = i_rows + ny;
    i_cols = i_cols + ny;
    if it > maximum_lag+1
        i_cols_A = i_cols_A + ny;
    end
end

% Evaluation of the maximum residual at the initial guess (steady state for the endogenous variables).
err = max(abs(res));

if options_.debug
    fprintf('\nLargest absolute residual at iteration %d: %10.3f\n',1,err);
    if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y))
        fprintf('\nWARNING: NaN or Inf detected in the residuals or endogenous variables.\n');
    end
    if ~isreal(res) || ~isreal(Y)
        fprintf('\nWARNING: Imaginary parts detected in the residuals or endogenous variables.\n');
    end
    skipline()
end

iA = iA(1:m,:);
A = sparse(iA(:,1),iA(:,2),iA(:,3),periods*ny,periods*ny);

% Try to update the vector of endogenous variables.
try
    Y(i_upd) =  Y(i_upd) - A\res;
catch
    % Normally, because the model is linear, the solution of the perfect foresight model should
    % be obtained in one Newton step. This is not the case if the model is singular.
    oo_.deterministic_simulation.status = false;
    oo_.deterministic_simulation.error = NaN;
    oo_.deterministic_simulation.iterations = 1;
    if verbose
        skipline()
        disp('Singularity problem! The jacobian matrix of the stacked model cannot be inverted.')
    end
    return
end

i_cols = ipcn+(maximum_lag-1)*ny;
i_rows = (1:ny)';
for it = (maximum_lag+1):(maximum_lag+periods)
    z(jendo) = Y(i_cols);
    z(jexog) = transpose(exo_simul(it,:));
    m = m + nv;
    res(i_rows) = jacobian*z;
    i_rows = i_rows + ny;
    i_cols = i_cols + ny;
end

ERR = max(abs(res));

if verbose
    fprintf('Iter: %s,\t Initial err. = %s,\t err. = %s,\t time = %s\n',num2str(1),num2str(err),num2str(ERR), num2str(etime(clock,h2)));
    printline(80);
end

if any(isnan(res)) || any(isinf(res)) || any(isnan(Y)) || any(isinf(Y)) || ~isreal(res) || ~isreal(Y)
    oo_.deterministic_simulation.status = false;% NaN or Inf occurred
    oo_.deterministic_simulation.error = ERR;
    oo_.deterministic_simulation.iterations = 1;
    oo_.endo_simul = reshape(Y,ny,periods+maximum_lag+M_.maximum_lead);
    if verbose
        skipline()
        if ~isreal(res) || ~isreal(Y)
            disp('Simulation terminated with imaginary parts in the residuals or endogenous variables.')
        else
            disp('Simulation terminated with NaN or Inf in the residuals or endogenous variables.')
        end
        disp('There is most likely something wrong with your model. Try model_diagnostics or another simulation method.')
    end
else
    oo_.deterministic_simulation.status = true;% Convergency obtained.
    oo_.deterministic_simulation.error = ERR;
    oo_.deterministic_simulation.iterations = 1;
    oo_.endo_simul = bsxfun(@plus,reshape(Y,ny,periods+maximum_lag+M_.maximum_lead),y_steady_state);
end

if verbose
    skipline();
end
