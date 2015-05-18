function myoutput = random_walk_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% function myoutput = random_walk_metropolis_hastings_core(myinputs,fblck,nblck,whoiam, ThisMatlab)
% Contains the most computationally intensive portion of code in
% random_walk_metropolis_hastings (the 'for xxx = fblck:nblck' loop). The branches in  that 'for'
% cycle are completely independent to be suitable for parallel execution.
%
% INPUTS
%   o myimput            [struc]     The mandatory variables for local/remote
%                                    parallel computing obtained from random_walk_metropolis_hastings.m
%                                    function.
%   o fblck and nblck    [integer]   The Metropolis-Hastings chains.
%   o whoiam             [integer]   In concurrent programming a modality to refer to the different threads running in parallel is needed.
%                                    The integer whoaim is the integer that
%                                    allows us to distinguish between them. Then it is the index number of this CPU among all CPUs in the
%                                    cluster.
%   o ThisMatlab         [integer]   Allows us to distinguish between the
%                                    'main' Matlab, the slave Matlab worker, local Matlab, remote Matlab,
%                                     ... Then it is the index number of this slave machine in the cluster.
% OUTPUTS
%   o myoutput  [struc]
%               If executed without parallel, this is the original output of 'for b =
%               fblck:nblck'. Otherwise, it's a portion of it computed on a specific core or
%               remote machine. In this case:
%                               record;
%                               irun;
%                               NewFile;
%                               OutputFileName
%
% ALGORITHM
%   Portion of Metropolis-Hastings.
%
% SPECIAL REQUIREMENTS.
%   None.
% 
% PARALLEL CONTEXT
% See the comments in the random_walk_metropolis_hastings.m funtion.


% Copyright (C) 2006-2015 Dynare Team
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

if nargin<4,
    whoiam=0;
end

% reshape 'myinputs' for local computation.
% In order to avoid confusion in the name space, the instruction struct2local(myinputs) is replaced by:

TargetFun=myinputs.TargetFun;
ProposalFun=myinputs.ProposalFun;
xparam1=myinputs.xparam1;
mh_bounds=myinputs.mh_bounds;
last_draw=myinputs.ix2;
last_posterior=myinputs.ilogpo2;
fline=myinputs.fline;
npar=myinputs.npar;
nruns=myinputs.nruns;
NewFile=myinputs.NewFile;
MAX_nruns=myinputs.MAX_nruns;
d=myinputs.d;
InitSizeArray=myinputs.InitSizeArray;
record=myinputs.record;
dataset_ = myinputs.dataset_;
dataset_info = myinputs.dataset_info;
bayestopt_ = myinputs.bayestopt_;
estim_params_ = myinputs.estim_params_;
options_ = myinputs.options_;
M_ = myinputs.M_;
oo_ = myinputs.oo_;

% Necessary only for remote computing!
if whoiam
    % initialize persistent variables in priordens()
    priordens(xparam1,bayestopt_.pshape,bayestopt_.p6,bayestopt_.p7, bayestopt_.p3,bayestopt_.p4,1);
end

MetropolisFolder = CheckPath('metropolis',M_.dname);
ModelName = M_.fname;
BaseName = [MetropolisFolder filesep ModelName];

options_.lik_algo = 1;
OpenOldFile = ones(nblck,1);
if strcmpi(ProposalFun,'rand_multivariate_normal')
    n = npar;
elseif strcmpi(ProposalFun,'rand_multivariate_student')
    n = options_.student_degrees_of_freedom;
end

%
% Now I run the (nblck-fblck+1) Metropolis-Hastings chains
%

proposal_covariance_Cholesky_decomposition = d*diag(bayestopt_.jscale);

block_iter=0;

for curr_block = fblck:nblck,
    block_iter=block_iter+1;
    try
        % This will not work if the master uses a random number generator not
        % available in the slave (different Matlab version or
        % Matlab/Octave cluster). Therefore the trap.
        %
        % Set the random number generator type (the seed is useless but needed by the function)
        set_dynare_seed(options_.DynareRandomStreams.algo, options_.DynareRandomStreams.seed);
        % Set the state of the RNG
        set_dynare_random_generator_state(record.InitialSeeds(curr_block).Unifor, record.InitialSeeds(curr_block).Normal);
    catch
        % If the state set by master is incompatible with the slave, we only reseed 
        set_dynare_seed(options_.DynareRandomStreams.seed+curr_block);
    end
    if (options_.load_mh_file~=0) && (fline(curr_block)>1) && OpenOldFile(curr_block) %load previous draws and likelihood
        load([BaseName '_mh' int2str(NewFile(curr_block)) '_blck' int2str(curr_block) '.mat'])
        x2 = [x2;zeros(InitSizeArray(curr_block)-fline(curr_block)+1,npar)];
        logpo2 = [logpo2;zeros(InitSizeArray(curr_block)-fline(curr_block)+1,1)];
        OpenOldFile(curr_block) = 0;
    else
        x2 = zeros(InitSizeArray(curr_block),npar);
        logpo2 = zeros(InitSizeArray(curr_block),1);
    end
    %Prepare waiting bars
    if whoiam
        prc0=(curr_block-fblck)/(nblck-fblck+1)*(isoctave || options_.console_mode);
        hh = dyn_waitbar({prc0,whoiam,options_.parallel(ThisMatlab)},['MH (' int2str(curr_block) '/' int2str(options_.mh_nblck) ')...']);
    else
        hh = dyn_waitbar(0,['Metropolis-Hastings (' int2str(curr_block) '/' int2str(options_.mh_nblck) ')...']);
        set(hh,'Name','Metropolis-Hastings');
    end
    accepted_draws_this_chain = 0;
    accepted_draws_this_file = 0;
    draw_index_current_file = fline(curr_block); %get location of first draw in current block
    draw_iter = 1;
    while draw_iter <= nruns(curr_block)
        par = feval(ProposalFun, last_draw(curr_block,:), proposal_covariance_Cholesky_decomposition, n);
        if all( par(:) > mh_bounds.lb ) && all( par(:) < mh_bounds.ub )
            try
                logpost = - feval(TargetFun, par(:),dataset_,dataset_info,options_,M_,estim_params_,bayestopt_,mh_bounds,oo_);
            catch
                logpost = -inf;
            end
        else
            logpost = -inf;
        end
        if (logpost > -inf) && (log(rand) < logpost-last_posterior(curr_block))
            x2(draw_index_current_file,:) = par;
            last_draw(curr_block,:) = par;
            logpo2(draw_index_current_file) = logpost;
            last_posterior(curr_block) = logpost;
            accepted_draws_this_chain = accepted_draws_this_chain + 1;
            accepted_draws_this_file = accepted_draws_this_file + 1;
        else
            x2(draw_index_current_file,:) = last_draw(curr_block,:);
            logpo2(draw_index_current_file) = last_posterior(curr_block);
        end
        prtfrc = draw_iter/nruns(curr_block);
        if (mod(draw_iter, 3)==0 && ~whoiam) || (mod(draw_iter,50)==0 && whoiam)
            dyn_waitbar(prtfrc,hh,[ 'MH (' int2str(curr_block) '/' int2str(options_.mh_nblck) ') ' sprintf('Current acceptance ratio %4.3f', accepted_draws_this_chain/draw_iter)]);
        end
        if (draw_index_current_file == InitSizeArray(curr_block)) || (draw_iter == nruns(curr_block)) % Now I save the simulations, either because the current file is full or the chain is done
            save([BaseName '_mh' int2str(NewFile(curr_block)) '_blck' int2str(curr_block) '.mat'],'x2','logpo2');
            fidlog = fopen([MetropolisFolder '/metropolis.log'],'a');
            fprintf(fidlog,['\n']);
            fprintf(fidlog,['%% Mh' int2str(NewFile(curr_block)) 'Blck' int2str(curr_block) ' (' datestr(now,0) ')\n']);
            fprintf(fidlog,' \n');
            fprintf(fidlog,['  Number of simulations.: ' int2str(length(logpo2)) '\n']);
            fprintf(fidlog,['  Acceptance ratio......: ' num2str(accepted_draws_this_file/length(logpo2)) '\n']);
            fprintf(fidlog,['  Posterior mean........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(mean(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(mean(logpo2)) '\n']);
            fprintf(fidlog,['  Minimum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(min(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(min(logpo2)) '\n']);
            fprintf(fidlog,['  Maximum value.........:\n']);
            for i=1:length(x2(1,:))
                fprintf(fidlog,['    params:' int2str(i) ': ' num2str(max(x2(:,i))) '\n']);
            end
            fprintf(fidlog,['    log2po:' num2str(max(logpo2)) '\n']);
            fprintf(fidlog,' \n');
            fclose(fidlog);
            accepted_draws_this_file = 0;
            if draw_iter == nruns(curr_block) % I record the last draw...
                record.LastParameters(curr_block,:) = x2(end,:);
                record.LastLogPost(curr_block) = logpo2(end);
            end
            % size of next file in chain curr_block
            InitSizeArray(curr_block) = min(nruns(curr_block)-draw_iter,MAX_nruns);
            % initialization of next file if necessary
            if InitSizeArray(curr_block)
                x2 = zeros(InitSizeArray(curr_block),npar);
                logpo2 = zeros(InitSizeArray(curr_block),1);
                NewFile(curr_block) = NewFile(curr_block) + 1;
                draw_index_current_file = 0;
            end
        end
        draw_iter=draw_iter+1;
        draw_index_current_file = draw_index_current_file + 1;
    end% End of the simulations for one mh-block.
    record.AcceptanceRatio(curr_block) = accepted_draws_this_chain/draw_iter;
    dyn_waitbar_close(hh);
    [record.LastSeeds(curr_block).Unifor, record.LastSeeds(curr_block).Normal] = get_dynare_random_generator_state();
    OutputFileName(block_iter,:) = {[MetropolisFolder,filesep], [ModelName '_mh*_blck' int2str(curr_block) '.mat']};
end% End of the loop over the mh-blocks.


myoutput.record = record;
myoutput.irun = draw_index_current_file;
myoutput.NewFile = NewFile;
myoutput.OutputFileName = OutputFileName;