function [xparams, logpost]=metropolis_draw(init)
% function [xparams, logpost]=metropolis_draw(init) 
% Builds draws from metropolis
%
% INPUTS:
%   init:              scalar equal to 1 (first call) or 0
%
% OUTPUTS:
%   xparams:           vector of estimated parameters
%   logpost:           log of posterior density
%   
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2003-2011 Dynare Team
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

global options_ estim_params_ M_
persistent mh_nblck NumberOfDraws BaseName FirstLine FirstMhFile MAX_nruns

xparams = 0;
logpost = 0;

if init
    nvx  = estim_params_.nvx;
    nvn  = estim_params_.nvn;
    ncx  = estim_params_.ncx;
    ncn  = estim_params_.ncn;
    np   = estim_params_.np ;
    npar = nvx+nvn+ncx+ncn+np;
    MetropolisFolder = CheckPath('metropolis',M_.dname);
    FileName = M_.fname;
    BaseName = [MetropolisFolder filesep FileName];
    load_last_mh_history_file(MetropolisFolder, FileName);
    FirstMhFile = record.KeepedDraws.FirstMhFile;
    FirstLine = record.KeepedDraws.FirstLine; 
    TotalNumberOfMhFiles = sum(record.MhDraws(:,2)); 
    LastMhFile = TotalNumberOfMhFiles; 
    TotalNumberOfMhDraws = sum(record.MhDraws(:,1));
    NumberOfDraws = TotalNumberOfMhDraws-floor(options_.mh_drop*TotalNumberOfMhDraws);
    MAX_nruns = ceil(options_.MaxNumberOfBytes/(npar+2)/8);
    mh_nblck = options_.mh_nblck;
    % set sub_draws option if empty
    if isempty(options_.sub_draws)
        options_.sub_draws = min(options_.posterior_max_subsample_draws, round(.25*NumberOfDraws));
    else
        if options_.sub_draws>NumberOfDraws
            skipline()
            disp(['Estimation::mcmc: The value of option sub_draws (' num2str(options_.sub_draws) ') is greater than the number of available draws in the MCMC (' num2str(NumberOfDraws) ')!'])
            disp('Estimation::mcmc: You can either change the value of sub_draws, reduce the value of mh_drop, or run another mcmc (with the load_mh_file option).')
            skipline()
            xparams = 1; % xparams is interpreted as an error flag
        end
    end
    return
end

ChainNumber = ceil(rand*mh_nblck);
DrawNumber  = ceil(rand*NumberOfDraws);

if DrawNumber <= MAX_nruns-FirstLine+1
    MhFilNumber = FirstMhFile;
    MhLine = FirstLine+DrawNumber-1;
else
    DrawNumber  = DrawNumber-(MAX_nruns-FirstLine+1);
    MhFilNumber = FirstMhFile+ceil(DrawNumber/MAX_nruns); 
    MhLine = DrawNumber-(MhFilNumber-FirstMhFile-1)*MAX_nruns;
end

load( [ BaseName '_mh' int2str(MhFilNumber) '_blck' int2str(ChainNumber) '.mat' ],'x2','logpo2');
xparams = x2(MhLine,:);
logpost= logpo2(MhLine);