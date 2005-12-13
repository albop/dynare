function indx = subset();
% stephane.adjemian@ens.fr [11-30-2005]
global options_ estim_params_ lgx_

ExcludedParamNames = options_.ExcludedParams;
VarObs = options_.varobs;
VarExo = lgx_;
info = options_.ParamSubSet;


nvx = estim_params_.nvx;
nvn = estim_params_.nvn;
ncx = estim_params_.ncx;
ncn = estim_params_.ncn;
np  = estim_params_.np;
nx  = nvx+nvn+ncx+ncn+np;

if strcmpi(info,'All')
  indx = (1:nx)';
elseif stcmpi('DeepParameters')
  indx = (nvx+nvn+ncx+ncn+1:nx)';
elseif strcmpi('StructuralShocks')
  indx = [(1:nvx),nvx+nvn+1:nvx+nvn+ncx]';
elseif stcmpi('StructuralShocksWithoutCorrelations')
  indx = (1:nvx)';
elseif strcmpi('MeasurementErrors')
  indx = [(nvx+1:nvx+nvn),(nvx+nvn+ncx+1:nvx+nvn+ncx+ncn)]';
elseif strcmpi('MeasurementErrorsWithoutCorrelations')
  indx = (nvx+1:nvx+nvn)';
elseif strcmpi(info,'AllWithoutMeasurementErros')
  indx = [(1:nvx),nvx+nvn+1:nvx+nvn+ncx,nvx+nvn+ncx+ncn+1:nx]';
end


if ~isempty(ExcludedParamNames)
  tt = [];
  for i = 1:length(ParamNames)
    tmp = strmatch(deblank(ExcludedParamNames(i,:)),lgx_);
    if ~isempty(tmp)% The parameter the user wants to exclude is related to the size of the structural innovations.
      if ncx
	disp(['I cannot exclude some the structural variances if the'])
	disp(['structural innovations are correlated...'])
	error
      end
      tt = [tt;tmp];
    elseif isempty(tmp) & nvn
      tmp = strmatch(deblank(ExcludedParamNames(i,:)),options_.varobs);
      if ~isempty(tmp)% The parameter the user wants to exclude is related to the size of the measurement errors.
		      % variances:
	tmp = nvx+tmp;
	if ncn
	  disp(['I cannot exclude some the measurement error variances if the'])
	  disp(['measurement errors are correlated...'])
	  error
	end
	tt = [tt;tmp];
      end
    else% Excluded parameters are deep parameters...  
      tmp = strmatch(deblank(ExcludedParamNames(i,:)),estim_params_.param_names);
      if ~isempty(tmp)
	tt = [tt,nvx+nvn+ncx+ncn+tmp];
      else
	disp('The parameter you want to exclude is unknown!')
	error
      end
    end
  end
  for i=1:size(ParamNames,1)
    indx = indx(indx~=tt)
  end
end