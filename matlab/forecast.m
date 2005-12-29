% Copyright (C) 2005 Michel Juillard
%
function forecast(var_list)
  global options_ dr_ ys_ endo_nbr exo_nbr exo_det_nbr ykmin_ y_ ex_det_ 
  global lgy_ lgx_det_ oo_ exe_det_
  
  old_options = options_;
  options_ = set_default_option(options_,'periods',40);
  if options_.periods == 0
    options_.periods = 40;
  end
  
  if size(y_,2) < ykmin_
    y0 = repmat(ys_,ykmin_);
  else
    y0 = y_(:,1:ykmin_);
  end
  
  if exo_det_nbr == 0
    [yf,var_yf] = forcst(dr_,y0,options_.periods,var_list);
  else
    ex = zeros(options_.periods,exo_nbr);
    if options_.periods > size(ex_det_,1)
      ex_det_ = [ ex_det_; repmat(exe_det_',options_.periods- ...
				  size(ex_det_,1),1)];
    end
    [yf,int_width] = simultxdet(y0,dr_,ex,ex_det_,options_.order, ...
				  var_list);
  end
  
  for i=1:endo_nbr
    eval(['oo_.forecast.Mean.' lgy_(i,:) '= yf(' int2str(i) ',:)'';']);
    eval(['oo_.forecast.HPDinf.' lgy_(i,:) '= yf(' int2str(i) ',2:end)''-' ...
		    ' int_width(:,' int2str(i) ');']);
    eval(['oo_.forecast.HPDsup.' lgy_(i,:) '= yf(' int2str(i) ',2:end)''+' ...
		    ' int_width(:,' int2str(i) ');']);
  end

  for i=1:exo_det_nbr
    eval(['oo_.forecast.Exogenous.' lgx_det_(i,:) '= ex_det_(:,' int2str(i) ');']);
  end
  
  options_ = old_options;