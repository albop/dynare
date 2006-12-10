function [dr_,info]=osr1(params,weights)
  global xkmax_ xkmin_ ykmin_ ykmax_ ys_ iy_ exo_nbr endo_nbr fname_ ...
      dynatol_ options_ it_ dr_ lgy_

  it_ = ykmin_ + 1 ;

  set_default_option(options_,'dr_algo',0);
  
  if exo_nbr == 0
    exe_ = [] ;
  end

  check_model;
  
  np = size(params,1);
  t0 = zeros(np,1);
  for i=1:np
    t0(i)=evalin('base',[params(i,:) ';']);
  end
  
  [f,vx,info] = osr_obj(t0,params,weights);
  if info > 0
    disp('OSR: bad initial value for the parameters');
    return
  end
  options = optimset('fminsearch');
%  options = optimset(options,'display','iter');
  [p,f]=fminsearch(@osr_obj,t0,options,params,weights);

  [f,vx,info] = osr_obj(p,params,weights);
  if info > 0
    print_info(info);
    disp(['OSR ends on a pathological case, try different initial values' ...
	  ' for the parameters']);
    return
  else
    disp('')
    disp('OPTIMAL VALUE OF THE PARAMETERS:')
    disp('')
    for i=1:np
      disp(sprintf('%16s %16.6g\n',params(i,:),p(i)))
    end
    disp(sprintf('Objective function : %16.6g\n',f));
    disp(' ')
    dr_=resol(ys_,0);
    disp(' ')
    disp('Contributions to the objective function')
    disp('Variables           (Co)variance         Weight     Contribution')
    vx1(dr_.order_var,dr_.order_var) = vx;
    for i=1:endo_nbr
      for j=1:i
	if weights(i,j) > 0
	  if i==j
	    x1 = weights(i,i);
	    vname = deblank(lgy_(i,:));
	  else
	    x1 = 2*weights(i,j);
	    vname = [deblank(lgy_(i,:)) ', ' deblank(lgy_(j,:))];
	  end
	  st = sprintf('%15s %14.4f  %14.1f  %14.4f',...
		       vname,vx1(i,j),x1,x1*vx1(i,j));
	  disp(st)
	end
      end
    end
    disp(' ')
    disp(' ')
  end

  % 05/10/03 MJ modified to work with osr.m and give full report