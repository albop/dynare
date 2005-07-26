function dr_=osr1(params,weights)
  global xkmax_ xkmin_ ykmin_ ykmax_ ys_ iy_ exo_nbr endo_nbr fname_ ...
      dynatol_ options_ it_

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
  
  [f,info] = osr_obj(t0,params,weights);
  if info > 0
    disp('OSR: bad initial value for the parameters');
  end
  [p,f]=fminsearch(@osr_obj,t0,[],params,weights);

  disp('')
  disp('OPTIMAL VALUE OF THE PARAMETERS:')
  disp('')
  for i=1:np
    disp(sprintf('%16s %16.6g\n',params(i,:),p(i)))
  end
  disp(sprintf('Objective function : %16.6g\n',f));
  disp(' ')
  dr_=resol(ys_,0);

  % 05/10/03 MJ modified to work with osr.m and give full report