function dr_=osr1(params,weights)
  global xkmax_ xkmin_ ykmin_ ykmax_ ys_ iy_ exo_nbr endo_nbr fname_ ...
      dynatol_ options_ it_

  xlen = xkmax_ + xkmin_ + 1;
  klen = ykmin_ + ykmax_ + 1;
  iyv = iy_';
  iyv = iyv(:);
  iyr0 = find(iyv) ;
  it_ = ykmin_ + 1 ;


  if exo_nbr == 0
    exe_ = [] ;
  end

  if ~ iy_(ykmin_+1,:) > 0
    error ('OLR: Error in model specification: some variables don"t appear as current') ;
  end

  if ykmax_ == 0
    error ('Backward or static model: no point in using OLR') ;
  end

  if xlen > 1
    error (['OLR: stochastic exogenous variables must appear only at the' ...
	    ' current period. Use additional endogenous variables']) ;
  end

  % check if ys is steady state
  fh = str2func([fname_ '_fff']);
  if max(abs(feval(fh,ys_))) > dynatol_
    [dr_.ys, check] = dynare_solve([fname_ '_fff'],ys_);
    if check
      error('OLR: convergence problem in DYNARE_SOLVE')
    end
  else
    dr_.ys = ys_;
  end

  
  np = size(params,1);
  t0 = zeros(np,1);
  for i=1:np
    t0(i)=evalin('base',[params(i,:) ';']);
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
  dr_=resol(ys_,options_.dr_algo,options_.linear,options_.order);

  % 05/10/03 MJ modified to work with osr.m and give full report