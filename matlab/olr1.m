% Copyright (C) 2001 Michel Juillard
%
function dr = olr1(ys,algo,olr_inst,bet,obj_var,W)

global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr endo_nbr
global ex_ valf_ it_ exe_ xkmin_ xkmax_ 
global fname_ means_ stderrs_ lgy_ maxit_
global dynatol_

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
tempex = ex_;
ex_ = exe_';
fh = str2func([fname_ '_fff']);
if max(abs(feval(fh,ys))) > dynatol_
  [dr.ys, check] = dynare_solve([fname_ '_fff'],ys);
  if check
    error('OLR: convergence problem in DYNARE_SOLVE')
  end
else
  dr.ys = ys;
end
dr = olr2(dr,olr_inst,bet,obj_var,W);
ex_ = tempex;
tempex = [];

% 04/13/03 MJ



