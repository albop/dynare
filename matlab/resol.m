% Copyright (C) 2001 Michel Juillard
%
function dr=resol(ys,algo,linear,iorder)

global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr exo_det_nbr endo_nbr
global ex_ ex_det_ valf_ it_ exe_ exe_det_ xkmin_ xkmax_ 
global fname_ means_ stderrs_ lgy_ maxit_
global dynatol_ options_
global BlanchardKahn_


if linear == 1
  iorder =1;
end

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
  error ('RESOL: Error in model specification: some variables don"t appear as current') ;
end

if xlen > 1
  error (['RESOL: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end

% check if ys is steady state
tempex = ex_;
tempexdet = ex_det_;
ex_ = exe_';
if exo_det_nbr > 0 
  ex_det_ = ones(ykmin_+1,1)*exe_det_';
end
fh = str2func([fname_ '_fff']);
if max(abs(feval(fh,ys))) > dynatol_
  if exist([fname_ '_steadystate'])
    [dr.ys,check] = feval([fname_ '_steadystate'],ys);
  else
    [dr.ys,check] = dynare_solve([fname_ '_fff'],ys);
  end
  if check
    error('RESOL: convergence problem in DYNARE_SOLVE')
  end
else
  dr.ys = ys;
end

dr.fbias = zeros(endo_nbr,1);
dr = dr1(iorder,dr,0);

if algo == 1 & iorder > 1
  dr.ys = dynare_solve('dr2',ys,dr);
  dr.fbias = 2*feval([fname_ '_fff'],dr.ys);
  dr = dr1(iorder,dr,0);
end
ex_det_ = tempexdet;
ex_ = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for ex_