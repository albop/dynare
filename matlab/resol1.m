% Copyright (C) 2001 Michel Juillard
%
function dr=resol1(ys,algo,linear,iorder)

global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr endo_nbr
global ex_ valf_ it_ exe_ xkmin_ xkmax_ 
global fname_ means_ stderrs_ lgy_ maxit_
global dynatol_ dr1_test_ bayestopt_

dr1_test_ = zeros(2,1);

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
ex_ = exe_';
fh = str2func([fname_ '_fff']);

if max(abs(feval(fh,ys))) > dynatol_
  % dirty trick to call either user function or dynare_solve
  [dr.ys, check] = feval(bayestopt_.static_solve,[fname_ '_fff'],ys);
  if check
    dr1_test_(1) = 1; % dynare_solve did not converge to the steady state.  
    resid = feval([fname_ '_fff'],dr.ys);
    dr1_test_(2) = resid'*resid; % penalty...
    disp('dynare_solve is unable to find the steady state.')
    return
  end
else
  dr.ys = ys;
end

dr.fbias = zeros(endo_nbr,1);
dr = dr11(iorder,dr,0);

if algo == 1 & iorder > 1
  dr.ys = dynare_solve('dr2',ys,dr);
  dr.fbias = 2*feval([fname_ '_fff'],dr.ys);
  dr = dr11(iorder,dr,0);
end
ex_ = tempex;
tempex = [];

% 01/01/2003 MJ added dr_algo == 1
% 08/24/2001 MJ uses Schmitt-Grohe and Uribe (2001) constant correction
%               in dr.ghs2 
% 05/26/2003 MJ added temporary values for ex_