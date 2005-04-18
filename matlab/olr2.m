% Copyright (C) 2001 Michel Juillard
%
function dr=olr2(dr,olr_inst,bet,obj_var,W)

global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr endo_nbr
global ex_ valf_ it_ exe_ xkmin_ xkmax_ ys_ stdexo_
global fname_ means_ Sigma_e_ lgy_
global eigval options_

options_ = set_default_option(options_,'loglinear',0);

xlen = xkmax_ + xkmin_ + 1;
klen = ykmin_ + ykmax_ + 1;
iyv = iy_';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = ykmin_ + 1 ;

inst_nbr = size(olr_inst,1);
inst_i = zeros(inst_nbr,1);
for i=1:inst_nbr
  k = strmatch(olr_inst(i,:),lgy_,'exact');
  if isempty(k)
    error(sprintf('OLR_INST %s isn''t a declared variable'));
  else
    inst_i(i) = k;
  end
end

if ykmax_ == 0
  error ('OLR : No forward variable: no point in using OLR') ;
end

if find(any(iy_([1:ykmin_ ykmin_+2:ykmax_],inst_i),2))
  error('OLR: instruments can only appear at the current period');
end

non_inst_i = setdiff([1:endo_nbr],inst_i);
iy1_ = iy_(:,non_inst_i);
endo_nbr_1 = endo_nbr - inst_nbr;

if exo_nbr == 0
  exe_ = [] ;
end

if ~ all(iy_(ykmin_+1,:) > 0)
  error ('Error in model specification: some variables don"t appear as current') ;
end

if xlen > 1
  error (['SS: stochastic exogenous variables must appear only at the' ...
	  ' current period. Use additional endogenous variables']) ;
end


dr=set_state_space(dr);
kstate = dr.kstate;
kad = dr.kad;
kae = dr.kae;
nstatic = dr.nstatic;
nfwrd = dr.nfwrd;
nsfwrd = dr.nsfwrd;
npred = dr.npred;
nspred = dr.nspred;
nboth = dr.nboth;
order_var = dr.order_var;
nd = size(kstate,1);
stat_var = order_var(1:nstatic);
stat_var = setdiff(stat_var,inst_i);
% static variables in objective function
[stat_obj_var] = intersect(obj_var,stat_var);
n_stat_obj_var = length(stat_obj_var);
pred_var = find(any(iy_(1:ykmin_,:),1))';
pred_var = [stat_obj_var; pred_var];
npred = length(pred_var)-nboth;
nstatic = length(stat_var);
nstatic1 = nstatic+inst_nbr;

endo_nbr_1 = endo_nbr-inst_nbr;
sdyn = endo_nbr - nstatic1;

order_var = [ inst_i; setdiff(stat_var,inst_i); order_var(nstatic1+1:end)];

% building QQ, RR and UU
iq = [];
iu = [];
io1 = [];
io2 = [];
for i=1:length(obj_var)
  i1 = find(obj_var(i)==order_var);
  if i1 > inst_nbr
    iq = [iq; i1-nstatic1];
    io1 = [io1; obj_var(i)];
  else
    iu = [iu; i1];
    io2 = [io2; obj_var(i)];
  end
end
QQ = zeros(sdyn,sdyn);
QQ(iq,iq) = W(io1,io1);
RR = zeros(sdyn,inst_nbr);
RR(iq,iu) = W(io1,io2);
UU = zeros(inst_nbr,inst_nbr);
UU(iu,iu) = W(io2,io2);

tempex = ex_;

it_ = ykmin_ + 1;
z = repmat(dr.ys,1,klen);
z = z(iyr0) ;
jacobia_=real(jacob_a('ff1_',[z; exe_])) ;

ex_ = tempex ;
tempex = [];

nz = size(z,1);
k1 = iy_(find([1:klen] ~= ykmin_+1),:);
b = jacobia_(1:endo_nbr_1,iy_(ykmin_+1,order_var(inst_nbr+1:end)));
a = b\jacobia_(1:endo_nbr_1,nonzeros(k1')); 
if any(isinf(a(:)))
  error('OLR: the model doesn''t determine current variables uniquely')
end
if exo_nbr
  fu = b\jacobia_(1:endo_nbr_1,nz+1:end);
end
% instruments' effects
b = b\jacobia_(1:endo_nbr_1,iy_(ykmin_+1,inst_i));

% reordered incidence matrix
siy = iy_(:,order_var);

sdyn1 = endo_nbr-nstatic-inst_nbr;          

ilambda = cell(klen,1);
jlambda = cell(klen,1);
klambda = cell(klen,1);
% number of multipliers corresponding to predetermined variables
aa = a(nstatic+1:end,:);
i_cum1 = any(aa(:,nonzeros(siy(1,:))),2);
ilambda{1} = find(i_cum1);
jlambda{1} = ilambda{1}
klambda{1} = [1:length(ilambda{1})]';
for i = 2:ykmin_
  i1 = any(aa(:,nonzeros(siy(i,:))),2);
  ilambda{i} = find(i1);
  i_cum1 = i_cum1 | i1;
  jlambda{i} = find(i_cum1);
  klambda{i} = find(i1(jlambda{i}));
end
%number of multipliers corresponding to forward looking variables
i_cum2 = any(aa(:,nonzeros(siy(klen,:))-endo_nbr),2);
ilambda{klen} = find(i_cum2);
jlambda{klen} = ilambda{klen};
klambda{klen} = [1:length(ilambda{klen})]';
for i = klen-1:-1:ykmin_+2
  i1 = any(aa(:,nonzeros(siy(i,:))-endo_nbr),2);
  ilambda{i} = find(i1);
  i_cum2 = i_cum2 | i1;
  jlambda{i} = find(i_cum2);
  klambda{i} = find(i1(jlambda{i}));
end
% the entries in ykmin_+1 are only used to set 1s corresponding to 
% lambda(t) in d matrix
i1 = ones(sdyn1,1);
j1 = find(i_cum1);
i1(j1) = zeros(length(j1),1);
j1 = find(i_cum2);
i2 = i1;
i2(j1) = ones(length(j1),1);
jlambda{ykmin_+1} = find(i2);
ilambda{ykmin_+1} = find(i1);
klambda{ykmin_+1} = find(i1(jlambda{ykmin_+1}));
nslambda = zeros(klen,1);
nslambdap = 0;
nslambdaf = 0;
for i=1:klen
  nlambda(i) = length(jlambda{i});
  if i <= ykmin_
    nslambdap = nslambdap + nlambda(i);
    % skiping ykmin_+1
  elseif i > ykmin_+1
    nslambdaf = nslambdaf + nlambda(i);
  end
end
nnslambda = nslambdap+nslambdaf;

% buildind D and E
nd1 = nd+nnslambda+inst_nbr;
d = zeros(nd1,nd1) ;
e = d ;

% variables order:
% z(+1) = [y(+r:+1)' lambda(+r:+1)' u' y(0:-s+1)' lambda(0:-s+1)']
%model dynamics
% future values of forward looking variables
k = find(kstate(:,2) >= ykmin_+2 & kstate(:,3));
d(1:sdyn1,k) = a(nstatic+1:end,kstate(k,3)) ;
% forward looking variables in period t
k1 = find(kstate(:,2) == ykmin_+2);
a1 = eye(sdyn1);
e(1:sdyn1,k1) =  -a1(:,kstate(k1,1)-nstatic1);
% previous values of predetermined variables
k2 = find(kstate(:,2) <= ykmin_+1 & kstate(:,4));
e(1:sdyn1,k2+nslambdap+inst_nbr) = -a(nstatic+1:end,kstate(k2,4)) ;
% purely predetermined variables in current period
k3 = find(kstate(:,2) == ykmin_+1);
k3 = k3(~ismember(kstate(k3,1),kstate(k1,1)));
d(1:sdyn1,k3+nslambdap+inst_nbr) = a1(:,kstate(k3,1)-nstatic1);
% instruments in current period
k4 = sum(kstate(:,2) >= ykmin_+2);
e(1:sdyn1,k4+nslambdap+[1:inst_nbr]) = -b(nstatic+1:end,:);

%first order condition from Lagrangian with respect to y
offsetc = nsfwrd;
a1 = eye(sdyn1);
order_var1 = order_var(nstatic+inst_nbr+1:end);
for i = 1:ykmin_
  kk1 = find(iy_(i,order_var1));
  kk2 = nonzeros(iy_(i,order_var1));
  d(sdyn1+kk1,offsetc+klambda{i}) = bet^(ykmin_+1-i)*aa(ilambda{i},kk2)';
  if i == ykmin_
    e(sdyn1+[1:sdyn1],offsetc+klambda{i}) = -a1(klambda{i},:)';
  end
  offsetc = offsetc+nlambda(i);
end 
offsetc = nsfwrd+nslambdap;
e(sdyn1+[1:sdyn1],offsetc+[1:inst_nbr]) = -RR;
d(sdyn1+[1:sdyn1],nslambdap+inst_nbr+k3) = 2*QQ(:,kstate(k3,1)-nstatic1);
e(sdyn1+[1:sdyn1],k1) =  -2*QQ(:,kstate(k1,1)-nstatic1);
offsetc = nsfwrd+nslambdap+inst_nbr+nspred;
d(sdyn1+[1:sdyn1],offsetc+klambda{ykmin_+1}) = a1(ilambda{ykmin_+1},:)';
for i = ykmin_+2:klen
  kk1 = find(iy_(i,order_var1));
  kk2 = nonzeros(iy_(i,order_var));
  e(sdyn1+kk1,offsetc+klambda{i}) = -bet^(ykmin_+1-i)*aa(ilambda{i},kk2- ...
						  endo_nbr)';
  offsetc = offsetc+nlambda(1);
end

%first order condition from Lagrangian with respect to u
d(2*sdyn1+[1:inst_nbr],nsfwrd+nslambdap+inst_nbr+k3) = RR(kstate(k3,1)-nstatic1,:)';
e(2*sdyn1+[1:inst_nbr],k1) =  -RR(kstate(k1,1)-nstatic1,:)';
e(2*sdyn1+[1:inst_nbr],nsfwrd+nslambdap+[1:inst_nbr]) = -2*UU;
bb = b(nstatic+1:end,:);
d(2*sdyn1+[1:inst_nbr],nsfwrd+nslambdap+inst_nbr+nspred+klambda{ykmin_+1}) = ...
    bb(ilambda{ykmin_+1},:)';
e(2*sdyn1+[1:inst_nbr],nsfwrd+nslambdap-nlambda(ykmin_)+klambda{ykmin_}) = ...
    -bb(ilambda{ykmin_},:)';

%auxiliary equations
if ~isempty(kad)
  for j = 1:size(kad,1)
    if kstate(kad(j),2) < ykmin_+2
      offsetc1 = nslambdap+inst_nbr;
    else
      offsetc1 = 0;
    end
    if kstate(kae(j),2) < ykmin_+2
      offsetc2 = nslambdap+inst_nbr;
    else
      offsetc2 = 0;
    end
    d(2*sdyn1+inst_nbr+j,offsetc1+kad(j)) = 1 ;
    e(2*sdyn1+inst_nbr+j,offsetc2+kae(j)) = 1 ;
  end
end
offsetr = 2*sdyn1+inst_nbr+size(kad,1)+1;
offsetc = nsfwrd;
for i=1:ykmin_-1
  [junk,kk1,kk2] = intersect(jlambda{i},jlambda{i+1});
  for j=1:length(junk)
    d(offsetr,offsetc+nlambda(i)+kk2(j)) = 1;
    e(offsetr,offsetc+kk1(j)) = 1;
    offsetr = offsetr + 1;
  end
  offsetc = offsetc + nlambda(i);
end
[junk,kk1,kk2] = intersect(jlambda{ykmin_},jlambda{ykmin_+2});
for j=1:length(junk)
  d(offsetr,nsfwrd+nslambdap+inst_nbr+nspred+kk2(j)) = 1;
  e(offsetr,offsetc+kk1(j)) = 1;
  offsetr = offsetr + 1;
end
offsetc = nsfwrd + nslambdap + inst_nbr + nspred;
for i=ykmin_+2:klen-1
  [junk,kk1,kk2] = intersect(jlambda{i},jlambda{i+1});
  for j=1:length(junk)
    d(offsetr,offsetc+nlambda(i)+kk2(j)) = 1;
    e(offsetr,offsetc+kk1(j)) = 1;
    offsetr = offsetr + 1;
  end
  offsetc = offsetc + nlambda(i);
end
 

options_ = set_default_option(options_,'qz_criterium',1.000001);
if ~exist('mjdgges')
  % using Chris Sim's routines
  use_qzdiv = 1;
  [ss,tt,qq,w] = qz(e,d);
  [tt,ss,qq,w] = qzdiv(options_.qz_criterium,tt,ss,qq,w);
  ss1=diag(ss);
  tt1=diag(tt);
  warning_state = warning;
  warning off;
  eigval = ss1./tt1 ;
  warning warning_state;
  nba = nnz(abs(eigval) > options_.qz_criterium);
else
  use_qzdiv = 0;
  [ss,tt,w,sdim,eigval,info] = mjdgges(e,d,options_.qz_criterium);
  if info & info ~= nd1+2;
    error(['ERROR' info ' in MJDGGES.DLL']);
  end
  nba = nd1-sdim;
end

nyf = nsfwrd+nslambdap+inst_nbr;

if nba ~= nyf;
  disp('WARNING: Blanchard-Kahn conditions are not satisfied. Run CHECK to learn more!');
  disp('Press any key to continue');
  pause
end

np = nd1 - nyf;
n2 = np + 1;
n3 = nyf;
n4 = n3 + 1;
% derivatives with respect to dynamic state variables
% forward variables
gx = -w(1:n3,n2:nd1)'\w(n4:nd1,n2:nd1)';
% predetermined variables
hx = w(1:n3,1:np)'*gx+w(n4:nd1,1:np)';
hx = (tt(1:np,1:np)*hx)\(ss(1:np,1:np)*hx);

if use_qzdiv
  gx = real(gx);
  hx = real(hx);
end

% including Lagrange multipliers in lgy_, order_var and kstate
for i=1:sdyn1;
  temp = ['mult_' int2str(i)];
  lgy_ = strvcat(lgy_,temp);
end

% reordering multipliers, predetermined - both - forward
im = zeros(sdyn1,1);
j = endo_nbr-inst_nbr+1;
for i=setdiff(jlambda{ykmin_+1},jlambda{ykmin_})';
  im(i) = j;
  j = j + 1;
end
for i=intersect(jlambda{ykmin_+1},jlambda{ykmin_})';
  im(i) = j;
  j = j + 1;
end
for i=setdiff(jlambda{ykmin_},jlambda{ykmin_+1})';
  im(i) = j;
  j = j + 1;
end

% adding instruments and multipliers to kstate
kstate(:,1) = kstate(:,1)-inst_nbr;
kstate = [kstate(1:nsfwrd,:);zeros(nslambdap+inst_nbr,4);kstate(nsfwrd+1:end,:);zeros(nslambdaf,4)];
offsetr = nsfwrd;
for i=1:ykmin_
  kstate(offsetr+[1:nlambda(i)],1:2) = [im(jlambda{i}) (klen-i+2)* ...
		    ones(nlambda(i),1)];
  offsetr = offsetr+nlambda(i);
end
kstate(offsetr+[1:inst_nbr],1:2) = ...
    [endo_nbr+nlambda(ykmin_+1)-inst_nbr+[1:inst_nbr] (ykmin_+2)*ones(inst_nbr,1)];
offsetr = nsfwrd + nslambdap + inst_nbr +nspred;
m = ykmin_+1;
for i=ykmin_+2:klen
  kstate(offsetr+[1:nlambda(i)],1:2) = [im(jlambda{i}) m* ...
		    ones(nlambda(i),1)];
  offsetr = offsetr+nlambda(i);
  m = m - 1;
end
disp(kstate);
%lead variables actually present in the model
% derivatives with respect to exogenous variables
if exo_nbr
  n1 = find(kstate(:,2) > ykmin_+1);
  ghu = -(d*[zeros(nd1,nsfwrd+nslambdap+inst_nbr) [gx;eye(nspred+nslambdaf)]]-e* ...
	  [[eye(nsfwrd+nslambdap+inst_nbr); zeros(nspred+nslambdaf, ...
						   nsfwrd+nslambdap+inst_nbr)] zeros(nd1,nspred+nslambdaf)])\[fu(nstatic+1:end,:); zeros(size(d,1)-size(fu,1)+nstatic,exo_nbr)];
end

nrgx = size(gx,1);
k1 = find((kstate(:,2) == ykmin_+1));
k2 = find((kstate(:,2) == ykmin_+2));
[junk,k3] = setdiff(kstate(k2,1),kstate(k1,1));
dr.ghx = [hx(k1-nrgx,:); gx(k2(k3),:)]; 
dr.ghu = ghu([k1; k2(k3)],:); 

% static variables
if nstatic > 0
  j3 = nonzeros(kstate(:,3));
  j4  = find(kstate(:,3));
  temp = -a(1:nstatic,j3)*gx(j4,:)*hx;
  temp = temp + b(1:nstatic,:)*gx(nsfwrd+[1:inst_nbr],:);
  j5 = find(kstate(n4:nd1,4));
  temp(:,j5) = temp(:,j5)-a(1:nstatic,nonzeros(kstate(:,4)));
  dr.ghx = [temp; dr.ghx];
  temp = -a(1:nstatic,j3)*gx(j4,:)*ghu(nsfwrd+nslambdap+inst_nbr+[1:nspred+nslambdaf],:);
  temp = temp + b(1:nstatic,:)*ghu(nsfwrd+[1:inst_nbr],:);
  temp = temp + fu(1:nstatic,:);
  dr.ghu = [temp; dr.ghu];
  temp = [];
end

dr.ys = [dr.ys; zeros(length(jlambda{ykmin_+1}),1)];
dr.nstatic = nstatic;
dr.npred = npred+nboth+length(jlambda{ykmin_+1});
dr.kstate = kstate;
dr.order_var = [order_var(inst_nbr+[1:nstatic+npred+nboth]);...
	     endo_nbr+[1:nlambda(ykmin_+1)]';...
	     order_var(nstatic+npred+nboth+[1:nfwrd]);...
	     order_var(1:inst_nbr)];
endo_nbr = endo_nbr+nlambda(ykmin_+1);
