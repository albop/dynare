% Copyright (C) 2001 Michel Juillard
%
function dr=dr1(iorder,dr,check)

global jacobia_ iy_ ykmin_ ykmax_ gstep_ exo_nbr exo_det_nbr endo_nbr
global ex_ ex_det_ valf_ it_ exe_ exe_det_ xkmin_ xkmax_ ys_ stdexo_
global fname_ means_ Sigma_e_ lgy_
global eigval options_ M_

options_ = set_default_option(options_,'loglinear',0);

xlen = xkmax_ + xkmin_ + 1;
klen = ykmin_ + ykmax_ + 1;
iyv = iy_';
iyv = iyv(:);
iyr0 = find(iyv) ;
it_ = ykmin_ + 1 ;


if exo_nbr == 0
  exe_ = [] ;
end

if ~ all(iy_(ykmin_+1,:) > 0)
  error ('Error in model specification: some variables don"t appear as current') ;
end

if ~check
  if xlen > 1
    error (['SS: stochastic exogenous variables must appear only at the' ...
	    ' current period. Use additional endogenous variables']) ;
  end
end

if (exo_det_nbr > 0) & (ykmin_ > 1 | ykmax_ > 1)
  error(['Exogenous deterministic variables are currently only allowed in' ...
	 ' models with leads and lags on only one period'])
end

dr=set_state_space(dr);
kstate = dr.kstate;
kad = dr.kad;
kae = dr.kae;
nstatic = dr.nstatic;
nfwrd = dr.nfwrd;
npred = dr.npred;
nboth = dr.nboth;
order_var = dr.order_var;
nd = size(kstate,1);

sdyn = endo_nbr - nstatic;


tempex = ex_;

it_ = ykmin_ + 1;
z = repmat(dr.ys,1,klen);
z = z(iyr0) ;
%jacobia_=real(diffext('ff1_',[z; exe_])) ;
jacobia_=real(jacob_a('ff1_',[z(:); exe_])) ;

ex_ = tempex ;
tempex = [];

nz = size(z,1);
k1 = iy_(find([1:klen] ~= ykmin_+1),:);
b = jacobia_(:,iy_(ykmin_+1,order_var));
a = b\jacobia_(:,nonzeros(k1')); 
if any(isinf(a(:)))
  error('DR1: the model doesn''t determine current variables uniquely')
end
if exo_nbr
  fu = b\jacobia_(:,nz+1:end);
end

if ykmax_ == 0;  % backward model
  dr.ghx = -a;
  if exo_nbr
    dr.ghu = -fu;
  end
  dr.eigval = eig(transition_matrix(dr));
  dr.rank = 0;
  return;
end

% buildind D and E
d = zeros(nd,nd) ;
e = d ;

k = find(kstate(:,2) >= ykmin_+2 & kstate(:,3));
d(1:sdyn,k) = a(nstatic+1:end,kstate(k,3)) ;
k1 = find(kstate(:,2) == ykmin_+2);
a1 = eye(sdyn);
e(1:sdyn,k1) =  -a1(:,kstate(k1,1)-nstatic);
k = find(kstate(:,2) <= ykmin_+1 & kstate(:,4));
e(1:sdyn,k) = -a(nstatic+1:end,kstate(k,4)) ;
k2 = find(kstate(:,2) == ykmin_+1);
k2 = k2(~ismember(kstate(k2,1),kstate(k1,1)));
d(1:sdyn,k2) = a1(:,kstate(k2,1)-nstatic);

if ~isempty(kad)
  for j = 1:size(kad,1)
    d(sdyn+j,kad(j)) = 1 ;
    e(sdyn+j,kae(j)) = 1 ;
  end
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
  if info & info ~= nd+2;
%    error(['ERROR' info ' in MJDGGES.DLL']);
  end
  nba = nd-sdim;
end

nyf = sum(kstate(:,2) > ykmin_+1);

if check
  dr.rank = rank(w(1:nyf,nd-nyf+1:end));
  dr.eigval = eig(e,d);
  return
end

if nba ~= nyf;
  disp('WARNING: Blanchard-Kahn conditions are not satisfied. Run CHECK to learn more!');
  disp('Press any key to continue');
  pause
end

np = nd - nyf;
n2 = np + 1;
n3 = nyf;
n4 = n3 + 1;
% derivatives with respect to dynamic state variables
% forward variables
gx = -w(1:n3,n2:nd)'\w(n4:nd,n2:nd)';
% predetermined variables
hx = w(1:n3,1:np)'*gx+w(n4:nd,1:np)';
hx = (tt(1:np,1:np)*hx)\(ss(1:np,1:np)*hx);

k1 = find(kstate(n4:nd,2) == ykmin_+1);
k2 = find(kstate(1:n3,2) == ykmin_+2);
dr.ghx = [hx(k1,:); gx(k2(nboth+1:end),:)];

%lead variables actually present in the model
j3 = nonzeros(kstate(:,3));
j4  = find(kstate(:,3));
% derivatives with respect to exogenous variables
if exo_nbr
  a1 = eye(endo_nbr);
  aa1 = [];
  if nstatic > 0
    aa1 = a1(:,1:nstatic);
  end
  dr.ghu = -[aa1 a(:,j3)*gx(j4,1:npred)+a1(:,nstatic+1:nstatic+ ...
						  npred) a1(:,nstatic+npred+1:end)]\fu;
end

% static variables
if nstatic > 0
  temp = -a(1:nstatic,j3)*gx(j4,:)*hx;
  j5 = find(kstate(n4:nd,4));
  temp(:,j5) = temp(:,j5)-a(1:nstatic,nonzeros(kstate(:,4)));
  dr.ghx = [temp; dr.ghx];
  temp = [];
end

if options_.loglinear == 1
    k = find(dr.kstate(:,2) <= ykmin_+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var;

    dr.ghx = repmat(1./dr.ys(k1),1,size(dr.ghx,2)).*dr.ghx.* ...
	     repmat(dr.ys(k1(klag(:,1)))',size(dr.ghx,1),1);
    dr.ghu = repmat(1./dr.ys(k1),1,size(dr.ghu,2)).*dr.ghu;
end

%necessary when using Sims' routines
if use_qzdiv
  gx = real(gx);
  hx = real(hx);
  dr.ghx = real(dr.ghx);
  dr.ghu = real(dr.ghu);
end

%exogenous deterministic variables
if exo_det_nbr > 1
  f1 = sparse(jacobia_(:,nonzeros(iy_(ykmin_+2:end,order_var))));
  f0 = sparse(jacobia_(:,nonzeros(iy_(ykmin_+1,order_var))));
  fudet = sparse(jacobia_(:,nz+exo_nbr+1:end));
  M1 = inv(f0+[zeros(endo_nbr,nstatic) f1*gx zeros(endo_nbr,nyf-nboth)]);
  M2 = M1*f1;
  dr.ghud = cell(M_.ex_det_length,1);
  dr.ghud{1} = -M1*fudet;
  for i = 2:M_.ex_det_length
    dr.ghud{i} = -M2*dr.ghud{i-1}(end-nyf+1:end,:);
  end
end

if iorder == 1
  return
end

% Second order
tempex = ex_ ;

%hessian = real(hessext('ff1_',[z; exe_]))' ;
kk = flipud(cumsum(flipud(iy_(ykmin_+1:end,order_var)),1));
if ykmin_ > 0
  kk = [cumsum(iy_(1:ykmin_,order_var),1); kk];
end
kk = kk';
kk = find(kk(:));
nk = size(kk,1)+exo_nbr;
k1 = iy_(:,order_var);
k1 = k1';
k1 = k1(:);
k1 = k1(kk);
k2 = find(k1);
kk1(k1(k2)) = k2;
kk1 = [kk1 length(k1)+1:length(k1)+exo_nbr];
kk = reshape([1:nk^2],nk,nk);
kk1 = kk(kk1,kk1);
hessian = zeros(endo_nbr,nk^2);
hessian(:,kk1(:)) = real(hessian_sparse('ff1_',[z; exe_])) ;

ex_ = tempex ;
clear tempex

n1 = 0;
n2 = np;
zx = zeros(np,np);
for i=2:ykmin_+1
  k1 = sum(kstate(:,2) == i);
  zx(n1+1:n1+k1,n2-k1+1:n2)=eye(k1);
  n1 = n1+k1;
  n2 = n2-k1;
end
kk = flipud(cumsum(flipud(iy_(ykmin_+1:end,order_var)),1));
k0 = [1:endo_nbr];
gx1 = dr.ghx;
zx = [zx; gx1];
for i=1:ykmax_
  k1 = find(kk(i+1,k0) > 0);
  gx1 = gx1(k1,:)*hx;
  zx = [zx; gx1];
  kk = kk(:,k0);
  k0 = k1;
end
zx=[zx; zeros(exo_nbr,np)];
[n1,n2] = size(zx);
if n1*n1*n2*n2 > 1e7
  rhs = zeros(endo_nbr,n2*n2);
  k1 = 1;
  for i1 = 1:n2
      for i2 = 1:n2
	rhs(:,k1) = hessian*kron(zx(:,i1),zx(:,i2));
	k1 = k1 + 1; 
      end
  end
else
  rhs = hessian*kron(zx,zx);
end
rhs = -rhs;

%lhs
n = endo_nbr+sum(kstate(:,2) > ykmin_+1 & kstate(:,2) < ykmin_+ykmax_+1);
A = zeros(n,n);
B = zeros(n,n);
kp = find(kstate(:,2) == ykmin_+ykmax_+1);
gx2 = gx(kp,:);
% variables with the highest lead
k1 = kstate(find(kstate(:,2) == ykmin_+ykmax_+1),3)+endo_nbr;
% Jacobian with respect to the variables with the highest lead
B1 = jacobia_(:,k1);
A(1:endo_nbr,1:endo_nbr) = jacobia_(:,iy_(ykmin_+1,order_var));
A(1:endo_nbr,nstatic+1:nstatic+npred)=...
    A(1:endo_nbr,nstatic+1:nstatic+npred)+B1*gx2(:,1:npred);
B(1:endo_nbr,end-length(k1)+1:end) = B1;
offset = endo_nbr;
k0 = find(kstate(:,2) == ykmin_+1);
for i=1:ykmax_-1
  k1 = find(kstate(:,2) == ykmin_+i+1);
  A(offset+[1:length(k1)],1:npred) = -gx(k1,1:npred); 
  [k2,junk,k3] = find(kstate(k1,3));
  A(1:endo_nbr,offset+k2) = jacobia_(:,k3);
  n1 = length(k1);
  A(offset+1:offset+n1,offset+1:offset+n1) = eye(n1);
  n0 = length(k0);
  E = eye(n0);
%  if i == 1
%    [junk,junk,k4]=intersect(kstate(k1,1),[1:endo_nbr]);
%  else
    [junk,junk,k4]=intersect(kstate(k1,1),kstate(k0,1));
%  end
  i1 = offset-n0+n1;
  B(offset+1:offset+n1,i1:offset) = -E(k4,:);
  k0 = k1;
  offset = offset + length(k1);
end
C = kron(hx,hx);
%C = hx;
D = [rhs; zeros(n-endo_nbr,size(rhs,2))];
x0 = sylvester3(A,B,C,D);
dr.ghxx = sylvester3a(x0,A,B,C,D);
%dr.ghxx = gensylv(2,A,B,C,D);

%ghxu
%rhs
hu = dr.ghu(nstatic+1:nstatic+npred,:);
zu=[zeros(np,exo_nbr);dr.ghu;gx(:,1:npred)*hu;eye(exo_nbr)];
%kk = reshape([1:np*np],np,np);
%kk = kk(1:npred,1:npred);
%rhs = -hessian*kron(zx,zu)-f1*dr.ghxx(end-nyf+1:end,kk(:))*kron(hx(1:npred,:),hu(1:npred,:));
if n1*n1*n2*exo_nbr > 1e7
  rhs = zeros(endo_nbr,n2*exo_nbr);
  k1 = 1;
  for i1 = 1:n2
      for i2 = 1:exo_nbr
	rhs(:,k1) = hessian*kron(zx(:,i1),zu(:,i2));
	k1 = k1 + 1; 
      end
  end
else
  rhs = hessian*kron(zx,zu);
end
nyf1 = sum(kstate(:,2) == ykmin_+2);
hu1 = [hu;zeros(np-npred,exo_nbr)];
rhs = -[rhs; zeros(n-endo_nbr,size(rhs,2))]-B*dr.ghxx*kron(hx,hu1);

%lhs
dr.ghxu = A\rhs;

%ghuu
%rhs
kk = reshape([1:np*np],np,np);
kk = kk(1:npred,1:npred);
if n1*n1*exo_nbr*exo_nbr > 1e7
  rhs = zeros(endo_nbr,exo_nbr*exo_nbr);
  k1 = 1;
  for i1 = 1:exo_nbr
      for i2 = 1:exo_nbr
	rhs(:,k1) = hessian*kron(zu(:,i1),zu(:,i2));
	k1 = k1 + 1; 
      end
  end
else
  rhs = hessian*kron(zu,zu);
end

rhs = -[rhs; zeros(n-endo_nbr,size(rhs,2))]-B*dr.ghxx*kron(hu1,hu1);

%lhs
dr.ghuu = A\rhs;

dr.ghxx = dr.ghxx(1:endo_nbr,:);
dr.ghxu = dr.ghxu(1:endo_nbr,:);
dr.ghuu = dr.ghuu(1:endo_nbr,:);


% dr.ghs2
% derivatives of F with respect to forward variables
% reordering predetermined variables in diminishing lag order
o_pred_dimin = [];
k0 = find(kstate(:,2) <= ykmin_+1);
for i=1:ykmin_
  o_pred_dimin = [o_pred_dimin; find(kstate(k0,2) == i+1)];
end
%reordering transition matrix
hx1 = hx(o_pred_dimin,o_pred_dimin);
%reordering second order transition matrix
nstate = size(hx,1);
k1 = reshape([1:nstate^2],nstate,nstate);
k1 = k1(o_pred_dimin,o_pred_dimin);
ghxx = dr.ghxx(:,k1(:));
hxx = dr.ghxx(nstatic+[1:npred],k1(:));
% computing selection index for hessian: k1
kh = reshape([1:nk^2],nk,nk);
ghx = dr.ghx(:,o_pred_dimin);
guu = dr.ghuu(end-nyf1+1:end,:);
gu = dr.ghu(end-nyf1+1:end,:);
[n1,n2] = size(gu);
kp = sum(kstate(:,2) <= ykmin_+1);
E1 = zeros(kp,npred);
E1(end-npred+1:end,:) = eye(npred);
H = E1;
%offset = size(Hs,1)-npred;
kk = find(kstate(:,2) == ykmin_+1);
for i=1:ykmin_-1
  kk1 = find(kstate(:,2) == ykmin_-i+1);
  [junk,kk2] = intersect(kstate(kk,1),kstate(kk1,1));
  E = eye(length(kk));
%  Hs(offset-length(kk1)+1:offset,offset+[1:length(kk)]) = E(kk2,:);
end
O1 = zeros(endo_nbr,nstatic);
O2 = zeros(endo_nbr,endo_nbr-nstatic-npred);
LHS = jacobia_(:,iy_(ykmin_+1,order_var));
RHS = zeros(endo_nbr,exo_nbr^2);
kk = find(kstate(:,2) == ykmin_+2);
gu = dr.ghu; 
guu = dr.ghuu; 
kk = find(kstate(:,2) == ykmin_+1);
Gu = zeros(size(hx1,1),exo_nbr);
Gu(end-length(kk)+1:end,:) = dr.ghu(kstate(kk,1),:);
Guu = zeros(size(hx1,1),exo_nbr^2);
Guu(end-length(kk)+1:end,:) = dr.ghuu(kstate(kk,1),:);
E = eye(endo_nbr);
iy_ordered = iy_(:,order_var)';
iy_ordered = iy_ordered(:);
k1 = find(iy_ordered);
iy_ordered(k1) = [1:length(k1)]';
iy_ordered =reshape(iy_ordered,endo_nbr,ykmin_+ykmax_+1)';
for i=1:ykmax_
  for j=i:ykmax_
    k2 = iy_(ykmin_+i+1,order_var);
    [junk,k2a,k2] = find(k2);
    k2b = nonzeros(iy_ordered(ykmin_+j+1,:));
    RHS = RHS + jacobia_(:,k2)*guu(k2a,:)+hessian(:,kh(k2b,k2b))* ...
	  kron(gu(k2a,:),gu(k2a,:));
  end

  % LHS
  k2 = iy_(ykmin_+i+1,order_var);
  [junk,k2a,k2] = find(k2);
  LHS = LHS + jacobia_(:,k2)*(E(k2a,:)+[O1(k2a,:) ghx(k2a,:)*H O2(k2a,:)]);
  
  if i == ykmax_ 
    break
  end
  
  kk = find(kstate(:,2) == ykmin_+i+1);
  gu = ghx*Gu;
  GuGu = kron(Gu,Gu);
  guu = ghx*Guu+ghxx*GuGu;
  Gu = hx1*Gu;
  Guu = hx1*Guu;
  Guu(end-npred+1:end,:) = Guu(end-npred+1:end,:) + hxx*GuGu;

  H = E1 + hx1*H;
end
RHS = RHS*Sigma_e_(:);
dr.fuu = RHS;
RHS = -RHS-dr.fbias;
dr.ghs2 = LHS\RHS;

% deterministic exogenous variables
if exo_det_nbr > 0
  hud = dr.ghud{1}(nstatic+1:nstatic+npred,:);
  zud=[zeros(np,exo_det_nbr);dr.ghud{1};gx(:,1:npred)*hud;zeros(exo_nbr,exo_det_nbr);eye(exo_det_nbr)];
  R1 = hessian*kron(zx,zud);
  dr.ghxud = cell(M_.ex_det_length,1);
  kf = [endo_nbr-nyf+1:endo_nbr];
  kp = nstatic+[1:npred];
  dr.ghxud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{1}(kp,:)));
  Eud = eye(exo_det_nbr);
  for i = 2:M_.ex_det_length
    hudi = dr.ghud{i}(kp,:);
    zudi=[zeros(np,exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi;zeros(exo_nbr+exo_det_nbr,exo_det_nbr)];
    R2 = hessian*kron(zx,zudi);
    dr.ghxud{i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hx,Eud)+dr.ghxx(kf,:)*kron(dr.ghx(kp,:),dr.ghud{i}(kp,:)))-M1*R2;
  end
  R1 = hessian*kron(zu,zud);
  dr.ghudud = cell(M_.ex_det_length,1);
  kf = [endo_nbr-nyf+1:endo_nbr];

  dr.ghuud{1} = -M1*(R1+f1*dr.ghxx(kf,:)*kron(dr.ghu(kp,:),dr.ghud{1}(kp,:)));
  Eud = eye(exo_det_nbr);
  for i = 2:M_.ex_det_length
    hudi = dr.ghud{i}(kp,:);
    zudi=[zeros(np,exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi;zeros(exo_nbr+exo_det_nbr,exo_det_nbr)];
    R2 = hessian*kron(zu,zudi);
    dr.ghuud{i} = -M2*dr.ghxud{i-1}(kf,:)*kron(hu,Eud)-M1*R2;
  end
  R1 = hessian*kron(zud,zud);
  dr.ghudud = cell(M_.ex_det_length,M_.ex_det_length);
  dr.ghudud{1,1} = -M1*R1-M2*dr.ghxx(kf,:)*kron(hud,hud);
  for i = 2:M_.ex_det_length
    hudi = dr.ghud{i}(nstatic+1:nstatic+npred,:);
    zudi=[zeros(np,exo_det_nbr);dr.ghud{i};gx(:,1:npred)*hudi+dr.ghud{i-1}(kf,:);zeros(exo_nbr+exo_det_nbr,exo_det_nbr)];
    R2 = hessian*kron(zudi,zudi);
    dr.ghudud{i,i} = -M2*(dr.ghudud{i-1,i-1}(kf,:)+...
			  2*dr.ghxud{i-1}(kf,:)*kron(hudi,Eud) ...
			  +dr.ghxx(kf,:)*kron(hudi,hudi))-M1*R2;
    R2 = hessian*kron(zud,zudi);
    dr.ghudud{1,i} = -M2*(dr.ghxud{i-1}(kf,:)*kron(hud,Eud)+...
			  dr.ghxx(kf,:)*kron(hud,hudi))...
	-M1*R2;
    for j=2:i-1
      hudj = dr.ghud{j}(kp,:);
      zudj=[zeros(np,exo_det_nbr);dr.ghud{j};gx(:,1:npred)*hudj;zeros(exo_nbr+exo_det_nbr,exo_det_nbr)];
      R2 = hessian*kron(zudj,zudi);
      dr.ghudud{j,i} = -M2*(dr.ghudud{j-1,i-1}(kf,:)+dr.ghxud{j-1}(kf,:)* ...
			    kron(hudi,Eud)+dr.ghxud{i-1}(kf,:)* ...
			    kron(hudj,Eud)+dr.ghxx(kf,:)*kron(hudj,hudi))-M1*R2;
    end
    
  end
end
% 01/08/2001 MJ put return on iorder == 1 after defining dr.kstate and dr.kdyn
% 01/17/2001 MJ added dr.delta_s: correction factor for order = 2
% 01/21/2001 FC correction of delta_s for more than 1 shock
% 01/23/2001 MJ suppression of redundant sum() in delta_s formula
% 02/22/2001 MJ stderr_ replaced by Sigma_e_
% 08/24/2001 MJ changed the order of the variables, separates static
%               variables and handles only one instance of variables both
%               in lead and lag
% 08/24/2001 MJ added sigma to state variables as in Schmitt-Grohe and
%               Uribe (2001)
% 10/20/2002 MJ corrected lags on several periods bug
% 10/30/2002 MJ corrected lags on several periods bug on static when some
%               intermediary lags are missing
% 12/08/2002 MJ uses sylvester3 to solve for dr.ghxx
% 01/01/2003 MJ added dr.fbias for iterative for dr_algo == 1
% 02/21/2003 MJ corrected bug for models without lagged variables
% 03/02/2003 MJ fixed second order for lag on several periods
% 05/21/2003 MJ add check call argument and make computation for CHECK
% 06/01/2003 MJ added a test for xkmax_ > 1 and order > 1
% 08/28/2003 MJ corrected bug in computation of 2nd order (ordering of
%               forward variable in LHS)
% 08/29/2003 MJ use Sims routine if mjdgges isn't available

   





