% the beginning and the end of this function may be adapted by the userx
function [z,vx]=osr_obj(x,params,weights);
  global ys_ Sigma_e_ endo_nbr exo_nbr optimal_Q_ it_ ykmin_
  
  % set parameters of the policiy rule
  np = size(params,1);
  for i=1:np
    assignin('base',deblank(params(i,:)),x(i))
  end
  
  % don't change below until the part where the loss function is computed
  it_ = ykmin_+1;
  dr_ = resol(ys_,1,0,1);
  nstatic = dr_.nstatic;
  npred = dr_.npred;
  ghx = dr_.ghx;
  ghu = dr_.ghu;
  order = dr_.order_var;
  k=[nstatic+1:nstatic+npred]';
  vx1 = ghu(k,:)*Sigma_e_*ghu(k,:)';
  
  % compute variance of predetermined variables
  vx = (eye(npred*npred)-kron(ghx(k,:),dr_.ghx(k,:)))\vx1(:);
  vx=reshape(vx,npred,npred);

  % compute variance of all variables
  if endo_nbr > npred
    qx = eye(npred);
    qu = zeros(npred,exo_nbr);
    if nstatic > 0
      qx = [ghx(1:nstatic,:);qx];
      qu = [ghu(1:nstatic,:);qu];
    end
    if endo_nbr > nstatic+npred
      qx = [qx;ghx(nstatic+npred+1:end,:)];
      qu = [qu;ghu(nstatic+npred+1:end,:)];
    end
    vx = qx*vx*qx'+qu*Sigma_e_*qu';
  end
  % end of the non touch region of the program
  
  % computes the loss function
  weights = weights(order,order);
  z = weights(:)'*vx(:);
  














