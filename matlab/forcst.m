function [yf,var_yf]=forcst(dr,y0,k,m)
  global endo_nbr exo_nbr ykmin_ Sigma_e_ ex_ options_
  
  options_.periods = k;
  make_ex_;
  yf = simult_(y0,dr,ex_(1:k,:),1);

  nstatic = dr.nstatic;
  npred = dr.npred;
  j = find(kstate(dr.kae,2) <= ykmin_+1);
  kae = dr.kae(j);
  nh = size(dr.ghx,2);
  hx = dr.ghx(nstatic+1:nstatic+npred,:);
  hu = dr.ghu(nstatic+1:nstatic+npred,:);
  if ~isempty(kae)
    n = length(kae);
    tmp = sparse([1:n]',kae-sum(dr.kstate(:,2)>ykmin_+1),ones(n,1),n,nh);
    hx = [hx; tmp];
    hu = [hu; zeros(n,exo_nbr)];
  end
  gx  = [];
  k2  = [];
  if nstatic > 0
    gx = dr.ghx(1:nstatic,:);
    k2 = [1:nstatic]';
  end
  if size(dr.ghx,1) > nstatic+npred
    gx = [gx; dr.ghx(nstatic+npred+1:end,:)];
    k2 = [k2; [nstatic+npred+1:size(dr.ghx,1)]'];
  end
  
  k1 = dr.order_var([nstatic+1:nstatic+npred]);
  k2 = dr.order_var(k2);
  
  sigma_u = hu*Sigma_e_*hu';
  sigma_y1 = 0;
  var_yf = zeros(k,endo_nbr);
  

  if isempty(k2)
    for i=1:k
      sigma_y1 = sigma_y1+sigma_u;
      var_yf(i,k1) = diag(sigma_y1(1:npred,1:npred))';
      if i == k
	break
      end
      sigma_u = hx*sigma_u*hx';
    end
  else
    for i=1:k
      sigma_y1 = sigma_y1+sigma_u;
      var_yf(i,k1) = diag(sigma_y1(1:npred,1:npred))';
      sigma_y2 = gx*sigma_y1*gx';
      var_yf(i,k2) = diag(sigma_y2)';
      if i == k
	  break
      end
      sigma_u = hx*sigma_u*hx';
    end
  end
  












