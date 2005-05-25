function yf=forcst2(ym,horizon,dr,n)
  global Sigma_e_ endo_nbr options_
  
  yf = zeros(horizon,endo_nbr,n);
  [A,B] = kalman_transition_matrix(dr);
  
  v0 = B*Sigma_e_*B';
  v = v0;
  for i=1:horizon
eig(v)
[c,p]=chol(v)
    %    c = chol(v);
%    yf(i,:,:) = randn(endo_nbr,n)*c+repmat(ym,1,n);
    if i < horizon;
      v = A*v*A'+v0;
    end
  end
  pause