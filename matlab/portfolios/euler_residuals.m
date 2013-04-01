function [resids,resids_dyn] = euler_residuals( drn , lind1 , lind2)

    % ind1 : list of indices of returns
    % ind2 : list of indices of stochastic discount factors

    global M_ options_ oo_
    
    
    keep_deltax = options_.portfolios.keep_deltax;
    only_covariances = options_.portfolios.only_covariances;


    [PI,VI,EI] = read_model_indices();

    n_pf_eqs = length(lind1);

    %lli = 
    %state_variables = find(M_.lead_lag_incidence(1,:));
    state_variables = find(M_.portfolios.lead_lag_incidence(1,:));
    
    n_state_variables = length(state_variables);
    pf_s_vars = M_.portfolios.pf_s_vars;
    n_pf_s_vars = length(pf_s_vars);


    resids = zeros(n_pf_eqs,1);
    resids_dyn = zeros(n_pf_eqs, n_state_variables);
    
    Me = M_.Sigma_e(:);

    for i = 1:length(lind1)

        ind1 = lind1(i);
        ind2 = lind2(i);
        
        res = 0;
        %res =  drn.ys(ind1) .* drn.ys(ind2);
        if keep_deltax == 1
            res = res + drn.ys(ind1) .* (drn.ghuu(ind2,:) * M_.Sigma_e(:));
            res = res + drn.ys(ind2) .* (drn.ghuu(ind1,:) * M_.Sigma_e(:));
        end
        
        res = res + diag( drn.ghu(ind1,:) * M_.Sigma_e * drn.ghu(ind2,:)');

        resids(i) = res;

        if (nargout >= 1) && options_.order ==2
          R2 = drn.ghu(ind1,:);
          E2 = drn.ghu(ind2,:);
  
          tR5 = drn.ghxu(ind1,:);
          R5 = reshape(tR5,M_.exo_nbr,n_state_variables);
    
          tE5 = drn.ghxu(ind2,:);
          E5 = reshape(tE5,M_.exo_nbr,n_state_variables);
          
  %         oo_.portfolios.R2 = R2;
  %         oo_.portfolios.E2 = E2;
  %         oo_.portfolios.R5 = R5;
  %         oo_.portfolios.E5 = E5;
  
          % variation in expected returns
          d_rx = drn.ys(ind1) * drn.ghx(ind2,:) + drn.ys(ind2) * drn.ghx(ind1,:);
          % variation in covariances
          d_cov = R2 * M_.Sigma_e * E5 + E2 * M_.Sigma_e * R5;
  
          if keep_deltax == 1
              % second order central term
              d_central_term = drn.ghuu(ind1,:) * Me * drn.ghx(ind2,:) +  drn.ghuu(ind2,:) * Me * drn.ghx(ind1,:);
          else
              d_central_term = 0;
          end

          if only_covariances == 1
              resids_dyn(i,:) =  d_cov;
          else
              resids_dyn(i,:) = d_rx + d_central_term + d_cov;
          end
        
        end
    end

          %n_pf_s_indices = M_.lead_lag_incidence(1,pf_s_vars);
          %resids_dyn = resids_dyn(:,n_pf_s_indices);

%        if length(options_.portfolios.pf_vars) == 1
%            ind1 = lind1(1);
%            ind2 = lind2(1);
            
            %%%%%% uncommment this part !!!
%            xis = options_.portfolios.xis;
            %a = drn.ghu(ind1,:) * M_.Sigma_e * drn.ghu(ind2,:)';
%            d = drn.ghu(ind2,xis) * (  drn.ghu(ind1,:) * M_.Sigma_e * drn.ghu(ind1,:)' );
%            d = d - drn.ghu(ind1,xis) * (  drn.ghu(ind2,:) * M_.Sigma_e * drn.ghu(ind2,:)' );
%              d = - d;
            %ds_der = -d;            
            %ds_der_dyn =  ( -d * ones(n_pf_eqs, n_state_variables) ) ;
%        end

    resids = resids  * 1000;
    resids_dyn = resids_dyn * 1000;
