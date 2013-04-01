function [dr] = dr_for_given_portfolios(alphas,dalphas)


    global M_ oo_ options_


    %% preamble
    
    [PI,VI,EI] = read_model_indices();
    
    pf_eqs = M_.portfolios.pf_eqs;
    pf_vars = M_.portfolios.pf_vars;
    n_pf_vars = length(pf_vars);

    non_pf_vars = M_.portfolios.non_pf_vars;
    non_pf_eqs = M_.portfolios.non_pf_eqs;

    aux_vars = M_.portfolios.aux_vars;
    aux_eqs = M_.portfolios.aux_eqs;
        
    s_vars = find(M_.lead_lag_incidence(1,:));
    n_s_vars = length(s_vars);

    pf_s_vars = M_.portfolios.pf_s_vars;
    n_pf_s_vars = length(pf_s_vars);
    

    %% replace dalphas by zeros if not supplied
    %if nargin == 1
    %  dalphas = zeros(n_pf_vars,n_pf_s_vars);
    %  
    %end

    if  length(M_.portfolios.aux_vars) == 0
        o_vars = non_pf_vars;
        o_eqs = non_pf_eqs;
    else
        o_vars = M_.portfolios.aux_vars;
        o_eqs = M_.portfolios.aux_eqs;
    end
            
    %% compute derivatives for given portfolio values
    ys = oo_.steady_state;
    ys( pf_vars ) = alphas;
    
    if options_.order == 1
      [res,jac] = compute_derivatives(ys);
    elseif options_.order ==2
      [res,jac,hes] = compute_derivatives(ys);
    end


    lli = M_.portfolios.lead_lag_incidence;
    order = M_.portfolios.dyn_var_order;

    jac = jac(: ,order);


    pf_s_vars = find(lli(1,:));

    jac(pf_eqs,:) = jac(pf_eqs,:) * 0 ;
    jac(pf_eqs,lli(1,pf_s_vars)) = - dalphas;
    jac(pf_eqs,lli(2,pf_vars)) =  jac(pf_eqs,lli(2,pf_vars)) + eye(n_pf_vars);

    if options_.order == 2
        uhes(pf_eqs,:) = hes(pf_eqs,:) * 0 ;
        hes_mask = reshape(1:size(hes,2),size(jac,2),size(jac,2));
        hes_mask = hes_mask(order,order);
        hes_mask = hes_mask(:);
        hes = hes(:,hes_mask);
    end

%    M_new = M_;
    M_new = struct;
    M_new.lead_lag_incidence = lli;

    M_new.exo_nbr = M_.exo_nbr; 
    M_new.endo_nbr = size(lli,2); 
    M_new.maximum_endo_lead = 1; 
    M_new.maximum_endo_lag = 1;
    M_new.exo_det_nbr = M_.exo_det_nbr;

    fwrd_var = find(any(lli(M_new.maximum_endo_lag+2:end,:),1));
    pred_var = find(any(lli(1:M_.maximum_endo_lag,:),1))';
    both_var = intersect(pred_var,fwrd_var);
    pred_var = setdiff(pred_var,both_var);
    fwrd_var = setdiff(fwrd_var,both_var);
    stat_var = setdiff([1:M_new.endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
    
    nboth = length(both_var);
    npred = length(pred_var);
    nfwrd = length(fwrd_var);
    nstatic = length(stat_var);
    order_var = [ stat_var(:); pred_var(:); both_var(:); fwrd_var(:)];
    inv_order_var(order_var) = (1:M_.endo_nbr);

    M_new.nboth = nboth;
    M_new.npred = npred;
    M_new.nfwrd = nfwrd;
    M_new.nstatic = nstatic;
    M_new.nfwrd = nfwrd;
    M_new.ndynamic = npred + nboth + nfwrd; 
    M_new.Sigma_e = M_.Sigma_e;
    
    dr = struct;
    dr = set_state_space(dr, M_new, options_);

    % number of forward variables in the state vector
    M_new.nsfwrd = sum(dr.kstate(:,2) > M_new.maximum_endo_lag+1);
    % umber of predetermined variables in the state vector
    M_new.nspred = sum(dr.kstate(:,2) <= M_new.maximum_endo_lag+1);



    dr.ys = ys;
    dr.jacobia_ = jac;
    oo_.portfolios.jac = jac;

    dr = dyn_first_order_solver(jac, M_new, dr, options_, 0);

    if options_.order == 2
        dr.hessian = hes;
        oo_.portfolios.hes = hes;
        dr = dyn_second_order_solver(jac, hes,  dr, M_new, options_, 0);
    end

   
end
