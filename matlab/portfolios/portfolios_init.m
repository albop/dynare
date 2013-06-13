function [] = portfolios_init()

    % reads the model to setup all model-dependent fields
    
    % all informations/indices are relative to the original model

    global oo_ M_ options_
    
    [res,jac,hes] = compute_derivatives(oo_.steady_state);
    
    if nnz(hes)==0
        error('dynare:portfolios','Second order model not available.')
        %throw(MException());
    end
    
    [res0] = compute_derivatives(oo_.steady_state*0);
    
    lli = M_.lead_lag_incidence;
    
    %% check type of portfolio model
    
    % is there first order indeterminacy ?
    r = rank(jac);
    s = size(jac,1);
    
    
    pf_eqs = select_from_table( M_.equations_tags , 'portfolio' );
    n_pf_eqs = length(pf_eqs);
    
    if n_pf_eqs ~= s - r
      error(['Only ' num2str(s -r)  ' moment equation(s) found. Expecting ' num2str(n_pf_eqs) ] );
    end
    
    %pf_hes = hes(pf_eqs,:);
    
    n_z = size(jac,2); 
    
    
    %% find variables appearing in the moment equations_tags
    
    mom_vars = zeros(n_pf_eqs,2);
    mom_var_names = cell(n_pf_eqs,2);
    
    for i = 1:n_pf_eqs
      inds = find( hes(pf_eqs(i),:) == 1);
      if length(inds) ~= 2
          error(['Equation ' num2str(pf_eqs(i)) ' was not recognized as a moment condition'])
      end  
      r = mod( inds(1) , n_z );
      q = (inds(1) - r) / n_z;
      
      i1 = find( lli(3,:) == q+1 ); % index of first variablepf_hes = hes(pf_eqs);
      i2 = find( lli(3,:) ==  r ); % index of second variable
      mom_vars(i,1) = i1;
      mom_vars(i,2) = i2;
      mom_var_names(i,1) = cellstr( M_.endo_names(i1,:) );
      mom_var_names(i,2) = cellstr( M_.endo_names(i2,:) );
    end
    
    %% Now look for portfolio variables
    
    pf_vars = zeros(n_pf_eqs,1);
    pf_var_names = cell(n_pf_eqs,1);
    
    for i = 1:n_pf_eqs
      temp =  M_.equations_tags( cell2mat(M_.equations_tags(:,1)) == pf_eqs(i) ,2:3 );
      s = temp( strmatch('portfolio',temp(:,1),'exact'), 2 );
      s = cell2mat(s);
      pf_vars(i) = strmatch( s, M_.endo_names, 'exact');
      pf_var_names(i) = cellstr(s);
    end

    
	%% construct result struct    
    
    M_.portfolios = struct;

    M_.portfolios.pf_var_names = pf_var_names;
    M_.portfolios.mom_var_names = mom_var_names;    
    M_.portfolios.mom_vars = mom_vars;
    
    % to be removed ...
    ll = M_.lead_lag_incidence(1,:);
    ll(pf_vars) = 0;
    s_vars = find(M_.lead_lag_incidence(1,:));
    pf_s_vars = find(ll);
    M_.portfolios.pf_s_vars = pf_s_vars;
    
    %portfolio_equations = zeros(M_.endo_nbr,1);
    %portfolio_equations(pf_eqs) = 1;

    M_.portfolios.pf_vars = pf_vars;
    M_.portfolios.pf_eqs = pf_eqs;
    
    non_pf_vars = 1:M_.endo_nbr;
    non_pf_vars(pf_vars) = 0;
    non_pf_vars = find(non_pf_vars);
    M_.portfolios.non_pf_vars = non_pf_vars;
    
    non_pf_eqs = 1:M_.endo_nbr;
    non_pf_eqs(pf_eqs) = 0;
    non_pf_eqs = find(non_pf_eqs);
    M_.portfolios.non_pf_eqs = non_pf_eqs;
    
    M_.portfolios.aux_vars = [];
    M_.portfolios.aux_eqs = [];
    
    M_.portfolios.pf_vars_mod = pf_vars; % these are the portfolio equations in the modfile
    M_.portfolios.pf_eqs_mod = pf_eqs;


        
    % experimental : change lead_lag_incidence matrix : not used yet
    lli_new = M_.lead_lag_incidence;
    lli_new(3,pf_vars) = lli_new(2,pf_vars);
    lli_new(2,pf_vars) = lli_new(1,pf_vars);
    lli_new(1,pf_vars) = 0;
    M_.portfolios.lead_lag_incidence_non_contiguous = lli_new;
    %M_.lead_lag_incidence = lli_new;
    
    lli_old = M_.lead_lag_incidence;
    
    lli = (M_.portfolios.lead_lag_incidence_non_contiguous)';
    order = 1:(size(jac,2));
    order(1:max(max(lli))) = [ lli(find(lli(:))) ];
    lli(find(lli)) = 1:length(find(lli));
    lli = lli';

    M_.portfolios.dyn_var_order = order;
    M_.portfolios.lead_lag_incidence = lli;


    M_.portfolios.pricing_kernels = mom_vars(:,1);   % ?
    M_.portfolios.returns = mom_vars(:,2);           % ?
