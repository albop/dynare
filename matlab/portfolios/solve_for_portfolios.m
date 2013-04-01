%function [alphas,dalphas,dr] = solve_for_portfolios(M_,oo_,options_)
function [M_,oo_,options_] = solve_for_portfolios(M_,oo_,options_)


    original_order = options_.order;

    ind1 = M_.portfolios.pricing_kernels;
    ind2 = M_.portfolios.returns;
    %[PI,VI,EI] = read_model_indices();

    pf_eqs = M_.portfolios.pf_eqs;
    pf_vars = M_.portfolios.pf_vars;
    
    
    n_pf_vars = length(pf_vars);
    n_s_vars = length(find(M_.lead_lag_incidence(1,:)));
    s_vars = find(M_.lead_lag_incidence(1,:));
    pf_s_vars = M_.portfolios.pf_s_vars;
    n_pf_s_vars = length(pf_s_vars);

    options = optimset();
    options.TolFun = 1e-10;
    options.Display = 'iter';
    options.MaxFunEvals = 10000;
    options.Jacobian = 'off';

    
%% static portfolios

    % initial value for portfolios : 

    alpha0 = zeros(n_pf_vars,1);
    dalpha0 = zeros(n_pf_vars, n_pf_s_vars);


    partial_alpha = @(xx) temp_fun_1(xx,dalpha0,ind1,ind2);   

    if isfield(options_.portfolios,'method') && strcmp( options_.portfolios.method , 'grad')
        disp('Using ''gradient method''')
        options_.order = 1;
        x0 = zeros( n_pf_vars,1);
        res0 = partial_alpha(x0);
        jac = jacobianest(partial_alpha,x0);
        alpha1 =  jac' \ (-res0) ;

        oo_.portfolios.method = 'grad';

    end

    if isfield(options_.portfolios,'method') && strcmp( options_.portfolios.method , 'devereux-sutherland')
%    if isfield(options_.portfolios,'method') && strcmp( options_.portfolios.method , 'devereux-sutherland')
        disp('Using ''Devereux and Sutherland'' method');
        disp(' ');

        if length(find(ind1==ind1(1))) ~= length(ind1)
            error('Devereux and Sutherland formulas have not been implemented for more than two agents')
        else
            Delta = ind1(1);
        end
        [PI,VI,EI] = read_model_indices();

        options_.order = 1;
        x0 = zeros( n_pf_vars,1);
        [res0,dr0] = partial_alpha(x0);
        drn = reorder_dr(dr0);


        R1 = drn.ghu(ind2,EI.xi);
        D1 = drn.ghu(Delta,EI.xi); 

        R2 = drn.ghu(ind2,:);
        D2 = drn.ghu(Delta,:);
        
        a = D2 * M_.Sigma_e * R2';
        lam_u = D1 * R2 * M_.Sigma_e * R2';
        %drn.ghu(Delta,EI.xi) * (  drn.ghu( ind2 ,:) * M_.Sigma_e * drn.ghu(ind2,:)' );
        %dd = drn.ghu(ind2,EI.xi) * diag(drn.ghu(Delta,:) * M_.Sigma_e * drn.ghu(ind2,:)') ;
        dd =  R2 * M_.Sigma_e * D2' * R1';
        lam_u = lam_u - dd;
        alpha1 = -lam_u\a';

        ds = struct;
        ds.R1 = R1;
        ds.R2 = R2;
        ds.D1 = D1;
        ds.D2 = D2;

        oo_.portfolios.ds = ds;
        oo_.portfolios.method = 'devereux-sutherland';

        % we check whether result is consistent with moment conditions
        [res,dr] = partial_alpha(alpha1);

        if max(abs(res)) > 1E-10
            error( ['DS method is not applicable. There are non zero residuals : ' num2str(res') ] )
        end
        oo_.portfolios.dr = dr;

    end
    

	disp( options_.portfolios.method );
    if strcmp( options_.portfolios.method , 'iterative')

        disp('Using ''iterative'' method');
        disp(' ');

        options_.order = 1;
        if options_.portfolios.solver == 'fsolve'
            alpha1 = fsolve( partial_alpha , alpha0 , options );
        else
            alpha1 = csolve( partial_alpha , alpha0 , [], options.TolFun, 500);
        end
    end

    disp('Steady-state portfolios : ');
    disp('');
    rownames = cellstr(M_.endo_names(pf_vars,:));
    colnames = {'---', 'constant'};
    r = [colnames; rownames, num2cell(alpha1)];
    disp(r)
    oo_.portfolios.alphas = alpha1;
  
%    if nargout <= 1
%        alphas = alpha1;
%        return;
%    end
    
%% dynamic portfolios

    options_.order = 2;
    
    partial_dalpha = @(xx) temp_fun_2(alpha1,xx,ind1,ind2,n_pf_vars, n_pf_s_vars);

    partial_alpha_all = @(xx) temp_fun_3(xx,ind1,ind2,n_pf_vars, n_pf_s_vars);


    if options_.portfolios.solver == 'fsolve'
        dalpha1 = fsolve(partial_dalpha,dalpha0(:),options);
    else
        dalpha1 = csolve(partial_dalpha,dalpha0(:), [], options.TolFun, 500);
    end
    dalpha1 = reshape(dalpha1,n_pf_vars, n_pf_s_vars);


    disp('Final portfolios : ')
    xx = [alpha1, dalpha1];
    if options_.portfolios.solver == 'fsolve'
	     a_da = fsolve(partial_alpha_all, xx(:), options);
    else
         a_da = csolve(partial_alpha_all, xx(:),  [], options.TolFun, 500);
    end
    a_da = reshape(a_da,n_pf_vars,n_pf_s_vars+1);
    alphas = a_da(:,1);
    dalphas = a_da(:, 2:end);
    rownames = cellstr(M_.endo_names(pf_vars,:));
    colnames = [{'---' 'constant'} cellstr(M_.endo_names(pf_s_vars,:))'];
    r = [colnames; rownames, num2cell([alphas,dalphas])];
    disp(r)   
    %end

	dr = dr_for_given_portfolios(alphas,dalphas);
    
    oo_.portfolios.dr = dr; 
        
    dr = oo_.portfolios.dr;
    
    
    options_.order = original_order;

end

%% auxiliary nested functions

function [res,dr] = temp_fun_1(xx,dalpha,ind1,ind2)
        dr = dr_for_given_portfolios(xx,dalpha);
        drn = reorder_dr( dr );
        [res,junk] = euler_residuals( drn , ind1 , ind2  );
        res = res(:,1);
        res = res;
        res = res;
		%ds_der = -eye(length(res)) * d;
end

function [dres] = temp_fun_2(alpha,xx,ind1,ind2,n_pf_vars, n_pf_s_vars)
        xx1 = reshape(xx,n_pf_vars, n_pf_s_vars);
        dr = dr_for_given_portfolios(alpha,xx1);
        drn = reorder_dr( dr );
        [res,dres1] = euler_residuals( drn, ind1, ind2  );
        dres = dres1(:);
        res = res;
        dres = dres;
        %dres_der = -eye(length(dres))*d;
end

function [out] = temp_fun_3(xx,ind1,ind2,n_pf_vars, n_pf_s_vars)
        xx = reshape(xx,n_pf_vars,n_pf_s_vars+1);
        a = xx(:,1);
        da = xx(:,2:end);
        dr = dr_for_given_portfolios(a,da);
        drn = reorder_dr( dr );
        [res,dres1] = euler_residuals( drn, ind1, ind2  );
        %[res,dres] = euler_residuals( drn, ind1, ind2  );        
        out = [res, dres1];
        out = out(:);
        %out_der = -eye(length(out)) * d;
        out = out;
        %out_der = out_der *1000;
end
