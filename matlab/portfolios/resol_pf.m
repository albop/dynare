function [dr,info]=resol_pf(ys,check_flag)
%function [alphas,dalphas,dr] = resol_pf()

    global M_ oo_ options_
    
    % ys and info are not used (yet)
    % info is set to 0 for now
    info = 0  % temporary

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

    solver_options = options_.solver_options;

    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% static portfolios %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    % initial value for portfolios : 

    alpha0 = zeros(n_pf_vars,1);
    dalpha0 = zeros(n_pf_vars, n_pf_s_vars);

    partial_alpha = @(xx) temp_fun_1(xx,dalpha0,ind1,ind2);   


    if strcmp( options_.portfolios.method , 'devereux-sutherland')
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

        dd =  R2 * M_.Sigma_e * D2' * R1';
        lam_u = lam_u - dd;
        alphas = -lam_u\a';

        ds = struct;
        ds.R1 = R1;
        ds.R2 = R2;
        ds.D1 = D1;
        ds.D2 = D2;

        oo_.portfolios.ds = ds;

        % we check whether result is consistent with moment conditions
        if options_.check_solution
            [res,dr] = partial_alpha(alphas);
            if max(abs(res)) > 1E-10
                disp(max(abs(res)))
                error( ['DS method failed. There are non zero residuals : ' num2str(res') ] )
            end
        end


    elseif strcmp( options_.portfolios.method , 'iterative')

        disp('Using ''iterative'' method');
        disp(' ');

        options_.order = 1;
        if options_.portfolios.solver == 'fsolve'
            alpha1 = fsolve( partial_alpha , alpha0 , solver_options );
        else
            alpha1 = csolve( partial_alpha , alpha0 , [], solver_options.TolFun, 500);
        end
        alphas = alpha1;
    else
        error(['Unknown method ' options_.portfolios.method])
    end

    oo_.portfolios.alphas = alphas;
  



    if ~strcmp( options_.portfolios.dynamic,'true')
        dr = dr_for_given_portfolios(alphas,dalpha0);
        % reorder before returning
        dr = rreorder_dr(dr, M_);
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    %% dynamic portfolios %%
    %%%%%%%%%%%%%%%%%%%%%%%%

    options_.order = 2;
    
    partial_dalpha = @(xx) temp_fun_2(alphas,xx,ind1,ind2,n_pf_vars, n_pf_s_vars);

    partial_alpha_all = @(xx) temp_fun_3(xx,ind1,ind2,n_pf_vars, n_pf_s_vars); % only used to check

    if strcmp( options_.portfolios.method , 'devereux-sutherland')
%        error('Dynamic Devereux-Sutherland solution not implemented')

        options_.order = 2;
        [res0,dr0] = partial_alpha(alphas); % use constant portfolios

        drn = reorder_dr(dr0);

        R1 = drn.ghu(ind2,EI.xi);
        D1 = drn.ghu(Delta,EI.xi); 

        R2 = drn.ghu(ind2,:);
        D2 = drn.ghu(Delta,:);
 
        n_state_variables = length(find(M_.portfolios.lead_lag_incidence(1,:)));

        tR5 = drn.ghxu(ind2,:);
        R5 = reshape(tR5,size(R2,1),M_.exo_nbr,n_state_variables);
        R5 = permute(R5,[1,3,2]);

        tD5 = drn.ghxu(Delta,:);        
        D5 = reshape(tD5,size(D2,1),M_.exo_nbr,n_state_variables);
        D5 = permute(D5,[1,3,2]);

        t1 = simple_tensor(simple_tensor(R5,M_.Sigma_e), D2' );
        t2 = simple_tensor(simple_tensor(D5,M_.Sigma_e), R2' );

        if length(size(t2)) == 3
            t2 = squeeze(t2)';
        end

        a = t1 + t2;

        lam_u = D1 * R2 * M_.Sigma_e * R2';
        dd =  R2 * M_.Sigma_e * D2' * R1';
        lam_u = lam_u - dd;

        dalphas = -lam_u\a;

        if options_.check_solution
            res = partial_dalpha(dalphas);
            if max(abs(res)) > 1E-10
                error( ['DS method failed. There are non zero residuals : ' num2str(res') ] )
            end
        end

        oo_.portfolios.dalphas = dalphas;


    elseif strcmp( options_.portfolios.method , 'iterative')
        if options_.portfolios.solver == 'fsolve'
            dalpha1 = fsolve(partial_dalpha,dalpha0(:),solver_options);
        else
            dalpha1 = csolve(partial_dalpha,dalpha0(:), [], solver_options.TolFun, 500);
        end
        dalpha1 = reshape(dalpha1,n_pf_vars, n_pf_s_vars);

    
        xx = [alpha1, dalpha1];
        if options_.portfolios.solver == 'fsolve'
    	     a_da = fsolve(partial_alpha_all, xx(:), solver_options);
        else
             a_da = csolve(partial_alpha_all, xx(:),  [], solver_options.TolFun, 500);
        end
        a_da = reshape(a_da,n_pf_vars,n_pf_s_vars+1);
        alphas = a_da(:,1);
        dalphas = a_da(:, 2:end);
    else
        error(['Unknown method ' options_.portfolios.method])
    end

        disp('Final portfolios : ')
        alpha1 = alphas;
        dalpha1 = dalphas;
        rownames = cellstr(M_.endo_names(pf_vars,:));
        colnames = [{'---' 'constant'} cellstr(M_.endo_names(pf_s_vars,:))'];
        r = [colnames; rownames, num2cell([alphas,dalphas])];
        disp(r)


	dr = dr_for_given_portfolios(alphas,dalphas);

    % reorder before returning
    
    %
    options_.order = 2;

    dr = rreorder_dr(dr, M_);

end

%% auxiliary nested functions

function [res,dr] = temp_fun_1(xx,dalpha,ind1,ind2)
        dr = dr_for_given_portfolios(xx,dalpha);
        drn = reorder_dr( dr );
        [res,junk] = euler_residuals( drn , ind1 , ind2  );
        res = res(:,1);
        res = res;
        res = res;
end

function [dres] = temp_fun_2(alpha,xx,ind1,ind2,n_pf_vars, n_pf_s_vars)
        xx1 = reshape(xx,n_pf_vars, n_pf_s_vars);
        dr = dr_for_given_portfolios(alpha,xx1);
        drn = reorder_dr( dr );
        [res,dres1] = euler_residuals( drn, ind1, ind2  );
        dres = dres1(:);
        res = res;
        dres = dres;
end

function [out] = temp_fun_3(xx,ind1,ind2,n_pf_vars, n_pf_s_vars)
        xx = reshape(xx,n_pf_vars,n_pf_s_vars+1);
        a = xx(:,1);
        da = xx(:,2:end);
        dr = dr_for_given_portfolios(a,da);
        drn = reorder_dr( dr );
        [res,dres1] = euler_residuals( drn, ind1, ind2  );
        out = [res, dres1];
        out = out(:);
end

function [C] = simple_tensor(A,B)
    % returns the tensor product of last dimension of A with first dimension of B
    s_A = size(A);
    s_B = size(B);
    C = zeros( [ s_A(1:end-1) s_B(2:end) ] );
    if s_A(end) ~= s_B(1)
        error('Dimensions mismatch');
    end
    t_A = s_A(1:(end-1));
    t_B = s_B(2:(end));
    AA = reshape(A, [prod(t_A),s_A(end)]);
    BB = reshape(B, [s_B(1),prod(t_B)]);
    CC = AA * BB;
    C = reshape( CC, [t_A,t_B] );
end
