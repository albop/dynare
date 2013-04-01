function [res,g1,g2,g3]  = compute_derivatives(ys)

    % Compute derivatives of a Dynare model at given steady-state

    global M_ oo_
    
    lli = M_.lead_lag_incidence;
    z = (lli' ~= 0) .* repmat((1:M_.endo_nbr)',1,size(lli,1));
    z = z(find(z));
    z = ys(z);

    ffname = [M_.fname '_dynamic'];
    
    if nargout == 1
        [res] = feval(ffname,z,oo_.exo_steady_state',M_.params,[],1);
    elseif nargout == 2
        [res,g1] = feval(ffname,z,oo_.exo_steady_state',M_.params,[],1);
    elseif nargout == 3
        [res,g1,g2] = feval(ffname,z,oo_.exo_steady_state',M_.params,[],1);
    elseif nargout == 4
        [res,g1,g2,g3] = feval(ffname,z,oo_.exo_steady_state',M_.params,[],1);
    %elseif nargout == 4
    %    [res,g1,g2,g3,g4] = feval(ffname,z,oo_.exo_steady_state',M_.params,1);        
    end
