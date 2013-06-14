function [dr] = convert_portfolio_solution(M_,oo_,options_)

lli = M_.lead_lag_incidence;
lli_pf = M_.portfolios.lead_lag_incidence;
pf_vars = M_.portfolios.pf_vars;

dr_pf = oo_.portfolios.dr;
% dr_pf.ghx(dr_pf.inv_order_var(pf_vars),:) = 0;
% dr_pf.ghu(dr_pf.inv_order_var(pf_vars),:) = 0;
% dr_pf.ghxx(dr_pf.inv_order_var(pf_vars),:) = 0;
% dr_pf.ghxu(dr_pf.inv_order_var(pf_vars),:) = 0;
% dr_pf.ghuu(dr_pf.inv_order_var(pf_vars),:) = 0;
% dr_pf.ghs2(dr_pf.inv_order_var(pf_vars),:) = 0;



dr = struct;
dr = set_state_space(dr, M_, options_);
reorder = dr_pf.inv_order_var(dr.order_var);


n_s = max(lli(1,:));
n_s_pf = max(lli_pf(1,:));

% next line compute indices of the states from
% the modified model in the original model
inds = lli_pf(1, lli(1,:)~=0);
inds = find(inds);


n_y = size( oo_.portfolios.dr.ghx, 1 );
n_u = size( oo_.portfolios.dr.ghu, 2 );

ghx = zeros( n_y, n_s );
ghxx = zeros( n_y, n_s, n_s );
ghxu = zeros( n_y, n_u, n_s );

ghx(:,inds) = dr_pf.ghx;
ghu = dr_pf.ghu;


ys = dr_pf.ys;
ys = ys(reorder);
ghx = ghx(reorder,:);
ghu = ghu(reorder,:);



dr.ys = ys;
dr.ghx = ghx;
dr.ghu = ghu;

if isfield(dr_pf, 'ghxu')
    
    
    ghxu(:,:,inds) = reshape( dr_pf.ghxu, n_y, n_u,  n_s_pf);
    ghxx(:,inds,inds) = reshape( dr_pf.ghxx, n_y, n_s_pf,  n_s_pf);
    ghuu = dr_pf.ghuu;
    
    ghxu = reshape(ghxu, n_y, n_s*n_u);
    ghxx = reshape(ghxx, n_y, n_s*n_s);
    
    ghxx = ghxx(reorder,:);
    ghxu = ghxu(reorder,:);
    ghuu = ghuu(reorder,:);
    ghs2 = dr_pf.ghs2(reorder);
    
    % deactivate second order
    dr.ghxx = ghxx*0;
    dr.ghxu = ghxu*0;
    dr.ghuu = ghuu*0;
    % dr.kstate = ?
    % fuu  ?
    dr.ghs2 = ghs2*0;
    
end

end
