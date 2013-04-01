function [new_dr] = rreorder_dr(old_dr,new_M)

    % old order in old_dr is supposed to correspond to old_lli
    
    old_row_order = old_dr.order_var;
%    old_states_order = old_dr.order_var( (old_dr.nstatic+1):(old_dr.nstatic+1+old_dr.nboth+old_dr.npred) );
    old_states_order = old_dr.order_var( (old_dr.nstatic+1):(old_dr.nstatic+old_dr.npred) );

    new_dr = struct;
    new_dr = set_state_space(new_dr,new_M);
    new_row_order = new_dr.order_var;
%    new_states_order = new_dr.order_var( (new_dr.nstatic+1):(new_dr.nstatic+1+new_dr.nboth+new_dr.npred) );
    new_states_order = new_dr.order_var( (new_dr.nstatic+1):(new_dr.nstatic+new_dr.npred) );

%    index(old_row_order,new_row_order)

    row_position = zeros(size(old_row_order));
    states_position = zeros(size(old_states_order));
    
    for i = 1:length( old_row_order );
        row_position(i) = find( new_row_order == old_row_order(i) );
    end

    for i = 1:length( old_states_order );
        states_position(i) = find( new_states_order == old_states_order(i) );
    end

    n_v = length(new_row_order);
    n_states = length(new_states_order);
    n_s = length(new_M.Sigma_e);
    
    % first order
    ghx = zeros(n_v,n_states);
    ghx(row_position,states_position) = old_dr.ghx;
    
    ghu = zeros(n_v,n_s);
    ghu(row_position,:) = old_dr.ghu;

    new_dr.ys = old_dr.ys;
    new_dr.ghx = ghx;
    new_dr.ghu = ghu;

    if isfield(old_dr,'ghuu')
        ghxx = zeros(n_v, n_states, n_states);
        ghxx(row_position,states_position,states_position) = reshape( old_dr.ghxx, n_v, length(old_states_order), length(old_states_order) );

        ghxu = zeros(n_v, n_states, n_s);
        ghxu(row_position,states_position,:) = reshape( old_dr.ghxu, n_v, length(old_states_order), n_s);

        ghuu = zeros(n_v, n_s, n_s);
        ghuu(row_position,:,:) = reshape( old_dr.ghuu, n_v, n_s, n_s );

        ghs2(row_position) = old_dr.ghs2;

        new_dr.ghxx = reshape(ghxx, n_v, n_states * n_states);
        new_dr.ghxu = reshape(ghxu, n_v, n_states * n_s);
        new_dr.ghuu = reshape(ghuu, n_v, n_s * n_s);
        new_dr.ghs2 = ghs2(:);

    end

end
