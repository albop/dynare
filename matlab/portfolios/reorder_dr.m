function [new_dr] = reorder_dr(dr)
% This returns a structure similar to Dynare's decision rule structure, where vectors are sorted according to declaration order (as in var declaration in the modfile).
% 
%
% INPUTS 
%
%   dr        [structure]    Dynare structure (for decision rule)
%    
% OUTPUTS 
%
%   dr        [structure]    Dynare structure (for reordered decision rule)


    global M_ options_

    new_dr = dr;
    new_dr.normal_order = 1;
    new_dr.order_var = 1:length(M_.endo_names);
    new_dr.inv_order_var = 1:length(M_.endo_names);

    nx =size(dr.ghx,2);
    nu =size(dr.ghu,2);
    k = find(dr.kstate(:,2) <= M_.maximum_lag+1);
    klag = dr.kstate(k,[1 2]);
    k1 = dr.order_var; % pour une variable en position i dans dr, la position dans M_ est donnée par k1(i)

    new_dr.nx = nx;
    new_dr.nu = nu;

    % ys :
    new_dr.ys = dr.ys;


    % ghx :
    % sur la colonne i il y a les réponses à la variable d'état k1(klag(i,1)),:)   (dr order)
    % in disp_dr lines chosen in var_list are printed in dr order

    state_variables = dr.order_var(klag(1:nx,1)); % indices for variables printed in each column
    [state_variables,perm] = sort( state_variables );
    
    new_dr.ghx = dr.ghx( dr.inv_order_var , perm );
    new_dr.state_variables = state_variables;
    % now lines of ghx correspond to each variable in M_ order
    % columns relate to state variables in M_ order

    % ghu :
    % shocks are sorted aaccording to declaration order
    new_dr.ghu = dr.ghu( dr.inv_order_var , : );
    
    
    if options_.order == 2
      % ghxu : 
      % should exist a nice way to vectorize it
      new_dr.ghxu = zeros(size(dr.ghxu));
      for i = 1:length(perm)
        orig = (M_.exo_nbr * (perm(i)-1) + 1):(M_.exo_nbr * (perm(i)));
        dest = (M_.exo_nbr * (i-1) + 1):(M_.exo_nbr * (i));
        new_dr.ghxu(:,dest) = dr.ghxu(:,orig);
      end
      new_dr.ghxu = new_dr.ghxu(dr.inv_order_var,:);
    
    
      % ghxu : 
      % should exist a nice way to vectorize it
      np = length(perm);
      new_dr.ghxx = zeros(size(dr.ghxx));
      for i = 1:length(perm)
        orig = np * (perm(i)-1) + perm;
        dest = (np * (i-1) + 1):(np * (i));
        new_dr.ghxx(:,dest) = dr.ghxx(:,orig);
      end
      new_dr.ghxx = new_dr.ghxx(dr.inv_order_var,:);
    
    
      % ghuu : 
      new_dr.ghuu = dr.ghuu( dr.inv_order_var , : );
      

      new_dr.ghs2 = dr.ghs2( dr.inv_order_var );
    
    end
