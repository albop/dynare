%
% Status : main Dynare file 
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

clear all
tic;
global M_ oo_ options_ ys0_ ex0_
options_ = [];
M_.fname = 'walsh1_old_ss';
%
% Some global variables initialization
%
global_initialization;
diary off;
logname_ = 'walsh1_old_ss.log';
if exist(logname_, 'file')
    delete(logname_)
end
diary(logname_)
M_.exo_names = 'e';
M_.exo_names_tex = 'e';
M_.exo_names = char(M_.exo_names, 'sigma');
M_.exo_names_tex = char(M_.exo_names_tex, 'sigma');
M_.endo_names = 'y';
M_.endo_names_tex = 'y';
M_.endo_names = char(M_.endo_names, 'c');
M_.endo_names_tex = char(M_.endo_names_tex, 'c');
M_.endo_names = char(M_.endo_names, 'k');
M_.endo_names_tex = char(M_.endo_names_tex, 'k');
M_.endo_names = char(M_.endo_names, 'm');
M_.endo_names_tex = char(M_.endo_names_tex, 'm');
M_.endo_names = char(M_.endo_names, 'n');
M_.endo_names_tex = char(M_.endo_names_tex, 'n');
M_.endo_names = char(M_.endo_names, 'R');
M_.endo_names_tex = char(M_.endo_names_tex, 'R');
M_.endo_names = char(M_.endo_names, 'pi');
M_.endo_names_tex = char(M_.endo_names_tex, 'pi');
M_.endo_names = char(M_.endo_names, 'z');
M_.endo_names_tex = char(M_.endo_names_tex, 'z');
M_.endo_names = char(M_.endo_names, 'u');
M_.endo_names_tex = char(M_.endo_names_tex, 'u');
M_.param_names = 'alpha';
M_.param_names_tex = 'alpha';
M_.param_names = char(M_.param_names, 'beta');
M_.param_names_tex = char(M_.param_names_tex, 'beta');
M_.param_names = char(M_.param_names, 'delta');
M_.param_names_tex = char(M_.param_names_tex, 'delta');
M_.param_names = char(M_.param_names, 'gamm');
M_.param_names_tex = char(M_.param_names_tex, 'gamm');
M_.param_names = char(M_.param_names, 'phi1');
M_.param_names_tex = char(M_.param_names_tex, 'phi1');
M_.param_names = char(M_.param_names, 'eta');
M_.param_names_tex = char(M_.param_names_tex, 'eta');
M_.param_names = char(M_.param_names, 'a');
M_.param_names_tex = char(M_.param_names_tex, 'a');
M_.param_names = char(M_.param_names, 'b');
M_.param_names_tex = char(M_.param_names_tex, 'b');
M_.param_names = char(M_.param_names, 'rho');
M_.param_names_tex = char(M_.param_names_tex, 'rho');
M_.param_names = char(M_.param_names, 'phi2');
M_.param_names_tex = char(M_.param_names_tex, 'phi2');
M_.param_names = char(M_.param_names, 'Psi');
M_.param_names_tex = char(M_.param_names_tex, 'Psi');
M_.param_names = char(M_.param_names, 'thetass');
M_.param_names_tex = char(M_.param_names_tex, 'thetass');
M_.exo_det_nbr = 0;
M_.exo_nbr = 2;
M_.endo_nbr = 9;
M_.param_nbr = 12;
M_.orig_endo_nbr = 9;
M_.aux_vars = [];
M_.Sigma_e = zeros(2, 2);
M_.H = 0;
options_.block=0;
options_.bytecode=0;
options_.use_dll=0;
erase_compiled_function('walsh1_old_ss_dynamic');
M_.lead_lag_incidence = [
 0 5 0;
 0 6 14;
 1 7 0;
 2 8 15;
 0 9 0;
 0 10 16;
 0 11 17;
 3 12 0;
 4 13 0;]';
M_.equations_tags = {
};
M_.exo_names_orig_ord = [1:2];
M_.maximum_lag = 1;
M_.maximum_lead = 1;
M_.maximum_endo_lag = 1;
M_.maximum_endo_lead = 1;
oo_.steady_state = zeros(9, 1);
M_.maximum_exo_lag = 0;
M_.maximum_exo_lead = 0;
oo_.exo_steady_state = zeros(2, 1);
M_.params = NaN(12, 1);
M_.NNZDerivatives = zeros(3, 1);
M_.NNZDerivatives(1) = 36;
M_.NNZDerivatives(2) = 66;
M_.NNZDerivatives(3) = -1;
M_.params( 1 ) = 0.36;
alpha = M_.params( 1 );
M_.params( 2 ) = 0.989;
beta = M_.params( 2 );
M_.params( 4 ) = 0.5;
gamm = M_.params( 4 );
M_.params( 3 ) = 0.019;
delta = M_.params( 3 );
M_.params( 5 ) = 2;
phi1 = M_.params( 5 );
M_.params( 10 ) = 0;
phi2 = M_.params( 10 );
M_.params( 6 ) = 1;
eta = M_.params( 6 );
M_.params( 7 ) = 0.95;
a = M_.params( 7 );
M_.params( 8 ) = 2.56;
b = M_.params( 8 );
M_.params( 9 ) = 0.95;
rho = M_.params( 9 );
M_.params( 11 ) = 1.47630583;
Psi = M_.params( 11 );
M_.params( 12 ) = 1.0125;
thetass = M_.params( 12 );
%
% SHOCKS instructions
%
make_ex_;
M_.exo_det_length = 0;
M_.Sigma_e(1, 1) = (0.007)^2;
M_.Sigma_e(2, 2) = (0.0089)^2;
M_.sigma_e_is_diagonal = 1;
steady;
save('walsh1_old_ss_results.mat', 'oo_', 'M_', 'options_');
diary off

disp(['Total computing time : ' dynsec2hms(toc) ]);
