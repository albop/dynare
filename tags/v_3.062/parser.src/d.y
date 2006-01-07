%{
#include <stdio.h>
#include <string.h>
#include "d.h"

  extern FILE *f_out;
  extern int yylineno, var_nbr, longname, dr_algo, simul_algo, drop, linear, order, replic, iter, nbr_tmpvar, ar, nocorr, nofunctions, nomoments, irf, hp_filter, hp_ngrid, simul, simul_seed;
  extern double shock_size;
  extern char * varexo;
  extern struct s_check check;
  int ms_flag, varlist_flag, tr_row, tr_col, permit_sigma_e_;
  struct queue *model;
  extern struct s_estim_params estim_params;
  char* optim_string;
  char* empty_string = "NaN";

%}

%union {
        char *string;
	struct queue *p_queue;
	struct token *p_tok;
        struct loop *p_loop;
      }

%token  END ENDVAL INITVAL HISTVAL MODEL PERIODS SHOCKS VALUES VAR VAREXO 
%token VAREXO_DET VARRECUR EQUAL
%token VARF VARP DEL SUM PROD TO MSHOCKS LONGNAMES PARAMETERS DYN2VEC RPLOT
%token DO TO ENDO PROD BY DOLLAR STEADY STDERR DSAMPLE
%token STOCH_SIMUL DR_ALGO SIMUL_ALGO SOLVE_ALGO 
%token DROP LINEAR ORDER REPLIC AR
%token SIGMA_E HP_FILTER HP_NGRID SIMUL_SEED
%token RESOL SIMUL IRF DISP_DR DISP_MOMENTS D_CORR SHOCK_SIZE OPTIM_WEIGHTS
%token OSR_PARAMS OSR CALIB_VAR CALIB AUTOCORR COVAR DYNATYPE DYNASAVE
%token OLR OLR_INST OLR_BETA CHECK
%token ESTIMATED_PARAMS GAMMA_PDF BETA_PDF NORMAL_PDF INV_GAMMA_PDF UNIFORM_PDF
%token INV_GAMMA1_PDF INV_GAMMA2_PDF
%token PREFILTER PRESAMPLE LIK_ALGO LIK_INIT CONF_SIG 
%token ESTIMATION DATAFILE NOBS FIRST_OBS VAROBS QZ_CRITERIUM MH_REPLIC MH_DROP
%token MH_JSCALE OPTIM MH_INIT_SCALE MODE_FILE MODE_COMPUTE MODE_CHECK
%token PRIOR_TRUNC MH_MODE MH_NBLOCKS LOAD_MH_FILE LOGLINEAR
%token UNIT_ROOT_VARS XTICK XTICKLABEL BAYESIAN_IRF RELATIVE_IRF
%token TEX FORECAST SMOOTHER MOMENTS_VARENDO FILTERED_VARS 
%token OBSERVATION_TRENDS ESTIMATED_PARAMS_INIT ESTIMATED_PARAMS_BOUNDS
%token KALMAN_ALGO KALMAN_TOL DIFFUSE_D NK
%token CORR MOMENTS FUNCTIONS DIAGNOSTIC PRINT GRAPH
%token NOCORR NOMOMENTS NOFUNCTIONS NODIAGNOSTIC NOPRINT NOGRAPH
%token MODEL_COMPARISON MODEL_COMPARISON_APPROXIMATION LAPLACE
%token MODIFIEDHARMONICMEAN SHOCKS_FILE IRF_TYPE
%token COMPILE_DEFINE COMPILE_IF COMPILE_ELSEIF COMPILE_ELSE COMPILE_ENDIF
%token <string> INUMBER DNUMBER NAME OPERATORS POUND EOL INDEX EQ NE LE GE
%token <p_tok> VAR_ID

%type <p_queue> equation_list equation other_inst p_expr expression elem_exp 
%type <p_queue> del_exp sum_exp prod_exp period_list value_list
%type <p_queue> value_list1 value_list2 expr1 do_loop
%type <p_queue> var_exp var_id_exp indexed_var_exp 
%type <p_queue> triangular_matrix triangular_row
%type <p_loop> loop_init
%type <string> value prior filename filename_elem compile_operator
%expect 5
%%

 statement_list: statement
      | statement_list statement 
      ;

 statement: periods 
      | var         
      | varexo      
      | varexo_det      
      | varrecur
      | model 
      | initval
      | endval 
      | histval
      | shocks
      | mshocks
      | obsolete
      | longnames    {longname=1;}
      | parameters
      | dyn2vec
      | rplot
      | steady
      | stoch_simul
      | resol
      | simul
      | dsample
      | check
      | irf
      | d_corr
      | disp_dr
      | disp_moments
      | optim_weights
      | osr_params
      | osr
      | calib_var
      | calib
      | dynatype
      | dynasave
      | sigma_e
      | olr
      | olr_inst
      | estimated_params
      | estimated_params_init
      | estimated_params_bounds
      | estimation
      | varobs
      | observation_trends
      | unit_root_vars
      | model_comparison
      | compile_statement
      | forecast
      ;
 
 longnames : LONGNAMES ';'
      ;

 periods : PERIODS INUMBER ';'{print_iter($2);}
      | PERIODS '=' INUMBER ';' {print_iter($3);}
      ;

var : VAR {varlist_flag=1;} varlist ';' {print_endo();}
      ;

varexo : VAREXO {varlist_flag=0;} varlist ';' {print_exo();}
      ;

varexo_det : VAREXO_DET {varlist_flag=5;} varlist ';' {print_exo_det();}
      ;

 varrecur : VARRECUR {varlist_flag=3;} varlist ';' {print_recur();}
      ;

 varlist : varlist NAME '-' NAME {add_var_range($2,$4,varlist_flag);}
      | varlist ',' NAME '-' NAME {add_var_range($3,$5,varlist_flag);}
      | varlist NAME              {add_var($2,varlist_flag);}
      | varlist ',' NAME          {add_var($3,varlist_flag);}
      | NAME '-' NAME              {add_var_range($1,$3,varlist_flag);}
      | NAME                       {add_var($1,varlist_flag);}
      ; 

 parameters : PARAMETERS {varlist_flag=4;} varlist ';' {print_param();}
      ;

 model : MODEL ';' equation_list END   {model=$3;print_model(model);}
      | MODEL '(' LINEAR ')' ';' equation_list END   {check.linear=1;model=$6;p_option("linear","1");print_model(model);}
      | MODEL ';' compile_statement equation_list END   {model=$4;print_model(model);}
      | MODEL '(' LINEAR ')' ';' compile_statement equation_list END   {check.linear=1;model=$7;p_option("linear","1");print_model(model);}
      ;

 initval : INITVAL ';' {p_initval();} initval_list END 
     {pe_initval();}
      ;

 endval : ENDVAL ';' {p_endval();} initval_list END
     {pe_endval();}
      ;

 histval : HISTVAL ';' {p_histval();} histval_list END
      ;

 shocks: SHOCKS ';' {ms_flag=0;p_i_shocks(ms_flag);} shock_list END {p_e_shocks();}
       | SHOCKS '(' options_shocks ')' ';' 
                    {ms_flag=0;p_i_shocks(ms_flag);} shock_list END {p_e_shocks();}
       ;

 mshocks: MSHOCKS ';' {ms_flag=1;p_i_shocks(ms_flag);} shock_list END {p_e_shocks();}
      | MSHOCKS '(' options_shocks ')' ';' 
                      {ms_flag=1;p_i_shocks(ms_flag);} shock_list END {p_e_shocks();}
      ;

 options_shocks: option_shocks
       | options_shocks ',' option_shocks
       ;

 option_shocks: o_shocks_file ; 

 equation_list : equation_list equation {$$=add_to_queue($1,$2);}
      | equation_list other_inst {$$=add_to_queue($1,$2);}
      | equation_list do_loop {$$=copy_queue($$,$2);}
      | equation {$$=create_queue($1);}
      | equation compile_statement {$$=$1;}
      | other_inst {$$=create_queue($1);}
      | do_loop
      ;

 equation : expression EQUAL expression ';' 
                 {$$=add_to_queue($1,token(" -(",-1));
		  $$=copy_queue($$,$3);
		  $$=add_to_queue($$,token(")",-1));}
      | expression ';'
      ;

 other_inst : POUND p_expr EOL {$$=$2; mark_pound($$);}
      ;

 do_loop : DO loop_init ';' equation_list ENDO ';' {$$=do_loop($2,$4);}
      ;

 loop_init : NAME EQUAL INUMBER TO INUMBER  {$$=initial_loop($1,$3,$5,"1");}
      | NAME EQUAL INUMBER TO INUMBER BY INUMBER {$$=initial_loop($1,$3,$5,$7);}
      | INDEX EQUAL INUMBER TO INUMBER  {$$=initial_loop($1,$3,$5,"1");}
      | INDEX EQUAL INUMBER TO INUMBER BY INUMBER {$$=initial_loop($1,$3,$5,$7);}
      ;

 p_expr : p_expr ';'  {$$ = $1;}
      |  p_expr ';' expression  {$$ = add_to_queue($1,token(";",-1)); $$=copy_queue($$,$3);}
      | expression {$$ = $1;}
      ;

 expression : expression elem_exp {$$=copy_queue($1,$2);}
      | elem_exp
      ;

 elem_exp : '(' expression ')' { $$ = create_queue(token("(",-1));
				 $$ = copy_queue($$,$2);
				 $$ = add_to_queue($$,token(")",-1));}
      | var_exp
      | expr1
      ;

 expr1: NAME {$$ = create_queue(token($1,-1));}
      | ',' {$$ = create_queue(token(",",-1));}
      | OPERATORS {$$ = create_queue(token($1,-1));}
      | ':' {$$ = create_queue(token(":",-1));}
      | INUMBER {$$ = create_queue(token($1,-1));}
      | DNUMBER {$$ = create_queue(token($1,-1));}
      | INDEX {$$ = create_queue(token($1,-2));}
      | del_exp
      | sum_exp
      | prod_exp
      ;

 var_exp : var_id_exp '(' INUMBER ')' {set_ll($1,$3,INUMBER);$$=$1;}
      | var_id_exp '(' OPERATORS INUMBER ')' {set_ll($1,$4,INUMBER);$$=$1;}
      | var_id_exp '(' INDEX ')' {set_ll($1,$3,INDEX);$$=$1;}
      | var_id_exp '(' OPERATORS INDEX ')' {set_ll($1,$4,INDEX);$$=$1;}
      | var_id_exp {set_ll($1,"0",INUMBER);$$=$1;}
      ;

 var_id_exp : VAR_ID {$$ = create_queue($1);}
      | indexed_var_exp
      ;

 indexed_var_exp : indexed_var_exp DOLLAR INDEX {$$=add_to_queue($1,token($3,-3));}
      | NAME DOLLAR INDEX {$$=create_queue(token($1,-3));$$=add_to_queue($$,token($3,-3));}
      ;

 del_exp : DEL '(' INUMBER ':' expression ')' {$$ = m_del($3,$5);}
      ;

 sum_exp : SUM '(' loop_init ':' expression ')'
      {$$ = operator_loop($3,$5,"+");}
      ;

 prod_exp : PROD '(' loop_init ':' expression ')'
      {$$ = operator_loop($3,$5,"*");}
      ;

 initval_list : initval_list initval_elem
      | initval_elem
      ;

 initval_elem : VAR_ID EQUAL expression ';' {p_init($1,$3);}
      | compile_statement
      ;

 histval_list : histval_list histval_elem
      | histval_elem
      | compile_statement
      ;

 histval_elem : VAR_ID '(' INUMBER ')' EQUAL expression ';' {p_hist($1,$3,$6);}
      ;

 shock_list : shock_list shock_elem
      | shock_elem
      ;

 shock_elem : VAR VAR_ID ';' PERIODS period_list ';' VALUES value_list ';'  {p_shocks($2,$5,$8,ms_flag);}
            | VAR VAR_ID ';' STDERR expression ';' {simul_algo=1; p_stderr($2,$5);}
            | VAR VAR_ID EQUAL expression ';' {simul_algo=1; p_variance($2,$2,$4);}
            | VAR VAR_ID ',' VAR_ID EQUAL expression ';' {simul_algo=1; p_variance($2,$4,$6);}
            | compile_statement
            ;

 period_list : period_list INUMBER  {$$=add_to_queue($1,periods($2,0));}
      | period_list ',' INUMBER  {$$=add_to_queue($1,periods($3,0));}
      | period_list INUMBER ':' INUMBER 
                                 {$$=add_to_queue($1,periods($2,$4));}
      | period_list ',' INUMBER ':' INUMBER 
                                 {$$=add_to_queue($1,periods($3,$5));}
      | INUMBER ':' INUMBER      {$$=create_queue(periods($1,$3));}
      | INUMBER                  {$$=create_queue(periods($1,0));}
      ;

 value_list: value_list1        
      | value_list2             
      ;

 value_list1: value_list1 DNUMBER  
                        {$$=add_to_queue($1,create_queue(token($2,-1)));}
      | value_list1 INUMBER  
                        {$$=add_to_queue($1,create_queue(token($2,-1)));}
      | DNUMBER         {$$=create_queue(create_queue(token($1,-1)));}
      | INUMBER         {$$=create_queue(create_queue(token($1,-1)));}
      ;

 value_list2: value_list2 '(' expression ')' 
                                 {$$=add_to_queue($1,$3);}
      | '(' expression ')'       {$$=create_queue($2);}
      ;

 obsolete : VARF varlist3 ';' {fputs("Warning: VARF is an obsolete feature\n",stderr);} 
      | VARP varlist3 ';' {fputs("Warning: VARP is an obsolete feature\n",stderr);}
      ;

 varlist3 : varlist3 NAME {;}
      | varlist3 ',' NAME {;}
      | NAME {;}
      ;

 dyn2vec : DYN2VEC ';' {dyn2vec(0);}
      | DYN2VEC {nbr_tmpvar = 0;} varlist4 ';' {dyn2vec(1);}
      ;

 varlist4 : varlist4 NAME {add_tmpvar($2,0);}
      | varlist4 NAME '=' NAME {add_tmpvar($2,$4);}
      | varlist4 ',' NAME {add_tmpvar($3,0);}
      | varlist4 ',' NAME '=' NAME {add_tmpvar($3,$5);}
      | NAME  {nbr_tmpvar = 0; add_tmpvar($1,0);}
      | NAME '=' NAME {nbr_tmpvar = 0; add_tmpvar($1,$3);}
      ;

 rplot : RPLOT {nbr_tmpvar = 0;} varlist4 ';' {print_rplot();}
      ; 

 steady : STEADY ';' {p_steady(model);}
      | STEADY '(' steady_options ')' ';' {p_steady_linear(model);check.linear=1;}
      ;

 steady_options: steady_options ',' steady_option
      | steady_option
      ;

 steady_option: o_solve_algo
      | o_print
      | o_noprint
      ;

 stoch_simul : STOCH_SIMUL ';' {p_stoch_simul(dr_algo,simul_algo,drop,linear,order,replic,ar,nocorr,nofunctions,nomoments,irf);}
      | STOCH_SIMUL '(' options_list1 ')' ';' 
                           {p_stoch_simul();}
      | STOCH_SIMUL varlist4 ';' {p_stoch_simul();}
      | STOCH_SIMUL '(' options_list1 ')' varlist4 ';' 
                           {p_stoch_simul();}
      ;

 o_dr_algo: DR_ALGO '=' INUMBER {p_option("dr_algo",$3);};
 o_solve_algo: SOLVE_ALGO '=' INUMBER {p_option("solve_algo",$3);};
 o_simul_algo: SIMUL_ALGO '=' INUMBER {p_option("simul_algo",$3);};
 o_linear: LINEAR {p_option("linear","1");};
 o_order: ORDER '=' INUMBER {p_option("order",$3);};
 o_replic: REPLIC '=' INUMBER {p_option("replic",$3);};
 o_drop: DROP '=' INUMBER {p_option("drop",$3);};
 o_ar: AR '=' INUMBER {p_option("ar",$3);};
 o_corr: CORR {p_option("nocorr","0");};
 o_nocorr: NOCORR {p_option("nocorr","1");};
 o_function: FUNCTIONS {p_option("nofunctions","0");};
 o_nofunction: NOFUNCTIONS {p_option("nofunctions","1");};
 o_moments: MOMENTS {p_option("nomoments","0");};
 o_nomoments: NOMOMENTS {p_option("nomoments","1");};
 o_irf: IRF '=' INUMBER {p_option("irf",$3);};
 o_hp_filter: HP_FILTER '=' INUMBER {p_option("hp_filter",$3);};
 o_hp_ngrid: HP_NGRID '=' INUMBER {p_option("hp_ngrid",$3);};
 o_periods: PERIODS '=' INUMBER {p_option("periods",$3);p_option("simul","1");};
 o_simul: SIMUL {p_option("simul","1");};
 o_simul_seed: SIMUL_SEED '=' INUMBER { p_option("simul_seed",$3)}
             | SIMUL_SEED '=' '(' expression ')' { p_option_e("simul_seed",$4)}
             ;
 o_qz_criterium: QZ_CRITERIUM '=' INUMBER { p_option("qz_criterium",$3)}
               | QZ_CRITERIUM '=' DNUMBER { p_option("qz_criterium",$3)}
               ;
 o_datafile: DATAFILE '=' NAME {p_s_option("datafile",$3);};   
 o_nobs: NOBS '=' INUMBER {p_option("nobs",$3);}
       | NOBS '=' '(' expression ')' {p_option_e("nobs",$4);}
       ;
/*
       | NOBS '=' '[' INUMBER ':' INUMBER ']' {p_interval_option("nobs",$4,$6);};
*/
 o_first_obs: FIRST_OBS '=' INUMBER {p_option("first_obs",$3);};
 o_prefilter: PREFILTER '=' INUMBER {p_option("prefilter",$3);};
 o_presample: PRESAMPLE '=' INUMBER {p_option("presample",$3);};
 o_lik_algo: LIK_ALGO '=' INUMBER {p_option("lik_algo",$3);}; 
 o_lik_init: LIK_INIT '=' INUMBER {p_option("lik_init",$3);}; 
 o_graph: GRAPH {p_option("nograph","0");}; 
 o_nograph: NOGRAPH {p_option("nograph","1");}; 
 o_print: PRINT {p_option("print","0");}; 
 o_noprint: NOPRINT {p_option("noprint","1");}; 
 o_conf_sig: CONF_SIG '=' DNUMBER {p_option("conf_sig",$3);}; 
 o_mh_replic: MH_REPLIC '=' INUMBER {p_option("mh_replic",$3);}; 
 o_mh_drop: MH_DROP '=' DNUMBER {p_option("mh_drop",$3);}; 
 o_mh_jscale: MH_JSCALE '=' value {p_option("mh_jscale",$3);}; 
 o_optim: OPTIM {p_optim_options("","",1);} '=' '(' optim_options ')'
{p_optim_options("","",3);};
 o_mh_init_scale :MH_INIT_SCALE '=' DNUMBER {p_option("mh_init_scale",$3);};
 o_mh_init_scale :MH_INIT_SCALE '=' INUMBER {p_option("mh_init_scale",$3);};
 o_mode_file : MODE_FILE '=' NAME {p_s_option("mode_file",$3);};
 o_mode_compute : MODE_COMPUTE '=' INUMBER {p_option("mode_compute",$3);};
 o_mode_check : MODE_CHECK {p_option("mode_check","1");};
 o_prior_trunc : PRIOR_TRUNC '=' DNUMBER {p_option("prior_trunc",$3);};
 o_mh_mode : MH_MODE '=' INUMBER {p_option("mh_mode",$3);};
 o_mh_nblcks : MH_NBLOCKS '=' INUMBER {p_option("mh_nblck",$3);};
 o_load_mh_file : LOAD_MH_FILE {p_option("load_mh_file","1");};
 o_loglinear : LOGLINEAR {p_option("loglinear","1");};
 o_diagnostic : DIAGNOSTIC {p_option("diagnostic","0");};
 o_nodiagnostic : NODIAGNOSTIC {p_option("nodiagnostic","1");};
 o_bayesian_irf : BAYESIAN_IRF {p_option("bayesian_irf","1");};
 o_tex : TEX {p_option("TeX","1");};
 o_forecast : FORECAST '=' INUMBER {p_option("forecast",$3);};
 o_smoother : SMOOTHER {p_option("smoother","1");};
 o_moments_varendo : MOMENTS_VARENDO {p_option("moments_varendo","1");};
 o_filtered_vars : FILTERED_VARS {p_option("filtered_vars","1");};
 o_relative_irf : RELATIVE_IRF {p_option("relative_irf","1");};
 o_kalman_algo : KALMAN_ALGO '=' INUMBER {p_option("kalman_algo",$3);};
 o_kalman_tol : KALMAN_TOL '=' value {p_option("kalman_tol",$3);};
 o_diffuse_d : DIFFUSE_D '=' INUMBER {p_option("diffuse_d",$3);};
 o_nk : NK '=' INUMBER {p_option("nk",$3);};
 o_model_comparison_approximation: MODEL_COMPARISON_APPROXIMATION '=' LAPLACE {p_s_option("model_comparison_approximation","Laplace");}
   | MODEL_COMPARISON_APPROXIMATION '=' MODIFIEDHARMONICMEAN {p_s_option("model_comparison_approximation","ModifiedHarmonicMean");}
 o_olr_beta : OLR_BETA '=' value {p_option("olr_beta",$3);};
 o_shocks_file : SHOCKS_FILE EQUAL NAME {p_s_option("shocks_file",$3);};


 optim_option1: '\'' NAME '\'' ',' '\'' NAME '\'' {p_optim_options($2,$6,2);}
              | '\'' NAME '\'' ',' value {p_optim_options($2,$5,2);}
              ;

 optim_options: optim_option1
              | optim_options ',' optim_option1;
              ;

 o_list1 : o_dr_algo
         | o_solve_algo
         | o_simul_algo
         | o_linear
         | o_order
         | o_replic
         | o_drop
         | o_ar
         | o_corr
         | o_function
         | o_moments
         | o_nocorr
         | o_nofunction
         | o_nomoments
         | o_irf
         | o_relative_irf
         | o_hp_filter
         | o_hp_ngrid
         | o_periods
         | o_simul
         | o_simul_seed
         | o_qz_criterium
         | o_print
         | o_noprint
         ;

 options_list1: options_list1 ',' o_list1
              | o_list1
              ;

 resol : RESOL ';' {p_resol(dr_algo,linear,order,model);}
      | RESOL '(' options_list2 ')' ';' {p_resol(dr_algo,linear,order,model);}
      ;

 options_list2 : options_list2 ',' DR_ALGO '=' INUMBER {dr_algo = atoi($5);}
      | options_list2 ',' LINEAR {linear = 1;}
      | options_list2 ',' ORDER '=' INUMBER {order = atoi($5);}
      | DR_ALGO '=' INUMBER {dr_algo = atoi($3);}
      | LINEAR {linear = 1;}
      | ORDER '=' INUMBER {order = atoi($3);}
      ;

 simul : SIMUL ';' {p_simul(0,linear,order,replic);}
       |SIMUL '(' options_list1 ')' ';' {p_simul(0,linear,order,replic);}
       ;

 dsample : DSAMPLE ';' {p_dsample(0);}
         | DSAMPLE INUMBER ';' {p_dsample(1,$2);}
         | DSAMPLE INUMBER INUMBER ';' {p_dsample(2,$2,$3);}
         ; 

 check : CHECK ';' {p_check();}
       | CHECK '(' check_options ')' ';' {p_check();} 
         ;

 check_options: check_options ',' check_option
       | check_option
       ;

 check_option: o_print
       | o_noprint
       | o_qz_criterium
       | o_solve_algo
       ;
        
 irf : IRF ';' {p_irf(varexo,shock_size,iter,drop,replic,order);}
     | IRF '(' option_list_irf ')' ';' {p_irf(varexo,shock_size,iter,drop,replic,order);}
     | IRF varlist4 ';' {p_irf(varexo,shock_size,iter,drop,replic,order);}
     | IRF '(' option_list_irf ')' varlist4 ';' {p_irf(varexo,shock_size,iter,drop,replic,order);}
     ;

 option_list_irf : option_list_irf ',' VAREXO '=' NAME {varexo = $5;}
                 | option_list_irf ',' SHOCK_SIZE '=' DNUMBER {shock_size = atof($5);} 
                 | option_list_irf ',' PERIODS '=' INUMBER {iter = atoi($5);}
                 | option_list_irf ',' DROP '=' INUMBER {drop = atoi($5);}
                 | option_list_irf ',' REPLIC '=' INUMBER {replic = atoi($5);}
                 | option_list_irf ',' ORDER '=' INUMBER {order = atoi($5);}
                 | VAREXO '=' NAME {varexo = $3;}
                 | SHOCK_SIZE '=' DNUMBER {shock_size = atof($3);}
                 | PERIODS '=' INUMBER {iter = atoi($3);}       
                 | DROP '=' INUMBER {drop = atoi($3);}
                 | REPLIC '=' INUMBER {replic = atoi($3);}         
                 | ORDER '=' INUMBER {order = atoi($3);}
                 ;

 d_corr : D_CORR ';' {p_d_corr();}
        | D_CORR varlist4  ';' {p_d_corr();}
        ;

 disp_dr : DISP_DR ';' {p_disp_dr(order);}
         | DISP_DR '(' ORDER '=' INUMBER ')' ';' {p_disp_dr(atoi($5));}
         | DISP_DR varlist4 ';' {p_disp_dr(order);}
         | DISP_DR '(' ORDER '=' INUMBER ')' varlist4 ';' {p_disp_dr(atoi($5));}
         ;

 disp_moments : DISP_MOMENTS ';' {p_disp_moments(drop,order);}
              | DISP_MOMENTS '(' option_list_disp_moments ')' ';' {p_disp_moments(drop,order);}
              | DISP_MOMENTS varlist4 ';' {p_disp_moments(drop,order);}
              | DISP_MOMENTS '(' option_list_disp_moments ')' varlist4 ';' {p_disp_moments(drop,order);}
              ;

 option_list_disp_moments : option_list_disp_moments ',' DROP '=' INUMBER {drop = atoi($5);}  
                          | option_list_disp_moments ',' ORDER '=' INUMBER {order = atoi($5);}
                          | DROP '=' INUMBER {drop = atoi($3);}   
                          | ORDER '=' INUMBER {order = atoi($3);}
                          ;

 optim_weights : OPTIM_WEIGHTS ';' {p_optim_weights_init();} optim_weights_list END
               ;

 optim_weights_list : optim_weights_list VAR_ID expression ';' {p_optim_weights1($2, $3);}
                    | optim_weights_list VAR_ID ',' VAR_ID expression ';' {p_optim_weights2($2, $4, $5);}
                    | VAR_ID expression ';' {p_optim_weights1($1, $2);}
                    | VAR_ID ',' VAR_ID expression ';' {p_optim_weights2($1, $3, $4);}
                    | compile_statement
                    ;

 osr_params : OSR_PARAMS osr_params_list ';'
            ;

 osr_params_list : osr_params_list ',' NAME {p_osr_params($3);}
                 | osr_params_list NAME {p_osr_params($2);}
                 | NAME {p_osr_params($1);}
                 | compile_statement   
                 ;

 osr : OSR ';' {p_osr();}
     | OSR '(' olr_options ')' ';' {p_osr();}
     | OSR varlist4 ';' {p_osr();}
     | OSR '(' olr_options ')' varlist4 ';' {p_osr();}
     ;
 
 olr : OLR ';' {p_olr();}
     | OLR '(' olr_options ')' ';' {p_olr();}
     | OLR varlist4 ';' {p_olr();}
     | OLR '(' olr_options ')' varlist4 ';' {p_olr();}
     ;
 
 olr_option : o_olr_beta
     | options_list1
     ;
 
 olr_options : olr_option
     | olr_options ',' olr_option
     ;

 olr_inst : OLR_INST {p_olr_inst_init();} olr_inst_list ';'
          ;

 olr_inst_list : olr_inst_list ',' NAME {p_olr_inst($3);}
              | olr_inst_list NAME {p_olr_inst($2);}
              | NAME {p_olr_inst($1);}
              ;

 calib_var : CALIB_VAR ';' {p_calib_init();} calib_var_list END
           ;

 calib_var_list : calib_var_list VAR_ID EQUAL expression ';' {p_calib_var($2, $4,"1");}
                | calib_var_list VAR_ID '(' DNUMBER ')' EQUAL expression ';' {p_calib_var($2,$7,$4);}
                | calib_var_list VAR_ID '(' INUMBER ')' EQUAL  expression ';' {p_calib_var($2,$7,$4);}
                |  calib_var_list VAR_ID ',' VAR_ID EQUAL expression ';' {p_calib_covar($2,$4,$6,"1");}
                | calib_var_list VAR_ID ',' VAR_ID '(' DNUMBER ')' EQUAL expression ';' {p_calib_covar($2,$4,$9,$6);}
                | calib_var_list VAR_ID ',' VAR_ID '(' INUMBER ')' EQUAL  expression ';' {p_calib_covar($2,$4,$9,$6);}
                | calib_var_list AUTOCORR VAR_ID '(' INUMBER ')' EQUAL expression ';' {p_calib_ac($3,$5,$8,"1");}
                | calib_var_list AUTOCORR VAR_ID '(' INUMBER ')' '(' DNUMBER ')' EQUAL expression ';' {p_calib_ac($3,$5,$11,$8);}
                | calib_var_list AUTOCORR VAR_ID '(' INUMBER ')' '(' INUMBER ')' EQUAL  expression ';' {p_calib_ac($3,$5,$11,$8);}
                | VAR_ID '(' DNUMBER ')' EQUAL expression ';' {p_calib_var($1,$6,$3);}
                | VAR_ID '(' INUMBER ')' EQUAL expression ';' {p_calib_var($1,$6,$3);}
                | VAR_ID EQUAL expression ';' {p_calib_var($1, $3, "1");}
                | VAR_ID ',' VAR_ID '(' DNUMBER ')' EQUAL expression ';' {p_calib_covar($1,$3,$8,$5);}
                | VAR_ID ',' VAR_ID '(' INUMBER ')' EQUAL expression ';' {p_calib_covar($1,$3,$8,$5);}
                | VAR_ID ',' VAR_ID EQUAL expression ';' {p_calib_covar($1,$3,$5,"1");}
                | AUTOCORR VAR_ID '(' INUMBER ')' '(' DNUMBER ')' EQUAL expression ';' {p_calib_ac($2,$4,$10,$7);}
                | AUTOCORR VAR_ID '(' INUMBER ')' '(' INUMBER ')' EQUAL expression ';' {p_calib_ac($2,$4,$10,$7);}
                | AUTOCORR VAR_ID '(' INUMBER ')' EQUAL expression ';' {p_calib_ac($2,$4,$7,"1");}
                | compile_statement
                ;

 calib : CALIB ';' {p_calib(0);}
       | CALIB '(' COVAR ')' ';' {p_calib(1);}
       ;

 dynatype : DYNATYPE '(' NAME ')'';' {p_dynatype($3,"");}
          | DYNATYPE '(' NAME ')' varlist4 ';' {p_dynatype($3,"");}
          | DYNATYPE NAME ';' {p_dynatype($2,"");}
          | DYNATYPE '(' NAME '.' NAME ')'';' {p_dynatype($3,$5);}
          | DYNATYPE '(' NAME '.' NAME ')' varlist4 ';' {p_dynatype($3,$5);}
          | DYNATYPE NAME '.' NAME ';' {p_dynatype($2,$4);};

 dynasave : DYNASAVE '(' NAME ')'';' {p_dynasave($3,"");}
          | DYNASAVE '(' NAME ')' varlist4 ';' {p_dynasave($3,"");}
          | DYNASAVE NAME ';' {p_dynasave($2,"");}
          | DYNASAVE '(' NAME '.' NAME ')'';' {p_dynasave($3,$5);}
          | DYNASAVE '(' NAME '.' NAME ')' varlist4 ';' {p_dynasave($3,$5);}
          | DYNASAVE NAME '.' NAME ';' {p_dynasave($2,$4);};

 sigma_e : SIGMA_E '=' '[' {tr_row = 0;} triangular_matrix ']' ';' {p_sigma_e($5);}
         ;

 triangular_matrix : triangular_matrix ';' {tr_col = 0;} triangular_row {++tr_row;tr_check(tr_row,tr_col);$$=add_to_queue($1,$4);}
                   | {tr_col = 0;} triangular_row {++tr_row;tr_check(tr_row,tr_col);$$=create_queue($2);}
                   ;

 triangular_row : triangular_row ',' '(' expression ')' {$$=add_to_queue($1,$4);++tr_col;}
                | triangular_row ',' DNUMBER {$$=add_to_queue($1,create_queue(token($3,-1)));++tr_col;}
                | triangular_row ',' INUMBER {$$=add_to_queue($1,create_queue(token($3,-1)));++tr_col;}
                | triangular_row '(' expression ')' {$$=add_to_queue($1,$3);++tr_col;}
                | triangular_row DNUMBER {$$=add_to_queue($1,create_queue(token($2,-1)));++tr_col;}
                | triangular_row INUMBER {$$=add_to_queue($1,create_queue(token($2,-1)));++tr_col;}
                | '(' expression ')' {$$=create_queue($2);++tr_col;}
                | DNUMBER {$$=create_queue(token($1,-1));$$=create_queue($$);++tr_col;}
                | INUMBER {$$=create_queue(token($1,-1));$$=create_queue($$);++tr_col;}
                ;

 estimated_params : ESTIMATED_PARAMS ';' {p_estimated_init();} estimated_list END
                  ;

 estimated_list : estimated_list {estim_params_init();} estimated_elem {p_estimated_elem();}
                | {estim_params_init();} estimated_elem {p_estimated_elem();}
                | compile_statement
                ;

 estimated_elem : estimated_elem1 ',' estimated_elem2 ';'
                ;

 estimated_elem1 : STDERR VAR_ID {estim_params.param_type = 1; estim_params.var_nbr = $2->nbr+1; estim_params.var_type = $2->endo_exo;}
                 | VAR_ID {estim_params.param_type = 2; estim_params.param_name = $1->name;}
                 | CORR VAR_ID ',' VAR_ID {estim_params.param_type = 3; estim_params.var_nbr = $2->nbr+1; estim_params.var_nbr2 = $4->nbr+1; estim_params.var_type = $2->endo_exo;}
                 ;

 estimated_elem2 : value {estim_params.init_val=strdup($1);}
                 | value ',' value ',' value {estim_params.init_val=strdup($1);estim_params.lb=strdup($3);estim_params.ub=strdup($5);}
                 | prior ',' estimated_elem3 {estim_params.prior=strdup($1);}
                 | value ',' prior ',' estimated_elem3 {estim_params.init_val=strdup($1);estim_params.prior=strdup($3);}
                 | value ',' value ',' value ',' prior ',' estimated_elem3 {estim_params.init_val=strdup($1);estim_params.lb=strdup($3);estim_params.ub=strdup($5);estim_params.prior=strdup($7);}
                 ;

 estimated_elem3 : value ',' value {estim_params.mean=strdup($1);estim_params.std=strdup($3);}
                 | value ',' value ',' value ',' value {estim_params.mean=strdup($1);estim_params.std=strdup($3);estim_params.p3=strdup($5);estim_params.p4=strdup($7);}
                 | value ',' value ',' value ',' value ',' value {estim_params.mean=strdup($1);estim_params.std=strdup($3);estim_params.p3=strdup($5);estim_params.p4=strdup($7);estim_params.jscale=strdup($9);}
                 ;

 estimated_params_init: ESTIMATED_PARAMS_INIT ';' estimated_init_list END
                      ;

 estimated_init_list : estimated_init_list estimated_init_elem
                     | estimated_init_elem
                     | compile_statement
                     ;

 estimated_init_elem : STDERR VAR_ID ',' value ';' {p_estimated_init_elem1($2,$4);}
                     | CORR VAR_ID ',' VAR_ID ',' value ';' {p_estimated_init_elem2($2,$4,$6);}
                     | VAR_ID ',' value ';' {p_estimated_init_elem3($1,$3);}
                     ;

 estimated_params_bounds: ESTIMATED_PARAMS_BOUNDS ';' estimated_bounds_list END
                      ;

 estimated_bounds_list : estimated_bounds_list estimated_bounds_elem
                     | estimated_bounds_elem
                     | compile_statement
                     ;

 estimated_bounds_elem : STDERR VAR_ID ',' value ',' value ';' {p_estimated_bounds_elem1($2,$4,$6);}
  | CORR VAR_ID ',' VAR_ID ',' value ',' value ';' {p_estimated_bounds_elem2($2,$4,$6,$8);}
  | VAR_ID ',' value ',' value ';' {p_estimated_bounds_elem3($1,$3,$5);}
  ;

 prior : BETA_PDF {$$ = "1";}
       | GAMMA_PDF {$$ = "2";}
       | NORMAL_PDF {$$ = "3";}
       | INV_GAMMA_PDF {$$ = "4";}
       | INV_GAMMA1_PDF {$$ = "4";}
       | UNIFORM_PDF {$$ = "5";}
       | INV_GAMMA2_PDF {$$ = "6";}
       ;

 value : {$$ = empty_string;} /* empty */ 
       | INUMBER
       | DNUMBER
       | NAME
       | OPERATORS NAME {$$ = my_strcat($1,$2);}
       ;

 estimation : ESTIMATION ';' {p_estimation();}
            | ESTIMATION '(' estimation_options ')' ';' {p_estimation();}
            | ESTIMATION varlist4 ';' {p_estimation();}
            | ESTIMATION '(' estimation_options ')' varlist4 ';' {p_estimation();}
            ;

 estimation_options : estimation_options ',' estimation_option
                    | estimation_option
                    | estimation_options ',' o_list1
                    | o_list1
                    ;

 estimation_option : o_datafile
                   | o_nobs
                   | o_first_obs
                   | o_prefilter
                   | o_presample
                   | o_lik_algo 
                   | o_lik_init 
                   | o_graph
                   | o_nograph
                   | o_conf_sig 
                   | o_mh_replic
                   | o_mh_drop
                   | o_mh_jscale
                   | o_optim
                   | o_mh_init_scale 
                   | o_mode_file 
                   | o_mode_compute 
                   | o_mode_check
                   | o_prior_trunc 
                   | o_mh_mode 
                   | o_mh_nblcks 
                   | o_load_mh_file 
                   | o_loglinear
                   | o_diagnostic
                   | o_nodiagnostic
                   | o_bayesian_irf
                   | o_tex
                   | o_forecast
                   | o_smoother
                   | o_moments_varendo
                   | o_filtered_vars
                   | o_kalman_algo
                   | o_kalman_tol
                   | o_diffuse_d
                   | o_nk
                   ;

 varobs : VAROBS {nbr_tmpvar = 0;} varlist4 ';' {print_varobs();}
        ;

 observation_trends : OBSERVATION_TRENDS ';' {p_trend_init();} trend_list END
                    ;

 trend_list : trend_list trend_element
            | trend_element
            | compile_statement
            ;

 trend_element :  VAR_ID '(' expression ')' ';' {p_trend_element($1,$3);}
               ;


 unit_root_vars : UNIT_ROOT_VARS {nbr_tmpvar = 0;} varlist4 ';' {print_unit_root_vars();}
                ;

 model_comparison : MODEL_COMPARISON {nbr_tmpvar=0;} model_comparison_args
 ;

 model_comparison_args : filename_list1 ';' {p_model_comparison(0);}
              | filename_list2 ';' {p_model_comparison(1);}
              | '(' model_comparison_options ')' filename_list1 ';' {p_model_comparison(0);}
              | '(' model_comparison_options ')' filename_list2 ';' {p_model_comparison(1);}
              ;

 model_comparison_options: model_comparison_options ',' model_comparison_option
              | model_comparison_option
              ;

 model_comparison_option : o_model_comparison_approximation
              | o_print
              | o_noprint
              ;

 filename_list1 : filename {add_tmpvar($1,"");}
        | filename_list1 ',' filename {add_tmpvar($3,"");}
        ;

 filename_list2 : filename '(' value ')' {add_tmpvar($1,$3);}
        | filename_list2 ',' filename '(' value ')' {add_tmpvar($3,$5);}
        ;

 filename : filename_elem {$$=$1;}
        | filename filename_elem {$$ = my_strcat($1,$2);}
        ;

 filename_elem : NAME
        | '\\' {$$ = "\\";}
        | '/' {$$ = "/";}
        | ':' {$$ = ":";}
        | '.' {$$ = ".";}
        ;

  compile_statement : compile_define_statement
        | compile_if_statement
        | compile_elseif_statement
        | COMPILE_ELSE {compile_else();}
        | COMPILE_ENDIF {compile_endif();}
        ;

  compile_define_statement:
            COMPILE_DEFINE NAME ';' {compile_define($2);}
          | COMPILE_DEFINE NAME INUMBER ';'{compile_define($2,$3);}
	  ;

  compile_operator: EQ 
          | '<' {$$=strdup("<");} 
          | '>' {$$=strdup(">");} 
          | NE
          | LE
          | GE
          ; 

  compile_if_statement: 
    COMPILE_IF NAME compile_operator INUMBER ';' {compile_if($2,$3,$4);};

  compile_elseif_statement: 
    COMPILE_ELSEIF NAME compile_operator INUMBER ';' {compile_elseif($2,$3,$4);};

  forecast: FORECAST ';' {p_forecast();}
          | FORECAST '(' forecast_options ')' ';' {p_forecast();}
          | FORECAST varlist4 ';' {p_forecast();}
          | FORECAST '(' forecast_options ')' varlist4 ';' {p_forecast();}
          ;

  forecast_options: forecast_option
          | forecast_options ',' forecast_option
          ;

  forecast_option: o_periods
          | o_conf_sig
          ;


%%
int yyerror (char *s)
     
{
  fprintf (stdout, "%s at line %d\n", s, yylineno);
  exit(1);
}

/*
03/15/03 MJ added varexo_st and option LINEAR for STEADY
01/01/03 MJ corrected dynatype for filename extension and added dynasave
04/06/02 MJ added optim_weights and optim_weights_list osr_params 
            osr_params_list osr
10/09/02 MJ added calib_var calib_var_list calib
05/26/03 MJ added option LINEAR to MODEL
08/28/03 MJ added option SOLVE_ALGO
*/
