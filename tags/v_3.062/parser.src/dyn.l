%{
#include <string.h>
#include "d_tab.h"
#include "d.h"
int previous_state;
int comment_level = 0;
int lineno = 1;
void yyerror(char *);
int name_var(char *);
%}



%option stack
%option yylineno
%x CMT
%x DYNARE1
%x DYNARE2
%x DYNARE3
%x DYNARE4
%x SKIP
%x THRU

%%



 /* literal keywords token */

_longnames {BEGIN DYNARE1; return LONGNAMES;}
initvalf {BEGIN THRU; ECHO;}
endval {BEGIN DYNARE2; return ENDVAL;}	       
initval {BEGIN DYNARE2; return INITVAL;}
histval {BEGIN DYNARE2; return HISTVAL;}
model {BEGIN DYNARE2; return MODEL;}
periods {BEGIN DYNARE1; return PERIODS;}
shocks {BEGIN DYNARE2; return SHOCKS;}
mshocks {BEGIN DYNARE2; return MSHOCKS;}
var {BEGIN DYNARE1; return VAR;}
varexo {BEGIN DYNARE1; return VAREXO;}
varexo_det {BEGIN DYNARE1; return VAREXO_DET;}
varrecur {BEGIN DYNARE1; return VARRECUR;}
parameters {BEGIN DYNARE1; return PARAMETERS;}
dyn2vec {BEGIN DYNARE1; return DYN2VEC;}
dyn2gau {BEGIN DYNARE1; return DYN2VEC;}
rplot {BEGIN DYNARE1; return RPLOT;}
steady {BEGIN DYNARE1; return STEADY;}
stoch_simul {BEGIN DYNARE1; return STOCH_SIMUL;}
resol {BEGIN DYNARE1; return RESOL;}
disp_dr {BEGIN DYNARE1; return DISP_DR;}
disp_moments {BEGIN DYNARE1; return DISP_MOMENTS;}
irf {BEGIN DYNARE1; return IRF;}
d_corr {BEGIN DYNARE1; return D_CORR;}
simul {BEGIN DYNARE1; return SIMUL;}
dsample {BEGIN DYNARE1; return DSAMPLE;}
check {BEGIN DYNARE1;  yylval.string=strdup(yytext); return CHECK;}
optim_weights {BEGIN DYNARE2; return OPTIM_WEIGHTS;}
osr_params {BEGIN DYNARE1; return OSR_PARAMS;}
osr {BEGIN DYNARE1; return OSR;}
calib_var {BEGIN DYNARE2; return CALIB_VAR;}
calib {BEGIN DYNARE1; return CALIB;}
dynatype {BEGIN DYNARE1; return DYNATYPE;}
dynasave {BEGIN DYNARE1; return DYNASAVE;}
sigma_e_ {BEGIN THRU; ECHO;}
sigma_e {BEGIN DYNARE3; return SIGMA_E;}
olr {BEGIN DYNARE1; return OLR;}
olr_inst {BEGIN DYNARE1; return OLR_INST;}
estimation {BEGIN DYNARE1; return ESTIMATION;}
varobs {BEGIN DYNARE1; return VAROBS;}
observation_trends {BEGIN DYNARE2; return OBSERVATION_TRENDS;}
estimated_params_init {BEGIN DYNARE2; return ESTIMATED_PARAMS_INIT;}
estimated_params_bounds {BEGIN DYNARE2; return ESTIMATED_PARAMS_BOUNDS;}
estimated_params {BEGIN DYNARE2; return ESTIMATED_PARAMS;}
unit_root_vars {BEGIN DYNARE1; return UNIT_ROOT_VARS;}
model_comparison {BEGIN DYNARE1; return MODEL_COMPARISON;}
forecast {BEGIN DYNARE1; return FORECAST;}

<*>@define {yy_push_state(DYNARE1); return COMPILE_DEFINE;}
<*>@if {yy_push_state(DYNARE1); return COMPILE_IF;}
<*>@elseif {yy_push_state(DYNARE1); return COMPILE_ELSEIF;}
<*>@else[ \t\n]*; {return COMPILE_ELSE;}
<*>@endif[ \t\n]*; {return COMPILE_ENDIF;}

<SKIP>. ;

 /* \r is always ignored */
<*>\r ;

 /* comments */
<*>[/][/].* ;
<*>"/*" {yy_push_state(CMT); ECHO;}
<CMT>[^*\n]* {ECHO;}
<CMT>"*"+[^*/\n]* {ECHO;}
<CMT>"*"+"/" {yy_pop_state(); ECHO;}
<CMT>[\n] {ECHO;}

 /* a semicolon in INITIAL state doesn't modify this state */
; ECHO;

 /* initial spaces or tabs are ignored */

[ \t] ;

 /* everything which is not a DYNARE command or a comment is ECHOED*/
[^ \t\n;] {BEGIN THRU;ECHO;}
<THRU>[^;\n] ECHO;
<THRU>[;\n] {BEGIN INITIAL; ECHO;}
 
<DYNARE1>; {BEGIN INITIAL; return yytext[0];}
<DYNARE2,DYNARE3,DYNARE4>; {return yytext[0];}
<DYNARE1,DYNARE2,DYNARE3,DYNARE4>[\(\)] {return yytext[0];}

 /* var, periods and values in shock statements */
<DYNARE2>var {return VAR;}
<DYNARE2>periods {return PERIODS;}
<DYNARE2>values {return VALUES;}
<DYNARE2>stderr {return STDERR;}

<DYNARE1>varexo {return VAREXO;}

<DYNARE1,DYNARE2,DYNARE3>end[ \t\n]*; {BEGIN INITIAL; return END;}

 /* lines starting with '#' */
<DYNARE2>^# {BEGIN DYNARE4; return POUND;}
<DYNARE4>[\n] {BEGIN DYNARE2; return EOL;}

 /* do loops in models */
<DYNARE2>do {return DO;}
<DYNARE2>to {return TO;}
<DYNARE2>by {return BY;}
<DYNARE2>endo {return ENDO;}
<DYNARE2>sum {return SUM;}
<DYNARE2>prod {return PROD;}

 /* options */
<DYNARE1>dr_algo {return DR_ALGO;}
<DYNARE1>simul_algo {return SIMUL_ALGO;}
<DYNARE1>solve_algo {return SOLVE_ALGO;}
<DYNARE1>drop {return DROP;}
<DYNARE1,DYNARE2>linear {return LINEAR;}
<DYNARE1>order {return ORDER;}
<DYNARE1>replic {return REPLIC;}
<DYNARE1>shock_size {return SHOCK_SIZE;}
<DYNARE1>moments {return MOMENTS;}
<DYNARE1>nomoments {return NOMOMENTS;}
<DYNARE1>corr {return CORR;}
<DYNARE1>nocorr {return NOCORR;}
<DYNARE1>ar {return AR;}
<DYNARE1>functions {return FUNCTIONS;}
<DYNARE1>nofunctions {return NOFUNCTIONS;}
<DYNARE1>irf {return IRF;}
<DYNARE1>irf_type {return IRF_TYPE;}
<DYNARE1>hp_filter {return HP_FILTER;}
<DYNARE1>hp_ngrid {return HP_NGRID;}
<DYNARE1>simul_seed {return SIMUL_SEED;}
<DYNARE1>simul {return SIMUL;}
<DYNARE1>periods {return PERIODS;}
<DYNARE1>qz_criterium {return QZ_CRITERIUM;}
<DYNARE1>optim {return OPTIM;}
<DYNARE2>shocks_file {return SHOCKS_FILE;}

 /* calib options */
<DYNARE2>autocorr {return AUTOCORR;}

 /* olr options */
<DYNARE1>olr_beta {return OLR_BETA;}

 /* estimation options */
<DYNARE1>datafile {return DATAFILE;}
<DYNARE1>nobs {return NOBS;}
<DYNARE1>first_obs {return FIRST_OBS;}
<DYNARE1>prefilter {return PREFILTER;} 
<DYNARE1>presample {return PRESAMPLE;} 
<DYNARE1>lik_algo  {return LIK_ALGO;}  
<DYNARE1>lik_init  {return LIK_INIT;}  
<DYNARE1>nograph   {return NOGRAPH;}  	  
<DYNARE1>graph   {return GRAPH;}  	  
<DYNARE1>noprint   {return NOPRINT;}  	  
<DYNARE1>print   {return PRINT;}  	  
<DYNARE1>conf_sig  {return CONF_SIG;}  
<DYNARE1>mh_replic {return MH_REPLIC;} 
<DYNARE1>mh_drop   {return MH_DROP;}   
<DYNARE1>mh_jscale   {return MH_JSCALE;}   
<DYNARE1>mh_init_scale {return MH_INIT_SCALE;}
<DYNARE1>mode_file {return MODE_FILE;}
<DYNARE1>mode_compute {return MODE_COMPUTE;}
<DYNARE1>mode_check {return MODE_CHECK;}
<DYNARE1>prior_trunc {return PRIOR_TRUNC;}
<DYNARE1>mh_mode {return MH_MODE;}
<DYNARE1>mh_nblocks {return MH_NBLOCKS;}
<DYNARE1>load_mh_file {return LOAD_MH_FILE;}
<DYNARE1>loglinear {return LOGLINEAR;}
<DYNARE1>diagnostic {return DIAGNOSTIC;}
<DYNARE1>nodiagnostic {return NODIAGNOSTIC;}
<DYNARE1>xtick {return XTICK;}
<DYNARE1>xticklabel {return XTICKLABEL;}
<DYNARE1>bayesian_irf {return BAYESIAN_IRF;}
<DYNARE1>tex {return TEX;}
<DYNARE1>forecast {return FORECAST;}
<DYNARE1>smoother {return SMOOTHER;}
<DYNARE1>moments_varendo {return MOMENTS_VARENDO;}
<DYNARE1>filtered_vars {return FILTERED_VARS;}
<DYNARE1>relative_irf {return RELATIVE_IRF;}
<DYNARE1>kalman_algo {return KALMAN_ALGO;}
<DYNARE1>kalman_tol {return KALMAN_TOL;}
<DYNARE1>diffuse_d {return DIFFUSE_D;}
<DYNARE1>nk {return NK;}
<DYNARE1>model_comparison_approximation {return MODEL_COMPARISON;}
<DYNARE1>laplace {return LAPLACE;}
<DYNARE1>modifiedharmonicmean {return MODIFIEDHARMONICMEAN;}

 /* estimated_params priors */
 <DYNARE2>gamma_pdf {return GAMMA_PDF;}
 <DYNARE2>beta_pdf {return BETA_PDF;}
 <DYNARE2>normal_pdf {return NORMAL_PDF;}
 <DYNARE2>inv_gamma_pdf {return INV_GAMMA_PDF;}
 <DYNARE2>inv_gamma1_pdf {return INV_GAMMA_PDF;}
 <DYNARE2>inv_gamma2_pdf {return INV_GAMMA_PDF;}
 <DYNARE2>uniform_pdf {return UNIFORM_PDF;}
 <DYNARE2>corr {return CORR;}

 /* spaces, tabs and EOL are ignored */ 
<DYNARE1,DYNARE2,DYNARE3>[ \t\n] ;

 /* names */
<DYNARE1,DYNARE3>[A-Za-z_][A-Za-z0-9_]* {
  yylval.string=strdup(yytext);
  return NAME;}

<DYNARE2,DYNARE4>[A-Za-z_][A-Za-z0-9_]* {
  return name_var(strdup(yytext));}

 /* floating-point numbers */
<DYNARE1,DYNARE2,DYNARE3,DYNARE4>-?((([0-9]*\.[0-9]+)|([0-9]+\.))([edED][-+]?[0-9]+)?)|([0-9]+[edED][-+]?[0-9]+) {
  yylval.string=strdup(yytext);
  return DNUMBER;}

 /* integer number */
<DYNARE1,DYNARE2,DYNARE3,DYNARE4>-?[0-9]+ {
  yylval.string=strdup(yytext); 
  return INUMBER;}
 
<DYNARE1,DYNARE2>: {return yytext[0];}
<DYNARE2>= {return EQUAL;}
<DYNARE1>"="+"=" {yylval.string=strdup(yytext); return EQ;}
<DYNARE1>"<"+"=" {yylval.string=strdup(yytext); return LE;}
<DYNARE1>">"+"=" {yylval.string=strdup(yytext); return GE;}
<DYNARE1>"!"+"=" {yylval.string=strdup(yytext); return NE;}
<DYNARE1,DYNARE3>= {return yytext[0];}
<DYNARE2>, {return yytext[0];}
<DYNARE2,DYNARE4>\$ {return DOLLAR;}
<DYNARE1>\$ {return yytext[0];}
<DYNARE3>] {BEGIN DYNARE1; return yytext[0];}
<DYNARE3>\[ {return yytext[0];}
<DYNARE1>, {return yytext[0];}
<DYNARE1>' {return yytext[0];}
<DYNARE1,DYNARE2,DYNARE3,DYNARE4>. {yylval.string=strdup(yytext); return OPERATORS;}
<*>. {return yytext[0];}

%%

int yywrap()
{
 return 1;
}
 
 
int name_var(char *s)
{
  struct var *p;

  p=var_search(s);
  if (p == NULL){
    yylval.string=s;
    return NAME;
  }
  else if (p->endo_exo == 2){
    yylval.string=s;
    return INDEX;
  }
  else{
    yylval.p_tok=create_var(p);
    return VAR_ID;
  }
}
 
void start_skipping(void)
{
  yy_push_state(SKIP);
}

void pop_flex_state(void)
{
  yy_pop_state();
}

/*
01/01/03 MJ added DYNASAVE
12/31/02 MJ corrected THRU mode to permit native lines without a final ;
         MJ put back CMT mode for Gauss style comments (/ *...* /)
04/06/02 MJ added OPTIM_WEIGHTS OSR_PARAMS OSR
*/

/*
*/




