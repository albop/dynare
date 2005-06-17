#ifdef MATLAB
#undef GAUSS
#undef SCILAB
#define PLATEFORM "MATLAB"
#define OPEN_COMMENTS "%"
#define CLOSE_COMMENTS ""
#define LEFT_PAREN "("
#define RIGHT_PAREN ")"
#define FILE_EXT "m"
#elif defined GAUSS
#undef MATLAB
#undef SCILAB
#define PLATEFORM "GAUSS"
#define OPEN_COMMENTS "/*"
#define CLOSE_COMMENTS "*/"
#define LEFT_PAREN "["
#define RIGHT_PAREN "]"
#define FILE_EXT "gau"
#elif defined SCILAB
#undef GAUSS
#undef MATLAB
#define PLATEFORM "SCILAB"
#define OPEN_COMMENTS "//"
#define CLOSE_COMMENTS ""
#define LEFT_PAREN "("
#define RIGHT_PAREN ")"
#define FILE_EXT "sci"
#endif


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "d.h"
#include "d_tab.h"

typedef struct var * t_var_p;
typedef struct queue * t_q_p;
typedef struct token * t_t_p;

int yydebug = 0;

/* global for parser program */
t_var_p var_list;
t_q_p do_loop_list;
int var_nbr, longname;
int endo_nbr, exo_nbr, exo_det_nbr, recur_nbr, y_max_lead, y_max_lag, x_max_lead, x_max_lag, z_max_lead, z_max_lag;
int *iy,*initval_check,*endval_check,*initval_check1;
char **tmpvar_list;
char * varexo;
int var_order_alpha = 0;
int stoch_shock_first_time = 1;
int ex_det_length = 0;
#if defined MATLAB || defined SCILAB
int param_nbr=0;
#endif
struct s_check check;
struct s_estim_params estim_params;
int inst_nbr = 0;

/* options, initialized in dynare_init */
int dr_algo, simul_algo, drop, linear, order, replic, iter, nbr_tmpvar;
int ar, irf, shock_size, nocorr, nofunctions, nomoments, hp_filter, hp_ngrid;
int simul, simul_seed;

int yyparse(void);
int var_comp1(const void *, const void *);
int var_comp2(const void *, const void *);
void str_output(char *);
void p_check(void);
void run_checks(void);
void set_options(int, char**,struct s_runtime_options*);
void make_empty(char *);
void p_model_comparison(int);

FILE *f_out;
char fname[200];

int main(int argc, char **argv)
{
  extern FILE *yyin,*yyout;
  char fn_in[200],fn_out[200],*cp,buffer[2000];
  struct s_runtime_options runtime_options;
  struct s_runtime_options* p_ro;
  p_ro = &runtime_options;
#ifdef GAUSS
  char fn_in2[200];
#endif
  if (argc > 1){
    strcpy(fname,argv[1]);
#ifdef GAUSS
    if ( (cp=strchr(fname,'.')) != NULL){
      if (strcasecmp(cp,".MDL")==0 || strcasecmp(cp,".MOD")==0 ){
	*cp='\0';
      }
      else{
	fprintf(stdout,
		"Error: %s is an invalid filename for DYNARE\n",
		argv[1]);
	exit(1);
      }
    }
    strcpy(fn_in,fname);
    strcpy(fn_out,fn_in);
    strcat(fn_in,".mod");
    strcat(fn_out,".gau");

    yyin = fopen(fn_in,"r");
    if (!yyin){
      strcpy(fn_in2,fname);
      strcat(fn_in2,".mdl");
      yyin = fopen(fn_in2,"r");
      if (!yyin){
	fprintf(stdout,"DYNARE: could not open %s or %s\n",fn_in,fn_in2);
	exit(1);
      }
    }
#elif defined MATLAB || defined SCILAB
    if ( (cp=strchr(fname,'.')) != NULL){
      if (strcasecmp(cp,".MOD")==0){
	*cp='\0';
      }
      else if (strcasecmp(cp,".DYN")==0){
	*cp='\0';
      }
      else{
	fprintf(stdout,
		"Error: %s is an invalid filename for DYNARE with %s\n",
		argv[1],PLATEFORM);
	exit(1);
      }
    }

    strcpy(fn_in,fname);
    strcat(fn_in,".dyn");
    yyin = fopen(fn_in,"r");
    if (!yyin)
      {
	strcpy(fn_in,fname);
	strcat(fn_in,".mod");
	yyin = fopen(fn_in,"r");
      }
    if (!yyin){
      fprintf(stdout,"DYNARE: could not open %s\n",fn_in);
      exit(1);
    }

    strcpy(fn_out,fname);
#ifdef MATLAB
    strcat(fn_out,".m");
#elif defined SCILAB
    strcat(fn_out,".sci");
#endif

#endif
    
    yyout = fopen(fn_out,"w");
    if (!yyout){
      fprintf(stdout,"DYNARE: could not open %s\n",fn_out);
      exit(1);
    }
    f_out=yyout;
    /* set runtime options */
    set_options(argc,argv,p_ro);
  }
  else{
    fprintf(stdout,"Error: DYNARE needs the name of the model file\n");
    exit(1);
  }
  if (runtime_options.debug == 1)
    {
      yydebug = 1;
    }
  dynare_init(fname,runtime_options);
  yyparse();
#if defined MATLAB
  sprintf(buffer,"save('%s_results','oo_','dr_');\n",fname);
  str_output(buffer);
  str_output("diary off\n");
#elif defined SCILAB
  str_output("mclose(fh_log);\n");
#elif defined GAUSS
  str_output("end;\n");
#endif
  run_checks();
  return 0;
}

void add_var(char *name,int endo_exo)
{
  var_list = (t_var_p)realloc(var_list,(var_nbr+1)*sizeof(struct var));
  var_list[var_nbr].name=strdup(name);
  var_list[var_nbr].endo_exo=endo_exo;
  var_nbr++;
}

void add_var_range(char *name1, char *name2, int endo_exo)
{
  int i1, i2, i3;
  char name[200],*p1,*p2;

  p1=name1;
  p2=name2;
  i1=0;
  while (*p1 == *p2)
    {
      p1++;
      p2++;
      i1++;
    }

  if (*p1 >= '0' && *p1 <= '9')
    {
      if (sscanf(p1,"%d",&i2) != 1)
	{
	  fprintf(stdout,"VAR, VAREXO or VARRECUR: error in name range expression\n");
	  return;
	}    
      if (sscanf(p2,"%d",&i3) != 1)
	{
	  fprintf(stdout,"VAR, VAREXO or VARRECUR: error in name range expression\n");
	  return;
	}    
      if (i3 < i2)
	{
	  fprintf(stdout,"VAR, VAREXO or VARRECUR: error in name range expression\n");
	  return;
	}
    
      var_list = (t_var_p)realloc(var_list,(var_nbr+i3-i2+1)*sizeof(struct var));
      strncpy(name,name1,i1);
      while (i2 <= i3)
	{
	  sprintf(name+i1,"%d\0",i2);
	  var_list[var_nbr].name=strdup(name);
	  var_list[var_nbr].endo_exo=endo_exo;
	  var_nbr++;
	  i2++;
	}
    }
  else
    {
      fprintf(stdout,"VAR, VAREXO or VARRECUR: error in name range expression\n");
      return;
    }
  return;
}

void sort_var()
{
  qsort(var_list,var_nbr,sizeof(struct var),var_comp1);
}

int var_comp1(const void *x1, const void *x2)
{
  return strcasecmp(((t_var_p)x1)->name,((t_var_p)x2)->name);
}

struct var *var_search(const char *key)
{
  return (t_var_p)bsearch(key,var_list,var_nbr,sizeof(struct var),var_comp2);
}

int var_comp2(const void *key, const void *element)
{
  return strcasecmp((char *)key,((t_var_p)element)->name);
}

struct queue *create_queue(void *p1)
{
  struct queue *p2;
  p2=(struct queue *)calloc(1,sizeof(struct queue));
  p2->ptr=(void **)calloc(1,sizeof(void *));
  (p2->ptr)[0] = p1;
  p2->imax = 1;
  p2->flag = 0;
  return p2;
}

struct queue *add_to_queue(struct queue *p, void *p1)
{
  p->ptr=(void **)realloc(p->ptr,(p->imax+1)*sizeof(void *));
  (p->ptr)[p->imax]=p1;
  p->imax++;
  return p;
}


void free_queue(struct queue *p)
{
  free(p->ptr);
  free(p);
}

struct queue *copy_queue(struct queue *p1, struct queue *p2)
{
  p1->ptr=(void **)realloc(p1->ptr,(p1->imax+p2->imax)*sizeof(void *));
  memcpy(p1->ptr+p1->imax,p2->ptr,p2->imax*sizeof(void *));
  p1->imax += p2->imax;
  free_queue(p2);
  return p1;
}

struct token *create_var(struct var *p)
{ 
    struct token *p1;
    p1 = (struct token *)calloc(1,sizeof(struct token));
    p1->endo_exo = p->endo_exo;
    p1->nbr = p->nbr;
    if (p1->endo_exo == 4){
      p1->name = p->name; /* It's a parameter */
    }
    else {
      p1->var_ptr = p; /* MJ 01/12/03 */
      p1->name = NULL;  /* It's a variable, not a string */
    }
    return p1;
}

void set_ll1(struct token *p_v, int i)
{
  int *temp;

  p_v->lead_lag=i;
  if (p_v->endo_exo == 1){
    if (i < -y_max_lag){
      temp=(int *)calloc(y_max_lead-i+1,endo_nbr*sizeof(int));
      memcpy(temp+(-i-y_max_lag)*endo_nbr,iy,(y_max_lag+y_max_lead+1)*endo_nbr*sizeof(int));
      free(iy);
      iy=temp;
      y_max_lag = -i;
    }
    else if (i > y_max_lead){
      iy=(int *)realloc(iy,(y_max_lag+i+1)*endo_nbr*sizeof(int));
      memset(iy+(y_max_lag+y_max_lead+1)*endo_nbr,0,(i-y_max_lead)*endo_nbr*sizeof(int));
      y_max_lead = i;
    }
    iy[p_v->nbr+(i+y_max_lag)*endo_nbr]=1;
  }
  else if (p_v->endo_exo == 0){
    if (i < -x_max_lag) x_max_lag = -i;
    else if (i > x_max_lead) x_max_lead = i;
    }
  else if (p_v->endo_exo == 3){
    if (i < -z_max_lag) z_max_lag = -i;
    else if (i > z_max_lead) z_max_lead = i;
    }
}

void set_ll(struct queue *p_q, char *s, int c)
{
  int i;
  struct token *p_v;

  p_v = *(struct token **)p_q->ptr;
  switch(c){
  case INUMBER:
    i=atoi(s);
    set_ll1(p_v, i);
    break;
  case INDEX:
    p_v->name=s;
  }
}

struct queue *m_del(char *s, struct queue *p_q)
{
;
}

struct token *token(char *s, int nbr)
{
  struct token *p;
  p=(struct token *)calloc(1,sizeof(struct token));
  p->name = s;
  p->nbr = nbr; /* > 0 : variable, -1: string, -2: loop index, -3: indexed variable */
  p->endo_exo = 0;
  return p;
}

void mark_pound(struct queue *p)
{
  p->flag=1;
}

int * periods(char *s1, char *s2)
{
  int *i;
  i=(int *)calloc(2,sizeof(int));
  i[0]=atoi(s1);
  if (s2 == 0) i[1]=0;
  else i[1]=atoi(s2);
  return i;
}

void p_shocks(struct token *var, struct queue *per, struct queue *val, int ms_flag)
{
  int i,j,i1,i2,i_par,flag;
  struct queue *p_q;
  struct token **p_t;
  char buffer[200];
  check.determ = 1;
  if (per->imax != val->imax){
    fprintf(stdout,"Error in SHOCKS: periods and values don't match for shocks on variable %s\n",var->var_ptr->name);
    return;
  }
  if (var->endo_exo != 0 && var->endo_exo != 5){
    fprintf(stdout,"Error in SHOCKS: %s isn't an exogenous variable\n",var->var_ptr->name); 
    return;
  }
  for (i=0; i < per->imax;i++){
    i_par=0;
    p_q = (struct queue *)val->ptr[i];
    if (var->endo_exo == 0)
      {
	i1= ((int *)(per->ptr[i]))[0]+x_max_lag;
	i2= ((int *)(per->ptr[i]))[1]+x_max_lag;
	flag = 0;
      }
    else
      {
	i1= ((int *)(per->ptr[i]))[0];
	i2= ((int *)(per->ptr[i]))[1];
	ex_det_length = (i1 > ex_det_length) ? i1 : ex_det_length; 
	ex_det_length = (i2 > ex_det_length) ? i2 : ex_det_length;
	i1 += y_max_lag;
	i2 += y_max_lag;
	flag = 2;
      }
    if (i2 == x_max_lag){
      if (ms_flag == 0)              /* mshocks */
	{
#ifdef GAUSS
	  sprintf(buffer,"_ex[%d,%d]=",i1,var->nbr+1);
#elif defined MATLAB || defined SCILAB
	  sprintf(buffer,"set_shocks(%d,%d,%d,",flag,i1,var->nbr+1);
	  i_par = 1;
#endif
	  str_output(buffer);
	}
      else
	{
#ifdef GAUSS
	  sprintf(buffer,"_ex[%d,%d]=_ex[%d,%d]*(",i1,var->nbr+1,i1,var->nbr+1);
	  i_par=1;
#elif defined MATLAB || defined SCILAB
	  sprintf(buffer,"set_shocks(%d,%d,%d,",flag+1,i1,var->nbr+1);
	  i_par = 1;
#endif
	  str_output(buffer);
	}
    }
    else{
      if (ms_flag == 0)              /* mshocks */
	{
#ifdef GAUSS
	  sprintf(buffer,"_ex[%d:%d,%d]=ones(%d,1).*",i1,i2+x_max_lag,var->nbr+1,i2-i1+x_max_lag+1);
#elif defined MATLAB || defined SCILAB
	  sprintf(buffer,"set_shocks(%d,[%d:%d],%d,",flag,i1,i2,var->nbr+1);
	  i_par = 1;
#endif
	  str_output(buffer);
	}
      else{
#ifdef GAUSS
	sprintf(buffer,"_ex[%d:%d,%d]=_ex[%d:%d,%d].*(",i1,i2+x_max_lag,var->nbr+1,i1,i2+x_max_lag,var->nbr+1);
#elif defined MATLAB || defined SCILAB
	sprintf(buffer,"set_shocks(%d,[%d:%d],%d,",flag+1,i1,i2,var->nbr+1);
#endif
	i_par=1;
	str_output(buffer);
      }
    }
    p_t=(struct token **) p_q->ptr;
    for (j=0;j < p_q->imax;j++){
      if ((*p_t)->nbr == -1)  str_output((*p_t)->name);
      else fprintf(stdout,"Reference to variable name isn't allowed in (m)shocks value statement");
      p_t++;
    }
    if (i_par) str_output(")");
    str_output(";\n");
  }
}

void p_stderr(struct token  *var, struct queue *expression)
{
  char buffer[200];
  check.stoch = 1;
#ifdef GAUSS
  if (stoch_shock_first_time)
    {
      sprintf(buffer,"_Sigma_e = zeros(%d,%d);\n",exo_nbr,exo_nbr);
      str_output(buffer);
      stoch_shock_first_time = 0;
    }
  sprintf(buffer,"_Sigma_e[%d,%d] = (",var->nbr+1,var->nbr+1);
  str_output(buffer);
  p_expression(expression);
  str_output(")^2;\n");
#elif defined MATLAB || defined SCILAB
  if (stoch_shock_first_time)
    {
      sprintf(buffer,"Sigma_e_ = zeros(%d,%d);\n",exo_nbr,exo_nbr);
      str_output(buffer);
      stoch_shock_first_time = 0;
    }
  sprintf(buffer,"Sigma_e_(%d,%d) = (",var->nbr+1,var->nbr+1);
  str_output(buffer);
  p_expression(expression);
  str_output(")^2;\n");
#endif
}

void p_variance(struct token  *var1, struct token  *var2, struct queue *expression)
{
  char buffer[200];
  check.stoch = 1;
#ifdef GAUSS
  if (stoch_shock_first_time)
    {
      sprintf(buffer,"_Sigma_e = zeros(%d,%d);\n",exo_nbr,exo_nbr);
      str_output(buffer);
      stoch_shock_first_time = 0;
    }
  sprintf(buffer,"_Sigma_e[%d,%d] = ",var1->nbr+1,var2->nbr+1);
  str_output(buffer);
  p_expression(expression);
  str_output(";\n");
  if (var1->nbr != var2->nbr)
    {
      sprintf(buffer,"_Sigma_e[%d,%d] = _Sigma_e[%d,%d];\n",var2->nbr+1,var1->nbr+1,var1->nbr+1,var2->nbr+1);
      str_output(buffer);
    }
#elif defined MATLAB || defined SCILAB
  if (stoch_shock_first_time)
    {
      sprintf(buffer,"Sigma_e_ = zeros(%d,%d);\n",exo_nbr,exo_nbr);
      str_output(buffer);
      stoch_shock_first_time = 0;
    }
  sprintf(buffer,"Sigma_e_(%d,%d) = ",var1->nbr+1,var2->nbr+1);
  str_output(buffer);
  p_expression(expression);
  str_output(";\n");
  if (var1->nbr != var2->nbr)
    {
      sprintf(buffer,"Sigma_e_(%d,%d) = Sigma_e_(%d,%d);\n",var2->nbr+1,var1->nbr+1,var1->nbr+1,var2->nbr+1);
      str_output(buffer);
    }
#endif
}

void p_init(struct token *v, struct queue *q)
{
  char buffer[200];
  int offset;

  offset = v->var_ptr-var_list;
  if (v->endo_exo == 1)
    {
      initval_check1[offset] = 1;
#ifdef GAUSS
      sprintf(buffer,"_ys[%d]=",v->nbr+1);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"ys_(%d)=",v->nbr+1);
#endif
    }
  else if (v->endo_exo == 0)
    {
      initval_check1[offset] = 1;
#ifdef GAUSS
      sprintf(buffer,"_exe[%d]=",v->nbr+1);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"exe_(%d)=",v->nbr+1);
#endif
    }
  else if (v->endo_exo == 3)
    {
#ifdef GAUSS
      sprintf(buffer,"_recurs[%d]=",v->nbr+1);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"recurs_(%d)=",v->nbr+1);
#endif
    }
    else if (v->endo_exo == 5)
    {
#ifdef GAUSS
      sprintf(buffer,"_exe_det[%d]=",v->nbr+1);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"exe_det_(%d)=",v->nbr+1);
#endif
    }
  str_output(buffer);

  p_expression(q);
  str_output(";\n");
}

void p_initval(void)
{
  char buffer[200];

  check.initval = 1;
  initval_check = (int *)calloc(var_nbr,sizeof(int));
  initval_check1 = initval_check;

#ifdef GAUSS
  str_output("/* INITVAL */\n");
  str_output("_valf = 0;\n_endval=0;\n");
#elif defined MATLAB || defined SCILAB
#ifdef MATLAB
  str_output("% ");
#else
  str_output("// ");
#endif
  str_output("INITVAL \n");
  str_output("valf_ = 0;\nendval_=0;\n");
#endif
  if(recur_nbr > 0){
#ifdef GAUSS
    sprintf(buffer,"_recurs = zeros(%d,1);\n",recur_nbr);
#elif defined MATLAB || defined SCILAB
    sprintf(buffer,"recurs_ = zeros(%d,1);\n",recur_nbr);
#endif
    str_output(buffer);
  }
#ifdef GAUSS
  str_output("_ys0 = 0;\n");
  str_output("_ex0 = 0;\n");
  str_output("_recurs0 = 0;\n");
#endif
}

void pe_initval(void)
{
#ifdef GAUSS
  fputs("if exo_nbr > 0;\n",f_out);
  fputs("  _ex=ones(_iter+_xkmin+_xkmax,1).*_exe';\n",f_out);
  fputs("endif;\n",f_out);
#elif defined MATLAB || defined SCILAB
  fputs("y_=[ys_*ones(1,ykmin_)];\n",f_out);
  fputs("if exo_nbr > 0;\n",f_out);
  fputs("  ex_=[ones(xkmin_,1)*exe_'];\n",f_out);
  fputs("end;\n",f_out);
  fputs("if exo_det_nbr > 0;\n",f_out);
  fputs("  ex_det_=[ones(ykmin_+1,1)*exe_det_'];\n",f_out);
  fputs("end;\n",f_out);
#endif
}

void p_endval(void)
{
  check.endval = 1;
  endval_check = (int *)calloc(var_nbr,sizeof(int));
  initval_check1 = endval_check;

#ifdef GAUSS  
  str_output("/* ENDVAL */\n");
  str_output("_ys0=_ys;\n_ex0=_exe;\n_recurs0=_recurs;\n_endval=1;\n");
#elif defined MATLAB || defined SCILAB
#ifdef MATLAB
  str_output("% ");
#else
  str_output("// ");
#endif
  str_output("ENDVAL \n");
  str_output("ys0_=ys_;\nex0_=exe_;\nex_det0_=exe_det_;\nrecurs0_=recurs_;\nendval_=1;\n");
#endif
}

void pe_endval(void)
{
#ifdef GAUSS
  fputs("if exo_nbr > 0;\n",f_out);
  fputs("  _ex[_xkmin+1:_iter+_xkmin+_xkmax,.]=ones(_iter+_xkmax,1).*_exe';\n",f_out);
  fputs("endif;\n",f_out);
#elif defined MATLAB || defined SCILAB
  fputs("y_=[y_ ys_*ones(1,iter_+ykmax_)];\n",f_out);
  fputs("if exo_nbr > 0;\n",f_out);
  fputs("  ex_=[ones(xkmin_,1)*ex0_';ones(iter_+xkmax_,1)*exe_'];\n",f_out);
  fputs("end;\n",f_out);
  fputs("if exo_det_nbr > 0;\n",f_out);
  fputs("  ex_det_=[ones(ykmin_+1,1)*ex_det0_';ones(iter_+ykmax_,1)*exe_det_'];\n",f_out);
  fputs("end;\n",f_out);
#endif
}

void p_hist(struct token *v, char *lag, struct queue *q)
{
  char buffer[200];
  int offset;

  offset = v->var_ptr-var_list;
  if (v->endo_exo == 1)
    {
#ifdef GAUSS
      sprintf(buffer,"_y[%d,_ykmin+(%s)]=",v->nbr+1,lag);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"y_(%d,ykmin_+(%s))=",v->nbr+1,lag);
#endif
    }
  else if (v->endo_exo == 0)
    {
      initval_check1[offset] = 1;
#ifdef GAUSS
      sprintf(buffer,"_ex[_xkmin+(%s),%d]=",lag,v->nbr+1);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"ex_(xkmin_+(%s),%d)=",lag,v->nbr+1);
#endif
    }
  else if (v->endo_exo == 5)
    {
      initval_check1[offset] = 1;
#ifdef GAUSS
      sprintf(buffer,"_ex_det[_xkmin+(%s),%d]=",lag,v->nbr+1);
#elif defined MATLAB || defined SCILAB
      sprintf(buffer,"ex_det_(ykmin_+(%s),%d)=",lag,v->nbr+1);
#endif
    }
/*   else if (v->endo_exo == 3) */
/*     { */
/* #ifdef GAUSS */
/*       sprintf(buffer,"_recurs[%d]=",v->nbr+1); */
/* #elif defined MATLAB || defined SCILAB */
/*       sprintf(buffer,"recurs_(%d)=",v->nbr+1); */
/* #endif */
/*     } */
  str_output(buffer);

  p_expression(q);
  str_output(";\n");
}

void p_histval(void)
{
  char buffer[200];

  check.histval = 1;

#ifdef GAUSS
  str_output("/* HISTVAL */\n");
#elif defined MATLAB || defined SCILAB
#ifdef MATLAB
  str_output("% ");
#else
  str_output("// ");
#endif
  str_output("HISTVAL \n");
#endif
}

void print_var()
{
  int i ;
  char buff[200];
  int n_exo, n_endo;
  n_exo = 0;
  n_endo = 0;
  for( i = 0; i < var_nbr; ++i)
    {
      if (var_list[i].endo_exo == 0)
	{
	  var_list[i].original_nbr = n_exo++;
	}
      else if (var_list[i].endo_exo == 1)
	{
	  var_list[i].original_nbr = n_endo++;
	}
    }
  sort_var();

#ifdef GAUSS  
  str_output("\nlet");
  if (longname) 
    str_output(" string");
  str_output(" _lgy = ");

  i=0;
  while (i < var_nbr){
    if (var_list[i].endo_exo == 1){
      var_list[i].nbr=endo_nbr;
      if (endo_nbr > 0) str_output(", ");
      strcpy(buff,"\"");
      strcat(buff,var_list[i].name);
      strcat(buff,"\"");
      str_output(buff);
      endo_nbr++;
    }
    i++;
  }
  str_output(";\n");

  i=0;
  while (i < var_nbr){
    if (var_list[i].endo_exo == 0){
      if (exo_nbr==0){
	str_output("\nlet");
	if (longname) str_output(" string");
	str_output(" _lgx = ");
      }
      var_list[i].nbr=exo_nbr;
      if (exo_nbr > 0) str_output(", ");
      strcpy(buff,"\"");
      strcat(buff,var_list[i].name);
      strcat(buff,"\"");
      str_output(buff);
      exo_nbr++;
    }
    i++;
  }
  if (exo_nbr > 0) str_output(";\n");

  i=0;
  while (i < var_nbr){
    if (var_list[i].endo_exo == 3){
      if (recur_nbr==0){
	str_output("\nlet");
	if (longname) str_output(" string");
	str_output(" _lgr = ");
      }
      var_list[i].nbr=recur_nbr;
      if (recur_nbr > 0) str_output(", ");
      strcpy(buff,"\"");
      strcat(buff,var_list[i].name);
      strcat(buff,"\"");
      str_output(buff);
      recur_nbr++;
    }
    i++;
  }
  if (recur_nbr > 0) str_output(";\n");

#elif defined MATLAB || defined SCILAB
  for (i=0; i < var_nbr; i++)
    {
      if ( var_list[i].endo_exo == 1)
	{
	  if ( endo_nbr == 0 )
	    {
	      sprintf(buff,"lgy_ = '%s';\n",var_list[i].name);
	    }
	  else
	    {
#ifdef MATLAB
	      sprintf(buff,"lgy_ = str2mat(lgy_,'%s');\n",var_list[i].name);
#else
	      sprintf(buff,"lgy_ = [lgy_;'%s'];\n",var_list[i].name);
#endif
	    }
	  str_output(buff) ;
	  var_list[i].nbr=endo_nbr ;
	  endo_nbr++ ;
	}
    }
  for (i=0; i < var_nbr; i++)
    {
      if ( var_list[i].endo_exo == 0)
	{
	  if ( exo_nbr == 0 )
	    {
	      sprintf(buff,"lgx_ = '%s';\n",var_list[i].name);
	    }
	  else
	    {
#ifdef MATLAB
	      sprintf(buff,"lgx_ = str2mat(lgx_,'%s');\n",var_list[i].name);
#else
	      sprintf(buff,"lgx_ = [lgx_;'%s'];\n",var_list[i].name);
#endif
	    }
	  str_output(buff);
	  var_list[i].nbr=exo_nbr ;
	  exo_nbr++ ;
	}
    }
  for (i=0; i < var_nbr; i++)
    {
      if ( var_list[i].endo_exo == 0)
	{
	  if ( var_list[i].nbr == 0 )
	    {
	      sprintf(buff,"lgx_orig_ord_ = [%d",var_list[i].original_nbr+1);
	    }
	  else
	    {
	      sprintf(buff," %d",var_list[i].original_nbr+1);
	    }
	  str_output(buff);
	}
    }
  if (exo_nbr > 0)
    {
      str_output("];\n");
    }
  for (i=0; i < var_nbr; i++)
    {
      if ( var_list[i].endo_exo == 3)
	{
	  if ( recur_nbr == 0 )
	    {
	      sprintf(buff,"lgr_ = '%s';\n",var_list[i].name);
	    }
	  else
	    {
#ifdef MATLAB
	      sprintf(buff,"lgr_ = str2mat(lgr_,'%s');\n",var_list[i].name);
#else
	      sprintf(buff,"lgr_ = [lgr_;'%s'];\n",var_list[i].name);
#endif
	    }
	  str_output(buff);
	  var_list[i].nbr=recur_nbr ;
	  recur_nbr++ ;
	}
    }
  for (i=0; i < var_nbr; i++)
    {
      if ( var_list[i].endo_exo == 5)
	{
	  if ( exo_det_nbr == 0 )
	    {
	      sprintf(buff,"lgx_det_ = '%s';\n",var_list[i].name);
	    }
	  else
	    {
#ifdef MATLAB
	      sprintf(buff,"lgx_det_ = str2mat(lgx_det_,'%s');\n",var_list[i].name);
#else
	      sprintf(buff,"lgx_det_ = [lgx_det_;'%s'];\n",var_list[i].name);
#endif
	    }
	  str_output(buff);
	  var_list[i].nbr=exo_det_nbr ;
	  exo_det_nbr++ ;
	}
    }
#endif

  sprintf(buff,"endo_nbr = %d;\n",endo_nbr);
  str_output(buff);
  sprintf(buff,"exo_nbr = %d;\n",exo_nbr);
  str_output(buff);
  sprintf(buff,"exo_det_nbr = %d;\n",exo_det_nbr);
  str_output(buff);
  sprintf(buff,"recur_nbr = %d;\n",recur_nbr);
  str_output(buff);
  iy=(int *)calloc(endo_nbr,sizeof(int));
}      
  
void str_output(char *str)
{
  int static column=0,p;
  char *ps;
  ps=str;
  p=0;
  while (*ps != '\n' && *ps != '\0') 
    {
      p++;
      ps++;
    }
#ifdef GAUSS
  if (column > 0 && (column+p) > 80)
    {
      fputc('\n',f_out);
      column=0;
    }
#elif defined MATLAB
  if (column > 0 && (column+p) > 76)
    {
      fputs(" ...\n",f_out);
      column=0;
    }
#endif
  fputs(str,f_out);
  if (*ps == '\n')
    {
      while (*ps != '\0')
	{
	  if (*ps == '\n') p = 0;
	  p++;
	  ps++;
	}
      column=p;
    }
  else column += p;
}

void i_output(int i)
{
  char s[20];
  sprintf(s,"%d",i);
  str_output(s);
}

void print_model(struct queue *p)
{
  int i,j,k,m;
  char buffer[200];

  i=0;
  j=1;
  for ( k=1; k <= (y_max_lag+y_max_lead+1); k++)
    {
      if (k == 1)
	{
#ifdef GAUSS
	  str_output("let _iy = ");
#elif defined MATLAB || defined SCILAB
	  str_output("iy_ = [");
#endif
	}
      else
	{
#ifdef GAUSS
	  str_output("let temp =");
#elif defined MATLAB || defined SCILAB
	  str_output("temp = [");
#endif
	}
      for( m=0 ; m < endo_nbr ; m++ )
	{
	  if (iy[i] == 1)
	    {
	      iy[i] = j;
	      str_output(" ");
	      i_output(j);
	      j++;
	    }
	  else 
	    {
	      str_output(" 0");
	    }
	  i++;
	}
#ifdef GAUSS
      str_output(";\n");
#elif defined MATLAB || defined SCILAB
      str_output("];\n");
#endif
      if (k==1) 
	{
#ifdef GAUSS
	  str_output("_iy = _iy';\n");
#endif
	}
      else
	{
#ifdef GAUSS
	  str_output("_iy = _iy | (temp');\n");
#elif defined MATLAB || defined SCILAB
	  str_output("iy_ = [ iy_ ; temp ];\n");
#endif
	}
    }

#ifdef GAUSS
  str_output("_ykmin = ");
  i_output(y_max_lag);
  str_output(";\n_ykmax = ");
  i_output(y_max_lead);
  str_output(";\n_xkmin = ");
  i_output(x_max_lag);
  str_output(";\n_xkmax = ");
  i_output(x_max_lead);
  str_output(";\n_zkmin = ");
  i_output(z_max_lag);
  str_output(";\n_zkmax = ");
  i_output(z_max_lead);
  str_output(";\n");
  sprintf(buffer,"_ys = zeros(%d,1);\n",endo_nbr);
  str_output(buffer);
  if (exo_nbr > 0)
    {
      sprintf(buffer,"_exe = zeros(%d,1);\n",exo_nbr);
      str_output(buffer);
    }
#elif defined MATLAB || defined SCILAB
  str_output("ykmin_ = ");
  i_output(y_max_lag);
  str_output(";\nykmax_ = ");
  i_output(y_max_lead);
  str_output(";\nxkmin_ = ");
  i_output(x_max_lag);
  str_output(";\nxkmax_ = ");
  i_output(x_max_lead);
  str_output(";\nzkmin_ = ");
  i_output(z_max_lag);
  str_output(";\nzkmax_ = ");
  i_output(z_max_lead);
  str_output(";\n");
  sprintf(buffer,"ys_ = zeros(%d,1);\n",endo_nbr);
  str_output(buffer);
  if (exo_nbr)
    {
      sprintf(buffer,"exe_ = zeros(%d,1);\n",exo_nbr);
      str_output(buffer);
    }
  if (exo_det_nbr)
    {
      sprintf(buffer,"exe_det_ = zeros(%d,1);\n",exo_det_nbr);
      str_output(buffer);
    }

#endif
print_model1(p,0);
print_model1(p,1);
}

void print_model1(struct queue *p, int flag_steady)
{
  struct queue **p_q;
  struct token **p_t;
  int i,j,k;

#ifdef GAUSS
  if (flag_steady == 0)
    {
      str_output("\nproc _ff(y);\n");
    }
  else
    {
      str_output("\nproc _fff(y);\n");
    }
  str_output("local z;\n");

  str_output("z=zeros(");
  i_output(endo_nbr);
  str_output(",1);\n");
#elif defined MATLAB || defined SCILAB
  FILE *f_temp;
  char fmod_name[200];
  
  f_temp = f_out;
  if (flag_steady == 0)
    {
      strcpy(fmod_name,fname);
      strcat(fmod_name,"_ff.");
      strcat(fmod_name,FILE_EXT);
      f_out = fopen(fmod_name,"w");
      strcpy(fmod_name,"\nfunction z=");
      strcat(fmod_name,fname);
      strcat(fmod_name,"_ff(y)\n");
      str_output(fmod_name);
    }
  else
    {
      strcpy(fmod_name,fname);
      strcat(fmod_name,"_fff.");
      strcat(fmod_name,FILE_EXT);
      f_out = fopen(fmod_name,"w");
      strcpy(fmod_name,"\nfunction z=");
      strcat(fmod_name,fname);
      strcat(fmod_name,"_fff(y)\n");
      str_output(fmod_name);
    }

  str_output("z=zeros(");
  i_output(endo_nbr);
  str_output(",1);\n");

#ifdef MATLAB
  str_output("global ex_ ex_det_ it_ recur_\n");
#endif

  param_nbr = 0;
  for ( i=0; i < var_nbr; i++ )
    {
      if (var_list[i].endo_exo == 4){
	var_list[i].nbr=param_nbr;
	if ((param_nbr % 10)== 0 ) str_output("\nglobal ");
	str_output(var_list[i].name);
	str_output(" ");
	param_nbr++;
      }
    }
  str_output("\n");
#endif

  p_q=(struct queue **)p->ptr;
  k=1;
  for (i=0;i < p->imax;i++){
    p_t=(struct token **)(*p_q)->ptr;
    if ((*p_q)->flag == 0){
#ifdef GAUSS
      str_output("z[");
      i_output(k);
      str_output("] = ");
#elif defined MATLAB || defined SCILAB
      str_output("z(");
      i_output(k);
      str_output(") = ");
#endif
      k++;
    }
    for (j=0;j < (*p_q)->imax;j++){
      if ((*p_t)->nbr == -1)  str_output((*p_t)->name);
      else var_output(*p_t,flag_steady);
      p_t++;
    }
    str_output(";\n");
    p_q++;
  }
  // check number of equations
  check.eq_nbr = k-1;
#ifdef GAUSS
  str_output("retp(z);\n");
  str_output("endp;\n");
#elif defined MATLAB || defined SCILAB
  if (flag_steady == 1)
    {
      str_output("if ~isreal(z)\n");
      str_output("  z = real(z)+imag(z).^2;\n");
      str_output("end\n");
    }
  fclose(f_out);
  f_out = f_temp;
#endif
}

void var_output(struct token * p_t, int flag_steady)
{
  char buffer[200];
  int i;
  if (p_t->endo_exo == 1){
#ifdef GAUSS
    if (flag_steady == 0)
      {      
	sprintf(buffer,"y[%d]",iy[(y_max_lag+p_t->lead_lag)*(endo_nbr)+(p_t->nbr)]);
      }
    else
      {
	sprintf(buffer,"y[%d]",p_t->nbr+1);
      }
#elif defined MATLAB || defined SCILAB
    if (flag_steady == 0)
      {      
	sprintf(buffer,"y(%d)",iy[(y_max_lag+p_t->lead_lag)*(endo_nbr)+(p_t->nbr)]);
      }
    else
      {
	sprintf(buffer,"y(%d)",p_t->nbr+1);
      }
#endif
    str_output(buffer);
  }
  else if (p_t->endo_exo == 0){
    i=p_t->lead_lag+x_max_lag-y_max_lag;
#ifdef GAUSS
    if (flag_steady == 1) sprintf(buffer,"_exe[%d]",p_t->nbr+1);
    else if (i == 0) sprintf(buffer,"_ex[_it,%d]",p_t->nbr+1);
    else sprintf(buffer,"_ex[_it%+d,%d]",i,p_t->nbr+1);
#elif defined MATLAB || defined SCILAB
    if (flag_steady == 1) sprintf(buffer,"exe_[%d]",p_t->nbr+1);
    if (i == 0) sprintf(buffer,"ex_(it_,%d)",p_t->nbr+1);
    else sprintf(buffer,"ex_(it_%+d,%d)",i,p_t->nbr+1);
#endif
    str_output(buffer);
  }
  else if (p_t->endo_exo == 3){
    i=p_t->lead_lag+z_max_lag-y_max_lag;
#ifdef GAUSS
    if (flag_steady == 1) sprintf(buffer,"_recurs[%d]",p_t->nbr+1);
    else if (i == 0) sprintf(buffer,"_recur[_it,%d]",p_t->nbr+1);
    else sprintf(buffer,"_recur[_it%+d,%d]",i,p_t->nbr+1);
#elif defined MATLAB || defined SCILAB
    if (flag_steady == 1) sprintf(buffer,"recurs_[%d]",p_t->nbr+1);
    else if (i == 0) sprintf(buffer,"recur_(it_,%d)",p_t->nbr+1);
    else sprintf(buffer,"recur_(it_%+d,%d)",i,p_t->nbr+1);
#endif
    str_output(buffer);
  }
  else if (p_t->endo_exo == 4){
    str_output(p_t->name);
  }
  else if (p_t->endo_exo == 5){
    i=p_t->lead_lag;
#ifdef GAUSS
    if (flag_steady == 1) sprintf(buffer,"_exe_det[%d]",p_t->nbr+1);
    else if (i == 0) sprintf(buffer,"_ex_det[_it,%d]",p_t->nbr+1);
    else sprintf(buffer,"_ex_det[_it%+d,%d]",i,p_t->nbr+1);
#elif defined MATLAB || defined SCILAB
    if (flag_steady == 1) sprintf(buffer,"exe_det_[%d]",p_t->nbr+1);
    if (i == 0) sprintf(buffer,"ex_det_(it_,%d)",p_t->nbr+1);
    else sprintf(buffer,"ex_det_(it_%+d,%d)",i,p_t->nbr+1);
#endif
    str_output(buffer);
  }
}

void p_i_shocks(void)
{
#ifdef GAUSS
  str_output("/* (M)SHOCKS */\n");
  str_output("if not _valf;\n");
  str_output("  _ex[_xkmin+1:_iter+_xkmin+_xkmax,.]=ones(_iter+_xkmax,1).*_exe';\n");
  str_output("endif;\n");
#elif defined MATLAB
  str_output("% (M)SHOCKS \n");
  str_output("make_ex_;\n");
#elif defined SCILAB
  str_output("// (M)SHOCKS \n");
  str_output("make_ex_();\n");
#endif
}

void p_e_shocks(void)
{
  char buffer[2000];
#if defined MATLAB || defined SCILAB
  sprintf(buffer,"M_.ex_det_length = %d;\n",ex_det_length);
  str_output(buffer);
#endif
}

void dynare_init(char* fname,struct s_runtime_options runtime_options)
{
  char buf[200];
#ifdef GAUSS
#ifdef OXGAUSS_GNU
  str_output("library pgraph_gnu;\n");
#else
  str_output("library pgraph;\n");
#endif
#ifdef OXGAUSS
  str_output("#include dynare1.src;\n");
  str_output("#include dynare2.src;\n");
  str_output("#include dynare3.src;\n");
#endif
  str_output("clear _scalv,_ex,_ys,_y,_exe,_recur,_recurs,_lgy,_lgx,_lgr;\n");
  str_output("_dsmpl=0;\n_dynatol=0.00001;\n_maxit=10;\n_slowc=1;\n");
  sprintf(buf,"_logname = \"%s.log\";\n",fname);
  str_output(buf);
  sprintf(buf,"output file = %s.log reset;\n",fname);
  str_output(buf);
#elif defined MATLAB
  if (runtime_options.clearall == 1)
    {
      str_output("clear all\n");
    }
  str_output("global scalv_ ex_ ex_det_ recur_ recurs_ ys_ y_ exe_ exe_det_ lgy_ lgx_ lgr_ dsmpl_ endval_\n");
  str_output("global endo_nbr exo_nbr exo_det_nbr iy_  ykmin_  ykmax_  xkmin_  xkmax_ zkmin_ zkmax_ iter_\n"); 
  str_output("global dynatol_ slowc_ maxit_ valf_ ys0_ recurs0_ ex0_ timing_ ct_ gstep_ Sigma_e_ fname_ lgx_orig_ord_ iter_ options_ dr_ oo_ trend_coeff_ eigenvalues_\n");
  str_output("global M_\n");
  str_output("M_.ex_det_length = 0;\n");
  str_output("M_.dname = '");
  str_output(fname);
  str_output("';\n");
  str_output("dsmpl_=0;\ndynatol_=0.00001;\nmaxit_=10;\nslowc_=1;\ntiming_=0;\nct_=0;\ngstep_=1e-2;\n");
  str_output("endval_=0;rplottype_=0;\nvalf_=0;\n");
  str_output("y_=[];\nex_=[];\nex_det_=[];\n");
  sprintf(buf,"fname_ = '%s';\n",fname);
  str_output(buf);
  sprintf(buf,"logname_ = '%s.log';\n",fname);
  str_output(buf);
  sprintf(buf,"diary off;\nwarning off;\ndelete %s.log;\nwarning on;\nwarning backtrace;\n",fname);
  str_output(buf);
  sprintf(buf,"diary %s.log;\n",fname);
  str_output(buf);
  str_output("options_ = [];\n");
#elif defined SCILAB
  str_output("clear\nclearglobal()\n");
  str_output("global ex_ ex_det_ y_ ys_ dsmpl_ i_exo_var_ iter_ options_ dr_\n");
  str_output("dsmpl_=0;\ndynatol_=0.00001;\nmaxit_=10;\nslowc_=1;\ntiming_=0;\nct_=0;\ngstep_=1e-2;\ndebug_=0;\nrecurs_=[];\nstart_simul=[];\n");
  str_output("endval_=0;rplottype_=0;\nvalf_=0;\n");
  sprintf(buf,"fname_ = '%s';\n",fname);
  str_output(buf);
  sprintf(buf,"logname_ = '%s.log';\n",fname);
  str_output(buf);
  str_output("fh_log=mopen(logname_,'w');\n");
  sprintf(buf,"getf('%s_ff.sci');\n",fname);
  str_output(buf);
  sprintf(buf,"getf('%s_fff.sci');\n",fname);
  str_output(buf);
  str_output("options_=tlist('opt');\n");
#endif

  varexo=strdup("");
  
}


void add_tmpvar(name1,name2)
char *name1,*name2;
{
  int static tmpvar_size = 0;

  if ( nbr_tmpvar+2 > tmpvar_size )
    {
      tmpvar_list = realloc(tmpvar_list, 200*sizeof(char *));
      tmpvar_size += 200;
    }

  tmpvar_list[nbr_tmpvar] = name1;
  nbr_tmpvar++;
  tmpvar_list[nbr_tmpvar] = name2;
  nbr_tmpvar++;
}

void dyn2vec(flag)
     int flag;                           /* flag: one or two names */
{
  int i;
  char buffer[200];
  if (flag == 0)
    {
#if defined GAUSS
      str_output("dyn2vec;\n");
#elif defined MATLAB || defined SCILAB
      str_output("global ");
      for(i=0;i < var_nbr; i++)
	{
	  if (var_list[i].endo_exo == 1)
	    {
	      str_output(var_list[i].name);
	      str_output(" ");
	    }
	}
      str_output("\ndyn2vec;\n");
#endif
    }
  else
    {
#ifdef GAUSS
      str_output("clear ");
      str_output(tmpvar_list[(tmpvar_list[1]?1:0)]);
      for( i = 2; i < nbr_tmpvar; i+=2)
	{
	  str_output(", ");
	  str_output(tmpvar_list[i+(tmpvar_list[i+1]?1:0)]);
	}
      str_output(";\n");
      for( i = 0; i < nbr_tmpvar; i+=2)
	{
	  if (tmpvar_list[i+1])
	    {
	      sprintf(buffer,"dyn2vec %s %s;\n",tmpvar_list[i],tmpvar_list[i+1]);
	    }
	  else
	    {
	      sprintf(buffer,"dyn2vec %s;\n",tmpvar_list[i]);
	    }
	  str_output(buffer);
	}
#elif defined MATLAB || defined SCILAB
      str_output("global ");
      for( i = 0; i < nbr_tmpvar; i+=2)
	{
	  str_output(tmpvar_list[i+(tmpvar_list[i+1]?1:0)]);
	  str_output(" ");
	}
      str_output("\n");
      for( i = 0; i < nbr_tmpvar; i+=2)
	{
	  if (tmpvar_list[i+1] > 0)
	    {
	      sprintf(buffer,"dyn2vec('%s', '%s');\n",tmpvar_list[i],tmpvar_list[i+1]);
	    }
	  else
	    {
	      sprintf(buffer,"dyn2vec('%s');\n",tmpvar_list[i]);
	    }
	  str_output(buffer);
	}
#endif
    }

}

void print_rplot(void)
{
  int i;
  char buffer[200];

#ifdef GAUSS
  strcpy(buffer,"rplot");
  for(i=0;i<nbr_tmpvar;i+=2)
    {
      strcat(buffer," ");
      strcat(buffer,tmpvar_list[i]);
    }
  strcat(buffer," ;\n");
#elif defined MATLAB || defined SCILAB
  strcpy(buffer,"plot_list_ = '");
  strcat(buffer,tmpvar_list[0]);
  strcat(buffer,"';\n");
  for(i=2;i<nbr_tmpvar;i+=2)
    {
#ifdef MATLAB
      strcat(buffer,"plot_list_ = str2mat(plot_list_,'");
      strcat(buffer,tmpvar_list[i]);
      strcat(buffer,"');\n");
#else
      strcat(buffer,"plot_list_ = [plot_list_;'");
      strcat(buffer,tmpvar_list[i]);
      strcat(buffer,"'];\n");
#endif
    }
  strcat(buffer,"rplot(plot_list_, rplottype_) ;\n");
#endif
  str_output(buffer);
}


void print_iter(char *s)
{
#ifdef GAUSS
  fprintf(f_out,"\n_iter = %s;\n",s);
#elif defined MATLAB || defined SCILAB
  fprintf(f_out,"\niter_ = %s;\n",s);
  fprintf(f_out,"options_.periods = %s;\n",s);
#endif
}

struct loop *initial_loop (char *name, char *initial, char *limit, char *increment)
{
  struct loop *loop;

  loop = (struct loop *)calloc(1,sizeof(struct loop));

  loop->loop_index = atoi(initial);
  loop->loop_limit = atoi(limit);
  loop->loop_increment = atoi(increment);
  loop->loop_index_name = strdup(name);

  add_var(name,2);
  var_list[var_nbr-1].nbr = 0;
  sort_var();

  return loop;
}

struct queue *do_loop(struct loop *loop, struct queue *q_eq1)
{
  struct queue *q_eq2,**p_q1,**p_q2;
  int i;

  q_eq2 = (struct queue *)calloc(1,sizeof(struct queue));

  q_eq2->ptr = (void **)calloc(1,sizeof(void *));
  q_eq2->imax = 0;
  while (loop->loop_index <= loop->loop_limit)
    {
      q_eq2->ptr=(void **)realloc(q_eq2->ptr,(q_eq2->imax+q_eq1->imax)*sizeof(void *));
      p_q1=(struct queue **)q_eq1->ptr;
      p_q2=(struct queue **)q_eq2->ptr+q_eq2->imax;
      q_eq2->imax += q_eq1->imax;
      for (i=0;i < q_eq1->imax;i++){
	*p_q2=(struct queue *)calloc(1,sizeof(struct queue));
	(*p_q2)->flag=(*p_q1)->flag;
	(*p_q2)->ptr=(void **)calloc((*p_q1)->imax,sizeof(struct token *));
	(*p_q2)->imax=0;
	copy_update_indexed_var((struct token **)(*p_q1)->ptr, *p_q2,(*p_q1)->imax,loop->loop_index_name,loop->loop_index);
	p_q1++;
	p_q2++;
      }
      loop->loop_index++;
    }

  return q_eq2;
}

void copy_update_indexed_var(struct token **p_t1, struct queue  *q2, int imax, char *index_name, int index)
{
  struct token **p_t3,**p_end;
  char buffer[200];
  int j, flag;
  struct var *p_var;

  p_end = p_t1+imax;
  for (j=0;j < imax;j++){
    if ((*p_t1)->nbr == -3){
      flag=1;
      for (p_t3=p_t1+1; p_t3 < p_end; p_t3++){
	if ((*p_t3)->nbr != -3){
	  break;
	}
	if ((*p_t3)->name != NULL ){
	  if (strcmp(index_name,(*p_t3)->name) == 0){
	    (*p_t3)->lead_lag = index;
	  }
	  else{
	    flag=0;
	  }
	}
      }
      p_t3=p_t1;
      if (flag == 1){
	strcpy(buffer,(*p_t3)->name);
	flag = 0;
	for (p_t3=p_t1+1; p_t3 < p_end; p_t3++){
	  if ((*p_t3)->nbr != -3){
	    break;
	  }
	  if ((*p_t3)->name == NULL || strcmp(index_name,(*p_t3)->name) == 0){
	    if (flag == 0){
	      sprintf(buffer+strlen(buffer),"%d\0",(*p_t3)->lead_lag);
	      flag = 1;
	    }
	    else{
	      sprintf(buffer+strlen(buffer),"_%d\0",(*p_t3)->lead_lag);
	    }
	  }
	  else{
	    strcat(buffer,(*p_t3)->name);
	  }
	}
	p_var=var_search(buffer);
	if (p_var == NULL){
	  fprintf(stdout,"%s is not an indexed variable\n",buffer);
	}
	else if (p_var->endo_exo == 3){
	  add_to_queue(q2,token(p_var->name,-1));
	}
	else{
	  add_to_queue(q2,create_var(p_var));
	}
	set_ll1((struct token *)(q2->ptr)[q2->imax-1],(*p_t1)->lead_lag);
      }
      else{
	add_to_queue(q2,token((*p_t3)->name,(*p_t3)->nbr));
	for (p_t3=p_t1+1; p_t3 < p_end; p_t3++){
	  if ((*p_t3)->nbr != -3){
	    break;
	  }
	  if ((*p_t3)->name == NULL || strcmp(index_name,(*p_t3)->name) == 0){
	    add_to_queue(q2,token(NULL,-3));
	    ((struct token *)(q2->ptr)[q2->imax-1])->lead_lag=(*p_t3)->lead_lag;
	  }
	  else{
	    add_to_queue(q2,token((*p_t3)->name,-3));
	  }
	}
      }      
      j += p_t3-p_t1-1;
      p_t1 = p_t3-1;
    }
    else if ((*p_t1)->nbr == -2 && strcmp(index_name,(*p_t1)->name) == 0){
      sprintf(buffer,"%d\0",index);
      add_to_queue(q2,token(strdup(buffer),-1));
      ((struct token *)(q2->ptr)[q2->imax-1])->lead_lag=(*p_t1)->lead_lag;
      ((struct token *)(q2->ptr)[q2->imax-1])->endo_exo=(*p_t1)->endo_exo;
    }
    else if ((*p_t1)->nbr == -1){
      add_to_queue(q2,token((*p_t1)->name,-1));
      ((struct token *)(q2->ptr)[q2->imax-1])->lead_lag=(*p_t1)->lead_lag;
      ((struct token *)(q2->ptr)[q2->imax-1])->endo_exo=(*p_t1)->endo_exo;
    }
    else {
      if ((*p_t1)->name == NULL || (strcmp(index_name,(*p_t1)->name) != 0)){
	add_to_queue(q2,token((*p_t1)->name,(*p_t1)->nbr));
	((struct token *)(q2->ptr)[q2->imax-1])->lead_lag=(*p_t1)->lead_lag;
      }
      else{
	add_to_queue(q2,token(NULL,(*p_t1)->nbr));
	set_ll1((struct token *)(q2->ptr)[q2->imax-1],index);
      }
      ((struct token *)(q2->ptr)[q2->imax-1])->lead_lag=(*p_t1)->lead_lag;
      ((struct token *)(q2->ptr)[q2->imax-1])->endo_exo=(*p_t1)->endo_exo;
    }
    p_t1++;
  }
}

struct queue *operator_loop(struct loop *loop, struct queue *q1, char *oper)
{
  struct queue *q2;

  q2=create_queue(token("((",-1));
  copy_update_indexed_var((struct token **)q1->ptr, q2,q1->imax,loop->loop_index_name,loop->loop_index);
  loop->loop_index++;
  add_to_queue(q2,token(")",-1));
  while (loop->loop_index <= loop->loop_limit){
      add_to_queue(q2,token(oper,-1));
      q2=add_to_queue(q2,token("(",-1));
      copy_update_indexed_var((struct token **)q1->ptr, q2,q1->imax,loop->loop_index_name,loop->loop_index);
      add_to_queue(q2,token(")",-1));
      loop->loop_index++;
  }
  add_to_queue(q2,token(")",-1));
     
  return q2;
}

void p_steady(struct queue *p)
{
  //  print_model1(p,1);
  str_output("steady(0);\n");
}

void p_steady_linear(struct queue *p)
{
  //  print_model1(p,1);
  str_output("steady(1);\n");
}

void print_param(void)
{
  int i;

#if defined MATLAB || defined SCILAB
  param_nbr = 0;
  for ( i=0; i < var_nbr; i++ )
    {
      if (var_list[i].endo_exo == 4){
	var_list[i].nbr=param_nbr;
	if ((param_nbr % 10)== 0 ) str_output("\nglobal ");
	str_output(var_list[i].name);
	str_output(" ");
	param_nbr++;
      }
    }
  str_output("\n");
#endif
}

void p_stoch_simul()
{
  char buffer[2000];
  int i;
  check.stoch_simul = 1;
#ifdef SCILAB
  str_output("global ");
  for(i=0;i < var_nbr; i++)
    {
      if (var_list[i].endo_exo == 1)
	{
	  str_output(var_list[i].name);
	  str_output(" ");
	}
    }
  str_output("\n");
#endif
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = [];\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  str_output("info=stoch_simul(var_list_);\n");
#endif
}

void p_resol(int dr_algo, int linear, int order, struct queue *p)
{
  char buffer[200];
  //  print_model1(p,1);
  sprintf(buffer,"dr_ = resol(ys_,%d,%d,%d);\n",dr_algo,linear,order);
  str_output(buffer);
}

void p_simul( void )
{
  check.simul = 1;
  str_output("simul(dr_);\n");
}

void p_dsample(int nargs, char *arg1, char *arg2)
{
#ifdef GAUSS
  str_output("dsample ");
  if (nargs > 0)
    {
      str_output(arg1);
    }
  if (nargs > 1)
    {
      str_output(" ");
      str_output(arg2);
    }
#elif defined MATLAB || defined SCILAB
  str_output("dsample");
  if (nargs > 0)
    {
      str_output("(");
      str_output(arg1);
      if (nargs > 1)
	{
	  str_output(",");
	  str_output(arg2);
	}
      str_output(")");
    }
#endif  
  str_output(";\n");
}

void p_check(void)
{
  check.check = 1;
#if defined MATLAB || defined GAUSS
  str_output("check;\n");
#else
  str_output("check();\n");
#endif
}

void p_irf(char * varexo, double shock_size, int iter, int drop, int replic, int order)
{
  char buffer[200];
  int i;
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = '';\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  sprintf(buffer,"y_=irf(dr_,'%s',%f,%d,%d,%d,%d,var_list_);\n",varexo,shock_size,iter+drop,drop,replic,order);
#endif
  str_output(buffer);
}

void p_d_corr()
{
  char buffer[2000];
  int i,flag;
  strcpy(buffer,"d_corr(");
  flag = 0;
  if (nbr_tmpvar == 0)
    {
      for (i=0; i < var_nbr; i++)
	{
	  if ( var_list[i].endo_exo == 1)
	    {
	      if (flag > 0)
		{
		  strcat(buffer,",");
		}
	      strcat(buffer,"'");
	      strcat(buffer,var_list[i].name);
	      strcat(buffer,"'");
	      flag = 1;
	    }
	}
    }
  else
    {
      for(i=0;i<nbr_tmpvar;i+=2)
	{
	  if (flag > 0)
	    {
	      strcat(buffer,",");
	    }
	  strcat(buffer,"'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'");
	  flag = 1;
	}
    }
  strcat(buffer,");\n");
  str_output(buffer);
}

void p_disp_dr(order)
{
  char buffer[200];
  int i;
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = '';\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  sprintf(buffer,"disp_dr(dr_,%d,var_list_);\n",order);
#endif
  str_output(buffer);
}

void p_disp_moments(order)
{
  char buffer[200];
  int i;
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = '';\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  sprintf(buffer,"disp_moments(y_,%d,var_list_);\n",order);
#endif
  str_output(buffer);
}

void p_optim_weights_init(void)
{
  str_output(OPEN_COMMENTS);
  str_output(" OPTIM_WEIGHTS\n");
  str_output(CLOSE_COMMENTS);

#if defined MATLAB
  str_output("optim_weights_ = sparse(endo_nbr,endo_nbr);\n");
  str_output("obj_var_ = [];\n\n");
#endif
}

void p_optim_weights1(struct token* p_t, struct queue* p_q)
{
  char buffer[200];

  if (p_t->endo_exo != 1 && p_t->endo_exo != 3)
    {
      fprintf(stdout,"OPTIM_WEIGHTS ERROR: only endogenous variables can have weights" );
      exit(1);
    }

  sprintf(buffer,"optim_weights_(%d,%d) = ",p_t->nbr+1,p_t->nbr+1);
  str_output(buffer);
  p_expression(p_q);
  str_output(";\n");
  sprintf(buffer,"obj_var_ = [obj_var_; %d];\n",p_t->nbr+1);
  str_output(buffer);
}

void p_optim_weights2(struct token* p_t1, struct token* p_t2, struct queue* p_q)
{
  char buffer[200];

  if ((p_t1->endo_exo != 1 && p_t1->endo_exo != 3) || (p_t2->endo_exo != 1 && p_t2->endo_exo != 3))
    {
      fprintf(stdout,"OPTIM_WEIGHTS ERROR: only endogenous variables can have weights" );
      exit(1);
    }

  sprintf(buffer,"optim_weights_(%d,%d) = ",p_t1->nbr+1,p_t2->nbr+1);
  str_output(buffer);
  p_expression(p_q);
  str_output(";\n");
  sprintf(buffer,"optim_weights_(%d,%d) = optim_weights_(%d,%d);\n",p_t2->nbr+1,p_t1->nbr+1,p_t1->nbr+1,p_t2->nbr+1);
  str_output(buffer);
  sprintf(buffer,"obj_var_ = [obj_var_; %d];\n",p_t1->nbr+1);
  sprintf(buffer,"obj_var_ = [obj_var_; %d];\n",p_t2->nbr+1);
}

void p_expression(struct queue* q)
{
  struct token** p_t;
  char buffer[200];
  int j;

  p_t=(struct token **)q->ptr;
  for (j=0;j < q->imax;j++){
    if ((*p_t)->nbr == -1)  str_output((*p_t)->name);
    else if ((*p_t)->endo_exo == 1)
      {
#ifdef GAUSS
	sprintf(buffer,"_ys[%d]",(*p_t)->nbr+1);
#elif defined MATLAB || defined SCILAB
	sprintf(buffer,"ys_(%d)",(*p_t)->nbr+1);
#endif
	str_output(buffer);
      }
    else if ((*p_t)->endo_exo == 0)
      {
#ifdef GAUSS
	sprintf(buffer,"_exe[%d]",(*p_t)->nbr+1);
#elif defined MATLAB || defined SCILAB
	sprintf(buffer,"exe_(%d)",(*p_t)->nbr+1);
#endif
 	str_output(buffer);
      }
    else if ((*p_t)->endo_exo == 3)
      {
#ifdef GAUSS
	sprintf(buffer,"_recurs[%d]",(*p_t)->nbr+1);
#elif defined MATLAB || defined SCILAB
	sprintf(buffer,"recurs_(%d)",(*p_t)->nbr+1);
#endif
 	str_output(buffer);
      }
    else if ((*p_t)->endo_exo == 4)
      {
	sprintf(buffer,"%s",(*p_t)->name);
 	str_output(buffer);
      }
    else if ((*p_t)->endo_exo == 5)
      {
#ifdef GAUSS
	sprintf(buffer,"_exe_det[%d]",(*p_t)->nbr+1);
#elif defined MATLAB || defined SCILAB
	sprintf(buffer,"exe_det_(%d)",(*p_t)->nbr+1);
#endif
 	str_output(buffer);
      }
    p_t++;
  }
}

void p_osr_params(char * name)
{
  static int flag = 0;
  char buffer[200];

#ifdef GAUSS
  if ( flag == 0)
    {
      sprintf(buffer,"let string _osr_params = \"%s\";\n",name);
      flag = 1;
    }
  else
    {
      sprintf(buffer,"osr_params = _osr_params | \"%s\";\n",name);
    }
#elif defined MATLAB
  if ( flag == 0)
    {
      sprintf(buffer,"osr_params_ = '%s';\n",name);
      flag = 1;
    }
  else
    {
      sprintf(buffer,"osr_params_ = str2mat(osr_params_,'%s');\n",name);
    }
#elif defined SCILAB
  if ( flag == 0)
    {
      sprintf(buffer,"osr_params_ = '%s';\n",name);
      flag = 1;
    }
  else
    {
      sprintf(buffer,"osr_params_ = [osr_params_;'%s'];\n",name);
    }
#endif
  str_output(buffer);
  str_output("\n");
}

void p_calib_init(void)
{
  int i;
  char buffer[200];

  str_output(OPEN_COMMENTS);
  str_output(" CALIB_VAR ");
  str_output(CLOSE_COMMENTS);
  str_output("\n\n");

#if defined MATLAB || defined SCILAB
  for(i=1;i<4;++i)
    {
      sprintf(buffer,"calib_var_index{%d} = [];\ncalib_targets{%d} = [];\ncalib_weights{%d}=[];\n",i,i,i);
      str_output(buffer);
    }
#endif
}

void p_calib_var(struct token* p_t, struct queue* p_q, char * weight)
{
  char buffer[200];
  if (p_t->endo_exo == 1 || p_t->endo_exo == 3)
    {
      sprintf(buffer,"calib_var_index{1} = [calib_var_index{1};%d];\n",p_t->nbr+1);
      str_output(buffer);
      str_output("calib_targets{1} =[calib_targets{1}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{1} = [calib_weights{1}; %s];\n",weight);
      str_output(buffer);
    }
  else if (p_t->endo_exo == 0)
    {
      sprintf(buffer,"calib_var_index{3} = [calib_var_index{3};%d %d];\n",p_t->nbr+1,p_t->nbr+1);
      str_output(buffer);
      str_output("calib_targets{3} =[calib_targets{3}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{3} = [calib_weights{3}; %s];\n",weight);
      str_output(buffer);
    }
  else
    {
      printf("ERROR in CALIB: one of the targets isn't legitimate\n");
    }
}

void p_calib_covar(struct token* p_t1, struct token* p_t2, struct queue* p_q, char * weight)
{
  char buffer[200];
  if ((p_t1->endo_exo == 0 && p_t2->endo_exo == 1)|| (p_t1->endo_exo == 1 && p_t2->endo_exo == 0))
    {
      printf("ERROR in CALIB: can't target correlation between an endogenous and an exogenous variable\n");
      exit(1);
    }
  else if (p_t1->endo_exo == 1 || p_t1->endo_exo == 3)
    {
      sprintf(buffer,"calib_var_index{2} = [calib_var_index{2};%d %d];\n",p_t1->nbr+1,p_t2->nbr+1);
      str_output(buffer);
      str_output("calib_targets{2} =[calib_targets{2}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{2} = [calib_weights{2}; %s];\n",weight);
      str_output(buffer);
    }
  else if (p_t1->endo_exo == 0)
    {
      sprintf(buffer,"calib_var_index{3} = [calib_var_index{2};%d %d];\n",p_t1->nbr+1,p_t2->nbr+1);
      str_output(buffer);
      str_output("calib_targets{3} =[calib_targets{3}; ");
      p_expression(p_q);
      str_output("];\n");
      sprintf(buffer,"calib_weights{3} = [calib_weights{3}; %s];\n",weight);
      str_output(buffer);
    }
  else
    {
      printf("ERROR in CALIB: one of the targets isn''t legitimate\n");
      exit(1);
    }
}
 
static int max_iar = 3;
void p_calib_ac(struct token* p_t, char * ar, struct queue* p_q, char * weight)
{
  char buffer[200];
  int i, iar;

  iar = atoi(ar)+3;
  if (iar > max_iar)
    {
      for(i=max_iar+1; i <= iar; ++i)
	{
	  sprintf(buffer,"calib_var_index{%d} = [];\ncalib_targets{%d} = [];\ncalib_weights{%d}=[];\n",i,i,i);
	  str_output(buffer);
	}
      max_iar = iar;
    }
  sprintf(buffer,"calib_var_index{%d} = [calib_var_index{%d};%d];\n",iar,iar,p_t->nbr+1);
  str_output(buffer);
  sprintf(buffer,"calib_targets{%d} =[calib_targets{%d}; ",iar,iar);
  str_output(buffer);
  p_expression(p_q);
  str_output("];\n");
  sprintf(buffer,"calib_weights{%d} = [calib_weights{%d}; %s];\n",iar,iar,weight);
  str_output(buffer);
}

void p_calib(int cova)
{
  char buffer[200];
  sprintf(buffer,"Sigma_e_=calib(calib_var_index,calib_targets,calib_weights,%d,%d,Sigma_e_);\n",max_iar-3,cova);
  str_output(buffer);
}

void p_dynatype(char *fname, char *ext)
{
  char buffer[2000];
  int i;
#ifdef SCILAB
  str_output("global ");
  for(i=0;i < var_nbr; i++)
    {
      if (var_list[i].endo_exo == 1)
	{
	  str_output(var_list[i].name);
	  str_output(" ");
	}
    }
  str_output("\n");
#endif
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = [];\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  if (strlen(ext) == 0)
    {
      sprintf(buffer,"dynatype('%s',var_list_);\n",fname);
    }
  else
    {
      sprintf(buffer,"dynatype('%s.%s',var_list_);\n",fname,ext);
    }
#else
  if (strlen(ext) == 0)
    {
      sprintf(buffer,"dynatype %s;\n",fname);
    }
  else
    {
      sprintf(buffer,"dynatype %s.%s;\n",fname,ext);
    }

#endif
  str_output(buffer);
}

void p_dynasave(char *fname, char *ext)
{
  char buffer[2000];
  int i;
#ifdef SCILAB
  str_output("global ");
  for(i=0;i < var_nbr; i++)
    {
      if (var_list[i].endo_exo == 1)
	{
	  str_output(var_list[i].name);
	  str_output(" ");
	}
    }
  str_output("\n");
#endif
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = [];\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  if (strlen(ext) == 0)
    {
      sprintf(buffer,"dynasave('%s',var_list_);\n",fname);
    }
  else
    {
      sprintf(buffer,"dynasave('%s.%s',var_list_);\n",fname,ext);
    }
#else
  if (strlen(ext) == 0)
    {
      sprintf(buffer,"dynasave %s;\n",fname);
    }
  else
    {
      sprintf(buffer,"dynasave %s.%s;\n",fname,ext);
    }

#endif
  str_output(buffer);
}

void p_sigma_e(struct queue *p_q)
{
  int up, i, j, i1, j1, p1, p2;
  int* order;
  int* po;
  struct queue **p_q1;
#if defined MATLAB || defined SCILAB
  order = (int*)calloc(exo_nbr,sizeof(int));
  po = order;
  for (i = 0; i < var_nbr; ++i)
    {
      if (var_list[i].endo_exo == 0)
	{
	  *po++ = var_list[i].original_nbr;
	}
    }
  str_output("Sigma_e_ = [\n");
  p_q1 = (struct queue **)(p_q->ptr);
  up = ((*p_q1)->imax > 1) ? 1 : 0;
  for(i1=0; i1 < exo_nbr; ++i1)
    {
      if (i1 > 0)
	{
	  str_output(";\n");
	}
      i = order[i1];
      for(j1=0; j1 < exo_nbr; ++j1)
	{
	  if (j1 > 0)
	    {
	      str_output(", ");
	    }
	  j = order[j1];
	  if (j <= i)
	    {
	      if (up)
		{
		  p1 = j;
		  p2 = i-j;
		}
	      else
		{
		  p1 = i;
		  p2 = j;
		}
	    }
	  else
	    {
	      if (up)
		{
		  p1 = i;
		  p2 = j-i;
		}
	      else
		{
		  p1 = j;
		  p2 = i;
		}
	    }
		  p_expression(*(struct queue**)((*(struct queue**)(p_q1+p1))->ptr+p2));
	}
    }
  str_output("\n];\n");
  free(order);
#endif
}

void tr_check(int row, int col)
{
  static int old_col = 0;
  static int up_down = 0;
  if (col == 0)
    {
      fprintf(stdout,"ERROR: empty row in Sigma_e!\n");
      exit(1);
    }
  if ( row > exo_nbr || col > exo_nbr)
    {
      fprintf(stdout,"ERROR: Sigma_e has more rows or columns than exogenous variables in the model!\n");
      exit(1);
    }
  if (old_col == 0)
    {
      old_col = col;
      up_down = (col > 1) ? -1 : 1;
    }
  else
    {
      old_col += up_down;
      if ( col != old_col)
	{
	  fprintf(stdout,"ERROR: Sigma_e isn't in triangular form!\n");
	  exit(1);
	}
    }
}

void p_option(char *name, char *value)
{
  char buffer[2000];
  sprintf(buffer,"options_.%s=%s;\n",name,value);
  str_output(buffer);
}

void p_option_e(char *name, struct queue* expression)
{
  char buffer[2000];
  sprintf(buffer,"options_.%s=",name);
  str_output(buffer);
  p_expression(expression);
  str_output(";\n");
}

void p_s_option(char *name, char *value)
{
  char buffer[2000];
  sprintf(buffer,"options_.%s='%s';\n",name,value);
  str_output(buffer);
}

void p_osr()
{
  char buffer[2000];
  int i;
  check.osr = 1;
#ifdef SCILAB
  str_output("global ");
  for(i=0;i < var_nbr; i++)
    {
      if (var_list[i].endo_exo == 1)
	{
	  str_output(var_list[i].name);
	  str_output(" ");
	}
    }
  str_output("\n");
#endif
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = [];\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  str_output("osr(var_list_,osr_params_,optim_weights_);\n");
#endif
}

void p_olr()
{
  char buffer[2000];
  int i;
  check.olr = 1;
#ifdef SCILAB
  str_output("global ");
  for(i=0;i < var_nbr; i++)
    {
      if (var_list[i].endo_exo == 1)
	{
	  str_output(var_list[i].name);
	  str_output(" ");
	}
    }
  str_output("\n");
#endif
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = [];\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  str_output("olr(var_list_,olr_inst_,obj_var_,optim_weights_);\n");
#endif
}

void p_olr_inst_init()
{
  str_output("%OLR_INST\n");
  str_output("olr_inst_ = [];\n");
}

void p_olr_inst(char* name)
{
  char buffer[2000];
  sprintf(buffer,"olr_inst_ = strvcat(olr_inst_,'%s');\n",name);
  str_output(buffer);
  inst_nbr++;
}

void run_checks(void)
 {
  int i,flag;
  if (check.olr && check.eq_nbr != endo_nbr - inst_nbr)
    {
      printf("ERROR: for OLR there must be %d equations for %d endogenous variables less %d instruments\n",endo_nbr-inst_nbr,endo_nbr,inst_nbr);
      exit(1);
    }
  else if (check.eq_nbr != endo_nbr - inst_nbr)
    {
      printf("ERROR: there are %d equations for %d endogenous variables\n",check.eq_nbr,endo_nbr);
      exit(1);
    }

  if (!check.linear && !check.initval)
    {
      printf("WARNING: there is no INITVAL, all variables initialized at 0\n");
    }

  if (!check.linear && check.initval)
    {
      flag = 1;
      for(i=0;i<var_nbr;++i)
	{
	  if (initval_check[i] == 0 && ((var_list[i].endo_exo == 0 && check.determ) || var_list[i].endo_exo == 1))
	    {
	      if (flag)
		{
		  printf("WARNING: the following variables aren't initialized in INITVAL and therefore set to 0\n");
		  flag =  0;
		}
	      printf("\t%s\n",var_list[i].name);
	    }
	}
    }             

  if (!check.linear && check.endval)
    {
      flag = 1;
      for(i=0;i<var_nbr;++i)
	{
	  if (endval_check[i] == 0 && ((var_list[i].endo_exo == 0 && check.determ) || var_list[i].endo_exo == 1))
	    {
	      if (flag)
		{
		  printf("WARNING: the following variables aren't initialized in ENDVAL and therefore set to 0\n");
		  flag =  0;
		}
	      printf("\t%s\n",var_list[i].name);
	    }
	}
              
    }
}

void p_estimated_init(void)
{
  str_output("global estim_params_\n");
  str_output("estim_params_.var_exo = [];\nestim_params_.var_endo = [];\nestim_params_.corrx = [];\nestim_params_.corrn = [];\nestim_params_.param_names = [];\nestim_params_.param_vals = [];\n");
}

void estim_params_init(void)
{
  estim_params.var_nbr = 0;
  estim_params.param_name=strdup("");
  estim_params.init_val=strdup("NaN");
  estim_params.mean=strdup("NaN");
  estim_params.std=strdup("NaN");
  estim_params.prior=strdup("NaN");
  estim_params.lb=strdup("-Inf");
  estim_params.ub=strdup("Inf");
  estim_params.p3=strdup("NaN");
  estim_params.p4=strdup("NaN");
  estim_params.jscale=strdup("NaN");
}

void p_estimated_elem(void)
{
      /* missing tests for arguments */
  char buffer[2000];

  if (estim_params.param_type == 1 && estim_params.var_type == 0)
    {
      str_output("estim_params_.var_exo = [estim_params_.var_exo; ");
      sprintf(buffer," %d,",estim_params.var_nbr);
      str_output(buffer);
    }
  else if(estim_params.param_type == 1 && estim_params.var_type == 1)
    {
      str_output("estim_params_.var_endo = [estim_params_.var_endo; ");
      sprintf(buffer," %d,",estim_params.var_nbr);
      str_output(buffer);
    }
  else if(estim_params.param_type == 2)
    {
      sprintf(buffer,"estim_params_.param_names = strvcat(estim_params_.param_names,'%s');\n",estim_params.param_name);
      str_output(buffer);
      str_output("estim_params_.param_vals = [estim_params_.param_vals;");
    }
  else if (estim_params.param_type == 3 && estim_params.var_type == 0)
    {
      str_output("estim_params_.corrx = [estim_params_.corrx; ");
      sprintf(buffer," %d,",estim_params.var_nbr);
      str_output(buffer);
      sprintf(buffer," %d,",estim_params.var_nbr2);
      str_output(buffer);
    }
  else if(estim_params.param_type == 3 && estim_params.var_type == 1)
    {
      str_output("estim_params_.corrn = [estim_params_.corrn; ");
      sprintf(buffer," %d,",estim_params.var_nbr);
      str_output(buffer);
      sprintf(buffer," %d,",estim_params.var_nbr2);
      str_output(buffer);
    }
  sprintf(buffer," %s",estim_params.init_val);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.lb);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.ub);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.prior);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.mean);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.std);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.p3);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.p4);
  str_output(buffer);
  sprintf(buffer,", %s",estim_params.jscale);
  str_output(buffer);
  str_output("];\n");
}
  
void p_estimated_elem1(struct token* id_var, char* init_val, char* lo_bound, char* up_bound, char* prior, char* p1, char *p2, char* p3, char *p4, char* jscale)
{
  char buffer[2000];
  if (id_var->endo_exo == 0)
    {
      /* missing tests for arguments */
      sprintf(buffer,"estim_params_.var_exo = [estim_params_.var_exo; %d %s %s %s %s %s %s %s %s %s];\n",id_var->nbr+1,init_val, lo_bound, up_bound, prior, p1, p2, p3, p4, jscale);
    }
  else if(id_var->endo_exo == 1)
    {
      /* missing tests for arguments */
      sprintf(buffer,"estim_params_.var_endo = [estim_params_.var_endo; %d %s %s %s %s %s %s %s %s %s];\n",id_var->nbr+1,init_val, lo_bound, up_bound, prior, p1, p2, p3, p4, jscale);
    }
  str_output(buffer);
}
 
void p_estimated_elem2(struct token* id_var1, struct token* id_var2, char* init_val, char* lo_bound, char* up_bound, char* prior, char* p1, char *p2, char* p3, char *p4, char* jscale)
{
  char buffer[2000];
  /* missing tests for arguments */
  if (id_var1->endo_exo == 0 && id_var2->endo_exo == 0)
    {
      sprintf(buffer,"estim_params_.corrx = [estim_params_.corrx; %d %d %s %s %s %s %s %s %s %s %s];\n",id_var1->nbr+1,id_var2->nbr+1,init_val, lo_bound, up_bound, prior, p1, p2, p3, p4, jscale);
    }
  else if (id_var1->endo_exo == 1 && id_var2->endo_exo == 1)
    {
      sprintf(buffer,"estim_params_.corrn = [estim_params_.corrn; %d %d %s %s %s %s %s %s %s %s %s];\n",id_var1->nbr+1,id_var2->nbr+1,init_val, lo_bound, up_bound, prior, p1, p2, p3, p4, jscale);
    }
  else
    {
      printf("Shocks and measurement errors can't be correlated\n");
      exit(1);
    }
  str_output(buffer);
}
 
void p_estimated_elem3(struct token* id_var, char* init_val, char* lo_bound, char* up_bound, char* prior, char* p1, char *p2, char* p3, char *p4, char *jscale)
{
  char buffer[2000];
  sprintf(buffer,"estim_params_.param_names = strvcat(estim_params_.param_names,'%s');\n",id_var->name);
  str_output(buffer);
  /* missing tests for arguments */
  sprintf(buffer,"estim_params_.param_vals = [estim_params_.param_vals; %s %s %s %s %s %s %s %s %s];\n",init_val, lo_bound, up_bound, prior, p1, p2, p3, p4, jscale);
  str_output(buffer);
}
 
void p_estimated_elem1a(struct token* id_var, char* prior, char* p1, char *p2, char* p3, char *p4, char* jscale)
{
  char buffer[2000];
  if (id_var->endo_exo == 0)
    {
      /* missing tests for arguments */
      sprintf(buffer,"estim_params_.var_exo = [estim_params_.var_exo; %d %s -Inf Inf %s %s %s %s %s %s];\n",id_var->nbr+1, p1, prior, p1, p2, p3, p4, jscale);
    }
  else if(id_var->endo_exo == 1)
    {
      /* missing tests for arguments */
      sprintf(buffer,"estim_params_.var_endo = [estim_params_.var_endo; %d %s -Inf Inf %s %s %s %s %s %s];\n",id_var->nbr+1, p1, prior, p1, p2, p3, p4, jscale);
    }
  str_output(buffer);
}
 
void p_estimated_elem2a(struct token* id_var1, struct token* id_var2, char* prior, char* p1, char *p2, char* p3, char *p4, char* jscale)
{
  char buffer[2000];
  /* missing tests for arguments */
  if (id_var1->endo_exo == 0 && id_var2->endo_exo == 0)
    {
      sprintf(buffer,"estim_params_.corrx = [estim_params_.corrx; %d %d %s -Inf Inf %s %s %s %s %s %s];\n",id_var1->nbr+1,id_var2->nbr+1, p1, prior, p1, p2, p3, p4, jscale);
    }
  else if (id_var1->endo_exo == 1 && id_var2->endo_exo == 1)
    {
      sprintf(buffer,"estim_params_.corrn = [estim_params_.corrn; %d %d %s -Inf Inf %s %s %s %s %s %s];\n",id_var1->nbr+1,id_var2->nbr+1, p1, prior, p1, p2, p3, p4, jscale);
    }
  else
    {
      printf("Shocks and measurement errors can't be correlated\n");
      exit(1);
    }
  str_output(buffer);
}
 
void p_estimated_elem3a(struct token* id_var, char* prior, char* p1, char *p2, char* p3, char *p4, char *jscale)
{
  char buffer[2000];
  sprintf(buffer,"estim_params_.param_names = strvcat(estim_params_.param_names,'%s');\n",id_var->name);
  str_output(buffer);
  /* missing tests for arguments */
  sprintf(buffer,"estim_params_.param_vals = [estim_params_.param_vals; %s -Inf Inf %s %s %s %s %s %s];\n", p1, prior, p1, p2, p3, p4, jscale);
  str_output(buffer);
}

void p_estimated_init_elem1(struct token* id_var, char* value)
{
  char buffer[2000];
  if (id_var->endo_exo == 0)
    {
      /* missing tests for arguments */
      sprintf(buffer,"tmp1 = find(estim_params_.var_exo(:,1)==%d);\nestim_params_.var_exo(tmp1,2) = %s;\n",id_var->nbr+1,value);
    }
  else if(id_var->endo_exo == 1)
    {
      /* missing tests for arguments */
      sprintf(buffer,"tmp1 = find(estim_params_.var_endo(:,1)==%d);\nestim_params_.var_endo(tmp1,2) = %s;\n",id_var->nbr+1,value);
    }
  str_output(buffer);
}
 
void p_estimated_init_elem2(struct token* id_var1, struct token* id_var2, char* value)
{
  char buffer[2000];
  /* missing tests for arguments */
  if (id_var1->endo_exo == 0 && id_var2->endo_exo == 0)
    {
      sprintf(buffer,"tmp1 = find(estim_params_.corrx(:,1)==%d & estim_params.corrx(:,2)==%d);\nestim_params_.corrx(tmp1,3) = %s;\n",id_var1->nbr+1,id_var2->nbr+1,value);
    }
  else if (id_var1->endo_exo == 1 && id_var2->endo_exo == 1)
    {
      sprintf(buffer,"tmp1 = find(estim_params_.corrn(:,1)==%d & estim_params.corrn(:,2)==%d);\nestim_params_.corrn(tmp1,3) = %s;\n",id_var1->nbr+1,id_var2->nbr+1,value);
    }
  else
    {
      printf("Shocks and measurement errors can't be correlated\n");
      exit(1);
    }
  str_output(buffer);
}
 
void p_estimated_init_elem3(struct token* id_var, char* value)
{
  char buffer[2000];
  /* missing tests for arguments */
  sprintf(buffer,"tmp1 = strmatch('%s',estim_params_.param_names,'exact');\nestim_params_.param_vals(tmp1,1) = %s;\n",id_var->name,value);
  str_output(buffer);
}
 
void p_estimated_bounds_elem1(struct token* id_var, char* value1, char* value2)
{
  char buffer[2000];
  if (id_var->endo_exo == 0)
    {
      /* missing tests for arguments */
      sprintf(buffer,"tmp1 = find(estim_params_.var_exo(:,1)==%d);\nestim_params_.var_exo(tmp1,3) = %s;\nestim_params_.var_exo(tmp1,4) = %s;\n",id_var->nbr+1,value1,value2);
    }
  else if(id_var->endo_exo == 1)
    {
      /* missing tests for arguments */
      sprintf(buffer,"tmp1 = find(estim_params_.var_endo(:,1)==%d);\nestim_params_.var_endo(tmp1,3) = %s;\nestim_params_.var_endo(tmp1,4) = %s;\n",id_var->nbr+1,value1,value2);
    }
  str_output(buffer);
}
 
void p_estimated_bounds_elem2(struct token* id_var1, struct token* id_var2, char* value1, char* value2)
{
  char buffer[2000];
  /* missing tests for arguments */
  if (id_var1->endo_exo == 0 && id_var2->endo_exo == 0)
    {
      sprintf(buffer,"tmp1 = find(estim_params_.corrx(:,1)==%d & estim_params.corrx(:,2)==%d);\nestim_params_.corrx(tmp1,4) = %s;\nestim_params_.corrx(tmp1,5) = %s;\n",id_var1->nbr+1,id_var2->nbr+1,value1,value2);
    }
  else if (id_var1->endo_exo == 1 && id_var2->endo_exo == 1)
    {
      sprintf(buffer,"tmp1 = find(estim_params_.corrn(:,1)==%d & estim_params.corrn(:,2)==%d);\nestim_params_.corrn(tmp1,4) = %s;\nestim_params_.corrn(tmp1,5) = %s;\n",id_var1->nbr+1,id_var2->nbr+1,value1,value2);
    }
  else
    {
      printf("Shocks and measurement errors can't be correlated\n");
      exit(1);
    }
  str_output(buffer);
}
 
void p_estimated_bounds_elem3(struct token* id_var, char* value1, char* value2)
{
  char buffer[2000];
  /* missing tests for arguments */
  sprintf(buffer,"tmp1 = strmatch('%s',estim_params_.param_names,'exact');\nestim_params_.param_vals(tmp1,2) = %s;\nestim_params_.param_vals(tmp1,3) = %s;\n",id_var->name,value1,value2);
  str_output(buffer);
}
 
void p_trend_init()
{
  str_output("trend_coeff_ = {};\n");
}

void p_trend_element(struct token* id_var, struct queue* expression)
{
  char buffer[20000];
  struct token** p_t;
  int j;

  sprintf(buffer,"tmp1 = strmatch('%s',options_.varobs,'exact');\n",id_var->var_ptr->name);
  str_output(buffer);
  sprintf(buffer,"trend_coeff_{tmp1} = sprintf('");
  p_t=(struct token **)expression->ptr;
  for (j=0;j < expression->imax;j++)
    {
      strcat(buffer,(*p_t)->name);
      ++p_t;
    }
  strcat(buffer,"');");
  str_output(buffer);
}

void p_estimation(void)
{
  char buffer[2000];
  int i;
#ifdef SCILAB
  str_output("global ");
  for(i=0;i < var_nbr; i++)
    {
      if (var_list[i].endo_exo == 1)
	{
	  str_output(var_list[i].name);
	  str_output(" ");
	}
    }
  str_output("\n");
#endif
#if defined MATLAB || defined SCILAB
  if (nbr_tmpvar == 0)
    {
      str_output("var_list_ = [];\n");
    }
  else
    {
      strcpy(buffer,"var_list_ = '");
      strcat(buffer,tmpvar_list[0]);
      strcat(buffer,"';\n");
      for(i=2;i<nbr_tmpvar;i+=2)
	{
#ifdef MATLAB
	  strcat(buffer,"var_list_ = str2mat(var_list_,'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"');\n");
#else
	  strcat(buffer,"var_list_ = [var_list_;'");
	  strcat(buffer,tmpvar_list[i]);
	  strcat(buffer,"'];\n");
#endif
	}
      str_output(buffer);
    }
  str_output("dynare_estimation(var_list_);\n");
#endif
}

void print_varobs(void)
{
  int i;
  char buffer[2000];
  str_output("options_.varobs = [];\n");
  for(i=0; i < nbr_tmpvar; i +=2)
    {
      sprintf(buffer,"options_.varobs = strvcat(options_.varobs,'%s');\n",tmpvar_list[i]);
      str_output(buffer);
    }
}

void print_unit_root_vars(void)
{
  int i;
  char buffer[2000];
  str_output("options_.unit_root_vars = {");
  for(i=0; i < nbr_tmpvar; i +=2)
    {
      if (i != 0)
	{
	  str_output("; ");
	}
      sprintf(buffer,"'%s'",tmpvar_list[i]);
      str_output(buffer);
    }
  str_output("};\n");
}

void p_optim_options(char* str1, char* str2, int task)
{
  static char optim_string[2000];
  static int start;
  switch(task){
  case 1:
    strcpy(optim_string,"options_.optim_opt = '");
    start = 0;
    return;
  case 2:
    if (start > 0)
      {
	strcat(optim_string,",");
      }
    else
      {
	start = 1;
      }
    strcat(optim_string,"''");
    strcat(optim_string,str1);
    strcat(optim_string,"'',");
    if (str2[0] >= 'A' && str2[0] <= 'z')
      {
	strcat(optim_string,"''");
	strcat(optim_string,str2);
	strcat(optim_string,"''");
      }
    else
      {
	strcat(optim_string,str2);
      }
    return;
  case 3:
    strcat(optim_string,"';\n");
    str_output(optim_string);
  }
}

char* my_strcat(char* s1, char* s2)
{
  char* result;
  result = calloc(strlen(s1)+strlen(s2),1);
  strcpy(result,s1);
  strcat(result,s2);
  return result;
}

void set_options(int argc, char** argv,struct s_runtime_options* runtime_options)
{
  int i;
  
  /* default values */

  runtime_options->clearall = 1;
  runtime_options->debug = 0;


  /* runtime arguments */
  for (i=2; i < argc; ++i)
    {
      if (strcasecmp(argv[i],"NOCLEARALL") == 0)
	{
	  runtime_options->clearall = 0;
	}
      if (strcasecmp(argv[i],"DEBUG") == 0)
	{
	  runtime_options->debug = 1;
	}
    }
}

void p_model_comparison(int flag_model_prior)
{
  int i;
  char buffer[2000];
  str_output("ModelNames_ = {");
  for(i=0; i < nbr_tmpvar; i +=2)
    {
      if (i != 0)
	{
	  str_output("; ");
	}
      sprintf(buffer,"'%s'",tmpvar_list[i]);
      str_output(buffer);
    }
  str_output("};\n");
  str_output("ModelPriors_ = [");
  if (flag_model_prior)
    {
      for(i=1; i < nbr_tmpvar; i +=2)
	{
	  if (i != 0)
	    {
	      str_output("; ");
	  }
	  sprintf(buffer,"%s",tmpvar_list[i]);
	  str_output(buffer);
	}
    }
  str_output("];\n");
  str_output("model_comparison(ModelNames_,ModelPriors_);\n");
}

/*
02/22/01 MJ added test for nbr of equations != nbr endogenous variables
            replaced stderr_ by Sigma_e_
02/23/01 MJ added global statement for all endogenous variables before 
            stoch_simul in MATLAB version
	    added p_dsample
02/28/01 MJ added "clear all" on top of MATLAB file
            corrected expression for Sigma_e_ in p_shocks
09/25/01 MJ added Scilab code
09/26/01 MJ added p_keyword
04/05/02 MJ corrected all compile warnings except m_del and \0 in format
04/06/02 MJ added p_optim_weights() p_optim_weights_init() p_expression()
            p_osr_params() p_osr()
            LEFT_PAREN RIGHT_PAREN OPEN_COMMENTS CLOSE_COMMENTS
10/09/02 MJ added p_calib_init() p_calib_var() p_calib()
            removed ";\n" from p_expression, corrected calling functions
10/20/02 MJ added automatic diary main() dynare_init()
01/16/03 MJ corrected bug in dyn2vec without argument
04/02/03 MJ added p_steady_linear()
            now ys_,exe_,exe_stoch_ are set in print_model() instead of 
	    p_initval()
04/28/03 MJ added p_option() and changed handling of options
05/02/03 MJ added p_olr(), p_olr_inst(), p_olr_inst_init()
            modified p_optim_weights1() and added p_optim_weights2()
05/18/03 MJ removed p_keyword(), modified p_steady() and p_steady_linear() 
            added p_check(); 
05/26/03 MJ added struct check and run_checks()
05/29/03 MJ p_shocks write set_shocks(flag,periods,ivar,values)
06/23/03 MJ added estimated_params and estimation
07/25/03 MJ added test for no exogenous for GAUSS version
10/31/03 MJ added oo for output variables
*/
