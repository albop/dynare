#define MAX_INDEX_NBR 10
struct var {
  char *name;
  int endo_exo;   /* 0: exogenous; 1: endogenous; 2: loop index (numerical)*/
  int nbr;
  int original_nbr;
};

struct token {
  int endo_exo;
  int nbr;        /* >=0: variable, -1: string, -2: loop index, -3 indexed variable */
  int lead_lag;
  char *name;     /* string content or variable lead_lag name */
  struct var *var_ptr; /* only for variables */
};

struct queue {
  void **ptr;
  int flag;
  int imax;
};

struct loop {
  int loop_index;
  int loop_limit;
  int loop_increment;
  char *loop_index_name;
};

struct s_check{
  int eq_nbr;
  int determ;
  int stoch;
  int initval;
  int endval;
  int histval;
  int steady;
  int linear;
  int stoch_simul;
  int simul;
  int check;
  int olr;
  int osr;
};

struct s_runtime_options{
  int clearall;
  int debug;
};

struct s_estim_params{
  int param_type;
  int var_type;
  int var_nbr;
  int var_nbr2;
  char* param_name;
  char* init_val;
  char* lb;
  char* ub;
  char* prior;
  char* mean;
  char* std;
  char* p3;
  char* p4;
  char* jscale;
};
  
void add_var(char *, int), print_var(), print_model(struct queue *);
void set_ll(struct queue *, char *, int), mark_pound(struct queue *);
struct queue * create_queue(void *);
struct queue * add_to_queue(struct queue *, void *);
struct queue * copy_queue(struct queue *, struct queue *);
struct token * create_var(struct var *); 
struct token * token(char *, int nbr);
struct queue * m_del(char *, struct queue *);
struct var *var_search(const char *);
struct loop *initial_loop(char *, char *, char *, char *);
struct queue * do_loop(struct loop *, struct queue *);
struct queue *operator_loop(struct loop *, struct queue *, char *);
void var_output(struct token *, int);
int * periods(char *, char *);
void p_shocks(struct token *, struct queue *, struct queue *, int);
void p_init(struct token *, struct queue *);
void p_initval(void),p_endval(void),p_i_shocks(void);
void dynare_init(char *, struct s_runtime_options);
void print_iter(char *);
void pe_initval(void);
void pe_endval(void);
int nbr_tmpvar;
void print_rplot(void);
void copy_update_indexed_var(struct token **, struct queue  *, int , char *, int);
void print_model(struct queue *);
void print_model1(struct queue *, int);
void print_param(void);
void p_optim_weights_init(void);
void p_optim_weights(struct token *, struct queue *);
void p_expression(struct queue *);
void p_osr_params(char *);
void p_osr(void);
char* my_strcat(char*, char*);
void p_estimated_elem(void);
void estimated_elem_init(void);

/*
04/06/02 MJ added p_optim_weights_init and p_optim_weights p_expression p_osr_params p_osr
*/
