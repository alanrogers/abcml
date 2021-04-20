/** likelihood.h: Prototypes for functions in amoeba.c **/

/****************************************************************
parameters: p[0] = kappa
            p[1] = alpha
****************************************************************/
/*****************************************************
     |p[0]| = kappa
      p[1]  = beta
      p[2]  = alpha[0]
           ....
  p[nagent] = alpha[nagent-2]
******************************************************/
#define ATTRITION 0     /* turn attrition on or off */
#if ATTRITION
#define NPARAM 3
#define BETA   p[1]
#else
#define NPARAM 2
#define BETA   0
#endif
#define KAPPA fabs(p[0])
/****************************************************************
randint(n) returns a random integer between 0 and n-1.  It is now
implemented as a macro, which should improve speed.  The float-int
conversion rounds down, thus giving an int uniformly distributed on
0,1,...,(n-1).  We will never get n because uni() is uniformly
distributed on [0,1), not [0,1].
****************************************************************/
#define randint(n)  ( (int) (uni() * (n)))

void get_Ey(real *Ey, real kappa, real beta, real *alpha, real lastalpha);
void get_cov(real **cov, real kappa, real beta, real *alpha, real
	     lastalpha);
real lnL(real *p);
int check_singularity(void);
int symmetric(real **c, int dim);
real **lumpmat(int a, int b, int dim);
real **id_mat(int dim);
void xTax(real **x, real **a, real **b, int ndim, int kdim);
void vTx(real *v, real **x, real *b, int ndim, int kdim);
void xaxT(real **x, real **a, real **b, int kdim, int ndim);
void mat_times_vec(real **x, real *v, real *b, int kdim, int ndim);
int dolump_ij(int a, int b, int rd, int np);
int dolump(int np, real **Fmean);
int dopc(int np, real **Fmean);
void getFmean(real **Fmean, struct agent **agent, int dim, real *alpha,
	      int nalpha);
void pr_reducmat(char **lbl, int rows, int cols, int dataset);
void free_reducmat(void);
void noreduc(void);
