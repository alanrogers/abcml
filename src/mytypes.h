/** mytypes.h: Type definitions **/

#include <float.h>
#define FTOL 1.0e-7
#define FINITE(x) ((x) > -HUGE_VAL && (x) < HUGE_VAL)
#define MISSING (0.0/0.0)
#if 0                       /* single precision */
#define PRECISION 1
typedef float real;  
#define RFMT "%f"
#define MACHEPS  FLT_EPSILON
#else                       /* double precision */
#define PRECISION 2
typedef double real; 
#define RFMT "%lf"
#define MACHEPS  DBL_EPSILON
#endif
#define START_COMMENT '#'
/* properties of faunal elements */
struct bonedef {
  int npart;
  char **label;            /* labels for faunal elements */
  int *live;               /* counts per living animal */
  real *density;
  real *sensitivity;       /* proportional to 1/density */
};

/* weights of faunal elements */
struct wgt {
  int npart;
  char **label;            /* labels for faunal elements */
  real *value;             /* value of this part */
  real *weight;            /* weight of this part */
};

/* define an agent of deposition */
struct agent {
  int npart, nconfig;
  char **label;            /* labels for faunal elements */
  int **cfg;               /* npart X nconfig matrix of configurations */
  real *pr;                /* nconfig-vector of probabilities */
  real *m;                 /* mean vector */
  real **F;                /* mean squares and cross products */
};

/* hold archeological data */
struct counts {
  int npart;        /* number of faunal elements */
  int ndataset;     /* number of data sets */
  char **label;     /* labels for faunal elements */
  int **y;          /* y[i] is a vector of npart ints */
};

/* results */
struct results {
  real *p;     /* parameter vector */
  int mni;
  real lnL;
  real ChiSq;
  real **SCov;
  int SCov_err;
};

