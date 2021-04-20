/** misc.h: Prototypes for functions in misc.c **/
#include <stdlib.h>

#define NULL_ARG 0
#define FLAG_ARG 1
#define UNKNOWN_ARG 2
#define CFG_ARG 3
#define BDF_ARG 4
#define CNT_ARG 5
#define MAU_ARG 6
#define WGT_ARG 7

real quantile(real p, real *vec, int size);
int compar(real *x, real *y);
void error(const char *s);
void *mustalloc(size_t bytes);
void set_mnemonics(real **alpha, real **beta, real **kappa, real **p);
char *strlwr(char *u);
int classify(char *arg);
void **alloc2d(int rows, int cols, int elsize);
void free2d(void **m, int rows);








