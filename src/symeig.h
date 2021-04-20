/** symeig.h: Prototypes for functions in symeig.c **/

int symeig(real **m, real *eigenval, real *work, int dim);
void tred2(real **a, int n, real *d, real *e);
void tqli(real *d, real *e, int n, real **z);
