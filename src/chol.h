/****************************************************************
    chol.h: Prototypes for functions in chol.c
****************************************************************/
int   cholesky(real **A, int n);
void  matmultLLT(real **L, real **B, int n);
real  log_determinant(real **L, int n);
real  dotprod(real *x, real *y, int n);
int   col_L_solve(real **L, int n, real *b);
int   row_L_solve(real **L, int n, real *b);
real  mahal(real *x, real **L, int n);
void  cholinv(real **L, real **s, int n);
int   srt_diag(real **A, int *ndx, int n);
int   srt_diag_compar(int *i, int *j);
void  invert_permutation(int *permutation, int *inv, int n);
