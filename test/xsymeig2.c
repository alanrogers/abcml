/* Driver for routine TQLI */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "symeig.h"
#include "misc.h"
#include "io.h"

#define DIM 3
#define TINY 1.0e-6

void main(void)
{
  int i,j,k;
  real *eigenval, *work, *f, **mat, **c;

  eigenval = (real *) mustalloc(DIM * sizeof(real));
  work     = (real *) mustalloc(DIM * sizeof(real));
  f        = (real *) mustalloc(DIM * sizeof(real));
  c        = (real **) alloc2d(DIM, DIM, sizeof(real));
  mat      = (real **) alloc2d(DIM, DIM, sizeof(real));
  if(c == NULL || mat == NULL)
    error("main: alloc2d");
  for (i=0; i<DIM; i++)
    for (j=0; j<DIM; j++)
      c[i][j] = mat[i][j] = 1.0 / ( 1.0 + i + j);

  /* get eigenvalues and eigenvectors of a real symmetric matrix */
  symeig(mat, eigenval, work, DIM);

  prfmat(c, DIM, DIM, "Original matrix");
  prfmat(mat, DIM, DIM, "Eigenvectors");
  putchar('\n');
  printf("Eigenvalues:\n");
  for(i=0; i<DIM; i++)
      printf(" %8.4f", eigenval[i]);
  putchar('\n');

  printf("\nEigenvectors for a real symmetric matrix\n");
  for (i=0; i<DIM; i++)
  {
    for (j=0; j<DIM; j++)
    {
      f[j]=0.0;
      for (k=0; k<DIM; k++)
	f[j] += (c[k][j]*mat[i][k]);
    }
    printf("Eigenvalue %d = %10.6f\n", i, eigenval[i]);
    printf("%12s %14s %9s\n","vector","matrix*vect","ratio");
    for (j=0; j<DIM; j++)
    {
      if (fabs(mat[i][j]) < TINY)
	printf("%12.6f %14.6f %12s\n", mat[i][j], f[j],"div. by 0");
      else
	printf("%12.6f %14.6f %12.6f\n", mat[i][j], f[j], f[j]/mat[i][j]);
    }
    printf("Press ENTER to continue...\n");
    getchar();
  }
  free2d((void **) mat, DIM);
  free(f);
  free(work);
  free(eigenval);
  exit(0);
}
