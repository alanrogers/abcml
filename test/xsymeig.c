/* Driver for routine TQLI */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "symeig.h"
#include "misc.h"

#define NP 10
#define TINY 1.0e-6
void error(char *msg);

void main(void)
{
  int i,j,k;
  real *eigenval, *work, *f, **mat;
  static real c[NP][NP]=
  { {5.0, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0,-2.0,-3.0,-4.0},
    {4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0,-2.0,-3.0},
    {3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0,-2.0},
    {2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0,-1.0},
    {1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0, 0.0},
    {0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0},
    {-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0},
    {-2.0,-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0},
    {-3.0,-2.0,-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 4.0},
    {-4.0,-3.0,-2.0,-1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0}};

  eigenval = (real *) mustalloc(NP * sizeof(real));
  work = (real *) mustalloc(NP * sizeof(real));
  f = (real *) mustalloc(NP * sizeof(real));
  mat = (real **) alloc2d(NP, NP, sizeof(real));
  if(mat == NULL)
    error("main: alloc2d");
  for (i=0; i<NP; i++)
    for (j=0; j<NP; j++)
      mat[i][j]=c[i][j];

  /* get eigenvalues and eigenvectors of a real symmetric matrix */
  symeig(mat, eigenval, work, NP);

  printf("\nEigenvectors for a real symmetric matrix\n");
  for (i=0; i<NP; i++)
  {
    for (j=0; j<NP; j++)
    {
      f[j]=0.0;
      for (k=0; k<NP; k++)
	f[j] += (c[k][j]*mat[i][k]);
    }
    printf("Eigenvalue %d = %10.6f\n", i, eigenval[i]);
    printf("%12s %14s %9s\n","vector","matrix*vect","ratio");
    for (j=0; j<NP; j++)
    {
      if (fabs(mat[i][j]) < TINY)
	printf("%12.6f %14.6f %12s\n", mat[i][j], f[j],"div. by 0");
      else
	printf("%12.6f %14.6f %12.6f\n", mat[i][j], f[j], f[j]/mat[i][j]);
    }
    printf("Press ENTER to continue...\n");
    getchar();
  }
  free2d((void **) mat, NP);
  free(f);
  free(work);
  free(eigenval);
  exit(0);
}
