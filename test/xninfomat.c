#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "mytypes.h"
#include "misc.h"
#include "io.h"
#include "chol.h"
#include "ninfomat.h"

int nparam = 2;

real f(real *p)
{
  real ans=0.0, x, y;

  x = p[0];
  y = p[1];
  ans = -2.0*(x-1.0)*(x-1.0)
        -3.0*(y-2.0)*(y-2.0)
        -(x-1.0)*(y-2.0);
  return(ans);
}

int main(int argc, char **argv)
{
  real **C, **Cinv, *p;
  int i;

  printf("\nBeginning of xninfomat");

  C = (real **) alloc2d(nparam, nparam, sizeof(real));
  Cinv = (real **) alloc2d(nparam, nparam, sizeof(real));
  if(C==NULL || Cinv==NULL)
  {
    fprintf(stderr,"\nNo memory\n");
    exit(1);
  }
  p = (real *) mustalloc(nparam * sizeof(real));

  for(i=0; i<nparam; i++)
    p[i] = i+1.0;

  ninfomat(C, p, nparam, f, 1);
  prfmat(C, nparam, nparam, "Matrix produced by ninfomat");
  cholesky(C, nparam);
  cholinv(C, Cinv, nparam);
  prfmat(Cinv, nparam, nparam, "Inverse");

  printf("\nEnd of xninfomat");
  putchar('\n');
  return 0;
}
