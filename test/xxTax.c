#include <stdio.h>
#include <assert.h>
#include "mytypes.h"
#include "bone.h"
#include "alloc2d.h"
#include "misc.h"
#include "io.h"

#if 0
int npart = 3, npart_full, nagent, **lump;
real **m, ***F, *sensitivity;
int *y;
#endif

int attrition=1;
int npart=3;

void main(void)
{
  int **x, i, j;
  real **a, **b;

  x = lumpmat(0,1,npart);
  a = (real **) alloc2d(npart, npart, sizeof(real));
  b = (real **) alloc2d(npart-1, npart-1, sizeof(real));
  if(x==NULL || a==NULL || b==NULL)
    error("memory");

  for(i=0; i<npart; i++)
  {
    for(j=0; j<npart; j++)
      a[i][j] = 1.0/(1.0 + i+j);
    a[i][i] += 1.0;
  }

  prfmat(a, npart, npart, "a");
  primat(x, npart, npart-1, "x");
  xTax(x, a, b, npart, npart-1);
  prfmat(b, npart-1, npart-1, "x' a x");
  printf("\nb is%s symmetric",
   (symmetric(b, npart-1) ? "" : " not"));
  exit(0);
}
