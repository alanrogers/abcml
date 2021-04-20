/****************************************************************
    symeig: Eigenvalues and eigenvectors of a symmetric matrix

    This code was adapted from that in Numerical Recipes in C, 2nd
    edition (by WH Press, SA Teukolksy, and WT Vetterling., Cambridge
    University Press).
****************************************************************/
#include <math.h>
#include "mytypes.h"
#include "misc.h"
#include "symeig.h"
#define SIGN(a,b) ((b)<0 ? -fabs(a) : fabs(a))

/****************************************************************
symeig: Eigenvalues and eigenvectors of a symmetric matrix.  The job
is done in two steps:

  1. Function tred2 reduces the matrix to tridiagonal form using
  Householder transformations.  The original matrix is replaced with
  a transformation matrix.

  2. Function tqli gets eigenvalues and vectors of the tridiagonal
  matrix using the QL algorithm.  It also converts the tridiagonal
  eigenvectors into eigenvectors of the original problem.

On entry,
  m is a dim X dim matrix of reals containing a symmetric matrix.
  eigenval and work are arbitrary vectors of dim reals.

On return
  eigenval contains the eigenvalues.
  The ith row of m contains the ith eigenvector.
  work is destroyed.
  The function returns 0 on success, 1 on failure.
****************************************************************/
int symeig(real **m, real *eigenval, real *work, int dim)
{
  tred2(m, dim, eigenval, work);   /* householder reduction */
  tqli(eigenval, work, dim, m);    /* QL algorithm */
  return(0);
}

/****************************************************************
Householder reduction of a real, symmetric matrix a[1..n][1..n].  On
output, a is replaced by the orthogonal matrix Q effecting the
transformation.  d[1..n] returns the diagonal elements of the
tridiagonal matrix, and e[1..n] the off-diagonal elements, with
e[1]=0.  Several statements, as noted in comments, can be omitted if
only eigenvalues are to be found, in which case a contains no useful
information on output.  Otherwise they are to be included.

Changes:
* Arrays and matrices are [0..(n-1)] rather than [1..n].  Alan Rogers

* The matrix a is now the transpose of a in the Numerical Recipes
version.  Alan Rogers
****************************************************************/
void tred2(real **a, int n, real *d, real *e)
{
  int i, j, k, l; 
  real scale, hh, h, g, f;

  for (i=n-1; i>=1; i--)
  {
    l = i - 1;
    h = scale = 0.0;
    if (l > 0)
    {
      for (k = 0; k <= l; k++)
	scale += fabs(a[k][i]);

      if (scale == 0.0)
	e[i] = a[l][i];
      else
      {
	for (k = 0; k <= l; k++)
	{
	  a[k][i] /= scale;
	  h += a[k][i] * a[k][i];
	}
	f = a[l][i];
	g = f>0 ? -sqrt(h) : sqrt(h);
	e[i] = scale * g;
	h -= f * g;
	a[l][i] = f - g;
	f = 0.0;
	for (j = 0; j<=l; j++)
	{
	  /* Next statement can be omitted if eigenvectors not wanted */
	  a[i][j] = a[j][i] / h;
	  g = 0.0;
	  for (k = 0; k <= j; k++)
	    g += a[k][j] * a[k][i];
	  for (k = j+1; k<=l; k++)
	    g += a[j][k] * a[k][i];
	  e[j] = g / h;
	  f += e[j] * a[j][i];
	}
	hh = f / (h + h);
	for (j = 0; j<=l; j++)
	{
	  f = a[j][i];
	  e[j] = g = e[j] - hh * f;
	  for (k = 0; k<=j; k++)
	    a[k][j] -= (f * e[k] + g * a[k][i]);
	}
      }
    }else
      e[i] = a[l][i];
    d[i] = h;
  }
  /* Next statement can be omitted if eigenvectors not wanted */
  d[0] = 0.0;
  e[0] = 0.0;
  /* Contents of this loop can be omitted if eigenvectors not
     wanted except for statement d[i] = a[i][i]; */
  for (i = 0; i<n; i++)
  {
    l = i - 1;
    if (d[i])
    {
      for (j = 0; j<=l; j++)
      {
	g = 0.0;
	for (k = 0; k<=l; k++)
	  g += a[k][i] * a[j][k];
	for (k = 0; k<=l; k++)
	  a[j][k] -= g * a[i][k];
      }
    }
    d[i] = a[i][i];
    a[i][i] = 1.0;
    for(j = 0; j<=l; j++)
      a[i][j] = a[j][i] = 0.0;
  }
}

/****************************************************************
QL algorithm with implicit shifts, to determine the eigenvalues and
eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
symmetric matrix previously reduced by tred2.  On input, d[1..n]
contains the diagonal elements of the tridiagonal matrix.  On output,
it returns the eigenvalues.  The vector e[1..n] inputs the subdiagonal
elements of the tridiagonal matrix, with e[1] arbitrary.  On output, e
is destroyed.  When finding only the eigenvalues, several lines may be
omitted, as noted in the comments.  If the eigenvectors of a
tridiagonal matrix are desired, the matrix z[1..n][1..n] is input as
the identity matrix.  If the eigenvectors of a matrix that has been
reduced by tred2 are required, then z is input as the matrix output by
tred2.  In either case, the kth column of z returns the normalized
eigenvector corresponding to d[k].

From Numerical Recipes, 2nd edition.

Changes:
* Re-written to use base-0 arrays rather than base-1 arrays. Alan Rogers
* Re-written to make z the transpose of the z in Numerical Recipes.
This makes the eigenvectors come out in rows rather than in columns.
****************************************************************/

void tqli(real *d, real *e, int n, real **z)
{
  int iter;
  int i, k, l, m;       /* indices for base-0 arrays */
  real s,r,p,g,f,dd,c,b;

  /* it is convenient to renumber the elements */
  for (i=1; i<n ;i++)
    e[i-1] = e[i];
  e[n-1] = 0.0;
  
  for (l=0; l < n; l++)
  {
    iter = 0;
    do {
      /* test for convergence */
      for (m=l; m < n-1; m++)
      {
	dd = fabs(d[m])+fabs(d[m+1]);
	if (fabs(e[m])+dd == dd)
	  break;
      }
      /* if m==l then we have converged */
      if (m != l)
      {
	if (iter++ == 30)
	  error("Too many iterations in TQLI");
	g = (d[l+1]-d[l])/(2.0*e[l]);
	r = sqrt((g*g)+1.0);
	g = d[m]-d[l]+e[l]/(g+SIGN(r,g));
	s = c = 1.0;
	p = 0.0;
	for (i=m-1; i>=l; i--)
	{
	  f = s*e[i];
	  b = c*e[i];
	  if (fabs(f) >= fabs(g))
	  {
	    c = g/f;
	    r = sqrt((c*c)+1.0);
	    e[i+1] = f*r;
	    c *= (s = 1.0/r);
	  }else
	  {
	    s = f/g;
	    r = sqrt((s*s)+1.0);
	    e[i+1] = g*r;
	    s *= (c = 1.0/r);
	  }
	  g = d[i+1]-p;
	  r = (d[i]-g)*s+2.0*c*b;
	  p = s*r;
	  d[i+1] = g+p;
	  g = c*r-b;
	  /* Next loop can be omitted if eigenvectors not wanted */
	  for (k = 0; k < n; k++)
	  {
	    f = z[i+1][k];
	    z[i+1][k] = s * z[i][k] + c * f;
	    z[i][k] = c * z[i][k] - s * f;
	  }
	}
	d[l] -= p;
	e[l] = g;
	e[m] = 0.0;
      }
    } while (m != l);
  }
}
