/****************************************************************
    likelihood: Calculate the likelihood of an assemblage
    Copyright (C) 1999 Alan R. Rogers

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

    Alan R. Rogers, Department of Anthropology, University of Utah,
    Salt Lake City, UT 84112. rogers@anthro.utah.edu
****************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>
#include "mytypes.h"
#include "likelihood.h"
#include "misc.h"
#include "chol.h"
#include "io.h"
#include "symeig.h"

int     icmp(int *i, int *j);
#ifndef NDEBUG
int     symmetric(real ** c, int dim);
#endif

real  **reducmat = NULL;	/* principal components matrix */
real   *eigval = NULL;		/* vector of eigenvalues */
int     dim_reducmat = 0;	/* dimension of reducmat as allocated */
int     singularity = 0;

extern int npart, rdim, nagent;
extern int attrition;
extern real **m, ***F, *sensitivity;
extern int *y;
int     constrain_lnL = 0;	/* set = 1 to make lnL constrained */
real    lnL_ChiSq;
real  **lnL_cov = NULL;

#if 0
void    get_Ey(real * Ey, real kappa, real beta, real * alpha, real lastalpha)
{
    int     i, j;
    extern struct agent **agent;
    extern struct bonedef *bones;

    /* get vector E{y} */
    for (i = 0; i < npart; i++) {
	Ey[i] = 0.0;
	for (j = 0; j < nagent - 1; j++)
	    Ey[i] += alpha[j] * agent[j]->m[i];
	Ey[i] += lastalpha * agent[nagent - 1]->m[i];
	Ey[i] *= kappa;
	if (attrition)
	    Ey[i] *= exp(-beta * bones->sensitivity[i]);
    }
}
#endif

/* calculate covariance matrix */
void    get_cov(real ** cov, real kappa, real beta, real * alpha, real
		lastalpha)
{
    int     i, j, k;
    extern struct agent **agent;
    extern struct bonedef *bones;

#if 0    
    printf("\nget_cov: kappa=%g beta=%g alpha=",
	   kappa, (attrition ? beta : 0.0));
    for (i = 0; i < nagent - 1; i++)
	printf(" %g", alpha[i]);
    printf(" %g", lastalpha);
#endif    


    /* get covariance matrix */
    for (i = 0; i < npart; i++) {
	cov[i][i] = 0.0;
	for (k = 0; k < nagent - 1; k++)
	    cov[i][i] += alpha[k] * agent[k]->F[i][i];
	cov[i][i] += lastalpha * agent[nagent - 1]->F[i][i];
	cov[i][i] *= kappa;
	if (attrition)
	    cov[i][i] *= exp(-2.0 * beta * bones->sensitivity[i]);

	for (j = 0; j < i; j++) {
	    cov[i][j] = 0.0;
	    for (k = 0; k < nagent - 1; k++)
		cov[i][j] += alpha[k] * agent[k]->F[i][j];
	    cov[i][j] += lastalpha * agent[nagent - 1]->F[i][j];
	    cov[i][j] *= kappa;
	    if (attrition)
		cov[i][j] *= exp(-beta * (bones->sensitivity[i]
					  + bones->sensitivity[j]));
	    cov[j][i] = cov[i][j];
	}
    }
#ifndef NDEBUG
    if (!symmetric(lnL_cov, npart)) {
	printf("\nget_cov: matrix cov is not symmetric");
	prfmat(cov, npart, npart, "cov");
	exit(1);
    }
#endif
}

/****************************************************************
log likelihood is multivariate normal.

We want the maximum of lnL subject to the constraints:

  kappa >= 0
  beta   >= 0
  alpha[i] >= 0
  sum(alpha) = 1

When the input violates these constraints, lnL returns -HUGE_VAL.  Its
maximum should thus be a legal value.

Parameter vector is:
  [kappa, beta, alpha0, alpha1, ...]    if attrition>0
or
  [kappa, alpha0, alpha1, ...]          if attrition==0

If there are nagent agents of deposition, then the parameter value
will have nagent-1 alpha values.  The nagent'th value is the
complement of the others.  Thus the length of p is  nagent+1 (if
attrition) and nagent (if !attrition).
****************************************************************/
#undef PCMAHAL
real    lnL(real * p)
{
    int     i, j;
    real    lndet, rval, u;
    real   *alpha, *beta, *kappa;
    real    lastalpha;
    static real *eigenval = NULL, *work = NULL;
    extern int attrition;
    extern struct agent **agent;
    extern struct bonedef *bones;
    static real **c2 = NULL;
    static real *x = NULL;
    static real *x2 = NULL;
    static int x2_size = 0;

#if 0    
    printf("\nEnter lnL(");
    for(i=0; i < (attrition ? nagent+1 : nagent); ++i)
	printf(" %g", p[i]);
    printf(")");
#endif    

    /* things that get allocated just once */
    if (lnL_cov == NULL) {
	x = (real *) mustalloc(npart * sizeof(real));
	eigenval = (real *) mustalloc(npart * sizeof(real));
	work = (real *) mustalloc(npart * sizeof(real));
	lnL_cov = (real **) alloc2d(npart, npart, sizeof(real));
	if (lnL_cov == NULL)
	    error("lnL: alloc2d");
    }
    /* things that get reallocated when rdim grows */
    if (x2_size < rdim) {
	x2_size = rdim;
	if (x2 != NULL)
	    free(x2);
	if (c2 != NULL)
	    free(c2);
	x2 = (real *) mustalloc(rdim * sizeof(real));
	c2 = (real **) alloc2d(rdim, rdim, sizeof(real));
	if (c2 == NULL)
	    error("lnL: alloc2d");
    }
    set_mnemonics(&alpha, &beta, &kappa, &p);

    if (constrain_lnL) {
	/* check inequality constraints */
	lastalpha = 1.0;
	for (i = 0; i < nagent - 1; i++) {	/* alpha[i] >= 0 */
	    if (alpha[i] < 0.0) {
#if 0
		printf("\nlnL returning %g", -HUGE_VAL);
#endif
		return (-HUGE_VAL);
	    }
	    lastalpha -= alpha[i];
	}
	if (*kappa < 0.0 || (attrition && *beta < 0.0)
	    || lastalpha < 0.0) {
#if 0
	    printf("\nlnL returning %g", -HUGE_VAL);
#endif
	    return (-HUGE_VAL);
	}
    } else {
	/* set lastalpha */
	lastalpha = 1.0;
	for (i = 0; i < nagent - 1; i++)
	    lastalpha -= alpha[i];
    }

    /* get vector x = y - E{y} */
    for (i = 0; i < npart; i++) {
	u = 0.0;
	for (j = 0; j < nagent - 1; j++)
	    u += alpha[j] * agent[j]->m[i];
	u += lastalpha * agent[nagent - 1]->m[i];
	u *= *kappa;
	if (attrition)
	    u *= exp(-*beta * bones->sensitivity[i]);
	x[i] = y[i] - u;
    }

#if 0
    privec(y, npart, "y");
    prfvec(x, npart, "y - mu");
#endif

    /* get covariance matrix */
    if (beta != NULL)
	get_cov(lnL_cov, *kappa, *beta, alpha, lastalpha);
    else
	get_cov(lnL_cov, *kappa, 0.0, alpha, lastalpha);

#if 0
    prfmat(lnL_cov, npart, npart, "lnL_cov");
#endif

    /* reduce dimension */
    mat_times_vec(reducmat, x, x2, rdim, npart);
    xaxT(reducmat, lnL_cov, c2, rdim, npart);

#if 0
    prfmat(c2, rdim, rdim, "reduced C");
#endif

    /* calculate Mahalanobis Distance */
    if (cholesky(c2, rdim) != 0) {	/* cholesky decomposition */
	singularity = 1;
#if 0
	printf("\nlnL: reduced C is not PD");
	printf("\nkappa=%g", *kappa);
	if (beta != NULL)
	    printf(" beta=%g", *beta);
	printf(" alpha = (");
	printf(" %g", alpha[0]);
	for (i = 1; i < nagent - 1; i++)
	    printf(" %g", alpha[i]);
	printf(" %g", lastalpha);
	printf(")");
	printf("\nreturning -Inf");
	putchar('\n');
#endif
	rval = 0.0 / 0.0;
    } else {
	singularity = 0;
	lndet = log_determinant(c2, rdim);	/* determinant */
	lnL_ChiSq = mahal(x2, c2, rdim);	/* mahalanobis distance */
#if 0
	printf("\nlndet=%g chiSq=%g", lndet, lnL_ChiSq);
#endif
	rval = -0.5 * (lndet + lnL_ChiSq);
    }

#if 0
    printf("\nlnL returning %g", rval);
#endif

    return (rval);
}

/* has a singularity been detected ? */
int     check_singularity(void)
{
    return (singularity);
}

/** set mnemonics the same way in every function **/
void    set_mnemonics(real ** alpha, real ** beta, real ** kappa, real ** p)
{
    extern int attrition;

    *kappa = *p;
    if (attrition) {
	*beta = *p + 1;
	if (nagent > 1)
	    *alpha = *p + 2;
	else
	    *alpha = NULL;
    } else {
	*beta = (real *) NULL;
	if (nagent > 1)
	    *alpha = *p + 1;
	else
	    *alpha = NULL;
    }
}

#ifndef NDEBUG
/* test a square matrix of reals for symmetry */
int     symmetric(real ** c, int dim)
{
    int     i, j;
    real    x, eps;

    if (sizeof(real) == sizeof(double))
	        eps = DBL_EPSILON;
    else
	eps = FLT_EPSILON;

    for (i = 0; i < dim; i++)
	for (j = 0; j < i; j++) {
	    if (fabs(c[i][j]) >= fabs(c[j][i]))
		x = fabs(c[i][j]);
	    else
		x = fabs(c[j][i]);
	    if (fabs(c[i][j] - c[j][i]) > x * eps) {
		printf("\nasymmetric: m[%d][%d]=%.20g m[%d][%d]=%.20g",
		       i, j, c[i][j], j, i, c[j][i]);
		printf("\n   diff = %g, eps=%g\n", c[i][j] - c[j][i], eps);
		return (0);	/* asymmetric */
	    }
	}
    return (1);			/* symmetric */
}
#endif

/* make a lump matrix */
/* transposed */
real  **lumpmat(int a, int b, int dim)
{
    int     row, col;
    real  **mat;

    mat = (real **) alloc2d(dim - 1, dim, sizeof(real));
    if (mat == NULL)
	return (mat);

    /* make a < b */
    if (b < a) {
	col = b;
	b = a;
	a = col;
    }
    assert(a < b);

    /* cols 0..(b-1) are an identity matrix */
    for (row = 0; row < dim - 1; row++) {
	for (col = 0; col < b; col++)
	    mat[row][col] = (row == col ? 1.0 : 0.0);
    }

    /* column b has 1 in row a */
    for (row = 0; row < dim - 1; row++)
	mat[row][b] = (row == a ? 1.0 : 0.0);

    /* columns b+1..(dim-1) are an offset identity matrix */
    for (row = 0; row < dim - 1; row++)
	for (col = b + 1; col < dim; col++)
	    mat[row][col] = (row == col - 1 ? 1.0 : 0.0);

    return (mat);
}

/* make an identity matrix of ints */
real  **id_mat(int dim)
{
    int     row, col;
    real  **mat;

    mat = (real **) alloc2d(dim, dim, sizeof(real));
    if (mat == NULL)
	return (mat);

    for (row = 0; row < dim; row++)
	for (col = 0; col < dim; col++)
	    mat[row][col] = (row == col ? 1.0 : 0.0);

    return (mat);
}

/**************************************************************** 
Form b = x' a x

x : matrix of ndim X kdim reals
a : symmetric matrix of ndim X ndim reals
b : a kdim X kdim real matrix into which the answer will be put

****************************************************************/
void    xTax(real ** x, real ** a, real ** b, int ndim, int kdim)
{
    int     i, j, k, l;

    assert(b != NULL);

    for (i = 0; i < kdim; i++) {
	b[i][i] = 0.0;
	for (k = 0; k < ndim; k++)
	    for (l = 0; l < ndim; l++)
		b[i][i] += (x[k][i] * x[l][i]) * a[k][l];
	for (j = 0; j <= i; j++) {
	    b[i][j] = 0.0;
	    for (k = 0; k < ndim; k++)
		for (l = 0; l < ndim; l++)
		    b[i][j] += (x[k][i] * x[l][j]) * a[k][l];
	    b[j][i] = b[i][j];
	}
    }
}

/**************************************************************** 
Form b = x a x'

x : matrix of kdim X ndim reals
a : symmetric matrix of ndim X ndim reals
b : a kdim X kdim real matrix into which the answer will be put

****************************************************************/
void    xaxT(real ** x, real ** a, real ** b, int kdim, int ndim)
{
    int     i, j, k, l;

    assert(b != NULL);

    for (i = 0; i < kdim; i++) {
	for (j = 0; j <= i; j++) {
	    b[i][j] = 0.0;
	    for (k = 0; k < ndim; k++)
		for (l = 0; l < ndim; l++)
		    b[i][j] += (x[i][k] * x[j][l]) * a[k][l];
	    if (i != j)
		b[j][i] = b[i][j];
	}
    }
}

/**************************************************************** 
  Form b = v' * x

x : matrix of ndim X kdim reals
v : vector of ndim reals
b : resulting vector of kdim reals  
****************************************************************/
void    vTx(real * v, real ** x, real * b, int ndim, int kdim)
{
    int     i, j;

    for (i = 0; i < kdim; i++) {
	b[i] = 0.0;
	for (j = 0; j < ndim; j++)
	    b[i] += v[j] * x[j][i];
    }
}

/**************************************************************** 
  Form b = x * v

x : matrix of kdim X ndim reals
v : vector of ndim reals
b : resulting vector of kdim reals  
****************************************************************/
void    mat_times_vec(real ** x, real * v, real * b, int kdim, int ndim)
{
    int     i, j;

    for (i = 0; i < kdim; i++) {
	b[i] = 0.0;
	for (j = 0; j < ndim; j++)
	    b[i] += v[j] * x[i][j];
    }
}

/* lump two skeletal parts */
/* transposed */
int     dolump_ij(int a, int b, int rd, int np)
{
    int     i, j, k;
    real  **lmat, **lump2;

    lmat = lumpmat(a, b, rd);
    if (reducmat == NULL) {
	reducmat = lmat;
	dim_reducmat = rd - 1;
	return (rd);
    } else {
	/* allocate new lump matrix */
	lump2 = (real **) alloc2d(rd - 1, np, sizeof(real));
	if (lump2 == NULL) {
	    fprintf(stderr, "\ndolump_ij:alloc2d\n");
	    exit(1);
	}
	/* matrix multiplication: lump2 = reducmat * lmat */
	for (i = 0; i < rd - 1; i++)
	    for (j = 0; j < np; j++) {
		lump2[i][j] = 0;
		for (k = 0; k < rd; k++)
		    lump2[i][j] += reducmat[k][j] * lmat[i][k];
	    }
	free_reducmat();
	free2d((void **) lmat, rd - 1);
	reducmat = lump2;
	dim_reducmat = rd - 1;
    }
    rd -= 1;
    return (rd);
}

/**************************************************************** 
  np = npart, the number of parts in the data
  Fmean = mean of F matrices
  Function returns the number of dimensions in the reduced data,
  rdim.
****************************************************************/
int     dolump(int np, real ** Fmean)
{
    real  **c;
    static real **c2 = NULL;
    int     j, k, l, rd;

    if (c2 == NULL) {
	c2 = (real **) alloc2d(np, np, sizeof(real));
	if (c2 == NULL)
	    error("dolump: alloc2d");
    }
    rd = npart;
    reducmat = id_mat(npart);
    dim_reducmat = npart;
    do {
	l = rd;
	if (rd == npart)
	    c = Fmean;
	else {
	    xaxT(reducmat, Fmean, c2, rd, npart);	/* quadratic form */
	    c = c2;
	}
	for (j = 0; j < rd && rd == l; j++)
	    for (k = 0; k < j && rd == l; k++) {
		if (fabs(c[j][k]) > sqrt(c[j][j] * c[k][k]) * 0.91)
		    rd = dolump_ij(j, k, rd, npart);
	    }
    } while (rd < l);
    return (rd);
}

#define MAXCOND 750		/* maximal condition number */
/* reduce dimension by principal components */
int     dopc(int np, real ** Fmean)
{
    int     i, j, rd;
    real   x;
    static real *work = NULL, **vptr;
    static int *ndx;
    static real **evect = NULL;

    printf("\nEnter dopc");

    /* these only need to be allocated once */
    if (eigval == NULL) {
	eigval = (real *) mustalloc(np * sizeof(real));
	work = (real *) mustalloc(np * sizeof(real));
	ndx = (int *) mustalloc(np * sizeof(int));
	vptr = (real **) mustalloc(np * sizeof(real *));
	evect = (real **) alloc2d(np, np, sizeof(real));
	if (evect == NULL)
	    error("dopc: no memory");
    }
    /* freed and reallocated for each data set */
    if (reducmat == NULL) {
	reducmat = (real **) alloc2d(np, np, sizeof(real));
	dim_reducmat = np;
	if (reducmat == NULL) {
	    fprintf(stderr, "\ndopc: No memory.\n");
	    exit(1);
	}
    }
    for (i = 0; i < np; i++) {
	memcpy(evect[i], Fmean[i], np * sizeof(real));
	ndx[i] = i;
    }

    /* get eigenvectors and eigvalues */
    (void) symeig(evect, eigval, work, np);
#if 0
    prfvec(eigval, np, "Unsorted eigenvalues");
#endif

    /* order ndx so that eigval[ndx[i]] is in descending order */
    qsort(ndx, (unsigned) np, sizeof(int),
	          (int (*)(const void *, const void *)) icmp);

    /* find reduced dimension (rd), copy eigenvectors into reducmat */
    for (rd = 0; rd < np; rd++) {
	if (fabs(eigval[ndx[0]]) > MAXCOND * fabs(eigval[ndx[rd]]))
	    break;
	x = sqrt(eigval[ndx[rd]]);  /* rescale eigenvectors */
	for(j=0; j<np; ++j)
	    evect[ndx[rd]][j] /= x;
	memcpy(reducmat[rd], evect[ndx[rd]], np * sizeof(real));
    }

#ifndef NDEBUG
    for (i = 1; i < np; i++) {
	if (fabs(eigval[ndx[i]]) > fabs(eigval[ndx[i - 1]])) {
	    fflush(stdout);
	    fprintf(stderr, "\nError: |eigval[%d]| < |eigval[%d]|",
		    i - 1, i);
	    fprintf(stderr, "\n       eigval[%d]=%g eigval[%d]=%g",
		    i - 1, eigval[ndx[i - 1]], i, eigval[ndx[i]]);
	    fprintf(stderr, "\nFLT_EPSILON=%g DBL_EPSILON=%g",
		    FLT_EPSILON, DBL_EPSILON);
	    putc('\n', stderr);
	    exit(1);
	}
    }
#endif

#if 0
    for (i = 0; i < rd; i++) {
	printf("\nEigenvalue %d=%g", i, eigval[ndx[i]]);
	/*    printf("\nEigenvector %d", i); */
	/*    prfvec(reducmat[i], np, ""); */
    }
    fflush(stdout);
#endif

    return (rd);		/* return reduced dimension */
}

int     icmp(int *i, int *j)
{
    if (fabs(eigval[*i]) > fabs(eigval[*j]))
	return (-1);
    if (fabs(eigval[*i]) < fabs(eigval[*j]))
	return (1);
    return (0);
}

/****************************************************************
Calculate the mean F matrix.  nalpha is the number of agents, but only the
1st nalpha-1 entries are used.  The last element is assumed to equal
the complement of the others.
****************************************************************/
void    getFmean(real ** Fmean, struct agent **agent, int dim, real * alpha,
		 int nalpha)
{
    int     i, j, k;
    real    lastalpha = 1.0;

    nalpha -= 1;		/* last entry of alpha is unused */
    for (i = 0; i < nalpha; i++) {
	assert(alpha[i] >= 0.0);
	assert(alpha[i] <= 1.0);
	lastalpha -= alpha[i];
    }
    assert(lastalpha >= 0.0);
    assert(lastalpha <= 1.0);

#if 0    
    printf("\ngetFmean: kappa=1 beta=0 alpha =");
    for (i = 0; i < nalpha; i++)
	printf(" %g", alpha[i]);
    printf(" %g", lastalpha);
#endif    

    for (j = 0; j < dim; j++)
	for (k = 0; k < dim; k++) {
	    Fmean[j][k] = 0.0;
	    for (i = 0; i < nalpha; i++)
		Fmean[j][k] += alpha[i] * agent[i]->F[j][k];
	    Fmean[j][k] += lastalpha * agent[nalpha]->F[j][k];
	}
}

/* transposed */
#define FMTSIZ 10
void    pr_reducmat(char **lbl, int rows, int cols, int dataset)
{
    int     maxlbl, i, j;
    char    fmt[FMTSIZ];

    maxlbl = 0;
    for (i = 0; i < npart; i++) {
	j = strlen(lbl[i]);
	if (j > maxlbl)
	    maxlbl = j;
    }
    if (maxlbl > 50) {
	fprintf(stderr, "\nWarning: labels truncated to 50 chars");
	maxlbl = 50;
    }
    sprintf(fmt, "\n%%-%ds", maxlbl);
    if (strlen(fmt) >= FMTSIZ)
	error("pr_reducmat: format buffer overflow");
    printf("\n\ntranspose of Reducmat matrix for dataset %d:", dataset);
    for (i = 0; i < rows; i++) {
	printf(fmt, lbl[i]);
	for (j = 0; j < cols; j++)
	    printf(" %g", reducmat[j][i]);
    }
}
#undef FMTSIZ

void    free_reducmat(void)
{
    if (reducmat != NULL)
	free2d((void **) reducmat, dim_reducmat);
    reducmat = NULL;
    dim_reducmat = 0;
}

/* make reducmat an identity matrix so that no reduction takes place */
void    noreduc(void)
{
    if (reducmat != NULL)
	free_reducmat();
    reducmat = id_mat(npart);
    dim_reducmat = npart;
}
