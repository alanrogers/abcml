/****************************************************************
    io: Input-out routines
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
#include <string.h>
#include "mytypes.h"
#include "io.h"
#include "getcic.h"
#include "misc.h"
#include "zeroin.h"

/* externals */
static int dim;
static int *live;
static real *density;

/* allocate and initialize a bonedef structure */
struct bonedef *alloc_bonedef(void)
{
  struct bonedef *a;

  a = (struct bonedef *) malloc(sizeof(struct bonedef));
  if (a == NULL)
    return (NULL);
  a->npart = 0;
  a->label = NULL;
  a->live = NULL;
  a->density = NULL;
  a->sensitivity = NULL;
  return (a);
}

/* allocate and initialize a weight structure */
struct wgt *alloc_wgt(void)
{
  struct wgt *a;

  a = (struct wgt *) malloc(sizeof(struct wgt));
  if (a == NULL)
    return (NULL);
  a->npart = 0;
  a->label = NULL;
  a->value = NULL;
  a->weight = NULL;
  return (a);
}

/* allocate and initialize an agent structure */
struct agent *alloc_agent(void)
{
  struct agent *a;

  a = (struct agent *) malloc(sizeof(struct agent));
  if (a == NULL)
    return (NULL);
  a->npart = 0;
  a->label = NULL;
  a->pr = NULL;
  a->m = NULL;
  a->F = NULL;
  return (a);
}

/* allocate and initialize a structure of counts */
struct counts *alloc_counts(void)
{
  struct counts *c;

  c = (struct counts *) malloc(sizeof(struct counts));
  if (c == NULL)
    return (NULL);
  c->npart = 0;
  c->ndataset = 0;
  c->label = NULL;
  c->y = NULL;
  return (c);
}

/* read a .bdf file, create a bonedef structure */
#define BUFFSIZ 100
struct bonedef *read_bonedef(char *fname, real * scale_factor_ptr)
{
  struct bonedef *b;
  real    scale_factor;
  char    buff[BUFFSIZ];
  int     i;
  FILE   *fp;

  fp = mustopen(fname, "r");

  b = alloc_bonedef();

  /* read number of faunal elements */
  getintic(&(b->npart), buff, BUFFSIZ, fp);
  fflush(stdout);

  /* allocations */
  b->label = (char **) mustalloc(b->npart * sizeof(char *));
  b->live = (int *) mustalloc(b->npart * sizeof(int));
  b->density = (real *) mustalloc(b->npart * sizeof(real));
  b->sensitivity = (real *) mustalloc(b->npart * sizeof(real));

  for (i = 0; i < b->npart; i++) {
    if ((b->label[i] = getwordic(buff, BUFFSIZ, fp)) == NULL)
      error("read_bonedef: unexpected EOF");
    b->label[i] = copystring(buff, BUFFSIZ);
    if (getintic(b->live + i, buff, BUFFSIZ, fp) == EOF)
      error("read_bonedef: unexpected EOF");
    if (getrealic(b->density + i, buff, BUFFSIZ, fp) == EOF)
      error("read_bonedef: unexpected EOF");
  }
  if (getwordic(buff, BUFFSIZ, fp) != NULL) {
    fprintf(stderr, "\nError: read_bonedef: Did not exhaust input file %s",
	    fname);
    fprintf(stderr, "\n  Extra word: %s\n", buff);
#if 0
    exit(1);
#endif
  }
  /* set sensitivities */
  dim = b->npart;		/* set externals used by sfun */
  live = b->live;
  density = b->density;
  /****************************************************************
  If *scale_factor_ptr is > 0.0 on input, then use the value as given.  
  Otherwise set scale_factor so that when beta=1, half the bones in a whole
  animal survive.
  ****************************************************************/
  if (scale_factor_ptr == NULL)
    scale_factor = -1.0;
  else
    scale_factor = *scale_factor_ptr;
  if (scale_factor <= 0.0)
    scale_factor = zeroin(0.0, 10.0, sfun, 0.0001);
  if (scale_factor_ptr != NULL)
    *scale_factor_ptr = scale_factor;

  for (i = 0; i < b->npart; i++)
    b->sensitivity[i] = scale_factor / b->density[i];

  fclose(fp);
  return (b);
}

/* read a .wgt file, create a wgt structure */
struct wgt *read_wgt(char *fname)
{
  struct wgt *b;
  char    buff[BUFFSIZ];
  int     i;
  FILE   *fp;

  fp = mustopen(fname, "r");

  b = alloc_wgt();

  /* read number of faunal elements */
  getintic(&(b->npart), buff, BUFFSIZ, fp);
  fflush(stdout);

  /* allocations */
  b->label = (char **) mustalloc(b->npart * sizeof(char *));
  b->value = (real *) mustalloc(b->npart * sizeof(real));
  b->weight = (real *) mustalloc(b->npart * sizeof(real));

  for (i = 0; i < b->npart; i++) {
    if ((b->label[i] = getwordic(buff, BUFFSIZ, fp)) == NULL)
      error("read_wgt: unexpected EOF");
    b->label[i] = copystring(buff, BUFFSIZ);
    if (getrealic(b->value + i, buff, BUFFSIZ, fp) == EOF)
      error("read_wgt: unexpected EOF");
    if (getrealic(b->weight + i, buff, BUFFSIZ, fp) == EOF)
      error("read_wgt: unexpected EOF");
  }
  if (getwordic(buff, BUFFSIZ, fp) != NULL) {
    fprintf(stderr, "\nError: read_wgt: Did not exhaust input file %s",
	    fname);
    fprintf(stderr, "\n  Extra word: %s\n", buff);
  }
  fclose(fp);
  return (b);
}

/* print a bonedef structure */
void    pr_bonedef(struct bonedef *b)
{
  int     i;

  printf("\n%-17s %4s %7s %11s", "part", "live", "density", "sensitivity");
  for (i = 0; i < b->npart; i++)
    printf("\n%-17s %4d %7.3f %11.4f", b->label[i], b->live[i], b->density[i],
	   b->sensitivity[i]);
  putchar('\n');
}

/* print an agent structure */
void    pr_agent(struct agent *a)
{
  int     i, j;
  real    minp;

  printf("\n%6d  # number of parts", a->npart);
  printf("\n%6d  # number of configurations", a->nconfig);

  if (a->pr != NULL) {
    /* find minimum val of a->pr[i] */
    for (i = 0; i < a->nconfig && a->pr[i] == 0.0; i++) ;
    if (i < a->nconfig && a->pr[i] > 0.0)
      minp = a->pr[i];
    else
      minp = 1.0;
    while (i < a->nconfig) {
      if (a->pr[i] > 0.0 && minp > a->pr[i])
	minp = a->pr[i];
      i += 1;
    }

    /* print probabilities rel to minp */
    printf("\n#Probabilities of configurations are proportional to:");
    printf("\n%17s ", "");
    for (i = 0; i < a->nconfig; i++)
      printf(" %g", a->pr[i] / minp);
    printf("\n#label");
    for (i = 0; i < a->npart; i++) {
      printf("\n%-17s", a->label[i]);
      for (j = 0; j < a->nconfig; j++)
	printf(" %2d", a->cfg[i][j]);
    }
  }
  if (a->m != NULL) {
    printf("\n\nMean:");
    for (i = 0; i < a->npart; i++) {
      printf("\n%-17s", a->label[i]);
      printf(" %7.3f", a->m[i]);
    }
  }
  if (a->F != NULL)
    prfmat(a->F, a->npart, a->npart, "\nF");

  putchar('\n');
}

/* used by zeroin to find scale factor for sensitivities */
double  sfun(double x)
{
  int     i;
  double  s = 0.0, t = 0.0;

  for (i = 0; i < dim; i++) {
    s += exp(-x / density[i]) * live[i];
    t += live[i];
  }
  return (0.5 - s / t);
}

char   *copystring(char *a, int size)
{
  char   *b;
  int     len;

  for (len = 0; len < size && a[len] != '\0'; len++) ;
  b = (char *) malloc((len + 1) * sizeof(char));
  if (b == NULL) {
    fprintf(stderr, "\nERROR: no memory in copystring\n");
    exit(1);
  }
  memcpy(b, a, (unsigned) len);	/* copy a into b */
  b[len] = '\0';
  return (b);
}

/* read agent configurations */
struct agent *read_agent_cfg(char *fname)
{
  struct agent *a;
  char    buff[BUFFSIZ];
  int     i, j;
  real    sum;
  FILE   *fp;

  fp = mustopen(fname, "r");

  /* read number of faunal elements */
  a = alloc_agent();
  getintic(&(a->npart), buff, BUFFSIZ, fp);	/* number of faunal elements */
  getintic(&(a->nconfig), buff, BUFFSIZ, fp);	/* number of configurations */
  fflush(stdout);

  /* allocations */
  a->label = (char **) mustalloc(a->npart * sizeof(char *));
  a->pr = (real *) mustalloc(a->nconfig * sizeof(real));
  a->cfg = (int **) alloc2d(a->npart, a->nconfig, sizeof(int));
  if (a->cfg == NULL)
    error("read_agent_cfg: alloc2d");

  /* read and normalize configuration probabilities */
  sum = 0.0;
  for (i = 0; i < a->nconfig; i++) {
    if (getrealic(a->pr + i, buff, BUFFSIZ, fp) == EOF)
      error("read_agent_cfg:unexpected EOF");
    sum += a->pr[i];
  }
  if (sum > 0.0) {
    for (i = 0; i < a->nconfig; i++)
      a->pr[i] /= sum;
  }
  /* col1: label, col2: config 1, col3: config 2 ... */
  for (i = 0; i < a->npart; i++) {
    if ((a->label[i] = getwordic(buff, BUFFSIZ, fp)) == NULL)
      error("read_agent_cfg: expecting label; got EOF");
    a->label[i] = copystring(buff, BUFFSIZ);
    for (j = 0; j < a->nconfig; j++)
      if (getintic(a->cfg[i] + j, buff, BUFFSIZ, fp) == EOF)
	error("read_agent_cfg: expecting integer; got EOF");
  }
  if (getwordic(buff, BUFFSIZ, fp) != NULL) {
    fprintf(stderr, "\nError: read_agent_cfg: Did not exhaust input file %s",
	    fname);
    fprintf(stderr, "\n  Extra word: %s\n", buff);
#if 0
    exit(1);
#endif
  }
  fclose(fp);
  return (a);
}

/* read archeological data */
struct counts *read_counts(char *fname)
{
  struct counts *c;
  char    buff[BUFFSIZ];
  int     i, j;
  FILE   *fp;

  fp = mustopen(fname, "r");

  /* read number of faunal elements */
  c = alloc_counts();
  getintic(&(c->npart), buff, BUFFSIZ, fp);	/* number of faunal elements */
  getintic(&(c->ndataset), buff, BUFFSIZ, fp);	/* number of datasets */
  fflush(stdout);

  /* allocations */
  c->label = (char **) mustalloc(c->npart * sizeof(char *));
  c->y = (int **) alloc2d(c->ndataset, c->npart, sizeof(int));
  if (c->y == NULL)
    error("read_counts: alloc2d");

  /* read labels and data */
  for (i = 0; i < c->npart; i++) {
    if ((c->label[i] = getwordic(buff, BUFFSIZ, fp)) == NULL)
      error("read_counts: Expected label; got EOF");
    c->label[i] = copystring(buff, BUFFSIZ);
    for (j = 0; j < c->ndataset; j++)
      if (getintic(c->y[j] + i, buff, BUFFSIZ, fp) == EOF)
	error("read_counts: Expected integer; got EOF");
  }
  if (getwordic(buff, BUFFSIZ, fp) != NULL) {
    fprintf(stderr, "\nread_counts: Did not exhaust input file %s", fname);
    fprintf(stderr, "\n  Extra word: %s\n", buff);
#if 0
    exit(1);
#endif
  }
  fclose(fp);
  return (c);
}

void    pr_counts(struct counts *c)
{
  int     i, j;
  char    buff[20];

  printf("\n%3d #number of parts", c->npart);
  printf("\n%3d #number of data sets", c->ndataset);
  printf("\n#%-17s", "part");
  for (i = 0; i < c->ndataset; i++) {
    sprintf(buff, "DS%d", i + 1);
    printf(" %4s", buff);
  }

  for (i = 0; i < c->npart; i++) {
    printf("\n%-18s", c->label[i]);
    for (j = 0; j < c->ndataset; j++)
      printf(" %4d", c->y[j][i]);
  }
  putchar('\n');
}

/* print a matrix of reals */
void    prfmat(real ** A, int nrow, int ncol, const char *label)
{
  int     row, col;

  printf("\n%s:", label);
  for (row = 0; row < nrow; row++) {
    printf("\n#row %d\n", row);
    for (col = 0; col < ncol; col++) {
      if ((col + 1) % 10 == 0)
	putchar('\n');
      printf(" %9.6f ", A[row][col]);
    }
  }
}

/* print a matrix of ints */
void    primat(int **A, int nrow, int ncol, char *label)
{
  int     row, col;

  printf("\n%s:", label);
  for (row = 0; row < nrow; row++) {
    putchar('\n');
    for (col = 0; col < ncol; col++)
      printf("%d ", A[row][col]);
  }
}

/* print a vector of reals */
void    prfvec(real * x, int n, const char *label)
{
  int     j;

  printf("\n%s:", label);
  for (j = 0; j < n; j++)
    printf(" %f", x[j]);
}

/* print a vector of ints */
void    privec(int *x, int n, char *label)
{
  int     j;

  printf("\n%s:", label);
  for (j = 0; j < n; j++)
    printf(" %d", x[j]);
}

/* calculate mean and mscp from configuration data */
/* returns 0 on success; 1 on failure */
int     configure_agent(struct agent *a)
{
  int     i, j, k;

  if (a->npart <= 0 || a->nconfig <= 0)
    return (1);

  a->m = (real *) mustalloc(a->npart * sizeof(real));
  a->F = (real **) alloc2d(a->npart, a->npart, sizeof(real));
  if (a->F == NULL || a->m == NULL)
    error("no memory in configure_agent");

  /* calculate mean vector (m) and SSCP matrix (F) */
  for (i = 0; i < a->npart; i++) {
    a->m[i] = a->F[i][i] = 0.0;
    for (k = 0; k < a->nconfig; k++)
      a->m[i] += a->pr[k] * a->cfg[i][k];
    for (j = 0; j <= i; j++) {
      a->F[i][j] = 0.0;
      for (k = 0; k < a->nconfig; k++)
	a->F[i][j] += a->pr[k] * a->cfg[i][k] * a->cfg[j][k];
      if (i != j)
	a->F[j][i] = a->F[i][j];
    }
  }
  return (0);
}

/** open a file, abort if unsuccessful **/
FILE   *mustopen(char *name, const char *mode)
{
  FILE   *fp;

  fp = fopen(name, mode);
  if (fp == NULL) {
    fflush(stdout);
    fprintf(stderr, "\nCan't open file \"%s\" with mode \"%s\"\n",
	    name, mode);
    exit(1);
  }
  return (fp);
}

int     lbls_differ(char **lbl1, char **lbl2, int len)
{
  int     i, differ = 0;

  for (i = 0; i < len; i++) {
    if (strcmp(lbl1[i], lbl2[i])) {
      fprintf(stderr, "\n%s != %s", lbl1[i], lbl2[i]);
      differ = 1;
    }
  }
  return (differ);
}

/* check consistency of part numbers and labels */
void    check_labels(int nagent, struct agent **agent,
		     struct bonedef *bones,
		     struct counts *counts)
{
  int     i;

  fflush(stdout);
  if (agent != NULL) {
    if (nagent <= 0)
      error("check_labels: agent!=NULL but nagent <= 0");
    for (i = 0; i < nagent - 1; i++) {
      if (agent[i]->npart != agent[i + 1]->npart) {
	fprintf(stderr, "\nError: .cfg files have different numbers of parts");
	fprintf(stderr, "\nagent[%d]->npart=%d agent[%d]->npart=%d\n",
		i, agent[i]->npart, i + 1, agent[i + 1]->npart);
	exit(1);
      }
      if (lbls_differ(agent[i]->label, agent[i + 1]->label, agent[i]->npart)) {
	fprintf(stderr, "\nLabels disagree btw agents %d and %d.\n",
		i, i + 1);
	exit(1);
      }
    }
  }
  if (counts != NULL) {
    if (bones != NULL) {
      if (counts->npart != bones->npart) {
	fprintf(stderr,
	  "\nError:.cnt and .bdf files have different numbers of parts");
	fprintf(stderr, "\n.cnt: %d .bdf: %d\n", counts->npart, bones->npart);
	exit(1);
      }
      if (lbls_differ(bones->label, counts->label, bones->npart))
	error("parts are not labeled consistently in .cnt and .bdf files");
    }
    if (agent != NULL) {
      if (counts->npart != agent[0]->npart) {
	fprintf(stderr,
	  "\nError:.cnt and .cfg files have different numbers of parts");
	fprintf(stderr, "\n.cnt: %d .cfg: %d\n",
		counts->npart, agent[0]->npart);
	exit(1);
      }
      if (lbls_differ(counts->label, agent[0]->label, counts->npart))
	error("parts are not labeled consistently in .cnt and .cfg files");
    }
  } else {
    /* the only comparison left is between bones and agent[0] */
    if (bones != NULL && agent != NULL) {
      if (bones->npart != agent[0]->npart) {
	fprintf(stderr,
		"\nError:.bdf and .cfg files have different numbers of parts");
	fprintf(stderr, "\n.bdf: %d .cfg: %d\n",
		bones->npart, agent[0]->npart);
	exit(1);
      }
      if (lbls_differ(bones->label, agent[0]->label, bones->npart))
	error("parts are not labeled consistently in .bdf and .cfg files");
    }
  }
}
