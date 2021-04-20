/****************************************************************
    mau2cfg: Translate from .mau format to .cfg format
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
#include <time.h>
#include <assert.h>
#include "mytypes.h"
#include "misc.h"
#include "io.h"
#include "header.h"
#include "getcic.h"
#include "zeroin.h"
#define YES(x) ((x) == 0 ? "No" : "Yes")
#define ON(x) ((x) == 0 ? "Off" : "On")
#define BUFFSIZ 100

/* define an agent of deposition in mau format*/
struct mau {
  int npart, nconfig;
  char **label;            /* labels for faunal elements */
  real **mau;              /* npart X nconfig matrix of configurations */
  real *pr;                /* nconfig-vector of probabilities */
  real *m;                 /* mean vector */
  real **F;                /* mean squares and cross products */
};

/* function prototypes */
void usage(const char *msg);
struct mau *alloc_mau(void);
struct mau *read_mau(char *fname);
void pr_mau(struct mau *a, struct bonedef *b);

/* external variables */
int npart;                  /* number of parts */
int nagent=0;               /* number of agents */
int K=100;                  /* number of animals in each assemblage */
int ndataset=1;             /* number of simulated datasets */
real *alpha;                /* vector of alpha values */
real beta=0.0;              /* intensity of attrition */

void usage(const char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nError on command line: %s", msg);
  fputs("\n\nusage: mau2cfg xxx.bdf xxx.mau",
	stderr);
  fputs("\nwhere:",stderr);
  fputs("\n      .bdf file defines characteristics of bones", stderr);
  fputs("\n      .mau file is a .cfg file in mau format",
	stderr);
  putc('\n', stderr);
  exit(1);
}

/* allocate and initialize an agent structure */
struct mau *alloc_mau(void)
{
  struct mau *a;

  a = (struct mau *) malloc(sizeof(struct mau));
  if(a==NULL)
    return(NULL);
  a->npart = 0;
  a->label = NULL;
  a->pr = NULL;
  a->m = NULL;
  a->F = NULL;
  a->mau = NULL;
  return(a);
}

/* read .mau file */
struct mau *read_mau(char *fname)
{
  struct mau *a;
  char buff[BUFFSIZ];
  int i, j;
  real sum;
  FILE *fp;

  fp = mustopen(fname, "r");

  /* read number of faunal elements */
  a = alloc_mau();
  getintic(&(a->npart), buff, BUFFSIZ, fp); /* number of faunal elements */
  getintic(&(a->nconfig), buff, BUFFSIZ, fp); /* number of configurations */
  fflush(stdout);

  /* allocations */
  a->label = (char **) mustalloc(a->npart * sizeof(char *));
  a->pr = (real *) mustalloc(a->nconfig * sizeof(real));
  a->mau = (real **) alloc2d(a->npart, a->nconfig, sizeof(real));
  if(a->mau == NULL)
    error("read_mau: alloc2d");

  /* read and configuration probabilities */
  sum=0.0;
  for(i=0; i<a->nconfig; i++)
  {
    if(getrealic(a->pr + i, buff, BUFFSIZ, fp) == EOF)
      error("read_mau:unexpected EOF");
  }

  /* col1: label, col2: config 1, col3: config 2 ... */
  for(i=0; i<a->npart; i++)
  {
    if( (a->label[i] = getwordic(buff, BUFFSIZ, fp)) == NULL)
      error("read_mau: expecting label; got EOF");
    a->label[i] = copystring(buff, BUFFSIZ);
    for(j=0; j<a->nconfig; j++)
      if(getrealic(a->mau[i]+j, buff, BUFFSIZ, fp) == EOF)
	error("read_mau: expecting float; got EOF");
  }
  if(getwordic(buff, BUFFSIZ, fp) != NULL)
  {
    fprintf(stdout,"\nError: read_mau: Did not exhaust input file %s",
	    fname);
    fprintf(stdout,"\n  Extra word: %s\n", buff);
    exit(1);
  }

  fclose(fp);
  return(a);
}

/* print a mau structure */
void pr_mau(struct mau *a, struct bonedef *b)
{
  int i, j;

  printf("\n%6d  # number of parts", a->npart);
  printf("\n%6d  # number of configurations", a->nconfig);

  if(a->pr != NULL)
  {
    printf("\n#Probabilities of configurations:");
    printf("\n%16s", "");
    for(i=0; i<a->nconfig; i++)
      printf(" %g", a->pr[i]);
    printf("\n#label");
    for(i=0; i<a->npart; i++)
    {
      printf("\n%-16s", a->label[i]);
      for(j=0; j<a->nconfig; j++)
	printf(" %4d", (int) ((a->mau[i][j] * b->live[i]) + 0.5));
    }
  }
  if(a->m != NULL)
  {
    printf("\n\nMean:");
    for(i=0; i<a->npart; i++)
    {
      printf("\n%-16s", a->label[i]);
      printf(" %7.3f", a->m[i]);
    }
  }
  if(a->F != NULL)
    prfmat(a->F, a->npart, a->npart, "\nF");

  putchar('\n');
}

int main(int argc, char **argv)
{
  int i;
  struct mau *mau=NULL;
  struct bonedef *bones=NULL;

  /* header */
  header("mau2cfg", "Translate .mau data into .cfg format", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);

  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case MAU_ARG:
      if(mau != NULL)
	usage("Only one .mau file is allowed.");
      mau = read_mau(argv[i]);
      printf("\n#Read .mau data from file %s", argv[1]);
      break;
    case BDF_ARG:
      if(bones != NULL)
	usage("Only one .bdf file is allowed.");
      bones = read_bonedef(argv[i], NULL);
      break;
    default:
      usage(argv[i]);
    }
  }

  if(mau == NULL)
    usage("No .mau file was specified");

  if(bones == NULL)
    usage("No .bdf file was specified");
  
  /* make sure data are consistent */
  if(mau->npart != bones->npart)
    error("input files don't agree about the number of parts");
  if(lbls_differ(mau->label, bones->label, mau->npart))
    error("input files don't agree about labels");

  /* output */
  pr_mau(mau, bones);

  putchar('\n');
  return(0);
}
