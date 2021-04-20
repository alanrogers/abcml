/****************************************************************
    abcsim: Simulate bone counts
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
#include "unirand.h"
#include "header.h"
#define YES(x) ((x) == 0 ? "No" : "Yes")
#define ON(x) ((x) == 0 ? "Off" : "On")

/* function prototypes */
void usage(const char *msg);
int real_lookup(real value, real *vec, int length);

/* external variables */
int npart;                  /* number of parts */
int nagent=0;               /* number of agents */
int K=100;                  /* number of animals in each assemblage */
int ndataset=1;             /* number of simulated datasets */
real *alpha;                /* vector of alpha values */
real beta=0.0;              /* intensity of attrition */
real scale_factor = -1.0;   /* for re-scaling attrition */

void usage(const char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nError on command line: %s", msg);
  fputs("\n\nusage: abcsim [options] xxx.bdf yyy.cfg zzz.cfg [...]",
	stderr);
  fputs("\nwhere:",stderr);
  fputs("\n      .bdf file defines characteristics of bones", stderr);
  fputs("\n      .cfg files configure agents; at least 2 are needed",
	stderr);
  fputs("\nFlag options may include:", stderr);
  fprintf(stderr,"\n -a x,x,...,x");
  fprintf(stderr,"\n       Set alpha vector.  Def: alpha[i] = 1/nagent");
  fprintf(stderr,"\n -b x  Set beta (attrition parameter). Def: %g", beta);
  fprintf(stderr,"\n -D x  Set sensitivity to x/density. Def: automatic");
  fprintf(stderr,"\n -K x  Set number of animals in assemblage. Def: %d", K);
  fprintf(stderr,"\n -n x  Set number of simulated datasets. Def: %d",
	  ndataset);
  putc('\n', stderr);
  exit(1);
}

/**************************************************************** 
Look up index of given value in a vector of reals sorted in ascending 
order.  Returns smallest integer i such that vec[i] >= value.
****************************************************************/
int real_lookup(real value, real *vec, int length)
{
  int lo, mid, hi;

  lo = 0;
  hi = length-1;
  if(vec[hi] < vec[lo] || value > vec[hi])
    return(-1);
  if(hi==0)
    return(0);

  while(lo < hi)
  {
    mid = lo + (hi-lo)/2;
    if(mid == lo)
      break;
    if(vec[mid] < value)
      lo = mid;
    else
      hi = mid;
  }
  if(vec[lo] >= value)
    hi = lo;

  assert(hi>=0);
  assert(hi<length);
  assert(vec[hi] >= value);
  assert(hi==0 || vec[hi-1]<value);

  return(hi);
}

int main(int argc, char **argv)
{
  int i, j, k, iagent, idata, ianimal, ipart, iconfig, tbone;
  real sbone, sum, u;
  int *n_per_agent;
  char *s;
  struct agent **agent=NULL;  /* vector of pointers to agents */
  struct bonedef *bones=NULL; /* .bdf data */
  struct counts *counts;      /* archeological data */
  real **cum;  /* cum[i][j] = cumulative prob of config j for agent i */
  real *survive;  /* vector of survival probabilities */

  /* header */
  header("ABCSIM", "Simulate Bone Counts", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);

  /* count .cfg arguments to determine number of agents */
  nagent = 0;
  for(i=1; i<argc; i++)
  {
    if(classify(argv[i]) == CFG_ARG)
      nagent += 1;
  }

  if(nagent < 1)
    usage("Need at least one .cfg file");

  agent = (struct agent **) mustalloc( nagent * sizeof(struct agent *));
  alpha = (real *) mustalloc( nagent * sizeof(real));
  for(i=0; i<nagent; i++)
    alpha[i] = 1.0/nagent;    /* default for alpha */

  /* process command line arguments */
  j = 0;
  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case FLAG_ARG:
      switch(argv[i][1])
      {
      case 'a' :
      	if(++i >= argc)
	  usage("Missing arg for -a");
	/* parse string of reals in format x,x,...,x */
	k = 0;
	s = strtok(argv[i], ",");
	while(s!=NULL && k < nagent)
	{
	  alpha[k++] = strtod(s, NULL);
	  s = strtok(NULL, ",");
	}
	if(k < nagent)
	  usage("The -a string has too few values");
	break;
      case 'b':
      	if(++i >= argc)
	  usage("Missing arg for -b");
	beta = strtod(argv[i], NULL);
	break;
      case 'D':
      	if(++i >= argc)
	  usage("Missing arg for -D");
	if(bones != NULL)
	  error("-D must come before .bdf on command line");
	scale_factor = strtod(argv[i], NULL);
	break;
      case 'K':
      	if(++i >= argc)
	  usage("Missing arg for -K");
	K = strtol(argv[i], NULL, 10);
	break;
      case 'n':
      	if(++i >= argc)
	  usage("Missing arg for -K");
	ndataset = strtol(argv[i], NULL, 10);
	break;
      default:
	usage(argv[i]);
      }
      break;
    case CFG_ARG:
      agent[j] = read_agent_cfg(argv[i]);
      printf("\n#Configured agent %d from file %s",
	     j, argv[i]);
      j += 1;
      break;
    case BDF_ARG:
      if(bones != NULL)
	usage("Only one .bdf file is allowed.");
      bones = read_bonedef(argv[i], &scale_factor);
      printf("\n#Configured bones from file %s", argv[i]);
      break;
    case CNT_ARG:
      usage("No .cnt file is allowed.");
    default:
      usage(argv[i]);
    }
  }

  if(bones == NULL)
    usage("Missing .bdf file");

  /* make sure data are consistent */
  check_labels(nagent, agent, bones, NULL);

  /* make alpha sum to 1 */
  sum = 0.0;
  for(i=0; i<nagent; i++)
    sum += alpha[i];
  for(i=0; i<nagent; i++)
    alpha[i] /= sum;

  /* get n_per_agent */
  n_per_agent = (int *) mustalloc( nagent * sizeof(int) );
  k = 0;
  for(i=0; i<nagent; i++)
  {
    n_per_agent[i] = (alpha[i] * K) + 0.5;  /* round to nearest int */
    k += n_per_agent[i];
  }

  /* reset alpha and K to reflect rounding */
  K = k;
  for(i=0; i<nagent; i++)
    alpha[i] = ((real) n_per_agent[i]) / ((real) K);

  /* get cumulative configuration probabilities */
  cum = (real **) mustalloc(nagent * sizeof(real *));
  for(i=0; i<nagent; i++)
  {
    cum[i] = (real *) mustalloc(agent[i]->nconfig * sizeof(real));
    cum[i][0] = agent[i]->pr[0];
    for(j=1; j<agent[i]->nconfig; j++)
      cum[i][j] = cum[i][j-1] + agent[i]->pr[j];
  }

  /* initialize vector of survival probabilities */
  survive = (real *) mustalloc( bones->npart * sizeof(real) );
  if(beta > 0.0)
    for(i=0; i<bones->npart; i++)
      survive[i] = exp(-beta * bones->sensitivity[i]);
  else
    for(i=0; i<bones->npart; i++)
      survive[i] = 1.0;

  /* make a structure of type struct counts */
  counts = alloc_counts();
  counts->npart = bones->npart;
  counts->ndataset = ndataset;
  counts->label = bones->label;
  counts->y = (int **) alloc2d(ndataset, bones->npart, sizeof(int));

  /* initialize random number generator */
  initrand(0);

  /* loop over datasets */
  for(idata=0; idata < ndataset; idata++)
  {
    /* initialize data vector */
    for(ipart=0; ipart<counts->npart; ipart++)
      counts->y[idata][ipart] = 0;

    for(iagent=0; iagent<nagent; iagent++)
    {
      for(ianimal=0; ianimal<n_per_agent[iagent]; ianimal++)
      {
	u = uni();   /* uniform r.v. */
	iconfig = real_lookup(u, cum[iagent], agent[iagent]->nconfig);
	if(iconfig < 0)    /* lookup error */
	{
	  fflush(stdout);
	  fprintf(stderr, "\nCouldn't look up %g in:", u);
	  for(k=0; k<agent[iagent]->nconfig; k++)
	    fprintf(stderr," %g", cum[iagent][k]);
	  putc('\n', stderr);
	  exit(1);
	}
	for(ipart=0; ipart<agent[iagent]->npart; ipart++)
	{
	  for(i=0; i<agent[iagent]->cfg[ipart][iconfig]; i++)
	    if(uni() <= survive[ipart])
	      counts->y[idata][ipart] += 1;
	}
      }
    }
  } /* end dataset loop */

  printf("\n### Simulated data");
  printf("\n# Number of animals in each assemblage: %d", K);
  printf("\n# Sensitivity to attrition: %0.20g / density", scale_factor);
  printf("\n# alpha: ");
  for(i=0; i<nagent; i++)
    printf(" %g", alpha[i]);

  printf("\n# beta = %g:", (beta >= 0.0 ? beta : 0.0));
  if(beta > 0)
  {
    /* calculate expected survival fraction */
    tbone=0;
    sbone=0.0;
    for(i=0; i<bones->npart; i++)
    {
      tbone += bones->live[i];
      sbone += bones->live[i] * survive[i];
    }
    printf(" A fraction %g of a whole animal will survive attrition\n",
	   sbone/tbone);
  }else
    printf(" No attrition");
  pr_counts(counts);
  putchar('\n');
  return(0);
}
