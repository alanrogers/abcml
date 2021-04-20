/****************************************************************
    gnaw: Subject an assemblage to density-mediated attrition
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
real beta=1.0;              /* intensity of attrition */
real scale_factor = -1.0;   /* for re-scaling attrition */

void usage(const char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nError on command line: %s", msg);
  fputs("\n\nusage: gnaw xxx.bdf xxx.cnt [options]",
	stderr);
  fputs("\nwhere:",stderr);
  fputs("\n      .bdf file defines characteristics of bones", stderr);
  fputs("\n      .cnt file configures agent",
	stderr);
  fprintf(stderr,"\nOptions may include");
  fprintf(stderr,"\n -b x  Set beta.  Def: %g", beta);
  fprintf(stderr,"\n -D x  Set sensitivity to x/density. Def: automatic");
  putc('\n', stderr);
  exit(1);
}

int main(int argc, char **argv)
{
  int i, data, part, bone, n;
  int n_before=0, n_after=0;
  real p, u;
  struct bonedef *bones=NULL; /* .bdf data */
  struct counts *counts=NULL;

  /* header */
  header("GNAW", "Subject an assemblage to attrition.", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);
  printf("\n#beta = %g", beta);

  /* process command line arguments */
  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case CNT_ARG:
      if(counts != NULL)
	  usage("Only one .cnt file is allowed.");
      counts = read_counts(argv[i]);
      break;
    case BDF_ARG:
      if(bones != NULL)
	usage("Only one .bdf file is allowed.");
      bones = read_bonedef(argv[i], &scale_factor);
      printf("\n#Configured bones from file %s", argv[i]);
      break;
    case FLAG_ARG:
      switch(argv[i][1])
      {
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
      default:
	usage(argv[i]);
      }
      break;
    default:
      usage(argv[i]);
    }
  }

  if(bones == NULL)
    usage("Missing .bdf file");
  if(counts == NULL)
    usage("Missing .cnt file");

  /* make sure data are consistent */
  check_labels(0, NULL, bones, counts);

  /* initialize random number generator */
  (void) initrand(0);

  for(data=0; data < counts->ndataset; data++)
  {
    for(part=0; part < counts->npart; part++)
    {
      n = 0;
      p = exp(-beta * bones->sensitivity[part]);
      for(bone=0; bone < counts->y[data][part]; bone++)
      {
	n_before++;
	u = uni();
	if(u < p)
	{
	  n++;
	  n_after++;
	}
      }
      counts->y[data][part] = n;
    }
  }

  printf("\n# Sensitivity to attrition: %0.20g / density", scale_factor);
  printf("\n# %d / %d bones survived attrition", n_after, n_before);
  pr_counts(counts);
  putchar('\n');
  return(0);
}
