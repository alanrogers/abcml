/****************************************************************
    cplcfg: Create a .cfg file that is the complement of an existing
    .cfg file.
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

void usage(const char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nError on command line: %s", msg);
  fputs("\n\nusage: cplcfg xxx.bdf xxx.cfg",
	stderr);
  fputs("\nwhere:",stderr);
  fputs("\n      .bdf file defines characteristics of bones", stderr);
  fputs("\n      .cfg file configures agent",
	stderr);
  putc('\n', stderr);
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j;
  struct agent *agent=NULL;  /* pointer to agent */
  struct bonedef *bones=NULL; /* .bdf data */

  /* header */
  header("CPLCFG", "Compute Complement of a .cfg File", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);

  /* process command line arguments */
  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case CFG_ARG:
      if(agent != NULL)
	usage("Only one .cfg file is allowed.");
      agent = read_agent_cfg(argv[i]);
      printf("\n#Configured agent from file %s", argv[i]);
      break;
    case BDF_ARG:
      if(bones != NULL)
	usage("Only one .bdf file is allowed.");
      bones = read_bonedef(argv[i], NULL);
      printf("\n#Configured bones from file %s", argv[i]);
      break;
    default:
      usage(argv[i]);
    }
  }

  if(bones == NULL)
    usage("Missing .bdf file");
  if(agent == NULL)
    usage("Missing .cfg file");

  /* make sure data are consistent */
  check_labels(1, &agent, bones, NULL);

  /* complement configurations */
  for(i=0; i<agent->npart; i++)
    for(j=0; j<agent->nconfig; j++)
    {
      if(agent->cfg[i][j] > bones->live[i])
	fprintf(stderr,"\n%d = cfg[%d][%d] > live[%d] = %d",
		agent->cfg[i][j], i, j, i, bones->live[i]);
      agent->cfg[i][j] = bones->live[i] - agent->cfg[i][j];
    }

  pr_agent(agent);
  putchar('\n');
  return(0);
}
