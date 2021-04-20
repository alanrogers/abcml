/****************************************************************
    tabcfg: Tabulate configurations in a .cfg file, combining duplicates
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

/* external variables */
int npart;                  /* number of parts */
int nagent=0;               /* number of agents */

void usage(const char *msg)
{
  fflush(stdout);
  if(msg != NULL)
    fprintf(stderr,"\n\nError on command line: %s", msg);
  fputs("\n\nusage: tabcfg xxx.cfg ",
	stderr);
  fputs("\nwhere:",stderr);
  fputs("\n  xxx.cfg is an agent configuration file", stderr);
  putc('\n', stderr);
  exit(1);
}

int main(int argc, char **argv)
{
  int i, j, k, l, ncfg, same, match;
  struct agent *agent=NULL, *new;  /* pointer to agent */
  char *fname;

  /* header */
  header("TABCFG", "Tabulate Configurations", stdout);

  /* echo command line */
  printf("\n#Cmd line: %s", argv[0]);
  for(i=1; i<argc; i++)
    printf(" %s", argv[i]);

  /* process command line arguments */
  for(i=1; i<argc; i++)
  {
    switch(classify(argv[i]))
    {
    case FLAG_ARG:
      switch(argv[i][1])
      {
      default:
	usage(argv[i]);
      }
      break;
    case CFG_ARG:
      if(agent != NULL)
	usage("Only one .cfg file is allowed");
      fname = argv[i];
      agent = read_agent_cfg(fname);
      printf("\n#Configured agent from file %s", fname);
      break;
    default:
      usage(argv[i]);
    }
  }

  if(agent==NULL)
    usage("Must specify .cfg file on command line");

  new = alloc_agent();
  new->npart = agent->npart;
  new->label = agent->label;
  new->pr = (real *) mustalloc(agent->nconfig * sizeof(real));
  new->cfg = (int **) alloc2d(agent->npart, agent->nconfig, sizeof(int));
  if(new->cfg == NULL)
    error("main: no memory");
  new->nconfig = 0;

  /* cfg 0 is same for new and old */
  for(j=0; j<agent->npart; j++)
    new->cfg[j][0] = agent->cfg[j][0];
  new->pr[0] = agent->pr[0];

  ncfg = 1;
  for(i=1; i<agent->nconfig; i++)
  {
    match = -1;
    for(k=0; k<ncfg; k++)
    {
      /* is cfg[k] == cfg[i]? */
      same = 1;
      for(l=0; l<agent->npart; l++)
      {
	if(new->cfg[l][k] != agent->cfg[l][i])
	{
	  same = 0;
	  break;
	}
      }
      if(same)  /* cfg[i] is the same as new[k] */
      {
	match = k;
	break;
      }
    }
    /* if match >= 0 then cfg[i] = cfg[match] */
    if(match >= 0)
      new->pr[match] += agent->pr[i];
    else
    {
      new->pr[ncfg] = agent->pr[i];
      for(l=0; l<agent->npart; l++)
	new->cfg[l][ncfg] = agent->cfg[l][i];
      ncfg += 1;
    }
  }

  new->nconfig = ncfg;
  printf("\n#Input file: %s", fname);
  pr_agent(new);
  putchar('\n');
  return(0);
}
