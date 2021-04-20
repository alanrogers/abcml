/****************************************************************
    misc: Miscellaneous functions
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
#include <ctype.h>
#include <string.h>
#include "mytypes.h"
#include "misc.h"

/* get a quantile from a sorted array of data */
real quantile(real p, real *vec, int size)
{
  int i;

  i = p*size - 1;
  if(i<0)
    i = 0;
  return(vec[i]);
}

/* Compare reals */
int compar(real *x, real *y)
{

  if(*x < *y)
    return(-1);
  if(*x > *y)
    return(1);
  return(0);
}

/** print an error message and quit **/
void error(const char *s)
{
  fflush(stdout);
  fprintf(stderr,"\nERROR: %s\n", s);
  exit(1);
}

/** allocate memory using malloc; abort on failure **/
void *mustalloc(size_t bytes)
{
  void *p;

  p = malloc(bytes);
  if(p==NULL)
  {
    fflush(stdout);
    fprintf(stderr,"\nmustalloc: Can't allocate %d bytes of memory.\n",
	    bytes);
    exit(1);
  }
  return(p);
}

/* lower case string */
char *strlwr(char *u)
{
  char *v;

  if(u==NULL)
    return(NULL);
  
  for(v = u; *v != '\0'; v++)
    *v = (char) tolower(*v);

  return(u);
}

/* classify command-line arguments */
int classify(char *arg)
{
  char *dot;
  
  if(arg==NULL)
    return(NULL_ARG);

  if(*arg == '-')
    return(FLAG_ARG);

  dot = strrchr(arg, '.');  /* find last dot */
  if(dot == NULL)
    return(UNKNOWN_ARG);
  dot += 1;

  (void) strlwr(dot);
  
  if(strncmp(dot, "cfg", 3)==0)
    return(CFG_ARG);

  if(strncmp(dot, "bdf", 3)==0)
    return(BDF_ARG);

  if(strncmp(dot, "cnt", 3)==0)
    return(CNT_ARG);

  if(strncmp(dot, "mau", 3)==0)
    return(MAU_ARG);

  if(strncmp(dot, "wgt", 3)==0)
    return(WGT_ARG);

  return(UNKNOWN_ARG);
}

/****************************************************************
This version of alloc2d allocates each row with a separate call to
malloc.
 ****************************************************************/
void **alloc2d(int rows, int cols, int elsize)
{
  void **m;
  int i;

  m = (void **) malloc( rows * sizeof(void *));
  if(m==NULL)
  {
    fprintf(stderr,"\nalloc2d: memory");
    exit(1);
  }
  for(i=0; i<rows; i++)
  {
    m[i] = (void *) malloc( (unsigned) (cols * elsize) );
    if(m[i] == NULL)
    {
      fprintf(stderr,"\nalloc2d: memory");
      exit(1);
    }
  }
  return(m);
}

void free2d(void **m, int rows)
{
  int i;

  for(i=0; i<rows; i++)
    free(m[i]);
  free(m);
}
