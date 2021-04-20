/****************************************************************
    getcic: Input routines that ignore comments
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
#include "mytypes.h"
#include "getcic.h"


/****************************************************************
Modify getc() to ignore comments, which begin with START_COMMENT and
end with newline or EOF.  When a comment appears getcic() returns the
newline (or EOF) character only.
****************************************************************/
int getcic(FILE *fp)
{
  int c;

  c = getc(fp);
  if(c==START_COMMENT)
  {
    do{
      c = getc(fp);             /* read to end of line */
    }while(c != '\n' && c != EOF);
  }
  return(c);
}
/*
 * get a word, delimited by whitespace, from file ifp
 * Ignore comments delimited by START_COMMENT...\n
 */
char *getwordic(char *buff, int bufsiz, FILE * ifp)
{
  char   *bp;
  int     c;

  bp = buff;
  do
  {
    c = getcic(ifp);
  }while ( isspace(c) && c != EOF);

  if (c == EOF)
    return ( NULL );
  while (!isspace(c) && c != EOF)
  {
    if (--bufsiz < 1)
    {
      ungetc(c, ifp);
      break;
    }
    *bp++ = (char) c;  /* compress int to char */
    c = getcic(ifp);
  }
  if(isspace(c))
    ungetc(c, ifp);
  *bp = (char) '\0';
  return (buff);
}

/* read a floating point number and place its value into x */
int getrealic(real *x, char *buff, int bufsiz, FILE *ifp)
{
  if(getwordic(buff, bufsiz, ifp) == NULL)
    return(EOF);
  *x = strtod(buff, NULL);
  return(0);
}

/* read an int and place its value into i */
int getintic(int *i, char *buff, int bufsiz, FILE *ifp)
{
  if(getwordic(buff, bufsiz, ifp) == NULL)
    return(EOF);
  sscanf(buff,"%d", i);
  return(0);
}


