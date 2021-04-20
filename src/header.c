/****************************************************************
    header: print a header identifying a program
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
#include <ctype.h>
#include <string.h>
#include "version.h"
#include "header.h"
#define BUFFSIZE 80
void header(const char *program, const char *description, FILE *fp)
{
    char buff[BUFFSIZE], word[20];
    int i;

    centerline(program,fp);
    centerline(description,fp);
    centerline("by Alan R. Rogers",fp);
    sprintf(buff, "Version %s", VERSION);
    centerline(buff,fp);
    sprintf(buff, "Git version %s", GIT_VERSION);
    centerline(buff,fp);
    for(i=0; program[i] != '\0'; i++)
        word[i] = (char) tolower(program[i]);
    word[i] = '\0';
    sprintf(buff,"Type `%s -- ' for help", word);
    centerline(buff, fp);
    fflush(fp);
}

/** center a character string on output line */
void centerline(const char *s, FILE *fp)
{
  int pad;

  pad = 80 - strlen(s);
  if(pad < 0)
    pad = 0;
  else
    pad /= 2;
  putc('#', fp);   /* begin line with comment character */
  pad--;
  while(--pad > 0)
    putc(' ', fp);
  fputs(s, fp);
  putc('\n', fp);
}
