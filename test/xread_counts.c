#include <stdio.h>
#include <stdlib.h>
#include "mytypes.h"
#include "io.h"

int attrition=1;

void main(int argc, char **argv)
{
  struct counts *a;

  if(argc != 2)
  {
    fprintf(stderr,"usage: xread_counts inputfile.cnt\n");
    exit(1);
  }
  a = read_counts(argv[1]);

  pr_counts(a);
  putchar('\n');
  exit(0);
}
