#include <stdio.h>
#include <stdlib.h>
#include "mytypes.h"
#include "io.h"

int main(int argc, char **argv)
{
  struct bonedef *b;
  real scale = 1.0;

  if(argc != 2)
  {
    fprintf(stderr,"usage: xbonedef inputfile.bdf\n");
    exit(1);
  }

  b = read_bonedef(argv[1], &scale);

  pr_bonedef(b);
  return 0;
}
