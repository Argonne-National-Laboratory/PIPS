#include <stdio.h>
#include <stdlib.h>

#include "globals2.h"

void
exit2R (const char *msg)
{
  fprintf (stderr, "%s\n", msg);
  exit(1);
} // exit2R

void
msg2R (const char *msg)
{
  fprintf (stderr, "%s\n", msg);
  return;
}
