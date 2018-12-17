#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <ctype.h>
#include <errno.h>

#define GDX_MAIN
#include "gdxstatic.h"

static int objectCount = 0;

#define XCreate xcreate
#define XFree   x2free

#if defined(__cplusplus)
extern "C" {
#endif
void XCreate (gdxHandle_t *pgdx);
void XFree   (gdxHandle_t *pgdx);

/* gdxCreate: return true always */
int gdxCreate (gdxHandle_t *pgdx, char *msgBuf, int msgBufSize)
{
  msgBuf[0] = '\0';
  XCreate(pgdx);
  if (*pgdx == NULL) {
    strcpy(msgBuf,"Error while creating object");
    return 0;
  }
  objectCount++;
  return 1;                     /* return true on successful library load */
} /* gdxCreate */

int gdxFree (gdxHandle_t *pgdx)
{
  XFree(pgdx);
  objectCount--;
  return 1;
} /* gdxFree */
#if defined(__cplusplus)
}
#endif
