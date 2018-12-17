#include "p3io.h"
#include "idglobal_p3.h"

/**** C code included from idglobal_p3.pas(39:1): 6 lines ****/
#if defined(_WIN32)
#include <windows.h>
#else
#include <sys/time.h>
#include <time.h>
#endif

Function(SYSTEM_cardinal ) IDGLOBAL_P3_gettickcount(void)
{
  SYSTEM_cardinal result;

  /**** C code included from idglobal_p3.pas(54:1): 15 lines ****/
#if defined(_WIN32)
  result = GetTickCount();
#else
# if 0
  result = clock() / (CLOCKS_PER_SEC / 1000);
# else
{
  struct timeval tv;

  gettimeofday (&tv, NULL);
  result = tv.tv_sec; /* force tv_sec to take the type of the result */
  result = result * 1000 + tv.tv_usec / 1000;
}
# endif
#endif
  return result;
}  /* gettickcount */

Function(SYSTEM_cardinal ) IDGLOBAL_P3_gettickdiff(
  SYSTEM_cardinal aoldtickcount,
  SYSTEM_cardinal anewtickcount)
{
  SYSTEM_cardinal result;

  if (anewtickcount >= aoldtickcount) { 
    result = anewtickcount - aoldtickcount;
  } else {
    result = 2147483647;
    result = result + 2147483647 + 1 - aoldtickcount + anewtickcount;
  } 
  return result;
}  /* gettickdiff */

/* unit idglobal_p3 */
void _Init_Module_idglobal_p3(void)
{
} /* _Init_Module_idglobal_p3 */

void _Final_Module_idglobal_p3(void)
{
} /* _Final_Module_idglobal_p3 */

