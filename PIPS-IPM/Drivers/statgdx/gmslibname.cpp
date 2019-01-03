#include "p3io.h"
#include "gmslibname.h"


Function(SYSTEM_ansichar *) GMSLIBNAME_gamslibnameold(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  _P3strcpy(result,_len_ret,s);
  _P3strcat(result,_len_ret,result,_P3str1("\00264"));
  return result;
}  /* gamslibnameold */
/**** C code included from gmslibname.pas(31:1): 3 lines ****/
#if ! defined (_WIN32)
# include <sys/utsname.h>
#endif

Function(SYSTEM_ansichar *) GMSLIBNAME_gamslibnamep3(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  /**** C code included from gmslibname.pas(41:1): 79 lines ****/
#if defined(_WIN32)
# if ! defined(GMS_DLL_PREFIX)
#   define GMS_DLL_PREFIX ""
# endif
# if ! defined(GMS_DLL_EXTENSION)
#   define GMS_DLL_EXTENSION ".dll"
# endif
# if ! defined(GMS_DLL_SUFFIX)
#   if defined(_WIN64)
#     define GMS_DLL_SUFFIX "64"
#   else
#     define GMS_DLL_SUFFIX ""
#   endif
# endif

#else  /* start non-Windows */

# if ! defined(GMS_DLL_PREFIX)
#   define GMS_DLL_PREFIX "lib"
# endif
# if ! defined(GMS_DLL_EXTENSION)
#   if defined(__APPLE__)
#     define GMS_DLL_EXTENSION ".dylib"
#   else
#     define GMS_DLL_EXTENSION ".so"
#   endif
# endif
# if ! defined(GMS_DLL_SUFFIX)
#   if defined(__WORDSIZE)
#     if 64 == __WORDSIZE
#       define GMS_DLL_SUFFIX "64"
#     else
#       define GMS_DLL_SUFFIX ""
#     endif
#   elif defined(__SIZEOF_POINTER__)
#     if 4 == __SIZEOF_POINTER__
#       define GMS_DLL_SUFFIX ""
#     elif 8 == __SIZEOF_POINTER__
#       define GMS_DLL_SUFFIX "64"
#     endif
#   elif ( defined(__xlc__) || defined(__xlC__) )
#     if defined(__64BIT__)
#       define GMS_DLL_SUFFIX "64"
#     else
#       define GMS_DLL_SUFFIX ""
#     endif
#   elif defined(__sparcv9)
#     define GMS_DLL_SUFFIX "64"
#   elif defined(__sparc)
/*    check __sparc after __sparcv9, both are defined for 64-bit */
#     define GMS_DLL_SUFFIX ""
#   endif
# endif /* ! defined(GMS_DLL_SUFFIX) */
#endif

#if ! defined(GMS_DLL_PREFIX)
# error "GMS_DLL_PREFIX expected but not defined"
#endif
#if ! defined(GMS_DLL_EXTENSION)
# error "GMS_DLL_EXTENSION expected but not defined"
#endif
#if ! defined(GMS_DLL_SUFFIX)
# error "GMS_DLL_SUFFIX expected but not defined"
#endif
{
  char *p;
  const char *s1;
  int n;

  p = (char *) result+1;
  s1 = (const char *) s+1;
  *p = '\0';
  strncpy(p, GMS_DLL_PREFIX   , 255);
  n = 255-strlen(p);
  strncat(p, s1               , (*s > n) ? n : *s);
  strncat(p, GMS_DLL_SUFFIX   , 255-strlen(p));
  strncat(p, GMS_DLL_EXTENSION, 255-strlen(p));
  result[0] = strlen(p);
}
  return result;
}  /* gamslibnamep3 */

/* unit gmslibname */
void _Init_Module_gmslibname(void)
{
} /* _Init_Module_gmslibname */

void _Final_Module_gmslibname(void)
{
} /* _Final_Module_gmslibname */

