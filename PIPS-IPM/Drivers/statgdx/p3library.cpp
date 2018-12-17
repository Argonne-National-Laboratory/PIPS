#include "p3io.h"
#include "system_p3.h"
#include "p3private.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "p3platform.h"
#include "p3library.h"

/**** C code included from p3library.pas(47:1): 23 lines ****/
#if ! defined(_WIN32)
# include <dlfcn.h>
#endif

#if defined(_WIN32)
static char *
winLastErr2Buf (int errNum, char *buf, int buflen)
{
  FormatMessage(
                FORMAT_MESSAGE_FROM_SYSTEM |
                FORMAT_MESSAGE_IGNORE_INSERTS,
                NULL,
                errNum,
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
                buf,
                buflen-1,
                NULL
                );
  buf[buflen-1] = '\0';
  /* sprintf (buf, "%d\n", errNum); */
  return buf;
} /* winLastErr2Buf */
#endif

Function(P3LIBRARY_tlibhandle ) P3LIBRARY_p3loadlibrary(
  const SYSTEM_ansichar *lib,
  SYSTEM_ansichar *loadmsg)
{
  P3LIBRARY_tlibhandle result;
  SYSTEM_integer lasterr;
  P3PRIVATE_shortstrbuf libbuf;
  SYSTEM_P3_pansichar libptr;
  P3PRIVATE_shortstrbuf errmsgbuf;
  SYSTEM_P3_pansichar dlerrmsg;

  libptr = P3PRIVATE_strtostrbuf(lib,libbuf);
  dlerrmsg = NULL;
  /**** C code included from p3library.pas(120:1): 31 lines ****/
#if defined(_WIN32)
{
  UINT oldMode;

  oldMode = SetErrorMode(SEM_FAILCRITICALERRORS);
  result = (SYSTEM_pointer) LoadLibrary((char *)libptr);
  lasterr = GetLastError();
  SetErrorMode(oldMode);
  if (NULL == result) {
    if (ERROR_BAD_EXE_FORMAT == lasterr) {
      sprintf ((char *)errmsgbuf, "File is not a valid Win%s DLL",
               (8==sizeof(void *)) ? "64" : "32");
      dlerrmsg = errmsgbuf;
    }
    else
      dlerrmsg = (SYSTEM_P3_pansichar)
                 winLastErr2Buf(GetLastError(),(char *)errmsgbuf,256);
  }
}
#else
/* need RTLD_GLOBAL, so loaded library can load yet another library
 * dynamically, see #2472 */
result = (SYSTEM_pointer) dlopen((char *)libptr, RTLD_NOW
#ifndef AIX
 | RTLD_GLOBAL
#endif
);
if (NULL == result) {
  dlerrmsg = (SYSTEM_P3_pansichar) dlerror();
}
#endif /* if defined(_WIN32) .. else .. */
  if (result == NULL) { 
    if (dlerrmsg == NULL) { 
      _P3strcpy(loadmsg,255,_P3str1("\024No message available"));
    } else 
      P3PRIVATE_pchartostr(loadmsg,255,dlerrmsg);
  } else 
    _P3strclr(loadmsg);
  return result;
}  /* p3loadlibrary */

Function(SYSTEM_pointer ) P3LIBRARY_p3getprocaddress(
  P3LIBRARY_tlibhandle handle,
  const SYSTEM_ansichar *name)
{
  SYSTEM_pointer result;
  P3PRIVATE_shortstrbuf namebuf;
  SYSTEM_P3_pansichar nameptr;

  nameptr = P3PRIVATE_strtostrbuf(name,namebuf);
  /**** C code included from p3library.pas(183:1): 12 lines ****/
#if defined(_WIN32)
result = GetProcAddress((HMODULE)handle, (char *)nameptr);
#else
{
  char *errMsg;
  dlerror(); /* clear the error state, will not happen on success */
  result = dlsym(handle, (char *)nameptr);
  errMsg = dlerror();
  if (NULL != errMsg)
    result = NULL;
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3getprocaddress */

Function(SYSTEM_boolean ) P3LIBRARY_p3freelibrary(
  P3LIBRARY_tlibhandle handle)
{
  SYSTEM_boolean result;

  /**** C code included from p3library.pas(210:1): 17 lines ****/
#if defined(_WIN32)
  result = FreeLibrary((HMODULE)handle);
#else
# if defined(DAR)
  result = ! dlclose(handle);
  /* SPD 30 Mar 2005: I checked with the P3 users in an email on 25 Feb:
   * they sent no response but when asked said they do not load-unload-load
   * any DLLs so it is OK if re-initialization does not work
   * and they do not need finalization either */
  if (! result) {
    /* fprintf (stdout, "HACK HACK: dlclose return hacked in P3FreeLibrary\n");*/
    result = 1;
  }
# else
  result = ! dlclose(handle);
# endif
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3freelibrary */

Function(SYSTEM_boolean ) P3LIBRARY_p3libhandleisnil(
  P3LIBRARY_tlibhandle handle)
{
  SYSTEM_boolean result;

  result = handle == NULL;
  return result;
}  /* p3libhandleisnil */

Function(SYSTEM_ansichar *) P3LIBRARY_p3makelibname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *path,
  const SYSTEM_ansichar *base)
{
  if (_P3strcmpE(path,_P3str1("\000"))) { 
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      SYSTEM_shortstring _t3;

      _P3strcat(result,_len_ret,_P3strcat(_t2,255,
        P3LIBRARY_p3libraryprefix(_t1,255),base),
        P3LIBRARY_p3libraryext(_t3,255));
    }
  } else 
    {
      SYSTEM_shortstring _t1;
      _P3STR_3 _t2;
      _P3STR_255 _t3;
      SYSTEM_shortstring _t4;
      _P3STR_255 _t5;
      _P3STR_255 _t6;
      SYSTEM_shortstring _t7;

      _P3strcat(result,_len_ret,_P3strcat(_t6,255,_P3strcat(_t5,255,
        _P3strcat(_t3,255,SYSUTILS_P3_excludetrailingpathdelimiter(
        _t1,255,path),_P3ch2str(_t2,1,SYSUTILS_P3_pathdelim)),
        P3LIBRARY_p3libraryprefix(_t4,255)),base),
        P3LIBRARY_p3libraryext(_t7,255));
    }
  return result;
}  /* p3makelibname */

Function(P3LIBRARY_tlibhandle ) P3LIBRARY_p3nillibhandle(void)
{
  P3LIBRARY_tlibhandle result;

  result = NULL;
  return result;
}  /* p3nillibhandle */

Function(SYSTEM_ansichar *) P3LIBRARY_p3libraryext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,P3PLATFORM_osdllextension[
    P3PLATFORM_osplatform()]);
  return result;
}  /* p3libraryext */

Function(SYSTEM_ansichar *) P3LIBRARY_p3libraryprefix(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,P3PLATFORM_osdllprefix[
    P3PLATFORM_osplatform()]);
  return result;
}  /* p3libraryprefix */

/* unit p3library */
void _Init_Module_p3library(void)
{
} /* _Init_Module_p3library */

void _Final_Module_p3library(void)
{
} /* _Final_Module_p3library */

