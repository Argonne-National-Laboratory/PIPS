#include "p3io.h"
#include "system_p3.h"
#include "p3private.h"
#include "p3platform.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "math_p3.h"
#include "p3library.h"
#include "p3utils.h"

/**** C code included from p3utils.pas(177:1): 42 lines ****/
#if defined(_WIN32)
# include <winsock2.h>
# include <io.h>
# include <windows.h>

typedef BOOL (WINAPI * GetFileSizeEx_t) (HANDLE h, PLARGE_INTEGER fileSize);
GetFileSizeEx_t pGetFileSizeEx = NULL;
int triedGetFileSizeEx = 0;

typedef BOOL (WINAPI * SetFilePointerEx_t)
    (HANDLE h, LARGE_INTEGER distance,
     PLARGE_INTEGER newPointer, DWORD whence);
SetFilePointerEx_t pSetFilePointerEx = NULL;
int triedSetFilePointerEx = 0;

# if ! defined(_WIN64)

WINBASEAPI BOOL WINAPI
GetFileSizeEx (HANDLE h, PLARGE_INTEGER fileSize);
WINBASEAPI BOOL WINAPI
SetFilePointerEx (HANDLE h, LARGE_INTEGER distance,
                  PLARGE_INTEGER newPointer, DWORD whence);
# endif /* ! defined(_WIN64) */
#else
# include <fcntl.h>
# include <sys/types.h>
# include <sys/stat.h>

/* these next for socket commo */
# include <sys/socket.h>
# include <netinet/in.h>
# if (defined(__linux__) || defined(__APPLE__) || defined(__HOS_AIX__) || defined(__sparc) || defined(__sun__)) /* at least, maybe for others too */
#  include <arpa/inet.h>
# endif
# include <netdb.h>

# include <sys/utsname.h>
# include <pwd.h>

#endif

#include <locale.h>

Function(SYSTEM_integer ) P3UTILS_p3chmod(
  const SYSTEM_ansichar *path,
  SYSTEM_integer mode)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(239:1): 13 lines ****/
#if defined(_WIN32)
result = 0;
#else
{
  char filename[256];
  int len;
  /* */
  len = path[0];
  memcpy(filename, path + 1, len);
  filename[len] = '\0';
  result = chmod(filename, mode);
}
#endif
  return result;
}  /* p3chmod */

Function(SYSTEM_double ) P3UTILS_realtrunc(
  SYSTEM_double x)
{
  SYSTEM_double result;

  result = SYSTEM_int(x);
  return result;
}  /* realtrunc */

Function(SYSTEM_double ) P3UTILS_realround(
  SYSTEM_double x)
{
  SYSTEM_double result;

  if (x >= 0) { 
    result = SYSTEM_int(x + 0.5);
  } else 
    result = SYSTEM_int(x - 0.5);
  return result;
}  /* realround */

Function(SYSTEM_ansichar *) P3UTILS_floattoe(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double y,
  SYSTEM_integer decimals)
{
  SYSTEM_integer e, i, j, k, n;
  SYSTEM_double x;
  SYSTEM_shortstring s;

  x = SYSTEM_abs_r(y);
  n = 0;
  if (x != 0) {
    while (x >= 10.0) {
      n = n + 1;
      x = x /  10.0;
    
}
    while (x < 1.0) {
      n = n - 1;
      x = x * 10.0;
    
}
    x = MATH_P3_roundto(x,decimals);
    x = x * MATH_P3_intpower(10.0,n);
  } 
  _P3str_d0(x,s,255);
  k = SYSUTILS_P3_lastdelimiter(_P3str1("\002+-"),s);
  j = SYSTEM_pos(_P3str1("\001."),s);
  if (k - j - 2 < decimals) 
    decimals = k - j - 2;
  _P3strcpy(result,_len_ret,_P3str1("\002  "));
  if (y < 0) 
    result[2] = _P3char('-');
  {
    SYSTEM_shortstring _t1;
    _P3STR_255 _t2;
    _P3STR_255 _t3;
    _P3STR_3 _t4;

    _P3strcat(result,_len_ret,_P3strcat(_t3,255,_P3strcat(_t2,255,
      result,SYSTEM_copy(_t1,255,s,j - 1,decimals + 2)),_P3str1("\001E")),
      _P3ch2str(_t4,1,s[k]));
  }
  {
    SYSTEM_shortstring _t1;

    _P3val_i(SYSTEM_copy(_t1,255,s,k,5),e,&i);
  }
  if (e > 99) { 
    {
      SYSTEM_shortstring _t1;

      _P3strcat(result,_len_ret,result,SYSUTILS_P3_inttostr(_t1,255,
        e));
    }
  } else 
    {
      SYSTEM_shortstring _t1;

      _P3strcat(result,_len_ret,result,SYSTEM_copy(_t1,255,s,ValueCast(
        SYSTEM_int32,SYSTEM_length(s)) - 1,2));
    }
  return result;
}  /* floattoe */

Function(SYSTEM_ansichar *) P3UTILS_replacefileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *extension)
{
  SYSUTILS_P3_changefileext(result,_len_ret,filename,extension);
  return result;
}  /* replacefileext */

Function(SYSTEM_ansichar *) P3UTILS_completefileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *extension)
{
  {
    SYSTEM_shortstring _t1;

    if (_P3strcmpE(SYSUTILS_P3_extractfileext(_t1,255,filename),_P3str1("\000"))) { 
      SYSUTILS_P3_changefileext(result,_len_ret,filename,extension);
    } else 
      _P3strcpy(result,_len_ret,filename);
  }
  return result;
}  /* completefileext */

Function(SYSTEM_ansichar *) P3UTILS_paramstrzero(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  SYSTEM_P3_paramstr(result,_len_ret,0);
  return result;
}  /* paramstrzero */
/**** C code included from p3utils.pas(462:1): 101 lines ****/

#define _EINVALIDCAST_RAISE_E(msg) _P3_RAISE(ValueCast(EXCEPTIONS_einvalidcast,SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_einvalidcast_CD)),msg)));

#if ! defined(_WIN32)

# if ! (defined(DAR) || defined(DEG) || defined(DEI) || defined(DIG) || defined(DII))
/* all of this grabbed pretty much verbatim from p3io.c */
extern char **environ;

/* hackEnvSPD
 * routine to remove envName from the environment
 * Oh! this is ugly, but putenv does not do it "right" for:
 *   AIX
 * should not be necessary for:
 *   Solaris, Linux
 */
static int
hackEnvSPD (const char *envName)
{
  char **p, **x;
  int envNameLen;

  if (NULL == envName) {
    return 1;
  }
  envNameLen = strlen(envName);
  if (0 == envNameLen) {
    return 2;
  }
  for (p = environ;  *p != NULL;  p++) {
    if (0 == strncmp (envName, *p, envNameLen) &&
        ('\0' == (*p)[envNameLen] || '=' == (*p)[envNameLen])
        ) {
      free(*p);
      x = p;
      do {
        *x = *(x+1);
        x++;
      } while (*x);
    }
  }

  return 0;
} /* hackEnvSPD */
# endif /* ! (defined(DAR) || defined(DIG) || defined(DII)) */


/* is name allowed to be NULL or empty?  Not sure what happens if
 * that is the case    SPD, August 2001
 */
static SYSTEM_boolean
setEnvironmentVariable (SYSTEM_char *name,
                        SYSTEM_char *value)
{
  int rc;

  if (NULL == name || '\0' == *name) {
    return SYSTEM_false;
  }
  if (NULL == value) {          /* delete name from the environment */
#if defined(BGP) || defined(DAR) || defined(DEG) || defined(DEI) || defined(DIG) || defined(DII) || defined(LNX) || defined(LEG) || defined(LEI)
    /* unsetenv will do the job here */
    unsetenv ((char *)name);
    rc = 0;
#elif defined(SUN)
    /* putenv will do the job here */
    rc = putenv ((char *)name);
#else
    rc = hackEnvSPD ((char *)name);
#endif
    return rc ? SYSTEM_false : SYSTEM_true;
  }

#if defined(BGP) || defined(LNX) || defined(LEG) || defined(LEI)
  rc = setenv ((char *)name, (char *)value, 1);
#else
  {
    int nameLen, valueLen;
    char *s;
    int nChars;

    nameLen = strlen ((char *)name);
    valueLen = strlen ((char *)value);
    nChars = nameLen + 1 + valueLen + 1;
    s = (char *) malloc (nChars);
    if (NULL == s) {
      return SYSTEM_false;
    }
    strcpy (s, (char *)name);
    s[nameLen] = '=';
    strcpy (s+nameLen+1, (char *)value);
    rc = putenv (s);
    /* free (s); */
    /* SPD, Sep 2002: this is a guaranteed memory leak, but we cannot
     * do anything else; putenv expects to keep the space */
  }
#endif

  return rc ? SYSTEM_false : SYSTEM_true;
}                                       /* setEnvironmentVariable */
#endif /* ! defined(_WIN32) */

Function(SYSTEM_boolean ) P3UTILS_prefixpath(
  const SYSTEM_ansichar *s)
{
  SYSTEM_boolean result;
  static _P3STR_7 cpath = {5,'P','A','T','H','\000'};
  SYSTEM_integer slen, plen;
  SYSTEM_P3_pansichar tptr;

  slen = SYSTEM_length(s);
  if (slen == 0) { 
    result = SYSTEM_true;
  } else {
    tptr = NULL;
    /**** C code included from p3utils.pas(618:1): 42 lines ****/
{
    char *p;

#if defined(_WIN32)
    plen = GetEnvironmentVariable((char *)cpath+1,NULL,0);
    p = (char *) malloc (slen + 1 + plen);
    if (NULL == p) {
      return SYSTEM_false;
    }
    memcpy(p, (char *)s+1, slen);
    if (plen > 0) {
      int tlen;
      p[slen] = SYSUTILS_P3_pathsep;
      tlen = GetEnvironmentVariable((char *)cpath+1,p+slen+1,plen);
      assert(tlen == plen -1);
    }
    else {
      p[slen] = '\0';
    }
    result = SetEnvironmentVariable((char *)cpath+1,p);
    free(p);
#else
    tptr = (SYSTEM_char *)getenv((char *)cpath+1);
    plen = 0;
    if (NULL != tptr)
      plen = strlen((char *)tptr);
    p = (char *) malloc (slen + 1 + plen + 1);
    if (NULL == p) {
      return SYSTEM_false;
    }
    memcpy(p,(char *)s+1,slen);
    if (plen > 0) {
      p[slen] = SYSUTILS_P3_pathsep;
      memcpy(p+slen+1,(char *)tptr,plen);
      p[slen+1+plen] = '\0';
    }
    else
      p[slen] = '\0';
    result = setEnvironmentVariable(cpath+1,(SYSTEM_char *)p);
    free(p);
#endif /* #if defined(_WIN32) .. else .. */
}
  } 
  return result;
}  /* prefixpath */

Function(SYSTEM_boolean ) P3UTILS_prefixloadpath(
  const SYSTEM_ansichar *dir)
{
  SYSTEM_boolean result;
  SYSTEM_integer slen, plen;
  SYSTEM_P3_pansichar tptr;
  SYSTEM_shortstring s;
  SYSTEM_shortstring ldpath;

  result = SYSTEM_false;
  slen = SYSTEM_length(dir);
  if (slen == 0) {
    {
      SYSTEM_shortstring _t2;
      SYSTEM_shortstring _t3;

      SYSUTILS_P3_excludetrailingpathdelimiter(s,255,
        SYSUTILS_P3_extractfilepath(_t2,255,P3UTILS_paramstrzero(
        _t3,255)));
    }
    slen = SYSTEM_length(s);
  } else 
    _P3strcpy(s,255,dir);
  if (P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin) {
    result = SYSTEM_false;
    return result;
  } 
  tptr = NULL;
  switch (P3PLATFORM_osplatform()) {
    case P3PLATFORM_osaix: 
      _P3strcpy(ldpath,255,_P3str1("\007LIBPATH"));
      break;
    case P3PLATFORM_osdarwin_i386: 
    case P3PLATFORM_osdarwin: 
      _P3strcpy(ldpath,255,_P3str1("\021DYLD_LIBRARY_PATH"));
      break;
    case P3PLATFORM_oslinux: 
    case P3PLATFORM_oslinux86_64: 
    case P3PLATFORM_ossunos_sparc32: 
    case P3PLATFORM_ossunos_sparc64: 
      _P3strcpy(ldpath,255,_P3str1("\017LD_LIBRARY_PATH"));
      break;
    default:
      result = SYSTEM_false;
      return result;
  }
  _P3strcat(ldpath,255,ldpath,_P3str1("\001\000"));
  /**** C code included from p3utils.pas(733:1): 24 lines ****/
#if ! defined (_WIN32)
  {
    char *p;

    tptr = (SYSTEM_char *)getenv((char *)ldpath+1);
    plen = 0;
    if (NULL != tptr)
      plen = strlen((char *)tptr);
    p = (char *) malloc (slen + 1 + plen + 1);
    if (NULL == p) {
      return SYSTEM_false;
    }
    memcpy(p,(char *)s+1,slen);
    if (plen > 0) {
      p[slen] = SYSUTILS_P3_pathsep;
      memcpy(p+slen+1,(char *)tptr,plen);
      p[slen+1+plen] = '\0';
    }
    else
      p[slen] = '\0';
    result = setEnvironmentVariable(ldpath+1,(SYSTEM_char *)p);
    free(p);
  }
#endif /* if ! defined (_WIN32) */
  return result;
}  /* prefixloadpath */

Function(SYSTEM_boolean ) P3UTILS_p3setenv(
  const SYSTEM_ansichar *name,
  const SYSTEM_ansichar *val)
{
  SYSTEM_boolean result;
  P3PRIVATE_shortstrbuf namebuf, valbuf;
  SYSTEM_P3_pansichar nameptr, valptr;

  nameptr = P3PRIVATE_strtostrbuf(name,namebuf);
  valptr = P3PRIVATE_strtostrbuf(val,valbuf);
  /**** C code included from p3utils.pas(778:1): 5 lines ****/
#if defined(_WIN32)
  result = SetEnvironmentVariable((char *)nameptr,(char *)valptr);
#else
  result = setEnvironmentVariable(nameptr,valptr);
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3setenv */

Procedure P3UTILS_p3unsetenv(
  const SYSTEM_ansichar *name)
{
  P3PRIVATE_shortstrbuf namebuf;
  SYSTEM_P3_pansichar nameptr;

  nameptr = P3PRIVATE_strtostrbuf(name,namebuf);
  /**** C code included from p3utils.pas(803:1): 5 lines ****/
#if defined(_WIN32)
  (void) SetEnvironmentVariable((char *)nameptr,NULL);
#else
  (void) setEnvironmentVariable(nameptr,NULL);
#endif /* if defined(_WIN32) .. else .. */
}  /* p3unsetenv */

Function(SYSTEM_boolean ) P3UTILS_p3issetenv(
  const SYSTEM_ansichar *name)
{
  SYSTEM_boolean result;
  P3PRIVATE_shortstrbuf namebuf;
  SYSTEM_P3_pansichar nameptr;

  nameptr = P3PRIVATE_strtostrbuf(name,namebuf);
  /**** C code included from p3utils.pas(828:1): 5 lines ****/
#if defined(_WIN32)
  result = (GetEnvironmentVariable((char *)nameptr,NULL,0) != 0);
#else
  result = (getenv((char *)nameptr) != NULL);
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3issetenv */

Procedure P3UTILS_p3setconsoletitle(
  const SYSTEM_ansichar *s)
{
  P3PRIVATE_shortstrbuf namebuf;
  SYSTEM_P3_pansichar nameptr;

  nameptr = P3PRIVATE_strtostrbuf(s,namebuf);
  /**** C code included from p3utils.pas(849:1): 5 lines ****/
#if defined(_WIN32)
  SetConsoleTitle((char *)nameptr);
#else
  /* do nothing for now */
#endif /* if defined(_WIN32) .. else .. */
}  /* p3setconsoletitle */

Procedure P3UTILS_p3nopopups(void)
{
  /**** C code included from p3utils.pas(869:1): 5 lines ****/
#if defined(_WIN32)
  SetErrorMode(SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX);
#else
  /* do nothing */
#endif /* if defined(_WIN32) .. else .. */
}  /* p3nopopups */
/**** C code included from p3utils.pas(880:1): 62 lines ****/
#if defined(_WIN32)
/* map the Windows codes returned by GetLastError to libc codes
 * use when making Windows API calls with P3, since P3 ioresult-ish codes
 * are expected to be libc codes on all platforms
 */
static int
win2c (int rc)
{
  int result;

  switch (rc) {
    case ERROR_FILE_NOT_FOUND:
    case ERROR_PATH_NOT_FOUND:
      result = ENOENT;
      break;
    case ERROR_TOO_MANY_OPEN_FILES:
      result = EMFILE;
      break;
    case ERROR_ACCESS_DENIED:
      result = EACCES;
      break;
    case ERROR_INVALID_HANDLE:
      result = EBADF;
      break;
    case ERROR_NOT_ENOUGH_MEMORY:
      result = ENOMEM;
      break;
    case ERROR_INVALID_ACCESS:
      result = EACCES;
      break;
    case ERROR_NO_MORE_FILES:
      result = ENFILE;
      break;
    case ERROR_SEEK_ON_DEVICE:
      result = ESPIPE;
      break;
    case ERROR_INVALID_PARAMETER:
    case ERROR_NEGATIVE_SEEK:
      result = EINVAL;
      break;
    default:
      result = 0; /* no guessing */
  }              /* case */
  return result;
} /* win2c */

static const DWORD accessMode[3] = {
  GENERIC_READ,
  GENERIC_WRITE,
  GENERIC_READ | GENERIC_WRITE };
/* this works for GDX so we do it: it is kind of silly to use the
 * shareMode var then but why not? */
static const DWORD shareMode[3] = {
  FILE_SHARE_READ | FILE_SHARE_WRITE,
  FILE_SHARE_READ | FILE_SHARE_WRITE,
  FILE_SHARE_READ | FILE_SHARE_WRITE };
static const DWORD createHow[3] = {
  OPEN_EXISTING,
  CREATE_ALWAYS,
  OPEN_ALWAYS };

#endif /* if defined(_WIN32) */

static Function(SYSTEM_boolean ) P3UTILS_isvalidhandle(
  P3UTILS_tp3filehandle h)
{
  SYSTEM_boolean result;

  /**** C code included from p3utils.pas(955:1): 5 lines ****/
#if defined(_WIN32)
  result = h && (INVALID_HANDLE_VALUE != (HANDLE) h);
#else
  result = (SYSTEM_int64)h > 0;
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* isvalidhandle */

Function(SYSTEM_integer ) P3UTILS_p3fileopen(
  const SYSTEM_ansichar *fname,
  P3UTILS_tp3fileopenaction mode,
  P3UTILS_tp3filehandle *h)
{
  SYSTEM_integer result;
  P3PRIVATE_shortstrbuf namebuf;
  SYSTEM_P3_pansichar nameptr;

  nameptr = P3PRIVATE_strtostrbuf(fname,namebuf);
  /**** C code included from p3utils.pas(1072:1): 81 lines ****/
#if defined(_WIN32)
{
  DWORD lowMode;
  HANDLE hFile;
  int rc;

  lowMode = mode & 3;
  if (3 == lowMode) {
    *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
    return ERROR_INVALID_PARAMETER;
  }

  if (nameptr[0] == 0) {
     if (P3UTILS_p3openread == mode)
        hFile = GetStdHandle(STD_INPUT_HANDLE);
     else if (P3UTILS_p3openwrite == mode)
        hFile = GetStdHandle(STD_OUTPUT_HANDLE);
     else {
        *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
        return ERROR_INVALID_PARAMETER;
     }
  }
  else
    hFile = CreateFile ((char *)nameptr, accessMode[lowMode], shareMode[lowMode], NULL,
                      createHow[lowMode], FILE_ATTRIBUTE_NORMAL, NULL);
  if (INVALID_HANDLE_VALUE == hFile) {
    *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
    rc = GetLastError();
    result = win2c(rc);
    if (0 == result) { /* ouch: just pick a likely but non-specific code */
      result = EACCES;
    }
  }
  else {
    *h = (P3UTILS_tp3filehandle) hFile;
    result = 0;
  }
}
#else
{
  struct stat statBuf;
  int fd, flags, rc;

  if (nameptr[0] == 0) {
     if (P3UTILS_p3openread == mode)
        *h = (P3UTILS_tp3filehandle) STDIN_FILENO;
     else if (P3UTILS_p3openwrite == mode)
        *h = (P3UTILS_tp3filehandle) STDOUT_FILENO;
     else {
        *h = (P3UTILS_tp3filehandle) 0;
        return -1;
     }
     return 0;
  }
  flags = mode & 3;
  if (flags > 0)            /* write-only or read-write */
    flags |= O_CREAT;
  if (flags & 1)
    flags |= O_TRUNC;
  fd = open ((char *)nameptr, flags, 0666);
  if (-1 == fd) {
    *h = (P3UTILS_tp3filehandle) 0;
    return errno;
  }
  result = 0;
  /* before calling this a success, check for directory on read-only */
  if (P3UTILS_p3openread == mode) {
    rc = fstat (fd, &statBuf);
    if (rc)
      result = errno;
    else if (S_ISDIR(statBuf.st_mode))
      result = EISDIR;
  }
  if (result) {
    close(fd);
    return result;
  }

  *h = (P3UTILS_tp3filehandle) (SYSTEM_nativeuint) fd;
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3fileopen */

Function(SYSTEM_integer ) P3UTILS_p3fileclose(
  P3UTILS_tp3filehandle *h)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(1186:1): 28 lines ****/
if ( ! P3UTILS_isvalidhandle (h)) {
  return EBADF;
}
#if defined(_WIN32)
{
  if (CloseHandle((HANDLE)*h)) {
    result = 0; /* success */
  }
  else {
    int rc = GetLastError();
    result = win2c(rc);
    if (0 == result) { /* ouch: just pick a likely but non-specific code */
      result = EIO;
    }
  }
  *h = (P3UTILS_tp3filehandle) INVALID_HANDLE_VALUE;
}
#else
{
  SYSTEM_int64 h64;

  h64 = (SYSTEM_int64)*h;
  result = 0;
  if (close((int)h64))   /* error!! */
    result = errno;
  *h = (P3UTILS_tp3filehandle) 0;
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3fileclose */

Function(SYSTEM_integer ) P3UTILS_p3fileread(
  P3UTILS_tp3filehandle h,
  SYSTEM_untyped *buffer,
  SYSTEM_longword buflen,
  SYSTEM_longword *numread)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(1246:1): 31 lines ****/
#if defined(_WIN32)
{
  if (ReadFile((HANDLE)h, buffer, buflen, (LPDWORD)numread, NULL)) {
    result = 0; /* success */
  }
  else {
    int rc = GetLastError();
    result = win2c(rc);
    if (0 == result) { /* ouch: just pick a likely but non-specific code */
      result = EIO;
    }
  }
}
#else
{
  int rc;
  SYSTEM_int64 h64;

  h64 = (SYSTEM_int64)h;
  result = 0;
  rc = read((int)h64, buffer, buflen);
  if (rc < 0) {
    result = errno;
    *numread = 0;
  }
  else {
    result = 0;
    *numread = (SYSTEM_longword) rc;
  }
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3fileread */

Function(SYSTEM_integer ) P3UTILS_p3filewrite(
  P3UTILS_tp3filehandle h,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword buflen,
  SYSTEM_longword *numwritten)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(1309:1): 31 lines ****/
#if defined(_WIN32)
{
  if (WriteFile((HANDLE)h, buffer, buflen, (LPDWORD)numwritten, NULL)) {
    result = 0; /* success */
  }
  else {
    int rc = GetLastError();
    result = win2c(rc);
    if (0 == result) { /* ouch: just pick a likely but non-specific code */
      result = EIO;
    }
  }
}
#else
{
  int rc;
  SYSTEM_int64 h64;

  h64 = (SYSTEM_int64)h;
  result = 0;
  rc = write((int)h64, buffer, buflen);
  if (rc < 0) {
    result = errno;
    *numwritten = 0;
  }
  else {
    result = 0;
    *numwritten = (SYSTEM_longword) rc;
  }
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3filewrite */

Function(SYSTEM_integer ) P3UTILS_p3filegetsize(
  P3UTILS_tp3filehandle h,
  SYSTEM_int64 *filesize)
{
  SYSTEM_integer result;

  *filesize =  -1;
  /**** C code included from p3utils.pas(1386:1): 49 lines ****/
if ( ! P3UTILS_isvalidhandle (h)) {
  return EBADF;
}
#if defined(_WIN32)
{
  BOOL frc;

  if (! triedGetFileSizeEx) {
    pGetFileSizeEx = (GetFileSizeEx_t) GetProcAddress(
      GetModuleHandle("kernel32"),"GetFileSizeEx");
    triedGetFileSizeEx = 1;
  }
  if (pGetFileSizeEx) {
    frc = pGetFileSizeEx((HANDLE)h, (PLARGE_INTEGER) filesize);
  }
  else {
    DWORD tt;

    frc = GetFileSize((HANDLE)h, &tt);
    *filesize = tt;
  }

  if (frc)
    result = 0;
  else {
    int rc = GetLastError();
    result = win2c(rc);
    if (0 == result) { /* ouch: just pick a likely but non-specific code */
      result = EACCES;
    }
  }
}
#else
{
  struct stat statBuf;
  int rc;
  SYSTEM_int64 h64;

  h64 = (SYSTEM_int64)h;
  rc = fstat ((int)h64, &statBuf);
  if (rc) {
    result = errno;
  }
  else {
    result = 0;
    *filesize = statBuf.st_size;
  }
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3filegetsize */

Function(SYSTEM_integer ) P3UTILS_p3filesetpointer(
  P3UTILS_tp3filehandle h,
  SYSTEM_int64 distance,
  SYSTEM_int64 *newpointer,
  SYSTEM_longword whence)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(1491:1): 75 lines ****/
if ( ! P3UTILS_isvalidhandle (h)) {
  return EBADF;
}
#if defined(_WIN32)
{
  int rc;
  BOOL frc;
  LARGE_INTEGER d; /* declared as a union - compiler rejects a cast */

  if (! triedSetFilePointerEx) {
    pSetFilePointerEx = (SetFilePointerEx_t) GetProcAddress(
      GetModuleHandle("kernel32"),"SetFilePointerEx");
    triedSetFilePointerEx = 1;
  }
  d.QuadPart = distance;
  if (pSetFilePointerEx) {
    frc = pSetFilePointerEx((HANDLE)h, d, (PLARGE_INTEGER) newpointer, whence);
  }
  else {
    DWORD tt, trc;;

    trc = SetFilePointer((HANDLE)h, (int)distance, NULL, whence);
    if (INVALID_SET_FILE_POINTER == trc) {
      frc = 0;
      *newpointer = 0;
    }
    else {
      *newpointer = trc;
      frc = 1;
    }
  }
  if (frc)
    result = 0;
  else {
    rc = GetLastError();
    result = win2c(rc);
    if (0 == result) {
      result = EINVAL;
    }
  }
}
#else
{
  int w = -1;
  off_t newPos, offset;
  SYSTEM_int64 h64;

  switch (whence) {
    case P3UTILS_p3_file_begin:
      w = SEEK_SET;
      break;
    case P3UTILS_p3_file_current:
      w = SEEK_CUR;
      break;
    case P3UTILS_p3_file_end:
      w = SEEK_END;
      break;
    default:
      return EINVAL;
  } /* switch */

  /* check if conversion to off_t loses info */
  offset = (off_t)distance;
#if ! defined(__alpha) /* no overflow possible on OSF, so not defined */
  if (offset != distance)
    return EOVERFLOW;
#endif
  h64 = (SYSTEM_int64)h;
  newPos = lseek ((int)h64, offset, w);
  if ((off_t)-1 == newPos)
    return errno;
  *newpointer = newPos;
  result = 0;
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3filesetpointer */

Function(SYSTEM_integer ) P3UTILS_p3filegetpointer(
  P3UTILS_tp3filehandle h,
  SYSTEM_int64 *filepointer)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(1603:1): 56 lines ****/
if ( ! P3UTILS_isvalidhandle (h)) {
  return EBADF;
}
#if defined(_WIN32)
{
  int rc;
  BOOL frc;
  LARGE_INTEGER d; /* declared as a union - compiler rejects a cast */

  if (! triedSetFilePointerEx) {
    pSetFilePointerEx = (SetFilePointerEx_t) GetProcAddress(
      GetModuleHandle("kernel32"),"SetFilePointerEx");
    triedSetFilePointerEx = 1;
  }
  d.QuadPart = 0;
  if (pSetFilePointerEx) {
    frc = pSetFilePointerEx((HANDLE)h, d, (PLARGE_INTEGER) filepointer,
                            P3UTILS_p3_file_current);
  }
  else {
    DWORD tt, trc;;

    trc = SetFilePointer((HANDLE)h, 0, NULL, P3UTILS_p3_file_current);
    if (INVALID_SET_FILE_POINTER == trc) {
      frc = 0;
      *filepointer = 0;
    }
    else {
      *filepointer = trc;
      frc = 1;
    }
  }

  if (frc)
    result = 0;
  else {
    rc = GetLastError();
    result = win2c(rc);
    if (0 == result) {
      result = EINVAL;
    }
  }
}
#else
{
  off_t newPos;
  SYSTEM_int64 h64;

  h64 = (SYSTEM_int64)h;
  newPos = lseek ((int)h64, 0, SEEK_CUR);
  if ((off_t)-1 == newPos)
    return errno;
  *filepointer = newPos;
  result = 0;
}
#endif /* if defined(_WIN32) .. else .. */
  return result;
}  /* p3filegetpointer */

Function(SYSTEM_int64 ) P3UTILS_p3allocmemsize64(void)
{
  SYSTEM_int64 result;

  /**** C code included from p3utils.pas(1674:1): 1 lines ****/
  result = SYSTEM_allocmemsize64;
  return result;
}  /* p3allocmemsize64 */

Function(SYSTEM_pointer ) P3UTILS_p3allocmem64(
  SYSTEM_int64 size)
{
  SYSTEM_pointer result;

  result = NULL;
  if (size <= 0) 
    return result;
  /**** C code included from p3utils.pas(1695:1): 13 lines ****/
  if (8 == sizeof(result)) {
    _P3_new64 (&result, size);
    (void) memset(result, 0, (size_t) size);
  }
  else {
    if (size > SYSTEM_maxint) {
      _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
    }
    else {
      _P3_new (&result, (SYSTEM_longint) size);
      (void) memset(result, 0, (size_t) size);
    }
  }
  return result;
}  /* p3allocmem64 */

Procedure P3UTILS_p3fillchar64(
  SYSTEM_untyped *p,
  SYSTEM_int64 size,
  SYSTEM_byte fillvalue)
{
  /**** C code included from p3utils.pas(1740:1): 11 lines ****/
  if (size <= 0)
    return;
  if (8 == sizeof(p)) {
    (void) memset(p, (int) fillvalue, (size_t) size);
  }
  else {
    size_t sz = (size_t) size;
    if ((SYSTEM_int64)sz != size)
      _EINVALIDCAST_RAISE_E(_P3str1("\034p3FillChar64: size too large"));
    (void) memset(p, (int) fillvalue, sz);
  }
}  /* p3fillchar64 */

Procedure P3UTILS_p3freemem64(
  SYSTEM_pointer *p,
  SYSTEM_int64 size)
{
  /**** C code included from p3utils.pas(1771:1): 11 lines ****/
  if (8 == sizeof(p)) {
    _P3_free64 (*p, size);
  }
  else {
    if (size < 0 || size > SYSTEM_maxint) {
      _P3_free (*p, 0);
    }
    else {
      _P3_free (*p, (SYSTEM_longint) size);
    }
  }
}  /* p3freemem64 */

Procedure P3UTILS_p3getmem64(
  SYSTEM_pointer *p,
  SYSTEM_int64 size)
{
  if (size <= 0) {
    *p = NULL;
    return;
  } 
  /**** C code included from p3utils.pas(1803:1): 11 lines ****/
  if (8 == sizeof( *p )) {
    _P3_new64 (p, size);
  }
  else {
    if (size > SYSTEM_maxint) {
      _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
    }
    else {
      _P3_new (p, (SYSTEM_longint) size);
    }
  }
}  /* p3getmem64 */

Procedure P3UTILS_p3move64(
  const SYSTEM_untyped *source,
  SYSTEM_untyped *dest,
  SYSTEM_int64 sz)
{
  /**** C code included from p3utils.pas(1833:1): 9 lines ****/
  if (8 == sizeof( source )) {
    (void) memmove (dest, source, (size_t) sz);
  }
  else {
    size_t ssz = (size_t) sz;
    if (ssz != sz)
      _EINVALIDCAST_RAISE_E(_P3str1("\030p3Move64: size too large"));
    (void) memmove (dest, source, ssz);
  }
}  /* p3move64 */

Procedure P3UTILS_p3reallocmem64(
  SYSTEM_pointer *p,
  SYSTEM_int64 size)
{
  if (size <= 0) 
    if (*p == NULL) { 
      return;
    } else {
      _P3freemem0(*p);
      *p = NULL;
    } 
  /**** C code included from p3utils.pas(1869:1): 11 lines ****/
  if (8 == sizeof(p)) {
    SYSTEM_reallocmem64 (p, size);
  }
  else {
    if (size > SYSTEM_maxint) {
      _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
    }
    else {
      SYSTEM_reallocmem (p, (SYSTEM_longint)size);
    }
  }
}  /* p3reallocmem64 */

Procedure P3UTILS_p3getfromurl(
  const SYSTEM_ansichar *servername,
  const SYSTEM_ansichar *filename,
  SYSTEM_word port,
  P3UTILS_thavedatacb havedata,
  SYSTEM_pointer usermem,
  SYSTEM_ansichar *msg)
{
  cnstdef {maxbuf = 4096};
  typedef SYSTEM_uint16 _sub_1P3GETFROMURL;
  typedef SYSTEM_ansichar _arr_0P3GETFROMURL[4096];
  _arr_0P3GETFROMURL buffer;
  SYSTEM_integer len;

  _P3strclr(msg);
  _P3strcpy(msg,255,_P3str1("\041Not implemented for this platform"));
  /**** C code included from p3utils.pas(1988:1): 114 lines ****/
#if defined(_WIN32)
#define WSACLEANUP   WSACancelBlockingCall(); WSACleanup(); return
{
  WSADATA wsaData;
  SOCKET conn;
  struct hostent *host;
  char sbuf[300];
  unsigned int addr;
  struct sockaddr_in server;

  rc = WSAStartup (0x101, &wsaData);
  if (rc) {
    strcpy ((char *)msg, "\026Winsock startup failed");
    return;
  }
  conn = socket (AF_INET, SOCK_STREAM, IPPROTO_TCP);
  if (conn < 0) {
    strcpy ((char *)msg, "\017No socket found");
    WSACLEANUP;
  }
  len = servername[0];
  strncpy (sbuf, (char *)servername+1, len);
  sbuf[len] = '\0';
  addr = inet_addr(sbuf);
  if (INADDR_NONE == addr)
    host = gethostbyname (sbuf);
  else
    host = gethostbyaddr ((char *) &addr, sizeof(addr), PF_INET);
  if (NULL == host) {
    closesocket (conn);
    strcpy ((char *)msg, "\017Unknown address");
    WSACLEANUP;
  }
  memcpy (&server.sin_addr.s_addr, host->h_addr, host->h_length);
  server.sin_family=AF_INET;
  server.sin_port=htons(port);
  if (connect (conn, (struct sockaddr *) &server, sizeof(server))) {
    closesocket(conn);
    strcpy ((char *)msg, "\021Could not connect");
    WSACLEANUP;
  }
  sprintf(sbuf,"GET /%.*s HTTP/1.0\r\n\r\n", *filename, filename+1);
  len = send (conn, sbuf, sizeof(sbuf), 0);
  for ( ; ; ) {
    len = recv (conn, (char *)buffer, sizeof(buffer), 0);
    if (len < 0) {
      closesocket(conn);
      strcpy ((char *)msg, "\024Error receiving data");
      WSACLEANUP;
    }
    if (0 == len)
      break;
    if (! havedata(buffer, len, usermem)) {
      strcpy ((char *)msg, "\023Stopped by callback");
      WSACLEANUP;
    }
  }  /* receive loop */
  closesocket(conn);
  *msg = '\0';
}
#else
#if ! defined(INADDR_NONE)
# define INADDR_NONE 0xffffffff
#endif
{
  struct hostent *host;
  in_addr_t addr;
  struct sockaddr_in server;
  int conn;
  char sbuf[300];

  conn = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
  if (conn < 0) {
    strcpy ((char *)msg, "\017No socket found");
    return;
  }
  len = servername[0];
  strncpy (sbuf, (char *) servername+1, len);
  sbuf[len] = '\0';
  addr = inet_addr(sbuf);
  if (INADDR_NONE == addr)
    host = gethostbyname (sbuf);
  else
    host = gethostbyaddr ((void *)&addr, sizeof(addr), PF_INET);
  if (NULL == host) {
    strcpy ((char *)msg, "\017Unknown address");
    return;
  }
  memcpy (&server.sin_addr.s_addr, host->h_addr_list[0], host->h_length);
  server.sin_family=AF_INET;
  server.sin_port=htons(port);
  if (connect (conn, (struct sockaddr *) &server, sizeof(server))) {
    strcpy ((char *)msg, "\021Could not connect");
    return;
  }
  sprintf(sbuf,"GET /%.*s HTTP/1.0\r\n\r\n", *filename, filename+1);
  len = send (conn, sbuf, sizeof(sbuf), 0);
  for ( ; ; ) {
    len = recv (conn, buffer, sizeof(buffer), 0);
    if (len < 0) {
      strcpy ((char *)msg, "\024Error receiving data");
      return;
    }
    if (0 == len)
      break;
    if (! havedata(buffer, len, usermem)) {
      strcpy ((char *)msg, "\023Stopped by callback");
      return;
    }
  }  /* receive loop */
  (void) shutdown (conn, SHUT_RDWR);
  *msg = '\0';
}
#endif /* if defined(_WIN32) .. else .. */
}  /* p3getfromurl */

Function(SYSTEM_integer ) P3UTILS_p3getwindowsversion(void)
{
  SYSTEM_integer result;
  SYSTEM_pointer libhandle;
  SYSTEM_shortstring loadmsg;
  SYSTEM_pointer wine_get_version;

  result = P3UTILS_os_unknown;
  libhandle = P3LIBRARY_p3loadlibrary(_P3str1("\011ntdll.dll"),loadmsg);
  if (libhandle != NULL) {
    result = P3UTILS_os_winnt;
    wine_get_version = P3LIBRARY_p3getprocaddress(libhandle,_P3str1("\020wine_get_version"));
    if (wine_get_version != NULL) 
      result = P3UTILS_os_wine;
    P3LIBRARY_p3freelibrary(libhandle);
  } else {
    libhandle = P3LIBRARY_p3loadlibrary(_P3str1("\014kernel32.dll"),
      loadmsg);
    if (libhandle != NULL) {
      result = P3UTILS_os_win9x;
      wine_get_version = P3LIBRARY_p3getprocaddress(libhandle,_P3str1("\021GetDKrnl32Version"));
      if (wine_get_version != NULL) 
        result = P3UTILS_os_hx;
      P3LIBRARY_p3freelibrary(libhandle);
    } 
  } 
  return result;
}  /* p3getwindowsversion */

Function(SYSTEM_ansichar *) P3UTILS_p3getcomputername(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,_P3str1("\007unknown"));
  /**** C code included from p3utils.pas(2158:1): 27 lines ****/
#if defined(_WIN32)
   {
     char compName[256];
     DWORD n;

     n = (DWORD) sizeof(compName);
     if (GetComputerName(compName, &n)) {
       /* success, copy compName to result */
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), compName, n);
     }
   }
#else
   {
     int rc, n;
     struct utsname uts;

     rc = uname (&uts);
     if (rc >= 0) {
       n = strlen(uts.nodename);
       if (n > 255)
         n = 255;
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), uts.nodename, n);
     }
   }
#endif /* #if defined(_WIN32) .. else .. */
  return result;
}  /* p3getcomputername */

Function(SYSTEM_ansichar *) P3UTILS_p3getusername(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,_P3str1("\007unknown"));
  /**** C code included from p3utils.pas(2208:1): 44 lines ****/
#if defined(_WIN32)
   {
     char userName[256];
     DWORD n;

     n = (DWORD) sizeof(userName);
     if (GetUserName(userName, &n)) {
       /* success, copy userName to result */
       n--;  /* output n includes space for the null byte */
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), userName, n);
     }
   }
#else
   {
     int rc, n;
     char loginName[256];
     char *p = NULL;

# if (defined(__USE_XOPEN) || defined(_XOPEN_SOURCE) || defined(_POSIX_PTHREAD_SEMANTICS) || defined(__EXTENSIONS__)) && ! defined(__APPLE__)
     /* the cuserid routine is preferred */
     p = cuserid (loginName);
     if (p) {
       loginName[sizeof(loginName)-1] = '\0';
       n = strlen(loginName);
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), loginName, n);
     }
# else
#  if ((defined(SIG)||defined(SOL)) && ! defined(_POSIX_PTHREAD_SEMANTICS))
     p = getlogin_r (loginName, sizeof(loginName));
     rc = (NULL == p);
#  else
     rc = getlogin_r (loginName, sizeof(loginName));
#  endif
     if (0 == rc) { /* success */
       loginName[sizeof(loginName)-1] = '\0';
       n = strlen(loginName);
       *result = (SYSTEM_byte) n;
       (void) memcpy ((char *)(result+1), loginName, n);
     }
# endif
   }
#endif /* #if defined(_WIN32) .. else .. */
  return result;
}  /* p3getusername */

Function(SYSTEM_boolean ) P3UTILS_p3senddatamessage(
  SYSTEM_boolean broadcast,
  const SYSTEM_ansichar *_ftmp1,
  const SYSTEM_ansichar *_ftmp2)
{
  SYSTEM_shortstring wintitle;
  SYSTEM_shortstring data;
  SYSTEM_boolean result;

  _P3strcpy(wintitle,255,_ftmp1);
  _P3strcpy(data,255,_ftmp2);
  result = SYSTEM_false;
  /**** C code included from p3utils.pas(2294:1): 38 lines ****/
#if defined(_WIN32)
{
   HWND receiver;
   COPYDATASTRUCT cds;
   LRESULT r;
   char dBuf[256], tBuf[256];
   unsigned int n;

   n = data[0];
   if (n > 255)
      n = 255;
   (void) memcpy (dBuf, (char *)data+1, n);
   dBuf[n] = '\0';
   cds.dwData = 0;
   cds.lpData = dBuf;
   cds.cbData = n+1; /* include the null byte in this count */
   r = 0;
   if (broadcast) {
      r = SendMessage (HWND_BROADCAST, (UINT)WM_COPYDATA, (WPARAM)0,
                       (LPARAM)&cds);
   }
   else {
      n = wintitle[0];
      if (n > 255)
         n = 255;
      if (n > 0) {
         (void) memcpy (tBuf, (char *)wintitle+1, n);
         tBuf[n] = '\0';
         receiver = FindWindow (NULL, tBuf);
         if (receiver)
            r = SendMessage (receiver, (UINT)WM_COPYDATA, (WPARAM)0,
                             (LPARAM)&cds);
      }
   }
   result = (SYSTEM_boolean) r;
}
#else
#endif /* #if defined(_WIN32) .. else .. */
  return result;
}  /* p3senddatamessage */

Function(SYSTEM_ansichar *) P3UTILS_p3pushdeflocale(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,_P3str1("\001C"));
  /**** C code included from p3utils.pas(2341:1): 10 lines ****/
{
  char *s;
  s = setlocale(LC_NUMERIC, NULL);
  if (('C' != s[0]) || ('\0' != s[1])) { /* not what we want */
    /* so save the current locale and reset */
    (void) strcpy((char *)result+1, s);
    *result = strlen(s);
    (void) setlocale(LC_NUMERIC, "C");
  }
}
  return result;
}  /* p3pushdeflocale */

Procedure P3UTILS_p3popdeflocale(
  const SYSTEM_ansichar *s)
{
  if (_P3char('\001') == s[0] && _P3char('C') == s[1]) 
    return;
  /**** C code included from p3utils.pas(2364:1): 10 lines ****/
{
  char buf[32];
  int n = sizeof(buf)-1;

  if ( *s < n)
    n = *s;
  (void) memcpy (buf, s+1, n);
  buf[n] = '\0';
  (void) setlocale(LC_NUMERIC, buf);
}
}  /* p3popdeflocale */
/**** C code included from p3utils.pas(2382:1): 1 lines ****/
#include "p3Custom2.h"

Function(SYSTEM_integer ) P3UTILS_p3getexecname(
  SYSTEM_ansichar *execname,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;

  result = 9;
  _P3strclr(execname);
  _P3strcpy(msg,255,_P3str1("\027P3: not yet implemented"));
  /**** C code included from p3utils.pas(2394:1): 1 lines ****/
  result = xGetExecName (execname, msg);
  return result;
}  /* p3getexecname */

Function(SYSTEM_integer ) P3UTILS_p3getlibname(
  SYSTEM_ansichar *libname,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;

  if (!SYSTEM_islibrary()) {
    result = 2;
    _P3strclr(libname);
    _P3strcpy(msg,255,_P3str1("\031Not called from a library"));
    return result;
  } 
  result = 9;
  _P3strclr(libname);
  _P3strcpy(msg,255,_P3str1("\027P3: not yet implemented"));
  /**** C code included from p3utils.pas(2436:1): 1 lines ****/
  result = xGetLibName (libname, msg);
  return result;
}  /* p3getlibname */

Function(SYSTEM_integer ) P3UTILS_p3someioresult(void)
{
  SYSTEM_integer result;

  /**** C code included from p3utils.pas(2473:1): 1 lines ****/
result = EIO;
  return result;
}  /* p3someioresult */

/* unit p3utils */
void _Init_Module_p3utils(void)
{
} /* _Init_Module_p3utils */

void _Final_Module_p3utils(void)
{
} /* _Final_Module_p3utils */

