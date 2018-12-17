#include "p3io.h"
#include "system_p3.h"
#include "exceptions.h"
#include "p3platform.h"
#include "sysutils_p3.h"

SYSTEM_ansichar SYSUTILS_P3_pathdelim, SYSUTILS_P3_drivedelim, 
  SYSUTILS_P3_pathsep;
_arr_13SYSUTILS_P3 SYSUTILS_P3_monthdays = {{31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}, 
  {31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31}};
/**** C code included from sysutils_p3.pas(289:1): 77 lines ****/
#if   defined(P3UNIX)
#include <sys/time.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>  /* for NAME_MAX, at least */
#include <fnmatch.h>
#include <dirent.h>
#include <time.h>
#elif defined(P3DOS)
#include <windows.h>
#include <io.h>
#endif

static SYSTEM_char *
mkP3Msg (const char *txt, SYSTEM_shortstring msg)
{
  int n;

  n = snprintf ((char *)msg+1, 255, "%s", txt);
  /* stupid MS _snprintf */
  if ((n > 255) || (n < 0))
    n = 255;
  *msg = (SYSTEM_byte) n;
  return msg;
} /* mkP3Msg */

static SYSTEM_char *
mkP3Msg2 (const char *txt1, const char *txt2, SYSTEM_shortstring msg)
{
  int n;

  if ((NULL == txt2) || ('\0' == *txt2))
    return mkP3Msg (txt1, msg);
  n = snprintf ((char *)msg+1, 255, "%s: %s", txt1, txt2);
  /* stupid MS _snprintf */
  if ((n > 255) || (n < 0))
    n = 255;
  *msg = (SYSTEM_byte) n;
  return msg;
} /* mkP3Msg2 */

#if defined(_WIN32)
static char *
winErrMsg (int errNum, char *buf, int bufSiz)
{
  char *p;
  BOOL brc;

  *buf = '\0';
  if (0 == errNum)
    return buf;
  brc = FormatMessage(
                FORMAT_MESSAGE_FROM_SYSTEM,
                NULL,
                errNum,
                MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // Default language
                buf,
                bufSiz-1,
                NULL
                );
  if (! brc) {
    *buf = '\0';
    return buf;
  }
  buf[bufSiz-1] = '\0';

  /* Trim the end of the line and terminate it with a null */
  p = buf;
  while (( *p > 31) || (9 == *p))
    ++p;
  do {
    *p-- = 0;
     } while ( (p >= buf) && (('.' == *p) || (33 > *p))	);
  return buf;
} /* winErrMsg */
#endif

typedef SYSTEM_uint32 _sub_14SYSUTILS_P3;
typedef SYSTEM_ansichar SYSUTILS_P3_tchararray[1000001];
typedef SYSUTILS_P3_tchararray *SYSUTILS_P3_pchararray;
static _P3STR_3 SYSUTILS_P3_filestopper, SYSUTILS_P3_extstopper;

static Procedure SYSUTILS_P3_pcharconcatpchar(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  SYSTEM_P3_pansichar psrc)
{
  SYSTEM_integer k;

  if (psrc != NULL) {
    k = 0;
    while ((*ValueCast(SYSUTILS_P3_pchararray,psrc))[k] != _P3char('\000')) {
      (*ValueCast(SYSUTILS_P3_pchararray,pdest))[*w] = (*ValueCast(
        SYSUTILS_P3_pchararray,psrc))[k];
      _P3inc0(*w);
      _P3inc0(k);
    
}
    (*ValueCast(SYSUTILS_P3_pchararray,pdest))[*w] = _P3char('\000');
  } 
}  /* pcharconcatpchar */

static Function(SYSTEM_ansichar *) SYSUTILS_P3_pchartostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pansichar p)
{
  SYSTEM_integer k;

  if (p == NULL) { 
    _P3strclr(result);
  } else {
    k = 0;
    do {
      if ((*ValueCast(SYSUTILS_P3_pchararray,p))[k] == _P3char('\000')) {
        _P3setlength(result,k,255);
        SYSTEM_break(BRK_1);
      } 
      result[k + 1] = (*ValueCast(SYSUTILS_P3_pchararray,p))[k];
      k = k + 1;
      if (k >= 255) {
        _P3strcpy(result,_len_ret,_P3str1("\023PCharToStr Overflow"));
        SYSTEM_break(BRK_1);
      } 
    CNT_1:;
    } while (SYSTEM_true);
BRK_1:;
  } 
  return result;
}  /* pchartostr */

static Procedure SYSUTILS_P3_divmod(
  SYSTEM_integer dividend,
  SYSTEM_word divisor,
  SYSTEM_word *result,
  SYSTEM_word *remainder)
{
  SYSTEM_integer quotient;

  if (dividend < 0) 
    _P3_RAISE(ValueCast(EXCEPTIONS_eintoverflow,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&EXCEPTIONS_eintoverflow_CD)),_P3str1("\020Integer overflow"))));
  if (divisor == 0) 
    _P3_RAISE(ValueCast(EXCEPTIONS_edivbyzero,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&EXCEPTIONS_edivbyzero_CD)),_P3str1("\020Division by zero"))));
  quotient = dividend /  divisor;
  if (quotient > 65535) 
    _P3_RAISE(ValueCast(EXCEPTIONS_eintoverflow,
      SYSTEM_exception_DOT_create(ValueCast(SYSTEM_exception,
      _P3alloc_object(&EXCEPTIONS_eintoverflow_CD)),_P3str1("\020Integer overflow"))));
  *result = quotient;
  *remainder = dividend - quotient * divisor;
}  /* divmod */

Function(SYSTEM_pointer ) SYSUTILS_P3_allocmem(
  SYSTEM_cardinal sz)
{
  SYSTEM_pointer result;

  _P3getmem(result,sz);
  /**** C code included from sysutils_p3.pas(472:1): 1 lines ****/
  (void) memset(result, 0, (size_t)sz);
  return result;
}  /* allocmem */

Function(SYSTEM_ansichar *) SYSUTILS_P3_uppercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer l;

  l = SYSTEM_length(s);
  _P3setlength(result,l,255);
  while (l > 0) {
    result[l] = SYSTEM_upcase(s[l]);
    _P3dec0(l);
  
}
  return result;
}  /* uppercase */

Function(SYSTEM_ansichar *) SYSUTILS_P3_lowercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer l;

  l = SYSTEM_length(s);
  _P3setlength(result,l,255);
  while (l > 0) {
    result[l] = s[l];
    if (result[l] >= _P3char('A') && result[l] <= _P3char('Z')) 
      _P3inc1(result[l],32);
    _P3dec0(l);
  
}
  return result;
}  /* lowercase */

Function(SYSTEM_integer ) SYSUTILS_P3_comparestr(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_integer result;
  SYSTEM_integer i;
  SYSTEM_integer m;

  m = SYSTEM_length(s1);
  if (m > ValueCast(SYSTEM_int32,SYSTEM_length(s2))) 
    m = SYSTEM_length(s2);
  { register SYSTEM_int32 _stop = m;
    if ((i = 1) <=  _stop) do {
      result = SYSTEM_ord(s1[i]) - SYSTEM_ord(s2[i]);
      if (result != 0) 
        return result;
    
    } while (i++ !=  _stop);

  }
  result = ValueCast(SYSTEM_int32,SYSTEM_length(s1)) - SYSTEM_length(
    s2);
  return result;
}  /* comparestr */

Function(SYSTEM_integer ) SYSUTILS_P3_comparetext(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_integer result;
  SYSTEM_integer i;
  SYSTEM_integer m;

  m = SYSTEM_length(s1);
  if (m > ValueCast(SYSTEM_int32,SYSTEM_length(s2))) 
    m = SYSTEM_length(s2);
  { register SYSTEM_int32 _stop = m;
    if ((i = 1) <=  _stop) do {
      result = SYSTEM_ord(SYSTEM_upcase(s1[i])) - SYSTEM_ord(
        SYSTEM_upcase(s2[i]));
      if (result != 0) 
        return result;
    
    } while (i++ !=  _stop);

  }
  result = ValueCast(SYSTEM_int32,SYSTEM_length(s1)) - SYSTEM_length(
    s2);
  return result;
}  /* comparetext */

Function(SYSTEM_boolean ) SYSUTILS_P3_sametext(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_boolean result;
  SYSTEM_integer i;

  result = SYSTEM_length(s1) == SYSTEM_length(s2);
  if (result) 
    { register SYSTEM_int32 _stop = SYSTEM_length(s1);
      if ((i = 1) <=  _stop) do {
        result = SYSTEM_upcase(s1[i]) == SYSTEM_upcase(s2[i]);
        if (!result) 
          return result;
      
      } while (i++ !=  _stop);

    }
  return result;
}  /* sametext */

Function(SYSTEM_ansichar *) SYSUTILS_P3_trim(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer i, l;

  l = SYSTEM_length(s);
  i = 1;
  while (i <= l && s[i] <= _P3char(' ')) {

    _P3inc0(i);
}
  if (i > l) { 
    _P3strclr(result);
  } else {
    while (s[l] <= _P3char(' ')) {

      _P3dec0(l);
}
    SYSTEM_copy(result,_len_ret,s,i,l - i + 1);
  } 
  return result;
}  /* trim */

Function(SYSTEM_ansichar *) SYSUTILS_P3_trimleft(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer i, l;

  l = SYSTEM_length(s);
  i = 1;
  while (i <= l && s[i] <= _P3char(' ')) {

    _P3inc0(i);
}
  SYSTEM_copy(result,_len_ret,s,i,SYSTEM_maxint);
  return result;
}  /* trimleft */

Function(SYSTEM_ansichar *) SYSUTILS_P3_trimright(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer i;

  i = SYSTEM_length(s);
  while (i > 0 && s[i] <= _P3char(' ')) {

    _P3dec0(i);
}
  SYSTEM_copy(result,_len_ret,s,1,i);
  return result;
}  /* trimright */

Function(SYSTEM_ansichar *) SYSUTILS_P3_inttostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n)
{
  SYSTEM_int64 w, w2;

  if (n < 0) {
    result[1] = _P3char('-');
    n = -n;
    w2 = 1;
  } else 
    w2 = 0;
  w = 255;
  do {
    result[w] = ValueCast(SYSTEM_ansichar,n % 10 + 48);
    _P3dec0(w);
    n = n /  10;
  } while (!(n == 0));
  while (w < 255) {
    _P3inc0(w2);
    _P3inc0(w);
    result[w2] = result[w];
  
}
  _P3setlength(result,w2,255);
  return result;
}  /* inttostr */

Function(SYSTEM_ansichar *) SYSUTILS_P3_inttohex(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 v,
  SYSTEM_integer d)
{
  static SYSTEM_shortstring hex = {16,'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F'};
  cnstdef {len = 16};
  _P3STR_31 buf;
  SYSTEM_integer i;

  i = len;
  if (d > 16) 
    d = len;
  _P3setlength(buf,len,16);
  do {
    buf[i] = hex[(15 & v) + 1];
    v = ValueCast(SYSTEM_uint64,v) >> 4;
    i = i - 1;
  } while (!(v == 0 && d + i <= 16));
  SYSTEM_copy(result,_len_ret,buf,i + 1,len);
  return result;
}  /* inttohex */

Function(SYSTEM_integer ) SYSUTILS_P3_strtoint(
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_int64 res;

  res = SYSUTILS_P3_strtoint64(s);
  if (res <  SYSTEM_minint || res > 2147483647) { 
    result =  SYSTEM_minint;
  } else 
    result = res;
  return result;
}  /* strtoint */

Function(SYSTEM_int64 ) SYSUTILS_P3_strtoint64(
  const SYSTEM_ansichar *s)
{
  SYSTEM_int64 result;
  SYSTEM_int64 res;
  SYSTEM_integer i;
  SYSTEM_boolean error, negative, hex;

  res = 0;
  i = 1;
  error = SYSTEM_false;
  while (i <= ValueCast(SYSTEM_int32,SYSTEM_length(s)) && s[i] == _P3char(' ')) {

    _P3inc0(i);
}
  if (i <= ValueCast(SYSTEM_int32,SYSTEM_length(s)) && s[i] == _P3char('-')) {
    negative = SYSTEM_true;
    _P3inc0(i);
  } else 
    negative = SYSTEM_false;
  if (i <= ValueCast(SYSTEM_int32,SYSTEM_length(s)) && s[i] == _P3char('$')) {
    hex = SYSTEM_true;
    _P3inc0(i);
  } else 
    hex = SYSTEM_false;
  if (hex) { 
    while (i <= ValueCast(SYSTEM_int32,SYSTEM_length(s))) {
      if (_P3SET_in_3(s[i],_P3char('0'),_P3char('9'))) { 
        res = 16 * res + SYSTEM_ord(s[i]) - 48;
      } else 
        if (_P3SET_in_3(s[i],_P3char('A'),_P3char('F'))) { 
          res = 16 * res + SYSTEM_ord(s[i]) - 65 + 10;
        } else 
          if (_P3SET_in_3(s[i],_P3char('a'),_P3char('f'))) { 
            res = 16 * res + SYSTEM_ord(s[i]) - 97 + 10;
          } else 
            error = SYSTEM_true;
      _P3inc0(i);
    
}
  } else 
    while (i <= ValueCast(SYSTEM_int32,SYSTEM_length(s))) {
      if (_P3SET_in_3(s[i],_P3char('0'),_P3char('9'))) { 
        res = 10 * res + SYSTEM_ord(s[i]) - 48;
      } else 
        error = SYSTEM_true;
      _P3inc0(i);
    
}
  if (error) { 
    res = SYSTEM_minint64;
  } else 
    if (negative) 
      res = -res;
  result = res;
  return result;
}  /* strtoint64 */

Function(SYSTEM_integer ) SYSUTILS_P3_fileage(
  const SYSTEM_ansichar *filename)
{
  SYSTEM_integer result;
  typedef SYSTEM_uint8 _sub_1FILEAGE;
  typedef SYSTEM_ansichar _arr_0FILEAGE[256];
  _arr_0FILEAGE fname;
  SYSTEM_byte len;

  len = ValueCast(SYSTEM_uint8,filename[0]);
  SYSTEM_move(&filename[1],fname,len);
  fname[len] = _P3char('\000');
  /**** C code included from sysutils_p3.pas(755:1): 32 lines ****/
#if defined(_WIN32)
{
  HANDLE          handle;
  WIN32_FIND_DATA findData;
  FILETIME        localFileTime;
  WORD *wPtr;
  /* */
  handle = FindFirstFile((char *)fname, &findData);
  if (INVALID_HANDLE_VALUE != handle) {
    FindClose(handle);
    if (! (findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)) {
      FileTimeToLocalFileTime(&findData.ftLastWriteTime, &localFileTime);
      wPtr = (WORD *) &result;
      if (FileTimeToDosDateTime(&localFileTime, wPtr+1, wPtr)) {
        return result;
      }
    }
  }
  result = -1;
}
#else
{
  struct stat statBuf;
  /* */
  if (0 == stat((char *)fname, &statBuf)) {  /* successful stat */
    result = statBuf.st_mtime;
  }
  else {
    result = -1;
  }
}
#endif
  return result;
}  /* fileage */

Function(SYSTEM_boolean ) SYSUTILS_P3_fileexists(
  const SYSTEM_ansichar *filename)
{
  SYSTEM_boolean result;

  /**** C code included from sysutils_p3.pas(837:1): 12 lines ****/
char buf[256];
unsigned char len;
/* */
len = filename[0];
memcpy(buf, filename + 1, len);
buf[len] = '\0';
/* */
#if defined(_WIN32)
result = (0 == _access(buf, 0));
#else
result = (0 == access(buf, F_OK));
#endif
  return result;
}  /* fileexists */

Function(SYSTEM_boolean ) SYSUTILS_P3_directoryexists(
  const SYSTEM_ansichar *directory)
{
  SYSTEM_boolean result;

  /**** C code included from sysutils_p3.pas(860:1): 25 lines ****/
char dirBuf[256];
unsigned char len;
/* */
len = directory[0];
memcpy(dirBuf, directory + 1, len);
dirBuf[len] = '\0';
/* */
#if defined(_WIN32)
{
  int attribs;
  /* */
  attribs = GetFileAttributes(dirBuf);
  result  = (-1 != attribs) && (attribs & FILE_ATTRIBUTE_DIRECTORY);
}
#else
{
  struct stat statBuf;
  /* */
  if (0 == stat(dirBuf, &statBuf)) {  /* successful stat */
    result = S_ISDIR(statBuf.st_mode);
  }
  else
    result = SYSTEM_false;
}
#endif
  return result;
}  /* directoryexists */

static Function(SYSTEM_integer ) SYSUTILS_P3_findmatchingfile(
  SYSUTILS_P3_tsearchrec *f)
{
  SYSTEM_integer result;

  /**** C code included from sysutils_p3.pas(992:1): 129 lines ****/
#if defined(_WIN32)
{
  int len;
  FILETIME lastWriteTime;
  FILETIME localFileTime;
  WORD *wPtr;
  while (f->finddata.dwfileattributes & f->excludeattr) {
    if (! FindNextFile((HANDLE)f->findhandle, (PWIN32_FIND_DATA) &(f->finddata))) {
      result = GetLastError();
      return result;
    }
  }
  memcpy(&lastWriteTime, &(f->finddata.ftlastwritetime), sizeof(FILETIME));
  FileTimeToLocalFileTime(&lastWriteTime, &localFileTime);
  wPtr = (WORD *) &f->time;
  FileTimeToDosDateTime(&localFileTime, wPtr+1, wPtr);
  f->size = f->finddata.nfilesizelow;
  f->attr = f->finddata.dwfileattributes;
  len = strlen((char *)f->finddata.cfilename);
  if (len > 255) len = 255;
  strncpy ((char *)f->name + 1, (char *)f->finddata.cfilename, len);
  f->name[0] = len;
  result = 0;
  return result;
}
#else
{
#if ! defined(NAME_MAX)  /* go for something safe */
# define NAME_MAX 512
#endif
  int attr;
  /* the struct dirent contains a field d_name, astring;
   * on SGI, they declare char d_name[1] but write as much as they like to it,
   * limited only by NAME_MAX */
  char hackingbuf[sizeof(struct dirent)+NAME_MAX];
  struct dirent * const entryPtr = (struct dirent *) hackingbuf;
  struct dirent *readdirResult;
  DIR *dp;
  struct stat statbuf;
  struct stat linkstatbuf;
  int len, rc;
  char pattern[256];
  char fname[256];
  mode_t mode;
  /* */
  result = -1;
  readdirResult = NULL;
  dp = (DIR *) f->findhandle;
  /*  readdir_r(F.FindHandle, @Scratch, PtrDirEnt); */
#if defined(SOL)
  /* I have read that avoiding the reentrant (_r) calls in Solaris
   * will help maintain compatibility with older Solaris systems */
  readdirResult = readdir (dp);
  if (NULL != readdirResult) {
    len = readdirResult->d_reclen;
    if (len > sizeof(hackingbuf)) len = sizeof(hackingbuf); /* this is exceptional */
    memcpy (entryPtr, readdirResult, len);
  }
#else
  rc = readdir_r (dp, entryPtr, &readdirResult);
#endif
  if (readdirResult != NULL) {
    len = f->pattern[0];
    memcpy(pattern, f->pattern + 1, len);
    pattern[len] = '\0';
  }
  else {
  }
  while (readdirResult != NULL) {
    rc = fnmatch(pattern, entryPtr->d_name, 0);
    if (0 == rc) {
      /* F.PathOnly must include trailing backslash */
      /* FName := F.PathOnly + ShortString(PtrDirEnt.d_name) + #0; */
      len = f->pathonly[0];
      memcpy(fname, f->pathonly + 1, len);
      strcpy(fname + len, entryPtr->d_name);
      rc = lstat(fname, &statbuf);
      if (0 == rc) {
        attr = 0;
        mode = statbuf.st_mode;
        if (S_ISDIR(mode))
          attr |= SYSUTILS_P3_fadirectory;
        else if (! S_ISREG(mode)) {
          /* directories should not be treated as system files */
          if (S_ISLNK(mode)) {
            attr |= SYSUTILS_P3_fasymlink;
            if (0 == lstat(fname, &linkstatbuf)
                && S_ISDIR(linkstatbuf.st_mode))
              attr |= SYSUTILS_P3_fadirectory;
          }
          attr |= SYSUTILS_P3_fasysfile;
        }
        if (entryPtr->d_name[0] == '.' && entryPtr->d_name[1] != '\0')
          if (!  (entryPtr->d_name[1] == '.' && entryPtr->d_name[2] == '\0'))
            attr |= SYSUTILS_P3_fahidden;
        /* if (euidaccess(fname, W_OK) != 0) */
        if (access(fname, W_OK) != 0)
          attr |= SYSUTILS_P3_fareadonly;
        if (0 == (attr & f->excludeattr)) {
          f->size = statbuf.st_size;
          f->attr = attr;
          f->mode = statbuf.st_mode;
          /* f->name = entryPtr->d_name; */
          len = strlen(entryPtr->d_name);
          if (len > 255) len = 255;
          strncpy((char*) f->name + 1, entryPtr->d_name, len);
          f->name[0] = len;
          f->time = statbuf.st_mtime;
          result = 0;
          break;
        } /* matching file found */
      } /* lstat returns OK */
    } /* matches desired pattern */
# if defined(SOL)
    /* I have read that avoiding the reentrant (_r) calls in Solaris
     * will help maintain compatibility with older Solaris systems */
    readdirResult = readdir (dp);
    if (NULL != readdirResult) {
      len = readdirResult->d_reclen;
      if (len > sizeof(hackingbuf)) len = sizeof(hackingbuf); /* this is exceptional */
      memcpy (entryPtr, readdirResult, len);
    }
# else
    rc = readdir_r (dp, entryPtr, &readdirResult);
# endif
    result = -1;
  } /* readdir loop */
} /* end C#### block */
#endif /* defined(_WIN32) */
  return result;
}  /* findmatchingfile */

Function(SYSTEM_integer ) SYSUTILS_P3_findfirst(
  const SYSTEM_ansichar *path,
  SYSTEM_integer attr,
  SYSUTILS_P3_tsearchrec *f)
{
  SYSTEM_integer result;
  cnstdef {faspecial = 30};

  f->excludeattr = ~attr;
  f->excludeattr = f->excludeattr & faspecial;
  SYSUTILS_P3_extractfilepath(f->pathonly,255,path);
  SYSUTILS_P3_extractfilename(f->pattern,255,path);
  if (_P3strcmpE(f->pathonly,_P3str1("\000"))) 
    {
      SYSTEM_shortstring _t2;

      SYSUTILS_P3_includetrailingpathdelimiter(f->pathonly,255,
        SYSUTILS_P3_getcurrentdir(_t2,255));
    }
  /**** C code included from sysutils_p3.pas(1180:1): 39 lines ****/
#if defined(_WIN32)
{
  char pathbuf[256];
  int len;
  HANDLE fHandle;
  /* */
  len = path[0];
  memcpy(pathbuf, path+1, len);
  pathbuf[len] = '\0';
  f->findhandle = (SYSTEM_pointer) fHandle
                = FindFirstFile(pathbuf, (PWIN32_FIND_DATA) &(f->finddata));
  if (INVALID_HANDLE_VALUE != fHandle) {
    result = SYSUTILS_P3_findmatchingfile(f);
    if (result != 0)
      SYSUTILS_P3_findclose (f);
  }
  else {
    result = GetLastError();
  }
}
#else
{
  DIR *dp;
  char pathonly[256];
  int len;
  len = f->pathonly[0];
  memcpy(pathonly, f->pathonly + 1, len);
  pathonly[len] = '\0';
  f->findhandle = (SYSTEM_pointer) (dp = opendir(pathonly));
  if (NULL != dp) {
    result = SYSUTILS_P3_findmatchingfile (f);
    if (result != 0)
      SYSUTILS_P3_findclose (f);
  }
  else {
    result = errno; /* what should this be?? */
  }
}
#endif /* defined(_WIN32) */
  return result;
}  /* findfirst */

Function(SYSTEM_integer ) SYSUTILS_P3_findnext(
  SYSUTILS_P3_tsearchrec *f)
{
  SYSTEM_integer result;

  /**** C code included from sysutils_p3.pas(1235:1): 10 lines ****/
#if defined(_WIN32)
  if (FindNextFile((HANDLE)f->findhandle, (PWIN32_FIND_DATA) &(f->finddata))) {
    result = SYSUTILS_P3_findmatchingfile(f);
  }
  else {
    result = GetLastError();
  }
#else
  result = SYSUTILS_P3_findmatchingfile(f);
#endif /* defined(_WIN32) */
  return result;
}  /* findnext */

Procedure SYSUTILS_P3_findclose(
  SYSUTILS_P3_tsearchrec *f)
{
  /**** C code included from sysutils_p3.pas(1267:1): 13 lines ****/
#if defined(_WIN32)
  if (INVALID_HANDLE_VALUE != (HANDLE) f->findhandle) {
    FindClose((HANDLE) f->findhandle);
    f->findhandle = (SYSTEM_pointer) INVALID_HANDLE_VALUE;
  }
#elif defined(AIX) || defined(BGP) || defined(__APPLE__) || defined(__linux__) || defined(SIG) || defined(SOL)
  if (NULL != f->findhandle) {
    closedir((DIR *) f->findhandle);
    f->findhandle = NULL;
  }
#else
# error "This OS not yet implemented"
#endif
}  /* findclose */

Function(SYSTEM_boolean ) SYSUTILS_P3_deletefile(
  const SYSTEM_ansichar *filename)
{
  SYSTEM_boolean result;
  typedef SYSTEM_uint8 _sub_1DELETEFILE;
  typedef SYSTEM_ansichar _arr_0DELETEFILE[256];
  _arr_0DELETEFILE fname;
  SYSTEM_P3_pansichar p;

  p = ValueCast(SYSTEM_P3_pansichar,&fname[0]);
  SYSUTILS_P3_strpcopy(p,filename);
  /**** C code included from sysutils_p3.pas(1300:1): 5 lines ****/
#if defined(_WIN32)
  result = DeleteFile((char *)p);
#else
  result = (unlink((char *)p) != -1);
#endif
  return result;
}  /* deletefile */

Function(SYSTEM_boolean ) SYSUTILS_P3_renamefile(
  const SYSTEM_ansichar *oldname,
  const SYSTEM_ansichar *newname)
{
  SYSTEM_boolean result;
  typedef SYSTEM_uint8 _sub_1RENAMEFILE;
  typedef SYSTEM_ansichar _arr_0RENAMEFILE[256];
  _arr_0RENAMEFILE alt, neu;
  SYSTEM_P3_pansichar palt, pneu;

  palt = ValueCast(SYSTEM_P3_pansichar,&alt[0]);
  SYSUTILS_P3_strpcopy(palt,oldname);
  pneu = ValueCast(SYSTEM_P3_pansichar,&neu[0]);
  SYSUTILS_P3_strpcopy(pneu,newname);
  /**** C code included from sysutils_p3.pas(1328:1): 5 lines ****/
#if defined(_WIN32)
  result = MoveFileA((char *)palt, (char *)pneu);
#else
  result = (0 == rename((char *)palt, (char *)pneu));
#endif
  return result;
}  /* renamefile */

Function(SYSTEM_integer ) SYSUTILS_P3_lastdelimiter(
  const SYSTEM_ansichar *delimiters,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer i;

  result = SYSTEM_length(s);
  while (result > 0) {
    { register SYSTEM_int32 _stop = SYSTEM_length(delimiters);
      if ((i = 1) <=  _stop) do {
        if (s[result] == delimiters[i]) 
          return result;
      } while (i++ !=  _stop);

    }
    _P3dec0(result);
  
}
  return result;
}  /* lastdelimiter */

Function(SYSTEM_ansichar *) SYSUTILS_P3_changefileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *extension)
{
  SYSTEM_integer i;

  i = SYSUTILS_P3_lastdelimiter(SYSUTILS_P3_extstopper,filename);
  if (i == 0 || filename[i] != _P3char('.')) { 
    i = SYSTEM_length(filename);
  } else 
    i = i - 1;
  {
    SYSTEM_shortstring _t1;

    _P3strcat(result,_len_ret,SYSTEM_copy(_t1,255,filename,1,i),
      extension);
  }
  return result;
}  /* changefileext */

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfilepath(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename)
{
  SYSTEM_integer i;

  {
    _P3STR_3 _t1;
    _P3STR_3 _t2;
    _P3STR_3 _t3;

    i = SYSUTILS_P3_lastdelimiter(_P3strcat(_t3,2,_P3ch2str(_t1,1,
      SYSUTILS_P3_pathdelim),_P3ch2str(_t2,1,
      SYSUTILS_P3_drivedelim)),filename);
  }
  SYSTEM_copy(result,_len_ret,filename,1,i);
  return result;
}  /* extractfilepath */

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfiledir(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename)
{
  SYSTEM_integer i;

  {
    _P3STR_3 _t1;
    _P3STR_3 _t2;
    _P3STR_3 _t3;

    i = SYSUTILS_P3_lastdelimiter(_P3strcat(_t3,2,_P3ch2str(_t1,1,
      SYSUTILS_P3_pathdelim),_P3ch2str(_t2,1,
      SYSUTILS_P3_drivedelim)),filename);
  }
  if (i > 1 && filename[i] == SYSUTILS_P3_pathdelim) 
    if (filename[i - 1] != SYSUTILS_P3_pathdelim && filename[i - 1] != 
      SYSUTILS_P3_drivedelim) 
      _P3dec0(i);
  SYSTEM_copy(result,_len_ret,filename,1,i);
  return result;
}  /* extractfiledir */

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfilename(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename)
{
  SYSTEM_integer i;

  {
    _P3STR_3 _t1;
    _P3STR_3 _t2;
    _P3STR_3 _t3;

    i = SYSUTILS_P3_lastdelimiter(_P3strcat(_t3,2,_P3ch2str(_t1,1,
      SYSUTILS_P3_pathdelim),_P3ch2str(_t2,1,
      SYSUTILS_P3_drivedelim)),filename);
  }
  SYSTEM_copy(result,_len_ret,filename,i + 1,SYSTEM_maxint);
  return result;
}  /* extractfilename */

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename)
{
  SYSTEM_integer i;

  i = SYSUTILS_P3_lastdelimiter(SYSUTILS_P3_extstopper,filename);
  if (i > 0 && filename[i] == _P3char('.')) { 
    SYSTEM_copy(result,_len_ret,filename,i,SYSTEM_maxint);
  } else 
    _P3strclr(result);
  return result;
}  /* extractfileext */

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractshortpathname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename)
{
  /**** C code included from sysutils_p3.pas(1420:1): 17 lines ****/
#if defined(_WIN32)
  {
  char inbuf[256], outbuf[MAX_PATH];
  char *p;
  int nChars, len;

  nChars = *filename;
  memcpy(inbuf, filename+1, nChars);
  inbuf[nChars] = '\0';
  len = GetShortPathName(inbuf, outbuf, sizeof(outbuf));
  if (len > _len_ret)
    len = _len_ret;
  memcpy((char *)result+1, outbuf, *result = (SYSTEM_char)len);
  }
#else
  *result = '\0';
#endif
  return result;
}  /* extractshortpathname */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_filedatetodatetime(
  SYSTEM_integer filedate)
{
  SYSTEM_P3_tdatetime result;

  result = 0;
  if (P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin) {
    result = SYSUTILS_P3_encodedate(ValueCast(SYSTEM_int32,(VariableCast(
      SYSUTILS_P3_longrec,&filedate,SYSTEM_int32))._u._c1.hi >> 9) + 1980,ValueCast(
      SYSTEM_int32,(VariableCast(SYSUTILS_P3_longrec,&filedate,
      SYSTEM_int32))._u._c1.hi >> 5) & 15,ValueCast(
      SYSTEM_int32,(VariableCast(SYSUTILS_P3_longrec,&filedate,
      SYSTEM_int32))._u._c1.hi) & 31) + SYSUTILS_P3_encodetime((VariableCast(
      SYSUTILS_P3_longrec,&filedate,SYSTEM_int32))._u._c1.lo >> 11,ValueCast(
      SYSTEM_int32,(VariableCast(SYSUTILS_P3_longrec,&filedate,
      SYSTEM_int32))._u._c1.lo >> 5) & 63,ValueCast(
      SYSTEM_uint32,ValueCast(SYSTEM_int32,(VariableCast(
      SYSUTILS_P3_longrec,&filedate,SYSTEM_int32))._u._c1.lo) & 31) << 1,0);
    return result;
  } 
  /**** C code included from sysutils_p3.pas(1459:1): 11 lines ****/
{
#if ! defined(_WIN32)
  struct tm ut;
  time_t tim;
  /* */
  tim = filedate;
  localtime_r (&tim, &ut);
  result = SYSUTILS_P3_encodedate(ut.tm_year+1900, ut.tm_mon + 1, ut.tm_mday)
         + SYSUTILS_P3_encodetime(ut.tm_hour, ut.tm_min, ut.tm_sec, 0);
#endif
}
  return result;
}  /* filedatetodatetime */

Function(SYSTEM_integer ) SYSUTILS_P3_datetimetofiledate(
  SYSTEM_P3_tdatetime datetime)
{
  SYSTEM_integer result;
  SYSTEM_word year, month, day, hour, minu, sec, msec;

  SYSUTILS_P3_decodedate(datetime,&year,&month,&day);
  if (P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin) { 
    if (year < 1980 || year > 2107) { 
      result = 0;
    } else {
      SYSUTILS_P3_decodetime(datetime,&hour,&minu,&sec,&msec);
      (VariableCast(SYSUTILS_P3_longrec,&result,SYSTEM_int32)).
        _u._c1.lo = ValueCast(SYSTEM_int32,sec >> 1) | minu << 5 | 
        hour << 11;
      (VariableCast(SYSUTILS_P3_longrec,&result,SYSTEM_int32)).
        _u._c1.hi = ValueCast(SYSTEM_int32,day) | month << 5 | ValueCast(
        SYSTEM_uint32,ValueCast(SYSTEM_int32,year) - 1980) << 9;
    } 
  } else 
    if (year < 1970 || year > 2038) { 
      result = 0;
    } else {
      SYSUTILS_P3_decodetime(datetime,&hour,&minu,&sec,&msec);
      /**** C code included from sysutils_p3.pas(1521:1): 17 lines ****/
#if defined(_WIN32)
      result = -1;
#else
{
      struct tm tm;
      tm.tm_sec  = sec;
      tm.tm_min  = minu;
      tm.tm_hour = hour;
      tm.tm_mday = day;
      tm.tm_mon  = month - 1;
      tm.tm_year = year - 1900;
      tm.tm_wday = 0; /* ignored anyway */
      tm.tm_yday = 0; /* ignored anyway */
      tm.tm_isdst = -1;
      result = mktime(&tm);
}
#endif
    } 
  return result;
}  /* datetimetofiledate */

Function(SYSTEM_ansichar *) SYSUTILS_P3_getcurrentdir(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  SYSTEM_boolean isok;
  SYSTEM_shortstring emsg;

  isok = SYSTEM_true;
  _P3strclr(emsg);
  /**** C code included from sysutils_p3.pas(1598:1): 57 lines ****/
{
  /* this implementation assumes the return is a ShortString, */
  /* i.e. char buf[256] or less */
  char buf[256];
  unsigned char len;
#if defined(_WIN32)
  int rc;
  rc = GetCurrentDirectory(sizeof(buf),buf);
  if (0 == rc) {
    isok = 0;
    (void) winErrMsg (GetLastError(), buf, sizeof(buf));
    (void) mkP3Msg2 ("GetCurrentDir failed", buf, emsg);
  }
  else if (rc > sizeof(buf)) {
    isok = 0;
    (void) mkP3Msg ("GetCurrentDir failed: result too large for shortString", emsg);
  }
#else
  if (getcwd(buf, 256) == NULL) {
    if (ERANGE == errno) {
      isok = 0;
      (void) mkP3Msg ("GetCurrentDir failed: result too large for shortString", emsg);
    }
    else {
      char *p;
      isok = 0;
      p = strerror (errno);
      if (p)
        (void) mkP3Msg2 ("GetCurrentDir failed", p, emsg);
      else
        (void) mkP3Msg2 ("GetCurrentDir failed", "libc failure", emsg);
    }
  }
  else {
    /* getcwd OK, but check if we can do better */
    /* # if defined(__linux__) realpath() expected everywhere */
    const char *sym, *ss3;
    sym = getenv("PWD");
    char absp[4096];
    if (sym) {     /* got something, check if it is really same as getcwd */
      /* realpath(p,absp) converts the relative or symlink-ish path p
       * to an absolute physical path absp */
      ss3 = realpath(sym,absp);
      if (ss3 && (0 == strcmp(buf,absp)) && (strlen(sym) < 256)) {
        strcpy(buf,sym);
      }
    }
    /* # endif */
  }
#endif
  if (isok) {
    len = strlen(buf);
    if (len > _len_ret)
      len = _len_ret;
    memcpy((char *)result+1, buf, *result = (SYSTEM_char)len);
  }
}
  if (!isok) 
    _P3_RAISE(ValueCast(SYSTEM_exception,SYSTEM_exception_DOT_create(ValueCast(
      SYSTEM_exception,_P3alloc_object(&SYSTEM_exception_CD)),emsg)));
  return result;
}  /* getcurrentdir */

Function(SYSTEM_boolean ) SYSUTILS_P3_setcurrentdir(
  const SYSTEM_ansichar *dir)
{
  SYSTEM_boolean result;
  typedef SYSTEM_uint8 _sub_1SETCURRENTDIR;
  typedef SYSTEM_ansichar _arr_0SETCURRENTDIR[256];
  _arr_0SETCURRENTDIR _dirname;
  SYSTEM_P3_pansichar p;

  p = ValueCast(SYSTEM_P3_pansichar,&_dirname[0]);
  SYSUTILS_P3_strpcopy(p,dir);
  /**** C code included from sysutils_p3.pas(1684:1): 5 lines ****/
#if defined(_WIN32)
  result = SetCurrentDirectory((char *)p);
#else
  result = (chdir((char *)p) == 0);
#endif
  return result;
}  /* setcurrentdir */

Function(SYSTEM_boolean ) SYSUTILS_P3_createdir(
  const SYSTEM_ansichar *dir)
{
  SYSTEM_boolean result;

  _Iminus_bgn();
  SYSTEM_mkdir(dir);
  _Iminus_end();
  result = SYSTEM_ioresult() == 0;
  return result;
}  /* createdir */

Function(SYSTEM_boolean ) SYSUTILS_P3_removedir(
  const SYSTEM_ansichar *dir)
{
  SYSTEM_boolean result;
  typedef SYSTEM_uint8 _sub_1REMOVEDIR;
  typedef SYSTEM_ansichar _arr_0REMOVEDIR[256];
  _arr_0REMOVEDIR _dirname;
  SYSTEM_P3_pansichar p;

  p = ValueCast(SYSTEM_P3_pansichar,&_dirname[0]);
  SYSUTILS_P3_strpcopy(p,dir);
  /**** C code included from sysutils_p3.pas(1717:1): 5 lines ****/
#if defined(_WIN32)
  result = RemoveDirectory((char *)p);
#else
  result = (rmdir((char *)p) == 0);
#endif
  return result;
}  /* removedir */

Function(SYSTEM_cardinal ) SYSUTILS_P3_strlen(
  SYSTEM_P3_pansichar str)
{
  SYSTEM_cardinal result;

  /**** C code included from sysutils_p3.pas(1731:1): 1 lines ****/
  result = strlen((char *)str);
  return result;
}  /* strlen */

Function(SYSTEM_P3_pansichar ) SYSUTILS_P3_strpcopy(
  SYSTEM_P3_pansichar dest,
  const SYSTEM_ansichar *src)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer len;

  len = SYSTEM_length(src);
  SYSTEM_move(&src[1],dest,len);
  /**** C code included from sysutils_p3.pas(1747:1): 1 lines ****/
  dest[len] = '\0';
  result = dest;
  return result;
}  /* strpcopy */

Function(SYSTEM_ansichar *) SYSUTILS_P3_floattostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v)
{
  _P3STR_95 s;
  SYSTEM_integer i, j, k, e;

  if (v == 0.0) {
    _P3strcpy(result,_len_ret,_P3str1("\0010"));
    return result;
  } 
  _P3str_d0(v,s,64);
  if (v < 0.0) 
    v = -v;
  k = SYSUTILS_P3_lastdelimiter(_P3str1("\002+-"),s);
  j = SYSTEM_pos(_P3str1("\001."),s);
  if (v >= 1e-4 && v < 1e15) {
    {
      SYSTEM_shortstring _t1;

      _P3val_i(SYSTEM_copy(_t1,255,s,k,5),e,&i);
    }
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((i = k - 1) <=  _stop) do {
        s[i] = _P3char('0');
      } while (i++ !=  _stop);

    }
    if (e >= 0) {
      { register SYSTEM_int32 _stop = j + e;
        if ((i = j + 1) <=  _stop) do {
          s[i - 1] = s[i];
        } while (i++ !=  _stop);

      }
      s[j + e] = _P3char('.');
      { register SYSTEM_int32 _stop = j + e + 1;
        if ((i = SYSTEM_length(s)) >=  _stop) do {
          if (s[i] == _P3char('0')) {
            s[i] = _P3char(' ');
            if (i == j + e + 1) 
              s[j + e] = _P3char(' ');
          } else 
            SYSTEM_break(BRK_2);
CNT_2:;
        } while (i-- !=  _stop);
BRK_2:;

      }
    } else {
      s[j] = s[j - 1];
      s[j - 1] = _P3char('0');
      e = -e;
      { register SYSTEM_int32 _stop = j;
        if ((i = k - 2) >=  _stop) do {
          s[i + e] = s[i];
        } while (i-- !=  _stop);

      }
      { register SYSTEM_int32 _stop = j + e - 1;
        if ((i = j + 1) <=  _stop) do {
          s[i] = _P3char('0');
        } while (i++ !=  _stop);

      }
      s[j] = _P3char('.');
      _P3setlength(s,k + e - 2,64);
      { register SYSTEM_int32 _stop = j + e + 1;
        if ((i = SYSTEM_length(s)) >=  _stop) do {
          if (s[i] == _P3char('0')) { 
            s[i] = _P3char(' ');
          } else 
            SYSTEM_break(BRK_3);
CNT_3:;
        } while (i-- !=  _stop);
BRK_3:;

      }
    } 
  } else {
    if (s[k] == _P3char('+')) 
      s[k] = _P3char(' ');
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((i = k + 1) <=  _stop) do {
        if (s[i] == _P3char('0')) {
          s[i] = _P3char(' ');
          if (i == ValueCast(SYSTEM_int32,SYSTEM_length(s))) 
            s[k - 1] = _P3char(' ');
        } else 
          SYSTEM_break(BRK_4);
CNT_4:;
      } while (i++ !=  _stop);
BRK_4:;

    }
    { register SYSTEM_int32 _stop = j + 1;
      if ((i = k - 2) >=  _stop) do {
        if (s[i] == _P3char('0')) {
          s[i] = _P3char(' ');
          if (i == j + 1) 
            s[j] = _P3char(' ');
        } else 
          SYSTEM_break(BRK_5);
CNT_5:;
      } while (i-- !=  _stop);
BRK_5:;

    }
  } 
  j = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      if (s[i] != _P3char(' ')) {
        j = j + 1;
        result[j] = s[i];
      } 
    } while (i++ !=  _stop);

  }
  _P3setlength(result,j,255);
  return result;
}  /* floattostr */

Function(SYSUTILS_P3_ttimestamp *) SYSUTILS_P3_datetimetotimestamp(
  SYSUTILS_P3_ttimestamp *result,
  SYSTEM_P3_tdatetime datetime)
{
  SYSTEM_int64 i8;
  SYSTEM_double x;

  x = SYSTEM_abs_r(SYSTEM_frac(datetime));
  i8 = SYSTEM_round(x * SYSUTILS_P3_msecsperday);
  result->time = ValueCast(SYSTEM_int32,i8);
  i8 = SYSTEM_trunc(datetime) + SYSUTILS_P3_datedelta;
  result->date = ValueCast(SYSTEM_int32,i8);
  return result;
}  /* datetimetotimestamp */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_timestamptodatetime(
  const SYSUTILS_P3_ttimestamp *timestamp)
{
  SYSTEM_P3_tdatetime result;
  SYSTEM_P3_tdatetime t;

  t = timestamp->time;
  t = t /  SYSUTILS_P3_msecsperday;
  result = timestamp->date - SYSUTILS_P3_datedelta;
  if (result < 0) { 
    result = result - t;
  } else 
    result = result + t;
  return result;
}  /* timestamptodatetime */

Function(SYSTEM_boolean ) SYSUTILS_P3_tryencodetime(
  SYSTEM_word hour,
  SYSTEM_word _min,
  SYSTEM_word sec,
  SYSTEM_word msec,
  SYSTEM_P3_tdatetime *time)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  if (hour < SYSUTILS_P3_hoursperday && _min < SYSUTILS_P3_minsperhour && 
    sec < SYSUTILS_P3_secspermin && msec < SYSUTILS_P3_msecspersec) {
    *time = ValueCast(SYSTEM_double,hour * (3600000) + _min * (60000) + ValueCast(
      SYSTEM_int32,sec) * SYSUTILS_P3_msecspersec + msec) /  
      SYSUTILS_P3_msecsperday;
    result = SYSTEM_true;
  } 
  return result;
}  /* tryencodetime */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_encodetime(
  SYSTEM_word hour,
  SYSTEM_word _min,
  SYSTEM_word sec,
  SYSTEM_word msec)
{
  SYSTEM_P3_tdatetime result;
  SYSTEM_double mday;

  mday = 86400000;
  result = (hour * 3600000 + ValueCast(SYSTEM_int32,_min) * 60000 + ValueCast(
    SYSTEM_int32,sec) * 1000 + msec) /  mday;
  return result;
}  /* encodetime */

Procedure SYSUTILS_P3_decodetime(
  SYSTEM_P3_tdatetime datetime,
  SYSTEM_word *hour,
  SYSTEM_word *_min,
  SYSTEM_word *sec,
  SYSTEM_word *msec)
{
  SYSTEM_word mincount, mseccount;

  {
    SYSUTILS_P3_ttimestamp _t1;

    SYSUTILS_P3_divmod(SYSUTILS_P3_datetimetotimestamp(&_t1,datetime)->
      time,60000,&mincount,&mseccount);
  }
  SYSUTILS_P3_divmod(mincount,SYSUTILS_P3_minsperhour,hour,_min);
  SYSUTILS_P3_divmod(mseccount,SYSUTILS_P3_msecspersec,sec,msec);
}  /* decodetime */

Function(SYSTEM_boolean ) SYSUTILS_P3_isleapyear(
  SYSTEM_word year)
{
  SYSTEM_boolean result;

  result = ValueCast(SYSTEM_int32,year) % 4 == 0 && ValueCast(
    SYSTEM_int32,year) % 4000 != 0 && (ValueCast(SYSTEM_int32,
    year) % 100 != 0 || ValueCast(SYSTEM_int32,year) % 400 == 0);
  return result;
}  /* isleapyear */

Function(SYSTEM_boolean ) SYSUTILS_P3_tryencodedate(
  SYSTEM_word year,
  SYSTEM_word month,
  SYSTEM_word day,
  SYSTEM_P3_tdatetime *date)
{
  SYSTEM_boolean result;
  SYSTEM_integer i;
  SYSUTILS_P3_pdaytable daytable;

  result = SYSTEM_false;
  daytable = ValueCast(SYSUTILS_P3_pdaytable,SYSUTILS_P3_monthdays[
    SYSUTILS_P3_isleapyear(year)]);
  if (year >= 1 && year <= 9999 && month >= 1 && month <= 12 && 
    day >= 1 && day <= (*daytable)[month - 1]) {
    { register SYSTEM_int32 _stop = ValueCast(SYSTEM_int32,month) - 1;
      if ((i = 1) <=  _stop) do {
        _P3inc1(day,(*daytable)[i - 1]);
      } while (i++ !=  _stop);

    }
    i = ValueCast(SYSTEM_int32,year) - 1;
    *date = i * 365 + i /  4 - i /  100 + i /  400 + day - 
      SYSUTILS_P3_datedelta;
    result = SYSTEM_true;
  } 
  return result;
}  /* tryencodedate */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_encodedate(
  SYSTEM_word year,
  SYSTEM_word month,
  SYSTEM_word day)
{
  SYSTEM_P3_tdatetime result;
  SYSTEM_longint yr;

  if (year == 1600 && month < 3) { 
    if (month == 1) { 
      result = ValueCast(SYSTEM_int32,day) + 1;
    } else 
      result = ValueCast(SYSTEM_int32,day) + 30;
  } else {
    if (month > 2) { 
      month = ValueCast(SYSTEM_int32,month) - 3;
    } else {
      month = ValueCast(SYSTEM_int32,month) + 9;
      year = ValueCast(SYSTEM_int32,year) - 1;
    } 
    yr = ValueCast(SYSTEM_int32,year) - 1600;
    result = yr /  100 * 146097 /  4 + yr % 100 * 1461 /  4 + (153 * 
      month + 2) /  5 + day + 59 - 109572 + 1;
  } 
  return result;
}  /* encodedate */

Function(SYSTEM_boolean ) SYSUTILS_P3_decodedatefully(
  SYSTEM_P3_tdatetime datetime,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day,
  SYSTEM_word *dow)
{
  SYSTEM_boolean result;
  cnstdef {d1 = 365};
  cnstdef {d4 = 1461};
  cnstdef {d100 = 36524};
  cnstdef {d400 = 146097};
  SYSTEM_word y, m, d, i;
  SYSTEM_integer t;
  SYSUTILS_P3_pdaytable daytable;

  {
    SYSUTILS_P3_ttimestamp _t1;

    t = SYSUTILS_P3_datetimetotimestamp(&_t1,datetime)->date;
  }
  if (t <= 0) {
    *year = 0;
    *month = 0;
    *day = 0;
    *dow = 0;
    result = SYSTEM_false;
  } else {
    *dow = t % 7 + 1;
    _P3dec0(t);
    y = 1;
    while (t >= d400) {
      _P3dec1(t,d400);
      _P3inc1(y,400);
    
}
    SYSUTILS_P3_divmod(t,d100,&i,&d);
    if (i == 4) {
      _P3dec0(i);
      _P3inc1(d,d100);
    } 
    _P3inc1(y,ValueCast(SYSTEM_int32,i) * 100);
    SYSUTILS_P3_divmod(d,d4,&i,&d);
    _P3inc1(y,ValueCast(SYSTEM_int32,i) * 4);
    SYSUTILS_P3_divmod(d,d1,&i,&d);
    if (i == 4) {
      _P3dec0(i);
      _P3inc1(d,d1);
    } 
    _P3inc1(y,i);
    result = SYSUTILS_P3_isleapyear(y);
    daytable = ValueCast(SYSUTILS_P3_pdaytable,SYSUTILS_P3_monthdays[
      result]);
    m = 1;
    while (SYSTEM_true) {
      i = (*daytable)[m - 1];
      if (d < i) 
        SYSTEM_break(BRK_6);
      _P3dec1(d,i);
      _P3inc0(m);
    
CNT_6:;
    }
BRK_6:;
    *year = y;
    *month = m;
    *day = ValueCast(SYSTEM_int32,d) + 1;
  } 
  return result;
}  /* decodedatefully */

Procedure SYSUTILS_P3_decodedate(
  SYSTEM_P3_tdatetime datetime,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day)
{
  SYSTEM_word dummy;

  SYSUTILS_P3_decodedatefully(datetime,year,month,day,&dummy);
}  /* decodedate */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_date(void)
{
  SYSTEM_P3_tdatetime result;

  /**** C code included from sysutils_p3.pas(2051:1): 19 lines ****/
int rc;
#if defined(_WIN32)
SYSTEMTIME st;
GetLocalTime(&st);
rc  = SYSUTILS_P3_tryencodedate(st.wYear, st.wMonth, st.wDay, &result);
if (rc != 1)
  result = 0;
#else
time_t t;
struct tm lt;
(void) time(&t);
if (NULL == localtime_r(&t, &lt)) {
  result = 0;
  return result;
}
rc = SYSUTILS_P3_tryencodedate (lt.tm_year+1900, lt.tm_mon+1, lt.tm_mday, &result);
if (rc != 1)
  result = 0;
#endif
  return result;
}  /* date */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_time(void)
{
  SYSTEM_P3_tdatetime result;

  /**** C code included from sysutils_p3.pas(2081:1): 21 lines ****/
int rc;
#if defined(_WIN32)
SYSTEMTIME st;
GetLocalTime(&st);
rc = SYSUTILS_P3_tryencodetime (st.wHour, st.wMinute, st.wSecond, st.wMilliseconds,
                                &result);
if (!rc) result = 0;
#else
struct timeval tv;
struct tm lt;
if (gettimeofday(&tv, NULL)) {
  result = 0;
  return result;
}
if (NULL == localtime_r(&tv.tv_sec, &lt)) {
  result = 0;
  return result;
}
rc = SYSUTILS_P3_tryencodetime(lt.tm_hour, lt.tm_min, lt.tm_sec, tv.tv_usec/1000, &result);
if (!rc) result = 0;
#endif
  return result;
}  /* time */

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_now(void)
{
  SYSTEM_P3_tdatetime result;

  /**** C code included from sysutils_p3.pas(2112:1): 30 lines ****/
int rc;
double dnow, tnow;
#if defined(_WIN32)
SYSTEMTIME st;
GetLocalTime(&st);
rc  = SYSUTILS_P3_tryencodedate (st.wYear, st.wMonth, st.wDay, &dnow);
rc += SYSUTILS_P3_tryencodetime (st.wHour, st.wMinute, st.wSecond, st.wMilliseconds,
                                 &tnow);
if (rc != 2)
  result = 0;
else
  result = dnow + tnow;
#else
struct timeval tv;
struct tm lt;
if (gettimeofday(&tv, NULL)) {
  result = 0;
  return result;
}
if (NULL == localtime_r(&tv.tv_sec, &lt)) {
  result = 0;
  return result;
}
rc  = SYSUTILS_P3_tryencodedate (lt.tm_year+1900, lt.tm_mon+1, lt.tm_mday, &dnow);
rc += SYSUTILS_P3_tryencodetime (lt.tm_hour, lt.tm_min, lt.tm_sec, tv.tv_usec/1000, &tnow);
if (rc != 2)
  result = 0;
else
  result = dnow + tnow;
#endif
  return result;
}  /* now */

Function(SYSTEM_ansichar *) SYSUTILS_P3_syserrormessage(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer errorcode)
{
  /**** C code included from sysutils_p3.pas(2153:1): 15 lines ****/
{
  int i;
  char *errMsg = strerror(errorcode);
  if (NULL == errMsg) {
    SYSTEM_shortstring msg, _t1;

    _P3strcat(msg,255,_P3str1("\016Unknown error "), SYSUTILS_P3_inttostr(_t1,255,errorcode));
    _P3strcpy(result,_len_ret,msg);
    return result;
  }
  for (i = 0; (i < _len_ret) && (errMsg[i] != '\0'); i++)
    result[i+1] = errMsg[i];
  *result = i;
  return result;
}
  return result;
}  /* syserrormessage */

Function(SYSTEM_ansichar *) SYSUTILS_P3_includetrailingpathdelimiter(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  if (SYSTEM_length(s) > 0 && SYSUTILS_P3_pathdelim == s[
    SYSTEM_length(s)]) { 
    _P3strcpy(result,_len_ret,s);
  } else 
    {
      _P3STR_3 _t1;

      _P3strcat(result,_len_ret,s,_P3ch2str(_t1,1,
        SYSUTILS_P3_pathdelim));
    }
  return result;
}  /* includetrailingpathdelimiter */

Function(SYSTEM_ansichar *) SYSUTILS_P3_excludetrailingpathdelimiter(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  _P3strcpy(result,_len_ret,s);
  if (SYSTEM_length(s) > 0 && SYSUTILS_P3_pathdelim == result[
    SYSTEM_length(result)]) 
    _P3setlength(result,ValueCast(SYSTEM_int32,SYSTEM_length(result)) - 1,255);
  return result;
}  /* excludetrailingpathdelimiter */

Procedure SYSUTILS_P3_sleep(
  SYSTEM_cardinal milliseconds)
{
  /**** C code included from sysutils_p3.pas(2203:1): 14 lines ****/
#if defined(_WIN32)
  Sleep(milliseconds);
#else
{
  long nano;
  struct timespec req, rem;

  req.tv_sec = milliseconds / 1000; /* whole seconds */
  nano = milliseconds % 1000;
  nano *= 1000000;
  req.tv_nsec = nano;
  (void) nanosleep (&req, &rem);
}
#endif
}  /* sleep */

Procedure SYSUTILS_P3_freeandnil(
  SYSTEM_untyped *obj)
{
  SYSTEM_tobject temp;

  temp = PointerCast(SYSTEM_tobject,obj);
  PointerCast(SYSTEM_pointer,obj) = NULL;
  SYSTEM_tobject_DOT_free(temp);
}  /* freeandnil */

Function(SYSTEM_ansichar *) SYSUTILS_P3_getenvironmentvariable(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *name)
{
  /**** C code included from sysutils_p3.pas(2238:1): 40 lines ****/
{
  char buf[256];
  char *s;
  int len;

  memcpy (buf,name+1,(size_t) *name);
  buf[*name] = '\0';
#if defined(_WIN32)
  len = GetEnvironmentVariable(buf,NULL,0);
  if (0 == len) {
    *result = '\0';
  }
  else {
    s = (char *) malloc(len);
    if (NULL == s)
      *result = '\0';
    else {
      (void) GetEnvironmentVariable(buf,s,len);
      len--;  /* we want the length of the string, not the sizeof the buffer */
      if (len > _len_ret) {             /* input buffer too short */
        len = _len_ret;
      }
      memcpy ((char*)result+1, s, *result = (SYSTEM_char) len);
      free(s);
    }
  }
#else
  s = getenv (buf);
  if (NULL == s) {
    *result = '\0';
  }
  else {
    len = strlen (s);
    if (len > _len_ret) {               /* input buffer too short */
      len = _len_ret;
    }
    memcpy ((char*)result+1, s, *result = (SYSTEM_char) len);
  }
#endif
}
  return result;
}  /* getenvironmentvariable */

/* unit sysutils_p3 */
void _Init_Module_sysutils_p3(void)
{
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      SYSUTILS_P3_pathdelim = _P3char('\\');
      SYSUTILS_P3_drivedelim = _P3char(':');
      SYSUTILS_P3_pathsep = _P3char(';');
      _P3strcpy(SYSUTILS_P3_filestopper,3,_P3str1("\002\\:"));
      _P3strcpy(SYSUTILS_P3_extstopper,3,_P3str1("\003\\:."));
      break;
    case P3PLATFORM_osfileunix: 
      SYSUTILS_P3_pathdelim = _P3char('/');
      SYSUTILS_P3_drivedelim = _P3char('\000');
      SYSUTILS_P3_pathsep = _P3char(':');
      _P3strcpy(SYSUTILS_P3_filestopper,3,_P3str1("\001/"));
      _P3strcpy(SYSUTILS_P3_extstopper,3,_P3str1("\002/."));
      break;
    default:
      SYSUTILS_P3_pathdelim = _P3char('?');
      SYSUTILS_P3_drivedelim = _P3char('?');
      SYSUTILS_P3_pathsep = _P3char('?');
      _P3strcpy(SYSUTILS_P3_filestopper,3,_P3str1("\001?"));
      _P3strcpy(SYSUTILS_P3_extstopper,3,_P3str1("\001?"));
  }
} /* _Init_Module_sysutils_p3 */

void _Final_Module_sysutils_p3(void)
{
} /* _Final_Module_sysutils_p3 */

