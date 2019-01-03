#include "p3io.h"
#include "p3platform.h"

_arr_0P3PLATFORM P3PLATFORM_osfiletypetext = {{3,'W','I','N'}, {4,'U','N','I','X'}, {4,'X','X','X','X'}};
_arr_1P3PLATFORM P3PLATFORM_osplatformtext = {{3,'D','O','S'}, {5,'W','i','n','9','5'}, {5,'W','i','n','N','T'}, {8,'W','i','n','6','4','E','M','T'}, {3,'A','I','X'}, {5,'H','P','-','U','X'}, {6,'I','R','I','X','6','4'}, {5,'L','i','n','u','x'}, {13,'S','u','n','O','S','-','s','p','a','r','c','3','2'}, {4,'O','S','F','1'}, {10,'D','a','r','w','i','n','-','p','p','c'}, {10,'L','i','n','u','x','8','6','_','6','4'}, {11,'S','u','n','O','S','-','i','8','6','p','c'}, {13,'S','u','n','O','S','-','s','p','a','r','c','6','4'}, {11,'D','a','r','w','i','n','-','i','3','8','6'}, {10,'D','a','r','w','i','n','-','x','6','4'}, {8,'B','l','u','e','G','e','n','e'}, {7,'M','i','s','s','i','n','g'}};
_arr_2P3PLATFORM P3PLATFORM_osdllextension = {{4,'.','d','l','l'}, {4,'.','d','l','l'}, {4,'.','d','l','l'}, {4,'.','d','l','l'}, {3,'.','s','o'}, {3,'.','s','l'}, {3,'.','s','o'}, {3,'.','s','o'}, {3,'.','s','o'}, {3,'.','s','o'}, {6,'.','d','y','l','i','b'}, {3,'.','s','o'}, {3,'.','s','o'}, {3,'.','s','o'}, {6,'.','d','y','l','i','b'}, {6,'.','d','y','l','i','b'}, {3,'.','s','o'}, {4,'.','X','X','X'}};
_arr_3P3PLATFORM P3PLATFORM_osdllprefix = {{0}, {0}, {0}, {0}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}, {3,'l','i','b'}};
/**** C code included from p3platform.pas(127:1): 24 lines ****/
#if   defined(P3UNIX)
# include <sys/utsname.h>
#elif defined(P3DOS)
# include <windows.h>
typedef BOOL (WINAPI *LPFN_ISWOW64PROCESS) (HANDLE, PBOOL);
/* return 1 if we are a Win32 process running under Win64
 *        0 otherwise
 */
int isWow64(void)
{
  BOOL bIsWow64;
  LPFN_ISWOW64PROCESS fnIsWow64Process;

  fnIsWow64Process = (LPFN_ISWOW64PROCESS)
    GetProcAddress(GetModuleHandle("kernel32"),"IsWow64Process");
  if (NULL != fnIsWow64Process) {
    if (!fnIsWow64Process(GetCurrentProcess(),&bIsWow64)) {
      return 0;
    }
    return bIsWow64;
  }
  return 0;
} /* isWow64 */
#endif
static P3PLATFORM_tosfiletype P3PLATFORM_localosfiletype;
static P3PLATFORM_tosplatform P3PLATFORM_localosplatform;
static SYSTEM_shortstring P3PLATFORM_localosnullfilename;
static SYSTEM_shortstring P3PLATFORM_localosconsolename;
static SYSTEM_shortstring P3PLATFORM_localoslanguagepascal;
static SYSTEM_shortstring P3PLATFORM_localoslanguagec;
static SYSTEM_boolean P3PLATFORM_localislittleendian;

Function(P3PLATFORM_tosfiletype ) P3PLATFORM_osfiletype(void)
{
  P3PLATFORM_tosfiletype result;

  result = P3PLATFORM_localosfiletype;
  return result;
}  /* osfiletype */

Function(P3PLATFORM_tosplatform ) P3PLATFORM_osplatform(void)
{
  P3PLATFORM_tosplatform result;

  result = P3PLATFORM_localosplatform;
  return result;
}  /* osplatform */

Function(SYSTEM_ansichar *) P3PLATFORM_osnullfilename(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,P3PLATFORM_localosnullfilename);
  return result;
}  /* osnullfilename */

Function(SYSTEM_ansichar *) P3PLATFORM_osconsolename(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,P3PLATFORM_localosconsolename);
  return result;
}  /* osconsolename */

Function(SYSTEM_ansichar *) P3PLATFORM_oslanguagepascal(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,P3PLATFORM_localoslanguagepascal);
  return result;
}  /* oslanguagepascal */

Function(SYSTEM_ansichar *) P3PLATFORM_oslanguagec(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  _P3strcpy(result,_len_ret,P3PLATFORM_localoslanguagec);
  return result;
}  /* oslanguagec */

Function(SYSTEM_boolean ) P3PLATFORM_nativeislittleendian(void)
{
  SYSTEM_boolean result;

  result = P3PLATFORM_localislittleendian;
  return result;
}  /* nativeislittleendian */
typedef SYSTEM_uint8 _sub_5P3PLATFORM;
typedef SYSTEM_byte _arr_4P3PLATFORM[4];
typedef struct P3PLATFORM_tint32rec_S {
  union{
    struct{
      SYSTEM_integer i;
    } _c1;
    struct{
      _arr_4P3PLATFORM bytes;
    } _c2;
  } _u;
} P3PLATFORM_tint32rec;

static P3PLATFORM_tint32rec P3PLATFORM_i32rec;

/* unit p3platform */
void _Init_Module_p3platform(void)
{
  _P3strcpy(P3PLATFORM_localoslanguagepascal,255,_P3str1("\004P3PC"));
  _P3strcpy(P3PLATFORM_localoslanguagec,255,_P3str1("\007Unknown"));
  /**** C code included from p3platform.pas(263:1): 17 lines ****/
#define localosfiletype P3PLATFORM_localosfiletype
#if   defined(P3DOS)
  localosfiletype = P3PLATFORM_osfilewin;
#elif defined(P3CMS)
  localosfiletype = P3PLATFORM_osfilecms;
#elif defined(P3TSO)
  localosfiletype = P3PLATFORM_osfiletso;
#elif defined(P3VMS)
  localosfiletype = P3PLATFORM_osfilevms;
#elif defined(P3UNIX) /* NOT: || defined(DJGPP) comes below */
  localosfiletype = P3PLATFORM_osfileunix;
#elif defined(P3MAC)
  localosfiletype = P3PLATFORM_osfilemac;
#else
# error "P3<OSTYPE> not defined for any recognized OS type"
  P3OS_TYPE_not_defined_for_recognized_OS_type__ERROR;
#endif
  P3PLATFORM_localosplatform = P3PLATFORM_osmissing;
  /**** C code included from p3platform.pas(285:1): 62 lines ****/
#define localosplatform P3PLATFORM_localosplatform

#if defined(P3DOS) /* windows */
{
# if defined(_WIN64)
  SYSTEM_INFO siSysInfo;
  GetSystemInfo(&siSysInfo);
  /* We assume no Itanium builds: what else could wProcessorArchitecture be? */
  if (PROCESSOR_ARCHITECTURE_AMD64==siSysInfo.wProcessorArchitecture)
    localosplatform  = P3PLATFORM_oswindows64emt;
  else
    localosplatform  = P3PLATFORM_osmissing;

# else
  DWORD winVersion = GetVersion();
  if ((winVersion >> 31) & 1) /* Win95 */
    localosplatform  = P3PLATFORM_oswindows95;
  else if (isWow64())
    localosplatform  = P3PLATFORM_oswindows64emt;
  else
    localosplatform  = P3PLATFORM_oswindowsnt;
# endif /* #if defined(_WIN64) */
}
#elif defined(P3UNIX)
{
  struct utsname uts;  int err, len;

  localosplatform  = P3PLATFORM_osmissing;
  if ((err = uname(&uts)) < 0)
    localosplatform  = P3PLATFORM_osmissing;
  else {
    if      (0==strcmp(uts.sysname,"AIX"))
      localosplatform  = P3PLATFORM_osaix;
    else if (0==strcmp(uts.sysname,"Linux")) {
      if (0==strcmp(uts.machine,"x86_64"))
        localosplatform  = P3PLATFORM_oslinux86_64;
      else if ((0==strcmp(uts.machine,"BGP"))
            || (0==strcmp(uts.machine,"ppc64")))
        localosplatform  = P3PLATFORM_osbluegene;
      else
        localosplatform  = P3PLATFORM_oslinux;
    }
    else if (0==strcmp(uts.sysname,"SunOS")) {
      if (0==strcmp(uts.machine,"i86pc"))
        localosplatform  = P3PLATFORM_ossunos_i86pc;
      else
        /* SPD, Oct 2014: we no longer build for sparc32 */
        localosplatform  = P3PLATFORM_ossunos_sparc64;
    }
    else if (0==strcmp(uts.sysname,"Darwin")) {
      if (0==strcmp(uts.machine,"i386"))
        localosplatform  = P3PLATFORM_osdarwin_i386;
      else
        localosplatform  = P3PLATFORM_osdarwin_x64;
    }
    else
      localosplatform  = P3PLATFORM_osmissing;
  }
}
#else
#error "localOSPlatform not yet defined"
#endif
  switch (P3PLATFORM_localosfiletype) {
    case P3PLATFORM_osfilewin: 
      _P3strcpy(P3PLATFORM_localosnullfilename,255,_P3str1("\003nul"));
      _P3strcpy(P3PLATFORM_localosconsolename,255,_P3str1("\003con"));
      break;
    case P3PLATFORM_osfileunix: 
      _P3strcpy(P3PLATFORM_localosnullfilename,255,_P3str1("\011/dev/null"));
      _P3strcpy(P3PLATFORM_localosconsolename,255,_P3str1("\010/dev/tty"));
      break;
    case P3PLATFORM_osfilemissing: 
      _P3strclr(P3PLATFORM_localosnullfilename);
      _P3strclr(P3PLATFORM_localosconsolename);
      break;
    default: break;
  }
  P3PLATFORM_i32rec._u._c1.i = 1;
  P3PLATFORM_localislittleendian = P3PLATFORM_i32rec._u._c2.bytes[0] == 1;
} /* _Init_Module_p3platform */

void _Final_Module_p3platform(void)
{
} /* _Final_Module_p3platform */

