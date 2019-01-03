#include "p3io.h"
#include "system_p3.h"
#include "p3platform.h"
#include "p3private.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "p3process.h"

#include "globals2.h"

void * const P3PROCESS_texecarglist_VT[] = {(void*)&
  P3PROCESS_texecarglist_DOT_destroy};

/* Class descriptor for 'texecarglist' */
const SYSTEM_classdescriptor_t P3PROCESS_texecarglist_CD = {
  _P3str1("\014texecarglist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(P3PROCESS_texecarglist_OD), P3PROCESS_texecarglist_VT, NULL};

static _P3STR_31 P3PROCESS_cmd_win7 = {27,'c',':','\\','w','i','n','d','o','w','s','\\','s','y','s','t','e','m','3','2','\\','c','m','d','.','e','x','e'};
static _P3STR_31 P3PROCESS_cmd_winnt = {25,'c',':','\\','w','i','n','n','t','\\','s','y','s','t','e','m','3','2','\\','c','m','d','.','e','x','e'};
/**** C code included from p3process.pas(180:1): 62 lines ****/
#if defined(_WIN32)
# include <tlhelp32.h>
/* turn off some bits when calling OpenProcess: XP cannot handle them being on
 * and we do not need them on in Windows versions newer than XP
 */
# define PROCESS_HACK_ACCESS (PROCESS_ALL_ACCESS & ~0x0000f000)
#else
# include <sys/types.h>
# include <sys/stat.h>
# include <sys/wait.h>
# include <signal.h>
# include <fcntl.h>
# if defined(__HOS_AIX__)
#  include <procinfo.h>
# elif defined(__APPLE__)
#  include <sys/sysctl.h>
# elif defined(__sparc)
/* use routines from p3Custom.c */
# elif (defined(__sun__) && defined(__amd64__))
#  include <procfs.h>
#  include <sys/proc.h>
# endif
#endif

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
     } while ((p >= buf) && (('.' == *p) || (33 > *p))
           );

  return buf;
} /* winErrMsg */
#endif


static Function(SYSTEM_integer ) asyncsystem4unix(
  SYSTEM_P3_pansichar cmdptr,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg);

static Function(SYSTEM_integer ) asyncsystem4win(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean newconsole,
  SYSTEM_boolean inheritedhandles,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg);
typedef SYSTEM_uint16 _sub_0P3PROCESS;
typedef SYSTEM_P3_pansichar P3PROCESS_targv[1001];
typedef P3PROCESS_targv *P3PROCESS_tpargv;
static SYSTEM_integer P3PROCESS_wshowwindow;

static Function(SYSTEM_ansichar *) P3PROCESS_whatquote(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *arg)
{
  SYSTEM_integer i;
  SYSTEM_shortstring targ;

  _P3strclr(result);
  SYSUTILS_P3_trim(targ,255,arg);
  if (SYSTEM_length(targ) > 1 && targ[1] == _P3char('\"') && 
    targ[SYSTEM_length(targ)] == _P3char('\"')) 
    return result;
  if (SYSTEM_length(arg) == 0) {
    _P3strcpy(result,_len_ret,_P3str1("\001\""));
    return result;
  } 
  { register SYSTEM_int32 _stop = SYSTEM_length(arg);
    if ((i = 1) <=  _stop) do {
      if (arg[i] <= _P3char(' ')) {
        _P3strcpy(result,_len_ret,_P3str1("\001\""));
        return result;
      } 
    } while (i++ !=  _stop);

  }
  return result;
}  /* whatquote */
/**** C code included from p3process.pas(433:1): 296 lines ****/
#define wshowwindow P3PROCESS_wshowwindow
int
Win32CreateProc(const char *exeName, char *cmdLine,
                int inheritedHandles, int *exeRC)
#if ! defined(_WIN32)
{ *exeRC = *exeName + *cmdLine;  return 1; } /*  bogus definition */
#else
{
  PROCESS_INFORMATION processinformation;
  STARTUPINFO         startupinfo;
  DWORD               exitcode;
  int                 result;
  BOOL brc;

  /* Initialise the startup information to be the same as that of the
   * calling application.  This is easier than initialising the many
   * individual startup information fields and should be fine in most
   * cases. */
  GetStartupInfo(&startupinfo);

  /* StartupInfo.wShowWindow determines whether the called application
   * will be initially displayed normal, maximises, minimised or some
   * other subtle variations */

  startupinfo.wShowWindow = wshowwindow;

  if (! CreateProcess(
    exeName,               /* ApplicationName */
    cmdLine,               /* lpCommandLine */
    NULL,                  /* lpProcessAttributes */
    NULL,                  /* lpThreadAttribute */
    inheritedHandles,      /* bInheritedHandles */
    0,                     /* dwCreationFlags */
    NULL,                  /* lpEnvironment */
    NULL,                  /* lpCurrentDirectory */
    &startupinfo,          /* lpStartupInfo */
    &processinformation    /* lpProcessInformation */
    )) {
    *exeRC = 0;
    result = GetLastError();  /* failed to execute */
  }
  else {
    WaitForSingleObject(processinformation.hProcess,INFINITE);
    brc = GetExitCodeProcess(processinformation.hProcess,&exitcode);
    CloseHandle(processinformation.hThread);
    CloseHandle(processinformation.hProcess);
    if (! brc         /* failed call to GetExitCodeProcess */
#if 0
       || (exitcode >> 31)
#endif
       || (255 == exitcode) ) {
      *exeRC = 0;
      return 1;
    }
    *exeRC = exitcode;
    result = 0;
  }
  return result;
} /* Win32CreateProc */
#endif /* if ! defined(_WIN32) .. else .. */

int
win32ASyncCreateProc (const char *exeName, char *cmdLine, int newConsole,
                      int inheritedHandles, P3PROCESS_pprocinfo procInfo)
#if ! defined(_WIN32)
{
  procInfo->pid = ~0;
  return 1; /* failure */
}
#else
{
  int brc, nc;
  PROCESS_INFORMATION processInformation;
  STARTUPINFO startupInfo;

  procInfo->pid = ~0;

  /* Initialize the startup information to be the same as that of the
   * parent.  This is easier than initializing the many
   * individual startup information fields and should be fine in most
   * cases. */
  GetStartupInfo (&startupInfo);

  if (newConsole)
  {

     /* This new settings allows us to send a CtrlC via the new "GAMS Message Interrupt" to a
        particular GAMS job without disturbing the parent and vice versa. Along this, we
        allowed the new async processes to gets their own stdin/out/err. */

     startupInfo.dwFlags |=  STARTF_USESHOWWINDOW;
     startupInfo.dwFlags &= ~STARTF_USESTDHANDLES;
     startupInfo.wShowWindow = SW_MINIMIZE;
     inheritedHandles = FALSE;
     nc = CREATE_NEW_CONSOLE;
  }
  else
  {
     startupInfo.wShowWindow = wshowwindow;    /* get the locally stored showWindow preference */
     nc = 0;
  }

  if (! CreateProcess(
    exeName,               /* ApplicationName */
    cmdLine,               /* lpCommandLine */
    NULL,                  /* lpProcessAttributes */
    NULL,                  /* lpThreadAttribute */
    inheritedHandles,      /* bInheritedHandles */
    nc,                    /* dwCreationFlags */
    NULL,                  /* lpEnvironment */
    NULL,                  /* lpCurrentDirectory */
    &startupInfo,          /* lpStartupInfo */
    &processInformation    /* lpProcessInformation */
    )) {
    return GetLastError();  /* failed to execute */
  }

  /* child is running now - just clean up and return the process info */
  procInfo->pid = (SYSTEM_cardinal) processInformation.dwProcessId;
  procInfo->tid = (SYSTEM_cardinal) processInformation.dwThreadId;
  procInfo->hprocess = (SYSTEM_nativeuint) processInformation.hProcess;

  CloseHandle (processInformation.hThread);
  /* CloseHandle (processInformation.hProcess); */
  return 0;
}
#endif /* if ! defined(_WIN32) .. else .. */

/* unixPidStatus returns:
 *   0: valid process but not a zombie
 *   1: zombie process
 *   2: process does not exist
 *   3: not implemented or other error
 */
int unixPidStatus (int p)
{
#if defined(_WIN32)
  return 3;     /* not defined for Windows */
#else

  if (p <= 0)
    return 2;

# if defined(__HOS_AIX__)
{
  int p2;
  int rc;
  struct procentry64 pe;

  p2 = p;
  rc = getprocs64 (&pe, sizeof(pe), NULL, 0, &p2, 1);
  if (-1 == rc)
    return 2;        /* not a process */
  else if (1 != rc)
    return 3;        /* error */
  else if (pe.pi_pid != p)
    return 2;        /* not a process */
  if (SZOMB == pe.pi_state)
    return 1;        /* zombie process */

  return 0;          /* valid, non-zombie process */
}

# elif defined(__linux__)
{
  char sbuf[1024];              /* buffer for content of stat file */
  char *tmp, *t2;
  char filename[80];
  int fd, p2;
  ssize_t numRead;
  char state;
  struct stat sb;               /* stat() buffer */

  sprintf (filename, "/proc/%d", p);
  if (-1 == stat(filename, &sb))
    return 2;

  sprintf (filename, "/proc/%d/stat", p);
  fd = open (filename, O_RDONLY, 0);
  if (-1 == fd)
    return 3;
  numRead = read (fd, sbuf, sizeof(sbuf) - 1);
  close (fd);
  if (numRead <= 0)
    return 3;

  /* start of file looks like
   * 2796 (firefox) S  where we have
   *  pid  cmdline  status */
  sbuf[numRead] = '\0';
  numRead = sscanf (sbuf, "%d", &p2);
  if (1 != numRead)
    return 3;
  tmp = strchr(sbuf, '(') + 1;
  t2 = strrchr(sbuf, ')');
  if (NULL == tmp || NULL == t2)
    return 3;
  *t2 = '\0';
  t2 = t2 + 2;                 // skip ") "
  numRead = sscanf (t2, "%c", &state);
  if (1 != numRead)
    return 3;
  // printf ("DEBUG: in isZombie: read %d (%s) %c\n", p2, tmp, state);
  switch (state) {
  case 'D':                     /* uninterruptible sleep */
  case 'R':                     /* running */
  case 'S':                     /* sleeping */
  case 'T':                     /* traced or stopped */
    return 0;                   /* valid, non-zombie process */
    break;
  case 'Z':                     /* zombie */
    return 1;
    break;
  default:
    return 3;
  }
}
# elif defined(__APPLE__)
{
  int mib[4], rc;
  size_t len;
  struct kinfo_proc kp;

  mib[0] = CTL_KERN;
  mib[1] = KERN_PROC;
  mib[2] = KERN_PROC_PID;
  mib[3] = p;

  len = (int) sizeof(kp);
  /* printf ("sizeof(kp) = %d\n", len); */
  rc = sysctl(mib, 4, &kp, &len, NULL, 0);
  if (0 != rc)
    return 3;                   /* error */
  if (0 == len)
    return 2;                   /* not a process */
  if ((int)sizeof(kp) != len)
    return 3;                   /* error */

  if (SZOMB == kp.kp_proc.p_stat)
    return 1;                   /* this is a zombie */

  return 0;          /* valid, non-zombie process */
} /* Sun Intel */

# elif defined(__sparc)
{
#  include "p3Custom.h"
  return sparcPidStatus ((pid_t) p);
}

# elif (defined(__sun__) && defined(__amd64__))
{
  psinfo_t psinfo;
  char filename[80];
  int fd;
  ssize_t numRead;
  struct stat sb;               /* stat() buffer */

  sprintf (filename, "/proc/%d", p);
  if (-1 == stat(filename, &sb))
    return 2;

  sprintf (filename, "/proc/%d/psinfo", p);
  fd = open (filename, O_RDONLY, 0);
  if (-1 == fd)
    return 3;
  numRead = read (fd, &psinfo, sizeof(psinfo));
  close (fd);
  if (numRead <= 0)
    return 3;

#  if 0
  printf ("DEBUG: pid = %d   uid = %d   gid = %d\n",
          (int) psinfo.pr_pid, (int) psinfo.pr_uid, (int) psinfo.pr_gid);
  printf ("DEBUG: in isZombie: read %d (%s) %d\n",
          psinfo.pr_pid, psinfo.pr_fname, psinfo.pr_flag);
  printf ("DEBUG: in isZombie: wstat = %d\n",
          psinfo.pr_wstat);
  printf ("DEBUG: in isZombie: %d %c\n",
          psinfo.pr_lwp.pr_state, psinfo.pr_lwp.pr_sname);
#  endif

  if (SZOMB == psinfo.pr_lwp.pr_state)
    return 1;                   /* this is a zombie */

  return 0;          /* valid, non-zombie process */
}
# else

# endif         /* various non-Windows flavors */

  return 3;     /* not implemented */

#endif          /* #if defined(_WIN32) .. #else .. */

}
/**** C code included from p3process.pas(786:1): 98 lines ****/
int
LibcForkExec(int argc, char *const argv[], int *exeRC)
#if defined(_WIN32)
{ *exeRC = argc + *(argv[0]);  return 1; } /* bogus definition */
#else
{
  int result = 1;
  int pid, pid2;
  int wstat;
  /* */
  pid = fork();
  if (pid < 0) {                /* could not fork */
    *exeRC = 0;
    result = 1;
  }
  else if (0 == pid) {          /* I am the child */
    execvp (argv[0], argv);
    execl("/bin/sh", "/bin/sh", "-c", "exit 255", NULL);
    /* if we are here, it is trouble */
#if 0 /* do not do this, exit() flushes stdio buffers of the parent */
    exit (255);                 /* -1 tells parent we could not exec */
#else
    /* _exit() is a more immediate termination, less likely to flush stdio */
    /* _exit (255); */                /* -1 tells parent we could not exec */
    exit2R("Failed exec after fork");
#endif
  }
  else {                        /* I am the parent */
    for ( ; ; ) {
      wstat = 0;
      pid2 = waitpid (pid, &wstat, 0);
      if (-1 == pid2)
        continue;
      if (pid != pid2) {
        *exeRC = 0;
        return 1;    /* failed wait so fail the entire fork/exec */
      }
      else
        break;
    }
    if (WIFEXITED(wstat)) {     /* normal exit from child */
      if (255 == WEXITSTATUS(wstat)) { /* because it couldn't exec */
        *exeRC = 0;
        result = 1;
      }
      else {
        *exeRC = WEXITSTATUS(wstat);
        result = 0;
      }
    }
    else {                      /* abnormal return from child */
      *exeRC = 0;
      result = 1;
    }
  } /* end parent code */
  return result;
} /* LibcForkExec */
#endif /* #if defined(_WIN32) .. else .. */

/* libcASyncForkExec does a fork/exec to start a process,
 * but it does not wait for it.  Instead it returns the PID.
 * also, it sets up a new process group for the child
 * result: 0 on success, ~0 on failure
 */
int
libcASyncForkExec (int argc, char *const argv[], SYSTEM_cardinal *pid)
#if defined(_WIN32)
{
  *pid = 0;
  return 1;
}
#else
{
  int result = 1;  /* failure */
  int lPid;        /* local pid, typed appropriately */
  int wstat;

  *pid = ~0;
  lPid = fork();
  if (lPid < 0) {               /* could not fork */
    result = 1;
  }
  else if (0 == lPid) {         /* I am the child */
    (void) setpgid (0,0);       /* make this process a new process group */
    execvp (argv[0], argv);

    /* if we are here, it is trouble */
    execl("/bin/sh", "/bin/sh", "-c", "exit 127", NULL);
    exit2R("Failed exec after fork");
    /* _exit (127); */                /* consistent with & usage in bash */
  }
  else {                        /* I am the parent */
    (void) setpgid (lPid,0);     /* make the child its own, new process group */
    /* we call setpgid for both parent and child to avoid a race condition */
    result = 0;
    *pid = (SYSTEM_cardinal) lPid;
  }
  return result;
} /* libcASyncForkExec */
#endif /* #if defined(_WIN32) .. else .. */

static Function(SYSTEM_P3_pansichar ) P3PROCESS_getparamshortstr(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *param);

static Procedure pushchar(
  SYSTEM_ansichar c,
  SYSTEM_P3_pansichar *_2r,
  SYSTEM_integer *_2len)
{
  if (*_2len < 255) {
    **_2r = c;
    _P3inc0(*_2len);
    _P3inc0(*_2r);
  } 
}  /* pushchar */

static Function(SYSTEM_P3_pansichar ) P3PROCESS_getparamshortstr(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *param)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer len;
  SYSTEM_P3_pansichar s;
  SYSTEM_P3_pansichar r;

  while (SYSTEM_true) {
    while (*p != _P3char('\000') && *p <= _P3char(' ')) {

      _P3inc0(p);
}
    s = p;
    _P3inc0(s);
    if (*p == _P3char('\"') && *s == _P3char('\"')) { 
      _P3inc1(p,2);
    } else 
      SYSTEM_break(BRK_1);
  
CNT_1:;
  }
BRK_1:;
  len = 0;
  r = ValueCast(SYSTEM_P3_pansichar,&param[1]);
  while (*p > _P3char(' ')) {

    if (*p == _P3char('\"')) {
      _P3inc0(p);
      while (*p != _P3char('\000') && *p != _P3char('\"')) {
        pushchar(*p,&r,&len);
        _P3inc0(p);
      
}
      if (*p != _P3char('\000')) 
        _P3inc0(p);
    } else {
      pushchar(*p,&r,&len);
      _P3inc0(p);
    } 
}
  _P3setlength(param,len,255);
  result = p;
  return result;
}  /* getparamshortstr */

Function(SYSTEM_integer ) P3PROCESS_p3execp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_integer argc, i;
  P3PROCESS_tpargv pargv;
  SYSTEM_P3_pansichar s;
  SYSTEM_shortstring param;

  result = 1;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      /**** C code included from p3process.pas(966:1): 1 lines ****/
      result = Win32CreateProc (NULL, (char *) cmdptr, 1, progrc);
      break;
    case P3PLATFORM_osfileunix: 
      argc =  -1;
      s = cmdptr;
      do {
        s = P3PROCESS_getparamshortstr(s,param);
        _P3inc0(argc);
      } while (!_P3strcmpE(param,_P3str1("\000")));
      if (argc == 0) {
        *progrc = 0;
        result = 1;
        return result;
      } 
      _P3getmem(pargv,(argc + 1) * sizeof(SYSTEM_pointer));
      s = cmdptr;
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 0) <=  _stop) do {
          s = P3PROCESS_getparamshortstr(s,param);
          (*pargv)[i] = P3PRIVATE_strtopchar(param);
          SYSTEM_assert(_P3strcmpN(param,_P3str1("\000")),_P3str1("\052cmd string should not be out of parameters"));
        
        } while (i++ !=  _stop);

      }
      P3PROCESS_getparamshortstr(s,param);
      SYSTEM_assert(_P3strcmpE(param,_P3str1("\000")),_P3str1("\036cmd string should be exhausted"));
      (*pargv)[argc] = NULL;
      /**** C code included from p3process.pas(1004:1): 1 lines ****/
      result =  LibcForkExec (argc, (char *const * )( *pargv), progrc);
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 0) <=  _stop) do {
          _P3freemem((*pargv)[i]);
        } while (i++ !=  _stop);

      }
      _P3freemem(pargv);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\044unimplemented P3ExecP for OSFileType"));
  }
  return result;
}  /* p3execp */

Function(SYSTEM_integer ) P3PROCESS_p3exec(
  const SYSTEM_ansichar *cmd,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_P3_pansichar cmdptr;
  P3PRIVATE_shortstrbuf cmdbuf;

  cmdptr = P3PRIVATE_strtostrbuf(cmd,cmdbuf);
  result = P3PROCESS_p3execp(cmdptr,progrc);
  return result;
}  /* p3exec */

Function(SYSTEM_integer ) P3PROCESS_p3exec2(
  const SYSTEM_ansichar *progname,
  const SYSTEM_ansichar *progparams,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_P3_pansichar prognameptr;
  P3PRIVATE_shortstrbuf prognamebuf;
  SYSTEM_P3_pansichar progparamsptr;
  P3PRIVATE_shortstrbuf progparamsbuf;
  SYSTEM_P3_pansichar cmdptr;
  SYSTEM_integer cmdlen;
  SYSTEM_integer argc, i, k;
  P3PROCESS_tpargv pargv;
  SYSTEM_P3_pansichar s;
  SYSTEM_shortstring param;
  SYSTEM_shortstring quote;

  result = 1;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      cmdlen = ValueCast(SYSTEM_int32,SYSTEM_length(progname)) + 3 + 
        SYSTEM_length(progparams) + 1;
      _P3getmem(cmdptr,cmdlen);
      k = 0;
      P3PROCESS_whatquote(quote,255,progname);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,progname);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,_P3str1("\001 "));
      P3PRIVATE_pcharconcatstr(cmdptr,&k,progparams);
      /**** C code included from p3process.pas(1072:1): 1 lines ****/
      result = Win32CreateProc (NULL, (char *) cmdptr, 1, progrc);
      _P3freemem(cmdptr);
      break;
    case P3PLATFORM_osfileunix: 
      prognameptr = P3PRIVATE_strtostrbuf(progname,prognamebuf);
      progparamsptr = P3PRIVATE_strtostrbuf(progparams,progparamsbuf);
      if (*prognameptr == _P3char('\000')) {
        *progrc = 0;
        result = 1;
        return result;
      } 
      argc = 0;
      s = progparamsptr;
      do {
        s = P3PROCESS_getparamshortstr(s,param);
        _P3inc0(argc);
      } while (!_P3strcmpE(param,_P3str1("\000")));
      _P3getmem(pargv,(argc + 1) * sizeof(SYSTEM_pointer));
      (*pargv)[0] = prognameptr;
      s = progparamsptr;
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 1) <=  _stop) do {
          s = P3PROCESS_getparamshortstr(s,param);
          (*pargv)[i] = P3PRIVATE_strtopchar(param);
          SYSTEM_assert(_P3strcmpN(param,_P3str1("\000")),_P3str1("\055params string should not be out of parameters"));
        
        } while (i++ !=  _stop);

      }
      P3PROCESS_getparamshortstr(s,param);
      SYSTEM_assert(_P3strcmpE(param,_P3str1("\000")),_P3str1("\041params string should be exhausted"));
      (*pargv)[argc] = NULL;
      /**** C code included from p3process.pas(1115:1): 1 lines ****/
      result =  LibcForkExec (argc, (char *const * )( *pargv), progrc);
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 1) <=  _stop) do {
          _P3freemem((*pargv)[i]);
        } while (i++ !=  _stop);

      }
      _P3freemem(pargv);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\044unimplemented P3Exec2 for OSFileType"));
  }
  return result;
}  /* p3exec2 */

Function(SYSTEM_integer ) P3PROCESS_p3execl(
  const SYSTEM_ansichar *progname,
  P3PROCESS_texecarglist progparams,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_boolean inheritedhandles;
  SYSTEM_integer i, k;
  P3PROCESS_tpargv pargv;
  SYSTEM_integer argc;
  SYSTEM_integer cmdlen;
  SYSTEM_P3_pansichar cmdptr;
  SYSTEM_shortstring quote;

  inheritedhandles = progparams->
    P3PROCESS_texecarglist_DOT_finheritedhandles;
  argc = 0;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      cmdlen = ValueCast(SYSTEM_int32,SYSTEM_length(progname)) + 3;
      { register SYSTEM_int32 _stop = progparams->
          P3PROCESS_texecarglist_DOT_fcount - 1;
        if ((i = 0) <=  _stop) do {
          {
            SYSTEM_shortstring _t1;

            cmdlen = cmdlen + SYSTEM_length(
              P3PROCESS_texecarglist_DOT_get(_t1,255,progparams,i)) + 3;
          }
        } while (i++ !=  _stop);

      }
      _P3getmem(cmdptr,cmdlen);
      k = 0;
      P3PROCESS_whatquote(quote,255,progname);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,progname);
      P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
      { register SYSTEM_int32 _stop = progparams->
          P3PROCESS_texecarglist_DOT_fcount - 1;
        if ((i = 0) <=  _stop) do {
          {
            SYSTEM_shortstring _t2;

            P3PROCESS_whatquote(quote,255,
              P3PROCESS_texecarglist_DOT_get(_t2,255,progparams,i));
          }
          P3PRIVATE_pcharconcatstr(cmdptr,&k,_P3str1("\001 "));
          P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
          {
            SYSTEM_shortstring _t1;

            P3PRIVATE_pcharconcatstr(cmdptr,&k,
              P3PROCESS_texecarglist_DOT_get(_t1,255,progparams,i));
          }
          P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
        
        } while (i++ !=  _stop);

      }
      break;
    case P3PLATFORM_osfileunix: 
      argc = 1 + progparams->P3PROCESS_texecarglist_DOT_fcount;
      _P3getmem(pargv,(argc + 1) * sizeof(SYSTEM_pointer));
      (*pargv)[0] = P3PRIVATE_strtopchar(progname);
      { register SYSTEM_int32 _stop = progparams->
          P3PROCESS_texecarglist_DOT_fcount - 1;
        if ((i = 0) <=  _stop) do {
          {
            SYSTEM_shortstring _t1;

            (*pargv)[i + 1] = P3PRIVATE_strtopchar(
              P3PROCESS_texecarglist_DOT_get(_t1,255,progparams,i));
          }
        } while (i++ !=  _stop);

      }
      (*pargv)[argc] = NULL;
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\044unimplemented P3Execl for OSFileType"));
  }
  /**** C code included from p3process.pas(1251:1): 5 lines ****/
#if defined(_WIN32)
  result = Win32CreateProc(NULL, (char *)cmdptr, inheritedhandles, progrc);
#else
  result = LibcForkExec(argc, (char *const * )( *pargv), progrc);
#endif /* #if defined(_WIN32) .. else .. */
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      _P3freemem(cmdptr);
      break;
    case P3PLATFORM_osfileunix: 
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 0) <=  _stop) do {
          _P3freemem((*pargv)[i]);
        } while (i++ !=  _stop);

      }
      _P3freemem(pargv);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\044unimplemented P3Execl for OSFileType"));
  }
  return result;
}  /* p3execl */

Function(SYSTEM_integer ) P3PROCESS_p3asyncexecp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean newconsole,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;
  SYSTEM_integer argc, i;
  P3PROCESS_tpargv pargv;
  SYSTEM_P3_pansichar s;
  SYSTEM_shortstring param;
  SYSTEM_cardinal pid;

  SYSTEM_P3_fillchar(procinfo,sizeof(P3PROCESS_tprocinfo),0);
  result = 1;
  _P3strclr(msg);
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      /**** C code included from p3process.pas(1305:1): 1 lines ****/
   result = win32ASyncCreateProc (NULL, (char *) cmdptr, newconsole, 1, procinfo);
      break;
    case P3PLATFORM_osfileunix: 
      argc =  -1;
      s = cmdptr;
      do {
        s = P3PROCESS_getparamshortstr(s,param);
        _P3inc0(argc);
      } while (!_P3strcmpE(param,_P3str1("\000")));
      if (argc == 0) 
        return result;
      _P3getmem(pargv,(argc + 1) * sizeof(SYSTEM_pointer));
      s = cmdptr;
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 0) <=  _stop) do {
          s = P3PROCESS_getparamshortstr(s,param);
          (*pargv)[i] = P3PRIVATE_strtopchar(param);
          SYSTEM_assert(_P3strcmpN(param,_P3str1("\000")),_P3str1("\052cmd string should not be out of parameters"));
        
        } while (i++ !=  _stop);

      }
      P3PROCESS_getparamshortstr(s,param);
      SYSTEM_assert(_P3strcmpE(param,_P3str1("\000")),_P3str1("\036cmd string should be exhausted"));
      (*pargv)[argc] = NULL;
      /**** C code included from p3process.pas(1333:1): 1 lines ****/
   result = libcASyncForkExec (argc, (char *const * )( *pargv), &pid);
      procinfo->pid = pid;
      { register SYSTEM_int32 _stop = argc - 1;
        if ((i = 0) <=  _stop) do {
          _P3freemem((*pargv)[i]);
        } while (i++ !=  _stop);

      }
      _P3freemem(pargv);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\051unimplemented p3ASyncExecP for OSFileType"));
  }
  return result;
}  /* p3asyncexecp */

Function(SYSTEM_integer ) P3PROCESS_p3asyncsystemp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean newconsole,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;

  _P3strclr(msg);
  result = 1;
  SYSTEM_P3_fillchar(procinfo,sizeof(P3PROCESS_tprocinfo),0);
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      result = asyncsystem4win(cmdptr,newconsole,SYSTEM_true,procinfo,
        msg);
      break;
    case P3PLATFORM_osfileunix: 
      result = asyncsystem4unix(cmdptr,procinfo,msg);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\053unimplemented p3ASyncSystemP for OSFileType"));
  }
  return result;
}  /* p3asyncsystemp */

Function(SYSTEM_integer ) P3PROCESS_p3asyncstatus(
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_integer *progrc,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;

  result = 0;
  _P3strclr(msg);
  if (0 == procinfo->pid) {
    _P3strcpy(msg,255,_P3str1("\013Invalid PID"));
    return result;
  } 
  /**** C code included from p3process.pas(1402:1): 101 lines ****/
{
#if defined(_WIN32)
  HANDLE h;
  DWORD p, rc, exitCode;
  char ebuf[256];

#if 0
  skip this stuff if we cannot get a thread ID from a process ID
  if ((0 == procinfo->tid) ^ (0 == procinfo->hprocess)) {
    _P3strcpy(msg,255,_P3str1("\031Corrupt or bogus procInfo"));
    return 0;
  }
#endif
  p = (DWORD) procinfo->pid;
  h = (HANDLE) procinfo->hprocess;
  if (NULL == h) {
    h = OpenProcess (PROCESS_HACK_ACCESS, FALSE, p);
    if (NULL == h) {
      rc = GetLastError();
      switch (rc) {
        case ERROR_INVALID_PARAMETER:
          /* system process but that is pid=0, we checked for that already */
          /* or expired or invalid PID */
          *msg = '\0';
          return 4; /* no such process */
          break;
        case ERROR_ACCESS_DENIED:
        default:
          (void) winErrMsg (rc, ebuf, sizeof(ebuf));
          _P3_pchar2str (msg, 255, (SYSTEM_char *)ebuf);
      } /* end switch */
      return 0;
    }
    procinfo->hprocess = (SYSTEM_nativeuint) h;
    /* we have no way to get a thread ID, given a process ID */
    /* procinfo->tid = 0; */
  }

  /* assume we have a PID and handle now */
  rc = WaitForSingleObject (h, 0);
  switch (rc) {
    case WAIT_OBJECT_0: /* signalled/completed */
      rc = GetExitCodeProcess (h, &exitCode);
      if (rc && (0xffffffff != exitCode)) {
        *progrc = (SYSTEM_integer) exitCode;
        result = 2;
      }
      else {
        result = 3;
      }
      CloseHandle (h);
      procinfo->hprocess = 0;
      return result;
      break;
    case WAIT_TIMEOUT: /* still running normally */
      return 1;
      break;
    default:
      _P3conp2str(msg,255,"Unexpected return from wait");
      return 0;
  } /* end switch */
  return result; /* should never get here */
#else
  pid_t pid, p2;
  int wstat;

  pid = (pid_t) procinfo->pid;
  if (pid <= 0) {                  /* PIDs are positive */
    _P3strcpy(msg,255,_P3str1("\013Invalid PID"));
    return 0;
  }
  if ((0 != procinfo->tid) || (0 != procinfo->hprocess)) {  /* we only use/set the pid on non-windows  */
    _P3strcpy(msg,255,_P3str1("\031Corrupt or bogus procInfo"));
    return 0;
  }
  p2 = waitpid (pid, &wstat, WNOHANG);
  if (pid == p2) { /* process p has changed state - assume it was to exit */
    /* consider using waitid() instead of waitpid() to get
     * "more precise control over which child state changes to wait for" */
    if (! WIFEXITED(wstat)) { /* no exit code is available */
      return 3;
    }
    *progrc = WEXITSTATUS(wstat);
    if (127 == *progrc) {   /* return for fork & failed exec */
      return 127;
    }
    else
      return 2;             /* we really have something to return */
  }
  else if (-1 == p2) {  /* error, e.g. no such process or not a child */
    _P3strcpy(msg,255,_P3str1("\036No such process or not a child"));
    return 4;
  }
  else if (0 == p2) {  /* child exists but has not exited */
    return 1;
  }

  _P3strcpy(msg,255,_P3str1("\033Unexpected return from wait"));
  return 0;
#endif
}
  return result;
}  /* p3asyncstatus */

Function(SYSTEM_cardinal ) P3PROCESS_p3getpid(void)
{
  SYSTEM_cardinal result;

  /**** C code included from p3process.pas(1578:1): 5 lines ****/
#if defined(_WIN32)
  result = (SYSTEM_cardinal) GetCurrentProcessId();
#else
  result = (SYSTEM_cardinal) getpid();
#endif
  return result;
}  /* p3getpid */

Function(SYSTEM_boolean ) P3PROCESS_p3ispidvalid(
  SYSTEM_cardinal pid)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  /**** C code included from p3process.pas(1600:1): 18 lines ****/
{
#if defined(_WIN32)
  HANDLE hProcess;
  DWORD p;

  p = (DWORD) pid;
  hProcess = OpenProcess (PROCESS_HACK_ACCESS, FALSE, p);
  if (NULL != hProcess) {
    CloseHandle (hProcess);
    result = SYSTEM_true;
  }
#else
  int rc;

  rc = unixPidStatus ((int) pid);
  result = (0 == rc) || (1 == rc);
#endif
}
  return result;
}  /* p3ispidvalid */

Function(SYSTEM_boolean ) P3PROCESS_p3ispidrunning(
  SYSTEM_cardinal pid)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  /**** C code included from p3process.pas(1649:1): 21 lines ****/
{
#if defined(_WIN32)
  HANDLE hProcess;
  DWORD p, rcWait;

  p = (DWORD) pid;
  hProcess = OpenProcess (PROCESS_HACK_ACCESS, FALSE, p);
  if (NULL != hProcess) {
    rcWait = WaitForSingleObject (hProcess, 0);
    /* should return either WAIT_OBJECT_0 for a zombie, or
     * WAIT_TIMEOUT for a runner */
    result = (WAIT_TIMEOUT == rcWait);
    CloseHandle (hProcess);
  }
#else
  int rc;

  rc = unixPidStatus ((int) pid);
  result = (0 == rc);
#endif
}
  return result;
}  /* p3ispidrunning */

Function(SYSTEM_boolean ) P3PROCESS_p3killprocess(
  const P3PROCESS_tprocinfo *procinfo,
  P3PROCESS_tkillhow how)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  /**** C code included from p3process.pas(1701:1): 52 lines ****/
{
#if defined(_WIN32)
  HANDLE hProcess;
  DWORD p, rc;
  const DWORD uExitCode = 0xffffffff;

  if (how == P3PROCESS_soft)
    return result;
  p = (DWORD) procinfo->pid;
  hProcess = (HANDLE) procinfo->hprocess;
  if (NULL == hProcess) {
    hProcess = OpenProcess (PROCESS_HACK_ACCESS, FALSE, p);
    if (NULL != hProcess) {
      result = TerminateProcess (hProcess, uExitCode);
      CloseHandle (hProcess);
    }
  }
  else {
    result = TerminateProcess (hProcess, uExitCode);
    CloseHandle (hProcess);
  }

#else
  int i, rc, wstat;
  pid_t p, p2;

  p = (pid_t) procinfo->pid;
  if (p > 0) {                  /* PIDs are positive */
    rc = kill (p, (P3PROCESS_soft == how) ? SIGINT: SIGKILL);
    result = (0 == rc);
    /* clean up the zombie.  If this is our child the PID is still valid
     * until we wait on the process
     * remaining issue: will this wait if this was NOT a child process? */
    if (0 == rc) {    /* signal sent successfully */
      for (i = 0;  i < 2;  i++) {
        rc = unixPidStatus (p);
        if (rc > 1)           /* not running */
          return result;
        else if (rc < 1) {      /* running, not a zombie */
          usleep (20000);     /* 20 millisecs */
          /* printf ("waiting for zombie\n"); */
        }
        else {                /* a zombie */
          p2 = waitpid (p, &wstat, 0);
          /* printf ("debug: after kill and waitpid: p2 = %d\n", p2); */
          return result;
        }
      } /* sleep loop */
    } /* if signal sent OK */
  }
#endif
}
  return result;
}  /* p3killprocess */
/**** C code included from p3process.pas(1790:1): 93 lines ****/
#if defined(_WIN32)
BOOL killProcessTree (DWORD myprocID)
{
  PROCESSENTRY32 pe;
  HANDLE hSnap;
  HANDLE hChildProc, hProc;

  memset (&pe, 0, sizeof(PROCESSENTRY32));
  pe.dwSize = sizeof(PROCESSENTRY32);

  hSnap = CreateToolhelp32Snapshot (TH32CS_SNAPPROCESS, 0);
  if (INVALID_HANDLE_VALUE == hSnap) {
    return FALSE;               /* failure */
  }
  if (! Process32First(hSnap, &pe)) {
    CloseHandle (hSnap);        /* clean up the snapshot object */
    return FALSE;               /* failure */
  }

  /* kill the main process */
  hProc = OpenProcess(PROCESS_HACK_ACCESS, FALSE, myprocID);
  if (hProc) {
    TerminateProcess (hProc, 1);
    CloseHandle (hProc);
  }

  for ( ; ; ) {
    if (pe.th32ParentProcessID == myprocID) {
      /* Recursion */
      killProcessTree (pe.th32ProcessID);

      hChildProc = OpenProcess(PROCESS_HACK_ACCESS, FALSE, pe.th32ProcessID);
      if (hChildProc) {
        TerminateProcess(hChildProc, 0);
        CloseHandle(hChildProc);
      }
    }
    if (! Process32Next(hSnap, &pe))
      break;
  }

  /* kill the main process */
  hProc = OpenProcess(PROCESS_HACK_ACCESS, FALSE, myprocID);
  if (hProc) {
    TerminateProcess (hProc, 1);
    CloseHandle (hProc);
  }
  return TRUE;                  /* success */
} /* killProcessTree */
#else
/* killProcGroupUnix
 * return:  true if the process is(was) running
 *              and our termination attempt succeeded
 *          false o/w
 */
SYSTEM_boolean killProcGroupUnix (pid_t p, P3PROCESS_tkillhow how)
{
  SYSTEM_boolean result;
  int i, rc, wstat;
  pid_t p2;

  result = SYSTEM_false;
  if (p > 0) {                  /* PIDs are positive */
# if defined(__APPLE__)
    p2 = getpgid (p);
    rc = killpg (p2, (P3PROCESS_soft == how) ? SIGINT: SIGKILL);
# else
    rc = kill   (-p, (P3PROCESS_soft == how) ? SIGINT: SIGKILL);
# endif
    result = (0 == rc);
    /* clean up the zombie.  If this is our child the PID is still valid
     * until we wait on the process
     * remaining issue: will this wait if this was NOT a child process? */
    if (0 == rc) {    /* signal sent successfully */
      for (i = 0;  i < 2;  i++) {
        rc = unixPidStatus (p);
        if (rc > 1)           /* not running */
          return result;
        else if (rc < 1) {      /* running, not a zombie */
          usleep (20000);     /* 20 millisecs */
          /* printf ("waiting for zombie\n"); */
        }
        else {                /* a zombie */
          p2 = waitpid (p, &wstat, 0);
          /* printf ("debug: after kill and waitpid: p2 = %d\n", p2); */
          return result;
        }
      } /* sleep loop */
    } /* if signal sent OK */
  }
  return result;
} /* killProcGroupUnix */
#endif  /* if defined(_WIN32) .. else .. */

Function(SYSTEM_boolean ) P3PROCESS_p3killprocgrouptp(
  const P3PROCESS_tprocinfo *procinfo,
  P3PROCESS_tkillhow how)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  /**** C code included from p3process.pas(1953:1): 18 lines ****/
{
#if defined(_WIN32)
  char cmdBuf[128];
  int rc, progRC;

  if (0 == procinfo->pid) {
    return 0;
  }

  result = killProcessTree (procinfo->pid);
  return result;

#else

  result = killProcGroupUnix ((pid_t) procinfo->pid, how);

#endif
}
  return result;
}  /* p3killprocgrouptp */

Function(SYSTEM_boolean ) P3PROCESS_p3killprocgrouptk(
  const P3PROCESS_tprocinfo *procinfo,
  P3PROCESS_tkillhow how)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  /**** C code included from p3process.pas(1989:1): 35 lines ****/
{
#if defined(_WIN32)
  HANDLE hProcess;
  char cmdBuf[128];
  int rc, progRC;

  if (0 == procinfo->pid) {
    return 0;
  }
  sprintf (cmdBuf, "taskkill /pid %u /t", (unsigned int) procinfo->pid);
  if (P3PROCESS_hard == how)
    strcat (cmdBuf, " /f");
#if 0
  strcat (cmdBuf, " > NUL");
  rc = P3PROCESS_p3systemp ((SYSTEM_P3_pansichar) cmdBuf, &progRC);
#else
  rc = P3PROCESS_p3execp ((SYSTEM_P3_pansichar) cmdBuf, &progRC);
#endif
  if (rc) {
    /* printf ("DEBUG p3KillProcGroupTK: taskkill not launched\n"); */
    return 0;
  }
  /* printf ("DEBUG p3KillProcGroupTK: taskkill launched OK: return = %d\n", progRC); */
  hProcess = (HANDLE) procinfo->hprocess;
  if (hProcess) {
    CloseHandle (hProcess);
  }
  return ! progRC;

#else

  result = killProcGroupUnix ((pid_t) procinfo->pid, how);

#endif
}
  return result;
}  /* p3killprocgrouptk */

static Function(SYSTEM_integer ) asyncsystem4unix(
  SYSTEM_P3_pansichar cmdptr,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;
  SYSTEM_cardinal pid;
  SYSTEM_integer argc, i;
  P3PROCESS_tpargv pargv;
  SYSTEM_P3_pansichar s;
  SYSTEM_shortstring param;

  if (P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin) {
    _P3strcpy(msg,255,_P3str1("\054asyncSystem4Unix not implemented for Windows"));
    result = 127;
    return result;
  } 
  _P3strclr(msg);
  s = P3PROCESS_getparamshortstr(cmdptr,param);
  if (_P3strcmpE(param,_P3str1("\000"))) {
    argc = 1;
    _P3getmem(pargv,(argc + 1) * sizeof(SYSTEM_pointer));
    (*pargv)[0] = P3PRIVATE_strtopchar(_P3str1("\007/bin/sh"));
  } else {
    argc = 3;
    _P3getmem(pargv,(argc + 1) * sizeof(SYSTEM_pointer));
    (*pargv)[0] = P3PRIVATE_strtopchar(_P3str1("\007/bin/sh"));
    (*pargv)[1] = P3PRIVATE_strtopchar(_P3str1("\002-c"));
    (*pargv)[2] = cmdptr;
  } 
  (*pargv)[argc] = NULL;
  /**** C code included from p3process.pas(2099:1): 1 lines ****/
  result = libcASyncForkExec (argc, (char *const * )( *pargv), &pid);
  procinfo->pid = pid;
  _P3freemem((*pargv)[0]);
  if (_P3strcmpN(param,_P3str1("\000"))) 
    _P3freemem((*pargv)[1]);
  _P3freemem(pargv);
  return result;
}  /* asyncsystem4unix */

static Function(SYSTEM_integer ) P3PROCESS_system4unix(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_integer rcode;
  SYSTEM_P3_pansichar newptr;

  /**** C code included from p3process.pas(2163:1): 31 lines ****/
#if defined(_WIN32)
  result = 127; /* should never happen */
#else
  if ('\0' == (char *) cmdptr) { /* special case, run default shell */
    newptr = (SYSTEM_char *) "sh";
  }
  else {
    newptr = cmdptr;
  }
  rcode = system((char *) newptr);
  if (WIFEXITED(rcode)) {       /* shell completed successfully */
    result = 0;
    *progrc = WEXITSTATUS(rcode);
    if (127 == *progrc) {       /* but cmd wasn't run (e.g file not found) */
      result = 127;
      *progrc = 0;
    }
    if (126 == *progrc) {       /* but cmd wasn't run (e.g permission denied not found) */
      result = 126;
      *progrc = 0;
    }
  }
  else if WIFSIGNALED(rcode) {  /* child stopped via a signal */
    result = 1;
    *progrc = WTERMSIG(rcode);
  }
  else {                        /* shell not completed successfully */
    result = 2;
    *progrc = 0;
  }
#endif /* #if defined(_WIN32) .. else .. */
  return result;
}  /* system4unix */

static Function(SYSTEM_integer ) asyncsystem4win(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean newconsole,
  SYSTEM_boolean inheritedhandles,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;
  SYSTEM_shortstring cs;
  SYSTEM_P3_pansichar csptr;
  P3PRIVATE_shortstrbuf csbuf;
  SYSTEM_P3_pansichar argptr;
  SYSTEM_integer arglen;

  _P3strclr(msg);
  SYSUTILS_P3_getenvironmentvariable(cs,255,_P3str1("\007COMSPEC"));
  if (_P3strcmpE(cs,_P3str1("\000"))) 
    if (SYSUTILS_P3_fileexists(P3PROCESS_cmd_win7)) { 
      _P3strcpy(cs,255,P3PROCESS_cmd_win7);
    } else 
      if (SYSUTILS_P3_fileexists(P3PROCESS_cmd_winnt)) { 
        _P3strcpy(cs,255,P3PROCESS_cmd_winnt);
      } else {
        result = 1;
        _P3strcpy(msg,255,_P3str1("\045COMSPEC not set and cmd.exe not found"));
        return result;
      } 
  csptr = P3PRIVATE_strtostrbuf(cs,csbuf);
  if (*cmdptr == _P3char('\000')) {
    arglen = SYSUTILS_P3_strlen(csptr) + 1;
    _P3getmem(argptr,arglen);
    arglen = 0;
    P3PRIVATE_pcharconcatstr(argptr,&arglen,cs);
  } else {
    arglen = SYSUTILS_P3_strlen(csptr) + 5 + SYSUTILS_P3_strlen(
      cmdptr);
    _P3getmem(argptr,arglen);
    arglen = 0;
    P3PRIVATE_pcharconcatstr(argptr,&arglen,cs);
    P3PRIVATE_pcharconcatstr(argptr,&arglen,_P3str1("\004 /C "));
    P3PRIVATE_pcharconcatpchar(argptr,&arglen,cmdptr);
  } 
  /**** C code included from p3process.pas(2249:1): 5 lines ****/
#if 0
  printf ("debug asyncSystem4Win: calling win32CreateProc (%s, %s, ...)\n",
          (char *) csptr, (char *) argptr);
#endif
  result = win32ASyncCreateProc ((char *) csptr, (char *) argptr, newconsole, inheritedhandles, procinfo);
  if (result != 0) 
    result = 2;
  if (NULL != argptr) 
    _P3freemem(argptr);
  return result;
}  /* asyncsystem4win */

static Function(SYSTEM_integer ) P3PROCESS_system4win(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean inheritedhandles,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_shortstring cs;
  SYSTEM_P3_pansichar csptr;
  P3PRIVATE_shortstrbuf csbuf;
  SYSTEM_P3_pansichar argptr;
  SYSTEM_integer arglen;

  SYSUTILS_P3_getenvironmentvariable(cs,255,_P3str1("\007COMSPEC"));
  if (_P3strcmpE(cs,_P3str1("\000"))) 
    if (SYSUTILS_P3_fileexists(P3PROCESS_cmd_win7)) { 
      _P3strcpy(cs,255,P3PROCESS_cmd_win7);
    } else 
      if (SYSUTILS_P3_fileexists(P3PROCESS_cmd_winnt)) { 
        _P3strcpy(cs,255,P3PROCESS_cmd_winnt);
      } else {
        result = 1;
        return result;
      } 
  csptr = P3PRIVATE_strtostrbuf(cs,csbuf);
  if (*cmdptr == _P3char('\000')) {
    arglen = SYSUTILS_P3_strlen(csptr) + 1;
    _P3getmem(argptr,arglen);
    arglen = 0;
    P3PRIVATE_pcharconcatstr(argptr,&arglen,cs);
  } else {
    arglen = SYSUTILS_P3_strlen(csptr) + 5 + SYSUTILS_P3_strlen(
      cmdptr);
    _P3getmem(argptr,arglen);
    arglen = 0;
    P3PRIVATE_pcharconcatstr(argptr,&arglen,cs);
    P3PRIVATE_pcharconcatstr(argptr,&arglen,_P3str1("\004 /C "));
    P3PRIVATE_pcharconcatpchar(argptr,&arglen,cmdptr);
  } 
  /**** C code included from p3process.pas(2314:1): 1 lines ****/
  result = Win32CreateProc ((char *) csptr, (char *) argptr, inheritedhandles, progrc);
  if (result != 0) 
    result = 2;
  if (NULL != argptr) 
    _P3freemem(argptr);
  return result;
}  /* system4win */

Function(SYSTEM_integer ) P3PROCESS_p3system(
  const SYSTEM_ansichar *cmd,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_P3_pansichar cmdptr;
  P3PRIVATE_shortstrbuf cmdbuf;

  cmdptr = P3PRIVATE_strtostrbuf(cmd,cmdbuf);
  result = 0;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      result = P3PROCESS_system4win(cmdptr,SYSTEM_true,progrc);
      break;
    case P3PLATFORM_osfileunix: 
      result = P3PROCESS_system4unix(cmdptr,progrc);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\045unimplemented P3system for OSFileType"));
  }
  return result;
}  /* p3system */

Function(SYSTEM_integer ) P3PROCESS_p3systemp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;

  result = 0;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      result = P3PROCESS_system4win(cmdptr,SYSTEM_true,progrc);
      break;
    case P3PLATFORM_osfileunix: 
      result = P3PROCESS_system4unix(cmdptr,progrc);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\046unimplemented P3SystemP for OSFileType"));
  }
  return result;
}  /* p3systemp */

Function(SYSTEM_integer ) P3PROCESS_p3system2(
  const SYSTEM_ansichar *progname,
  const SYSTEM_ansichar *progparams,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_P3_pansichar cmdptr;
  SYSTEM_integer cmdlen;
  SYSTEM_integer k;

  cmdlen = ValueCast(SYSTEM_int32,SYSTEM_length(progname)) + 1 + 
    SYSTEM_length(progparams);
  _P3getmem(cmdptr,cmdlen + 1);
  k = 0;
  P3PRIVATE_pcharconcatstr(cmdptr,&k,progname);
  P3PRIVATE_pcharconcatstr(cmdptr,&k,_P3str1("\001 "));
  P3PRIVATE_pcharconcatstr(cmdptr,&k,progparams);
  SYSTEM_assert(k == cmdlen,_P3str1("\040Strange result of PCharConcatStr"));
  result = 0;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      result = P3PROCESS_system4win(cmdptr,SYSTEM_true,progrc);
      break;
    case P3PLATFORM_osfileunix: 
      result = P3PROCESS_system4unix(cmdptr,progrc);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\045unimplemented P3system for OSFileType"));
  }
  _P3freemem(cmdptr);
  return result;
}  /* p3system2 */

Function(SYSTEM_integer ) P3PROCESS_p3systeml(
  const SYSTEM_ansichar *progname,
  P3PROCESS_texecarglist progparams,
  SYSTEM_integer *progrc)
{
  SYSTEM_integer result;
  SYSTEM_boolean inheritedhandles;
  SYSTEM_P3_pansichar cmdptr;
  SYSTEM_integer cmdlen;
  SYSTEM_integer k, i;
  SYSTEM_shortstring quote;

  inheritedhandles = progparams->
    P3PROCESS_texecarglist_DOT_finheritedhandles;
  cmdlen = ValueCast(SYSTEM_int32,SYSTEM_length(progname)) + 3;
  { register SYSTEM_int32 _stop = progparams->
      P3PROCESS_texecarglist_DOT_fcount - 1;
    if ((i = 0) <=  _stop) do {
      {
        SYSTEM_shortstring _t1;

        cmdlen = cmdlen + SYSTEM_length(P3PROCESS_texecarglist_DOT_get(
          _t1,255,progparams,i)) + 3;
      }
    } while (i++ !=  _stop);

  }
  _P3getmem(cmdptr,cmdlen);
  k = 0;
  P3PRIVATE_pcharconcatstr(cmdptr,&k,progname);
  { register SYSTEM_int32 _stop = progparams->
      P3PROCESS_texecarglist_DOT_fcount - 1;
    if ((i = 0) <=  _stop) do {
      {
        SYSTEM_shortstring _t2;

        P3PROCESS_whatquote(quote,255,
          P3PROCESS_texecarglist_DOT_get(_t2,255,progparams,i));
      }
      P3PRIVATE_pcharconcatstr(cmdptr,&k,_P3str1("\001 "));
      P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
      {
        SYSTEM_shortstring _t1;

        P3PRIVATE_pcharconcatstr(cmdptr,&k,
          P3PROCESS_texecarglist_DOT_get(_t1,255,progparams,i));
      }
      P3PRIVATE_pcharconcatstr(cmdptr,&k,quote);
    
    } while (i++ !=  _stop);

  }
  result = 0;
  switch (P3PLATFORM_osfiletype()) {
    case P3PLATFORM_osfilewin: 
      result = P3PROCESS_system4win(cmdptr,inheritedhandles,progrc);
      break;
    case P3PLATFORM_osfileunix: 
      result = P3PROCESS_system4unix(cmdptr,progrc);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\045unimplemented P3system for OSFileType"));
  }
  _P3freemem(cmdptr);
  return result;
}  /* p3systeml */
static SYSTEM_P3_pansichar P3PROCESS_unixcmdline;

static Function(SYSTEM_P3_pansichar ) P3PROCESS_unixgetcommandline(void)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer buflen, i, k;
  SYSTEM_shortstring s;

  if (NULL != P3PROCESS_unixcmdline) {
    result = P3PROCESS_unixcmdline;
    return result;
  } 
  buflen = 1;
  { register SYSTEM_int32 _stop = SYSTEM_P3_paramcount();
    if ((i = 0) <=  _stop) do {
      SYSTEM_P3_paramstr(s,255,i);
      buflen = buflen + SYSTEM_length(s) + 3;
    
    } while (i++ !=  _stop);

  }
  _P3getmem(result,buflen);
  k = 0;
  P3PRIVATE_pcharconcatstr(result,&k,_P3str1("\001\""));
  {
    SYSTEM_shortstring _t1;

    P3PRIVATE_pcharconcatstr(result,&k,SYSTEM_P3_paramstr(_t1,255,0));
  }
  P3PRIVATE_pcharconcatstr(result,&k,_P3str1("\002\" "));
  { register SYSTEM_int32 _stop = SYSTEM_P3_paramcount();
    if ((i = 1) <=  _stop) do {
      SYSTEM_P3_paramstr(s,255,i);
      P3PRIVATE_pcharconcatstr(result,&k,_P3str1("\001 "));
      {
        SYSTEM_shortstring _t1;

        P3PRIVATE_pcharconcatstr(result,&k,SYSTEM_P3_paramstr(_t1,255,
          i));
      }
    
    } while (i++ !=  _stop);

  }
  P3PROCESS_unixcmdline = result;
  return result;
}  /* unixgetcommandline */

Function(SYSTEM_P3_pansichar ) P3PROCESS_p3getcommandline(void)
{
  SYSTEM_P3_pansichar result;

  if (P3PLATFORM_osfilewin != P3PLATFORM_osfiletype()) { 
    result = P3PROCESS_unixgetcommandline();
  } else 
    /**** C code included from p3process.pas(2501:1): 5 lines ****/
#if defined(_WIN32)
    result = (SYSTEM_P3_pansichar) GetCommandLine();
#else
    result = NULL;
#endif
  return result;
}  /* p3getcommandline */
static P3PROCESS_tctrlhandler P3PROCESS_ctrlhandler;
/**** C code included from p3process.pas(2524:1): 6 lines ****/
#if ! defined(_WIN32)
static sigset_t sigSet;
static struct sigaction newAction;
static struct sigaction oldAction;
C_LINKAGE(static void ( * oldHandler) (int);)
#endif
/**** C code included from p3process.pas(2556:1): 28 lines ****/
#define ctrlhandler P3PROCESS_ctrlhandler

#if defined(_WIN32)
static BOOL WINAPI
P3Handler(DWORD s)
{
  BOOL result;
  /* */
  result = FALSE;
  if (CTRL_C_EVENT == s) {
    if (NULL != ctrlhandler) {
      ctrlhandler();
      result = TRUE;
    }
  }
  return result;
} /* P3Handler */
#else

C_LINKAGE(static void p3CtrlCHandler (int sig);)
static void
p3CtrlCHandler (int sig)
{
  if (NULL != ctrlhandler)
    ctrlhandler();
  return;
}
#endif

Function(SYSTEM_integer ) P3PROCESS_p3installctrlhandler(
  P3PROCESS_tctrlhandler newhandler)
{
  SYSTEM_integer result;
  SYSTEM_integer rc;

  if (NULL == ValueCast(SYSTEM_pointer,newhandler)) {
    result = P3PROCESS_p3uninstallctrlhandler();
    return result;
  } 
  if (ValueCast(SYSTEM_pointer,P3PROCESS_ctrlhandler) != NULL) {
    result = P3PROCESS_p3ctrlhandlerok;
    P3PROCESS_ctrlhandler = ValueCast(P3PROCESS_tctrlhandler,
      newhandler);
  } else 
    /**** C code included from p3process.pas(2647:1): 30 lines ****/
{
#if defined(_WIN32)
  ctrlhandler = newhandler;
  if (SetConsoleCtrlHandler(P3Handler,TRUE)) {
    result = P3PROCESS_p3ctrlhandlerok;
  }
  else {
    result = P3PROCESS_p3ctrlhandlersysfail;
    ctrlhandler = NULL;
  }
#else
  rc = sigemptyset(&sigSet);
  if (0 != rc) {
    return P3PROCESS_p3ctrlhandlersysfail;
  }
  newAction.sa_handler = p3CtrlCHandler;
  newAction.sa_mask    = sigSet;
  newAction.sa_flags   = 0;
  ctrlhandler = newhandler;
  rc = sigaction(SIGINT, &newAction, &oldAction);
  if (0 != rc) {
    result = P3PROCESS_p3ctrlhandlersysfail;
    ctrlhandler = NULL;
  }
  else {
    oldHandler = oldAction.sa_handler;
    result = P3PROCESS_p3ctrlhandlerok;
  }
#endif /* #if defined(_WIN32) .. else .. */
}
  return result;
}  /* p3installctrlhandler */

Function(SYSTEM_integer ) P3PROCESS_p3uninstallctrlhandler(void)
{
  SYSTEM_integer result;
  SYSTEM_integer rc;

  if (ValueCast(SYSTEM_pointer,P3PROCESS_ctrlhandler) == NULL) { 
    result = P3PROCESS_p3ctrlhandlerwasempty;
  } else {
    /**** C code included from p3process.pas(2719:1): 21 lines ****/
#if defined(_WIN32)
  if (SetConsoleCtrlHandler(P3Handler,FALSE))
    result = P3PROCESS_p3ctrlhandlerok;
  else
    result = P3PROCESS_p3ctrlhandlersysfail;
#else
  rc = sigemptyset(&sigSet);
  if (0 != rc) {
    return P3PROCESS_p3ctrlhandlersysfail;
  }

  /* newAction.sa_handler = SIG_DFL; */
  newAction.sa_handler = oldHandler;
  newAction.sa_mask    = sigSet;
  newAction.sa_flags   = 0;
  rc = sigaction(SIGINT, &newAction, &oldAction);
  if (rc != 0 || oldAction.sa_handler != p3CtrlCHandler)
    result = P3PROCESS_p3ctrlhandlersysfail;
  else
    result = P3PROCESS_p3ctrlhandlerok;
#endif /* #if defined(_WIN32) .. else .. */
    P3PROCESS_ctrlhandler = NULL;
  } 
  return result;
}  /* p3uninstallctrlhandler */

Function(P3PROCESS_tctrlhandler ) P3PROCESS_p3getctrlhandler(void)
{
  P3PROCESS_tctrlhandler result;

  result = P3PROCESS_ctrlhandler;
  return result;
}  /* p3getctrlhandler */

static Function(SYSTEM_ansichar *) P3PROCESS_getstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pshortstring ps)
{
  if (ps == NULL) { 
    _P3strclr(result);
  } else 
    _P3strcpy(result,_len_ret,*ps);
  return result;
}  /* getstring */

static Function(SYSTEM_P3_pshortstring ) P3PROCESS_newstring(
  const SYSTEM_ansichar *s)
{
  SYSTEM_P3_pshortstring result;

  if (_P3strcmpE(s,_P3str1("\000"))) { 
    result = NULL;
  } else {
    _P3getmem(result,ValueCast(SYSTEM_int32,SYSTEM_length(s)) + 1);
    _P3strcpy(*result,255,s);
  } 
  return result;
}  /* newstring */

static Procedure P3PROCESS_disposestring(
  SYSTEM_P3_pshortstring ps)
{
  if (ps != NULL) 
    _P3freemem2(ps,ValueCast(SYSTEM_int32,SYSTEM_length(*ps)) + 1);
}  /* disposestring */

Constructor(P3PROCESS_texecarglist ) P3PROCESS_texecarglist_DOT_create(
  P3PROCESS_texecarglist self)
{
  ValueCast(P3PROCESS_texecarglist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->P3PROCESS_texecarglist_DOT_finheritedhandles = SYSTEM_true;
  return self;
}  /* create */

Destructor(P3PROCESS_texecarglist ) P3PROCESS_texecarglist_DOT_destroy(
  P3PROCESS_texecarglist self)
{
  P3PROCESS_texecarglist_DOT_clear(self);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_ansichar *) P3PROCESS_texecarglist_DOT_getlast(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  P3PROCESS_texecarglist self)
{
  if (self->P3PROCESS_texecarglist_DOT_fcount <= 0) { 
    _P3strclr(result);
  } else 
    P3PROCESS_getstring(result,_len_ret,ValueCast(
      SYSTEM_P3_pshortstring,(*self->P3PROCESS_texecarglist_DOT_flist)[
      self->P3PROCESS_texecarglist_DOT_fcount - 1]));
  return result;
}  /* getlast */

Procedure P3PROCESS_texecarglist_DOT_clear(
  P3PROCESS_texecarglist self)
{
  SYSTEM_integer n;

  for (n = self->P3PROCESS_texecarglist_DOT_fcount - 1;n >= (
    SYSTEM_int32)0;--n) {
    P3PROCESS_texecarglist_DOT_freeitem(self,n);
  }
  self->P3PROCESS_texecarglist_DOT_fcount = 0;
  P3PROCESS_texecarglist_DOT_setcapacity(self,0);
}  /* clear */

Procedure P3PROCESS_texecarglist_DOT_delete(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index)
{
  P3PROCESS_texecarglist_DOT_freeitem(self,index);
  _P3dec0(self->P3PROCESS_texecarglist_DOT_fcount);
  if (index < self->P3PROCESS_texecarglist_DOT_fcount) 
    SYSTEM_move(&(*self->P3PROCESS_texecarglist_DOT_flist)[index + 1],&(*
      self->P3PROCESS_texecarglist_DOT_flist)[index],(self->
      P3PROCESS_texecarglist_DOT_fcount - index) * sizeof(
      SYSTEM_pointer));
}  /* delete */

Procedure P3PROCESS_texecarglist_DOT_grow(
  P3PROCESS_texecarglist self)
{
  SYSTEM_integer delta;

  if (self->P3PROCESS_texecarglist_DOT_fcapacity >= 1048576) { 
    delta = self->P3PROCESS_texecarglist_DOT_fcapacity /  4;
  } else 
    if (self->P3PROCESS_texecarglist_DOT_fcapacity == 0) { 
      delta = 16;
    } else 
      delta = 7 * self->P3PROCESS_texecarglist_DOT_fcapacity;
  P3PROCESS_texecarglist_DOT_setcapacity(self,self->
    P3PROCESS_texecarglist_DOT_fcapacity + delta);
}  /* grow */

Procedure P3PROCESS_texecarglist_DOT_setcapacity(
  P3PROCESS_texecarglist self,
  SYSTEM_integer newcapacity)
{
  if (newcapacity != self->P3PROCESS_texecarglist_DOT_fcapacity) {
    if (newcapacity < self->P3PROCESS_texecarglist_DOT_fcount) 
      newcapacity = self->P3PROCESS_texecarglist_DOT_fcount;
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      P3PROCESS_texecarglist_DOT_flist),newcapacity * sizeof(
      SYSTEM_pointer));
    self->P3PROCESS_texecarglist_DOT_fcapacity = newcapacity;
  } 
}  /* setcapacity */

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_splitappend(
  P3PROCESS_texecarglist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = P3PROCESS_texecarglist_DOT_split(self,1,s);
  return result;
}  /* splitappend */

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_splitprepend(
  P3PROCESS_texecarglist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = P3PROCESS_texecarglist_DOT_split(self,0,s);
  return result;
}  /* splitprepend */

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_split(
  P3PROCESS_texecarglist self,
  SYSTEM_integer append,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer count;
  SYSTEM_P3_pansichar sptr, p;
  P3PRIVATE_shortstrbuf sbuf;
  SYSTEM_shortstring param;

  sptr = P3PRIVATE_strtostrbuf(s,sbuf);
  p = sptr;
  count = 0;
  while (SYSTEM_true) {
    p = P3PROCESS_getparamshortstr(p,param);
    if (_P3strcmpE(param,_P3str1("\000"))) 
      SYSTEM_break(BRK_2);
    if (append != 0) { 
      P3PROCESS_texecarglist_DOT_add(self,param);
    } else 
      P3PROCESS_texecarglist_DOT_insert(self,count,param);
    _P3inc0(count);
  
CNT_2:;
  }
BRK_2:;
  result = count;
  return result;
}  /* split */

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_add(
  P3PROCESS_texecarglist self,
  const SYSTEM_ansichar *item)
{
  SYSTEM_integer result;

  result = self->P3PROCESS_texecarglist_DOT_fcount;
  if (result == self->P3PROCESS_texecarglist_DOT_fcapacity) 
    P3PROCESS_texecarglist_DOT_grow(self);
  (*self->P3PROCESS_texecarglist_DOT_flist)[result] = 
    P3PROCESS_newstring(item);
  _P3inc0(self->P3PROCESS_texecarglist_DOT_fcount);
  return result;
}  /* add */

Procedure P3PROCESS_texecarglist_DOT_freeitem(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index)
{
  P3PROCESS_disposestring(ValueCast(SYSTEM_P3_pshortstring,(*self->
    P3PROCESS_texecarglist_DOT_flist)[index]));
}  /* freeitem */

Function(SYSTEM_ansichar *) P3PROCESS_texecarglist_DOT_get(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  P3PROCESS_texecarglist self,
  SYSTEM_integer index)
{
  P3PROCESS_getstring(result,_len_ret,ValueCast(SYSTEM_P3_pshortstring,(*
    self->P3PROCESS_texecarglist_DOT_flist)[index]));
  return result;
}  /* get */

Procedure P3PROCESS_texecarglist_DOT_insert(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *item)
{
  if (self->P3PROCESS_texecarglist_DOT_fcount == self->
    P3PROCESS_texecarglist_DOT_fcapacity) 
    P3PROCESS_texecarglist_DOT_grow(self);
  if (index < self->P3PROCESS_texecarglist_DOT_fcount) 
    SYSTEM_move(&(*self->P3PROCESS_texecarglist_DOT_flist)[index],&(*
      self->P3PROCESS_texecarglist_DOT_flist)[index + 1],(self->
      P3PROCESS_texecarglist_DOT_fcount - index) * sizeof(
      SYSTEM_pointer));
  (*self->P3PROCESS_texecarglist_DOT_flist)[index] = 
    P3PROCESS_newstring(item);
  _P3inc0(self->P3PROCESS_texecarglist_DOT_fcount);
}  /* insert */

Procedure P3PROCESS_texecarglist_DOT_put(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring item;

  _P3strcpy(item,255,_ftmp1);
  P3PROCESS_texecarglist_DOT_freeitem(self,index);
  (*self->P3PROCESS_texecarglist_DOT_flist)[index] = 
    P3PROCESS_newstring(item);
}  /* put */

Procedure P3PROCESS_p3defaultshowwindow(void)
{
  /**** C code included from p3process.pas(2939:1): 5 lines ****/
#if defined(_WIN32)
  wshowwindow = SW_SHOWNA;
#else
  wshowwindow = 0;
#endif /* #if defined(_WIN32) .. else .. */
}  /* p3defaultshowwindow */

Procedure P3PROCESS_p3setshowwindow(
  SYSTEM_integer showwindow)
{
  P3PROCESS_wshowwindow = showwindow;
}  /* p3setshowwindow */

Function(SYSTEM_integer ) P3PROCESS_p3getshowwindow(void)
{
  SYSTEM_integer result;

  result = P3PROCESS_wshowwindow;
  return result;
}  /* p3getshowwindow */

Function(SYSTEM_integer ) P3PROCESS_p3getnumberofprocessors(void)
{
  SYSTEM_integer result;

  /**** C code included from p3process.pas(2969:1): 7 lines ****/
#if defined(_WIN32)
  SYSTEM_INFO siSysInfo;
  GetSystemInfo(&siSysInfo);
  result = (int) siSysInfo.dwNumberOfProcessors;
#else
  result = (int) sysconf(_SC_NPROCESSORS_ONLN);
#endif /* #if defined(_WIN32) .. else .. */
  return result;
}  /* p3getnumberofprocessors */

/* unit p3process */
void _Init_Module_p3process(void)
{
  P3PROCESS_ctrlhandler = NULL;
  P3PROCESS_unixcmdline = NULL;
  P3PROCESS_p3defaultshowwindow();
} /* _Init_Module_p3process */

void _Final_Module_p3process(void)
{
} /* _Final_Module_p3process */

