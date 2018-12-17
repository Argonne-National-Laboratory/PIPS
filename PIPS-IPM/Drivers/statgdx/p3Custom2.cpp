/* $Id: p3Custom2.cpp 52165 2015-05-19 13:33:37Z sdirkse $
 * This code is sensitive to the order of the #defines, so inlining is not OK
 */

#if defined(__linux)
# define _GNU_SOURCE  /* required for dladdr() but no desired globally */
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>

#if defined(_WIN32)
# define WIN32_LEAN_AND_MEAN   /* google it */
# include <windows.h>
# define snprintf _snprintf
#else
# include <unistd.h>
#endif

#if defined(__linux) || defined(__sun)
# include <dlfcn.h>
#endif
#if defined(__APPLE__)
# include <dlfcn.h>
# include <libproc.h>
#endif
#if defined(__HOS_AIX__)
# include <procinfo.h>
# include <sys/procfs.h>
# include <dlfcn.h>
# include <sys/types.h>
# include <sys/ldr.h>
#endif

#include "p3Custom2.h"

/* local use only: be sure to call with enough space for the sprintf */
static void myStrError (int n, char *buf, size_t bufSiz)
{
#if defined(_WIN32)
  if (strerror_s (buf, bufSiz, n))
    (void) sprintf (buf, "errno = %d", n);
#else
  if (strerror_r (n, buf, bufSiz))
    (void) sprintf (buf, "errno = %d", n);
#endif
} /* myStrError */

static void pchar2SS (unsigned char *ss, const char *p)
{
  char *d = (char *)ss+1;
  const char *s = p;
  const char *end = s + 255;

  while (*s && (s < end))
    *d++ = *s++;
  *ss = (unsigned char) (s - p);
} /* pchar2SS */


int xGetExecName (unsigned char *execName, unsigned char *msg)
{
  char execBuf[4096];
  char msgBuf[2048];
  char tmpBuf[2048];
  int rc;

  *msgBuf = '\0';
  rc = 8;

#if defined(__APPLE__)
  {
    int k;
    pid_t pid = getpid();
    k = proc_pidpath (pid, execBuf, sizeof(execBuf));
    if (k <= 0) {
      myStrError (errno, tmpBuf, sizeof(tmpBuf));
      (void) snprintf (msgBuf, sizeof(msgBuf),
                       "proc_pidpath(pid=%ld) failed: %s",
                       (long int) pid, tmpBuf);
      msgBuf[sizeof(msgBuf)-1] = '\0';
      *execBuf = '\0';
      rc = 4;
    }
    else
      rc = 0;
  }

#elif defined(__HOS_AIX__)
  {
    char procPath[128];
    struct stat statData;
    struct psinfo psinfoData;
    ssize_t ssz;
    pid_t pid;
    int k, fd;

    pid = getpid();
    sprintf (procPath, "/proc/%ld/psinfo", (long int) pid);
    k = stat (procPath, &statData);
    if (k) { /* error */
      myStrError (errno, tmpBuf, sizeof(tmpBuf));
      (void) snprintf (msgBuf, sizeof(msgBuf), "stat(%s,...) failure: %s",
                       procPath, tmpBuf);
      msgBuf[sizeof(msgBuf)-1] = '\0';
      *execBuf = '\0';
      rc = 4;
    }
    else {
      fd = open (procPath, O_RDONLY);
      if (fd <= 0) { /* error */
        myStrError (errno, tmpBuf, sizeof(tmpBuf));
        (void) snprintf (msgBuf, sizeof(msgBuf),
                         "open(%s, O_RDONLY) failure: %s", procPath, tmpBuf);
        msgBuf[sizeof(msgBuf)-1] = '\0';
        *execBuf = '\0';
        rc = 5;
      }
      else {
        ssz = read (fd, &psinfoData, sizeof(psinfoData));
        close (fd);
        if (ssz < 0) { /* error */
          myStrError (errno, tmpBuf, sizeof(tmpBuf));
          (void) snprintf (msgBuf, sizeof(msgBuf),
                           "reading %ld bytes of psinfoData from %s failed: %s",
                           sizeof(psinfoData), procPath, tmpBuf);
          msgBuf[sizeof(msgBuf)-1] = '\0';
          *execBuf = '\0';
          rc = 6;
        }
        else {
          const char *argv0;
          argv0 = ((const char ***) psinfoData.pr_argv)[0][0];
          if (realpath(argv0,execBuf)) /* success */
            rc = 0;
          else {
            myStrError (errno, tmpBuf, sizeof(tmpBuf));
            (void) snprintf (msgBuf, sizeof(msgBuf),
                             "realpath() failure: %s", tmpBuf);
            msgBuf[sizeof(msgBuf)-1] = '\0';
            *execBuf = '\0';
            rc = 7;
          }
        } /* read worked OK */
      }   /* open worked OK */
    }     /* stat worked OK */
  }

#elif defined(__linux)
  {
    ssize_t ssz;

    ssz = readlink ("/proc/self/exe", execBuf, sizeof(execBuf));
    if (ssz < 0) {
      myStrError (errno, tmpBuf, sizeof(tmpBuf));
      (void) snprintf (msgBuf, sizeof(msgBuf),
                       "readlink(/proc/self/exe,...) failure: %s", tmpBuf);
      msgBuf[sizeof(msgBuf)-1] = '\0';
      *execBuf = '\0';
      rc = 4;
    }
    else {
      if (ssz >= (ssize_t) sizeof(execBuf))
        ssz = (ssize_t) sizeof(execBuf) - 1;
      execBuf[ssz] = '\0';
      rc = 0;
    }
  }

#elif defined(__sun)
  {
    const char *execname;

    execname = getexecname();
    if (NULL == execname) {
      sprintf (msgBuf, "getexecname() failure");
      *execBuf = '\0';
      rc = 4;
    }
    else {
      if (realpath(execname,execBuf))
        rc = 0;
      else {
        myStrError (errno, tmpBuf, sizeof(tmpBuf));
        (void) snprintf (msgBuf, sizeof(msgBuf), "realpath() failure: %s",
                         tmpBuf);
        msgBuf[sizeof(msgBuf)-1] = '\0';
        *execBuf = '\0';
        rc = 5;
      }
    }
  }

#elif defined(_WIN32)
  {
    HMODULE h;
    int k = GetModuleFileName (NULL, execBuf, sizeof(execBuf));
    if (0 == k) {
      sprintf (msgBuf, "GetModuleFileName() failure: rc=%d", k);
      *execBuf = '\0';
      rc = 4;
    }
    else
      rc = 0;
  }

#else
  *execBuf = '\0';
  (void) strcpy (msgBuf, "not implemented for this platform");
  rc = 8;
#endif

  pchar2SS (execName, execBuf);
  pchar2SS (msg, msgBuf);
  if ((0 == rc) && (strlen(execBuf) > 255))
    rc = 1;
  return rc;
} /* xGetExecName */

int xGetLibName (unsigned char *libName, unsigned char *msg)
{
  char libBuf[4096];
  char msgBuf[2048];
  char tmpBuf[2048];
  int rc, k;

  *msgBuf = '\0';
  rc = 8;

#if defined(__linux) || defined(__sun) || defined(__APPLE__)
  {
    Dl_info dlInfo;

    k = dladdr((void *)(&xGetLibName), &dlInfo);
    if (k > 0) {
      strncpy (tmpBuf, dlInfo.dli_fname, sizeof(tmpBuf));
      tmpBuf[sizeof(tmpBuf)-1] = '\0';
      if (realpath(tmpBuf,libBuf))
        rc = 0;
      else {
        myStrError (errno, tmpBuf, sizeof(tmpBuf));
        sprintf (msgBuf, "realpath() failure: %s", tmpBuf);
        *libBuf = '\0';
        rc = 5;
      }
    }
    else {
      sprintf (msgBuf, "dladdr() failure");
      *libBuf = '\0';
      rc = 4;
    }
  }

#elif defined(__HOS_AIX__)
  {
    /* some truly ugly code to get the shared library name:
     * so ugly it's beautiful */
    const size_t bufSize = 8192;
    char buf[bufSize], *pBuf;
    struct ld_info *pldi = (struct ld_info *) buf;
    void *pMe, *pLo, *pUp;

    k = loadquery (L_GETINFO, pldi, bufSize);
    if (k) {
      myStrError (errno, tmpBuf, sizeof(tmpBuf));
      (void) snprintf (msgBuf, sizeof(msgBuf), "loadquery() failure: %s",
                       tmpBuf);
      msgBuf[sizeof(msgBuf)-1] = '\0';
      *libBuf = '\0';
      rc = 4;
    }
    else { /* loadquery OK */
      int done = 0;
      /* pMe = &xGetLibName; */
      /* on AIX &xGetLibName is a TOC (Table-Of-Contents) pointer, must dereference
       * For more on this, I found this URL helpful:
       * http://physinfo.ulb.ac.be/divers_html/powerpc_programming_info/intro_to_ppc/ppc4_runtime4.html
      */
      pMe = (void *) * (unsigned long *) (void *) &xGetLibName;
      pBuf = buf;
      while (! done) {
        pLo = pldi->ldinfo_textorg;
        pUp = (void *) ((char *)pLo + pldi->ldinfo_textsize);
        /* printf ("DEBUG AIX dllPath: rng= [%p %p] : %s  len=%ld\n",
                    pLo, pUp, pldi->ldinfo_filename, pldi->ldinfo_next); */
        if ((pLo <= pMe) && (pMe < pUp)) {
          if (realpath(pldi->ldinfo_filename,libBuf))
            rc = 0;
          else {
            myStrError (errno, tmpBuf, sizeof(tmpBuf));
            (void) snprintf (msgBuf, sizeof(msgBuf),
                             "realpath() failure: %s", tmpBuf);
            msgBuf[sizeof(msgBuf)-1] = '\0';
            *libBuf = '\0';
            rc = 5;
          }
          done = 1;
        }
        if (0 == pldi->ldinfo_next)
          done = 1;
        else {
          pBuf += pldi->ldinfo_next;
          pldi = (struct ld_info *) pBuf;
        }
      }    /* while not done */
    }      /* loadquery OK */
  }

#elif defined(_WIN32)
  {
    HMODULE h;
    k = GetModuleHandleEx (GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                           GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                           (LPCTSTR)&xGetLibName, &h);
    if (k) {  /* OK: got a handle */
      k = GetModuleFileName (h, libBuf, sizeof(libBuf));
      if (0 == k) {
        sprintf (msgBuf, "GetModuleFileName() failure: rc=%d", k);
        *libBuf = '\0';
        rc = 5;
      }
      else {
        rc = 0;
      }
    }
    else {
      sprintf (msgBuf, "GetModuleHandleEx() failure: rc=%d", k);
      *libBuf = '\0';
      rc = 4;
    }
  }

#else
  *libBuf = '\0';
  (void) strcpy (msgBuf, "not implemented for this platform");
  rc = 8;
#endif

  pchar2SS (libName, libBuf);
  pchar2SS (msg, msgBuf);
  if ((0 == rc) && (strlen(libBuf) > 255))
    rc = 1;
  return rc;
} /* xGetLibName */
