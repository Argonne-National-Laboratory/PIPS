/*******************************************************/
/*         P3    RUNTIME SUPPORT ROUTINES              */
/*               (c) 1991++ Soren Nielsen              */
/*                                                     */
/*******************************************************/

/************************** history *********************

25 May 06 SYSTEM_delete: Change args index,count from byte to integer
          SYSTEM_pos: Change function type from byte to integer
19 May 06 SPD: remove old spawn stuff - we don't use this anymore
          whacked GC_spawnvp, GC_spawnv, GC_shell, P3GC_gsspawnvp,
          p3gsexec, p3_spawn, p3gssystem, P3GC_gsspawn, P3GC_gsshell
19 May 06 SPD: do some other cleanup too
          get_environment_variable, set_environment_variable, hackEnvSPD
          Could still look at paramstr/moduleFileName stuff
1 Jan 06  SPD: modify _P3_new/_P3_free to handle size<=0 and nil
04/25/03  Added Steve's new routines for paramcount/paramstr for _WIN32
          and merged them with the old Unix stuff (which is fine as it is).
03/03/18  Get rid of ifdefs on EXCEPTION_MODEL, always use setjmp/longjmp.
03/03/18  Remove runtime_error; replace by _P3_Exception calls.
03/03/14  Add definitions for SYSTEM_exception and all its stuff.
          Keep all the old stuff for now (p3_signum etc.)
02/09/14  _P3strcmp: Initial test if args are identical returned 1! chg to 0.
02/09/06  Re-insert _P3setlength as procedure, safer than inlining.
02/08/22  Steve updated getting the error message for an errno to more general.
02/07/03  Check realloc's return value for NULL/runtime_error.
02/06/22  Discovered SIGIOT=SIGABRT on HP so don't use abort(), replace by exit(3);
02/06/22  _P3_Exception, cleanup, signal handlers, etc. etc.
02/06/20  Added Iplus and File, Line parms to: eof, eoln, seekeof/eoln, Filepos, Filesize.
02/05/31  In _P3block_read_write, further mod from Steve.
02/05/30  In _P3block_read_write, don't set _P3_errno if 4th parm was specified.
02/04/16  SYSTEM_delete: Change byte parameters to integers.
02/04/04  Chg BIND_DEFERRED to BIND_IMMEDIATE for HP dll loads
02/01/20  Chg _P3set_i so 2nd parm is _SYSTEM_longint (avoid trunc to byte)
02/01/16  Change _P3_Close to not close stdin or stdout. Add P3CLOSED.
02/01/08  SYSTEM_paramstr(0) now returns fully qualified name on all platf's
01/12/14  Change _P3fileopn again: Re. filemode do precisely as Delphi.
          Also cosmetic changes, adding _P3RESET, _P3REWRITE etc.
01/12/13  Change _P3fileopn behavior regarding filemode.
01/12/10  Follow Steve's suggestion and call GetModuleFileName in Windows
01/12/02: Remove tests on __alpha,use AXU instead (defining _init vs __init)
01/10/10: Change _P3_Assign back so it doesn check for empty name.
          Empty name indicates standard input/output!!!
01/09/17: Use __FILE__, __LINE__ info in P3assert.
01/09/17: New set/get_environment_variable routines from Steve.
01/09/16: Add __FILE__ and __LINE__ info to runtime I/O check error message
          (should generate .pas info in .c files...)
01/09/16: Remove the 't' or 'b' flags from fopen in Unix, leave in DOS
01/09/16: Chgd filemode default from 3 to 2: ReadWrite.
01/08/22: Add Steve's unixGetModuleFileName routine(s).
01/08/22: Change in _P3_Assign: Return error if length of file name is 0.
01/08/14: New integer types: (u)int8/16/32.
01/08/14: dispose/freemem: Allow freeing non-Lvalue, don't pass **.
01/07/13: New version of P3_Str_d2 from Erwin
01/07/13: Added _P3VariableCastError routine (unequal sizes in VariableCast)
01/07/12: Added void on all empty formal parm liste
01/07/11: Added checks on all scanf's return codes (0 or EOF).
01/07/11: Using Erwin's float->string routines (_P3_Str_d[012]),
          also use them for float printing (write_rx, write_ry).
01/02/01: Added support for DJGPP: gpp for DOS. This gpp works except:
          (1) no fork(), use system instead (*usually* get rc),
          (2) no -shared, i.e., can't do shared objects/dll's.
01/01/30: Replace P3CG_gsspawnv by P3CG_gsspawn to spawnv or spawnvp.
01/01/04: Cast argv to char** (from SYSTEM_char**); fix Alpha reallocmem.
00/12/17: Added SYSTEM_reallocmem; added include <malloc.h>
00/10/22: Added p3_spawnv which unpacks (parses) parameter list, then
          calls P3GC_gsspawnv.
00/10/18: Added p3gsexec, p3gssystem (copied/renamed from gamsext.c).
          Also needed gsshell, gsspawnv, gsspawnvp renamed p3...,
          from gcprocs.c. Changed #if defined(WAT)|| etc. etc. to
          #if define(P3DOS).
00/09/27: Various changes in integer output formats (%d to %ld or %lu):
          Add FMT_(routine) defines dependent upon integer size -> p3io.h
00/09/27: Take out section on Machine flags - never used anyway.
00/09/21: Added TObject.free procedure, plus some cosmetics.
00/09/14: SSN: Add _P3_abstract_call error handling fcts.
00/08/22: New code for Val_i from Erwin.
00/05/18: Code for writeln(x:a:b) changed back to %f from %g.
00/21/02: Various changes in write_rx, write_ry (closer to Delphi)
00/21/02: Make _P3_strcpy into a function so can use for Tobject.ClassName
00/19/01: Changed tms_cutime back to tms_utime.
               Removed inling of sqr and abs (compiler problems).
99/12/15: Started adding class stuff (TObject representation)
99/12/14: Cleaned up some of Erwin's Str/Val code (_P3_Cstr2Pstr)
99/12/13: Added Steve's bug correction to _P3read_c (int return code),
               remove a few i>=0 tests for i unsigned (gcc warnings).
99/11/16: Correct _P3_insert (EK bug).
99/11/15: Merge EK's str and val code
99/09/27: Change unix get_cpu_time so returns user+child times (cu_time)
99/09/17: Implement _P3SET_ic for constant i in "i in set"
99/08/29: On unix, remove trailing CarriageReturn (13) in read string
               and blank out CtrlZ's (added 12/99).
99/08/29: Correct bug in seekeoln/seekeof
99/08/13: Change all occurrences of memcpy to memmove.
99/08/10: Add Filepos function.
99/08/06: Change _P3copy (Erwin/Steve's change of test)
99/08/06: Make string constants False, True upper case as Delphi
99/07/30: Make P3str2pa macro in p3io.h; remove from p3io.c
99/06/27: Add block_size parameter to untyped reset/rewrite
               Add seek funktion.
99/06/26: Error checks after all I/O; set _P3_errno
99/06/21: Change prototype for _P3set_i + add check agnst len

99/06/12: Enter Erwins changes resulting compiling as C++
99/03/22: Merged Gams' and my p3io.c/h.

02/09/1999 PV: enabled/disabled IOchecking
               (the run-time I/O read/write routines do
                not set the error! decided to leave it
                because of the performance penalty
02/09/1999 SD: closed file when opening a directory
               changed P3errno and errno logic
02/03/1999 PV: removed system time from get_cpu_time
05/06/1998 SD: NOSTRERR stuff
04/20/1998 SD: fileopn fails for directory on UNIX
04/20/1998 SD: changed error output from stderr to stdout
11/20/98: SSN: Init stdin/stdout in _P3input/output in PGM_init
11/14/98: SSN: Minor changes in _P3fileopn and _P3_Close
11/14/98: SSN: Fixed block_read_write (wrong parm order and r/w flag
11/13/98: SSN: Changed test message in fileopn
10/24/98: SSN: Re-introduce P3str2pa w. any length combination
07/13/98: SSN: Change so ParamStr(0) returns prog. name
07/12/98: SSN: Goes with .h of same date.
11/10/97: SSN: Changed all prototypes to ANSI.
11/09/97: SSN: Added PC clock routine (Microsoft C++)
                 plus includes: stdlib.h, direct.h.
                          Removed mkdir's second parameter under DOS
03/13/97: SD: use HZ etc from include files; added system time
03/13/97: PV: removed size parameter in FreeMem
12/05/96: EK: added blockread/write/untyped files/filemode
 7/24/96: SPD: changed definition of macro _P3read_fu0 to _P3read_fi0
 7/16/96: SPD: added gsspecialwat
 6/26/96: SPD: added gsfreset
 5/24/96: SPD: added gsdostype
 5/22/96: SPD: added gssystem, gsexec, and hooks to gclib
 2/02/92: writing unsigned and pointers: changed parameters to u
11/13/91: Added _P3_seekeof, _P3_seekeoln.
10/27/91: Corrected _P3rangechk; didn't return anyhing!
             Same with get_env and get_cwd.

****************************** end history ***************************/


#include "p3io.h"
#include "exceptions_c.h"

#include "globals2.h"

#include <limits.h>
#include <sys/stat.h>

/* Would like to find better error codes for the next ones than just EIO */
#define READ_CONVERSION_ERROR  EIO  /* I/O error - good code? */
#define READ_AT_END_OF_FILE    EIO  /* EOF - can't find a better code. There's no EEOF */

_P3err_t _P3_err = { 0, 0, 0, ""};

static void
p3ErrZap (_P3err_t *e)
{
  (void) memset (e, 0, sizeof(_P3err_t));
} /* p3ErrZap */

/* set error info for operations with _P3file arg */
static void
p3ErrSet (_P3err_t *e, int n, const _P3file *p3fp,
          unsigned char verb, unsigned char notOpened)
{
  e->n = n;
  e->verb = verb;
  e->notOpened = notOpened;
  if (p3fp) {
    (void) memcpy(e->nam, p3fp->nam, ((size_t) p3fp->nam[0]) + 2);
  }
} /* p3ErrSet */

/* set error info for file operations with string arg */
static void
p3ErrCSet (_P3err_t *e, int n, const SYSTEM_char *s, unsigned char verb)
{
  e->n = n;
  e->verb = verb;
  e->notOpened = 0;
  if (s) {
    (void) memcpy(e->nam, s, ((size_t) s[0]) + 1); /* copy shortString */
    e->nam[s[0]+1] = '\0';                         /* and add null byte */
  }
} /* p3ErrCSet */


#if 0
# define _P3_IO_NOTOPEN(REC, NN, P3FP, VERB) REC.n = NN,  REC.p3fp = P3FP,  REC.verb = VERB, REC.notOpened = 1
# define _P3_IO_ERR(REC, NN, P3FP, VERB)     REC.n = NN,  REC.p3fp = P3FP,  REC.verb = VERB
#else
# define _P3_IO_NOTOPEN(REC, NN, P3FP, VERB) p3ErrSet (&REC, NN, P3FP, VERB, 1)
# define _P3_IO_ERR(REC, NN, P3FP, VERB)     p3ErrSet (&REC, NN, P3FP, VERB, 0)
#endif

static const char blankBuf[] = "                                                   ";
SYSTEM_byte    SYSTEM_filemode;

/* DLL stuff: */
SYSTEM_integer SYSTEM_exitcode = 0;      /* system.ExitCode */
_P3void_procT  SYSTEM_exitproc = NULL;   /* system.ExitProc */
SYSTEM_boolean _P3islibrary;             /* system.IsLibrary */
SYSTEM_integer SYSTEM_dll_refcount = 0;  /* Refcounting Handles to this DLL */

#if _P3_EXC_MODEL_ == 1
/* C++ try/catch */

static SYSTEM_tobject
_P3_Create_Exception2 (EXCEPTIONS_texceptions code, char *CMsg)
{
  SYSTEM_char buf[257];  /* to store Pascal string */
  int len = (int) strlen(CMsg);

  EXCTST(printf("P3IO.C: Creating Exception Object, Code %d, Message (%d) '%s'\n",
         (int)code, len, CMsg););
  memmove(buf+1, CMsg, *buf = len); /* One-liner to create Pascal string */
  return EXCEPTIONS_create_exception_by_code(code, buf);
} /* _P3_Create_Exception2 */

/* handle exceptions in the C++ try/catch model */
void _P3_Std_Exception_Handler (SYSTEM_tobject exception)
{
#if defined(EXC_DBG)
  printf("ENTERING _P3_Std_Exception_Handler...Exception = %d\n",
         (SYSTEM_nativeuint) exception);
#endif

  if (NULL == exception)
    return;

# if 0
  /* If it's an Exception object but not an EAbort then print its message: */
  if ( _P3_is(exception, &SYSTEM_exception_CD) &&
      !_P3_is(exception, &EXCEPTIONS_eabort_CD)) {
    /* Print the exception's built-in message: */
    printf("P3 Standard Exception Handler: ");
    _P3_write_s0(((SYSTEM_exception)exception)->SYSTEM_exception_DOT_message);
    printf("\n");
    fflush(stdout);
  }
  _P3_ABORT();
# else
  {
    char tbuf[300], *p;
    const unsigned char *m;

    strcpy (tbuf, "P3 Standard Exception Handler: ");
    for (p = tbuf;  *p;  p++)
      ;
    m = ((SYSTEM_exception)exception)->SYSTEM_exception_DOT_message;
    memcpy(p, m+1, *m);
    p[*m] = '\0';
    exit2R(tbuf);
  }
# endif

} /* _P3_Std_Exception_Handler */


#endif  /* if _P3_EXC_MODEL_ == 1 */


/* assume input s is long enough */
static char *
getNoun (int v, char *s)
{
  strcpy (s, "file");
  switch (v) {
  case _P3_ERR_VERB_RMDIR:
  case _P3_ERR_VERB_MKDIR:
  case _P3_ERR_VERB_CHDIR:
    strcpy (s, "directory");
    break;
  } /* switch */
  return s;
} /* getNoun */

/* assume input s is long enough */
static char *
getVerb (int v, char *s)
{
  strcpy (s, "unknown action");
  switch (v) {
  case _P3_ERR_VERB_READ:
    strcpy (s, "read");
    break;
  case _P3_ERR_VERB_WRITE:
    strcpy (s, "write");
    break;
  case _P3_ERR_VERB_FLUSH:
    strcpy (s, "flush");
    break;
  case _P3_ERR_VERB_SEEK:
    strcpy (s, "seek");
    break;
  case _P3_ERR_VERB_EOFCHK:
    strcpy (s, "eof-check");
    break;
  case _P3_ERR_VERB_SEEKEOF:
    strcpy (s, "seekeof");
    break;
  case _P3_ERR_VERB_EOLNCHK:
    strcpy (s, "eoln-check");
    break;
  case _P3_ERR_VERB_SEEKEOLN:
    strcpy (s, "seekeoln");
    break;
  case _P3_ERR_VERB_FPOS:
    strcpy (s, "filepos");
    break;
  case _P3_ERR_VERB_FSIZE:
    strcpy (s, "filesize");
    break;
  case _P3_ERR_VERB_CLOSE:
    strcpy (s, "close");
    break;
  case _P3_ERR_VERB_APPEND:
    strcpy (s, "append");
    break;
  case _P3_ERR_VERB_REWRITE:
    strcpy (s, "rewrite");
    break;
  case _P3_ERR_VERB_RESET:
    strcpy (s, "reset");
    break;
  case _P3_ERR_VERB_ERASE:
    strcpy (s, "erase");
    break;
  case _P3_ERR_VERB_RMDIR:
    strcpy (s, "rmdir");
    break;
  case _P3_ERR_VERB_MKDIR:
    strcpy (s, "mkdir");
    break;
  case _P3_ERR_VERB_CHDIR:
    strcpy (s, "chdir");
    break;
  case _P3_ERR_VERB_UNSPEC:
    strcpy (s, "unspecified op");
    break;
  } /* switch */
  return s;
} /* getVerb */


void _P3_Exception (int eCode, const char *s)
/* This routine is the central one that gets called for
 * every exception that is generated from within P3 (i.e. p3io.c)
 *
 * eCode indicates from which situation it was called,
 * s gives more information in some cases.
 *
 * The routine should have all the info it needs to generate
 * a reasonable Exception Object in the Delphi style.
 *
 *****************************************************************/
{
  char mess[1024] = "\0";
  EXCEPTIONS_texceptions ExcCode = EXCEPTIONS_cnoexception;

  EXCTST(printf("_P3_Exception called, eCode = %d\n", eCode); fflush(stdout);)


  switch (eCode) {
  case _P3_EXC_CODE_INVALIDCAST:
    mess[0] = '\0';             /* just use the message passed in */
    ExcCode = EXCEPTIONS_cinvalidcast;
    break;
  case _P3_EXC_CODE_INOUTERROR:
    /* sprintf(mess, "I/O Error. "); */
    mess[0] = '\0';
    ExcCode = EXCEPTIONS_cinouterror;
    break;
  case _P3_EXC_CODE_RANGEERROR:
    sprintf(mess, "Range check error:  ");
    ExcCode = EXCEPTIONS_crangeerror;
    break;
  case _P3_EXC_CODE_ASSERTIONFAILED:
    mess[0] = '\0';             /* just use the message passed in */
    ExcCode = EXCEPTIONS_cassertionfailed;
    break;
  case _P3_EXC_CODE_ADDRANGE:
    sprintf(mess, "Empty Set Error. "); /* better never happen*/
    ExcCode = EXCEPTIONS_cexception;
    break;
  case _P3_EXC_CODE_OUTOFMEMORY:
    sprintf(mess, "Out of memory");
    ExcCode = EXCEPTIONS_coutofmemory;
    break;
  case _P3_EXC_CODE_ABSTRACTERROR:
    mess[0] = '\0';             /* just use the message passed in */
    ExcCode = EXCEPTIONS_cabstracterror;
    break;
  default:
    sprintf(mess, "Unknown cause. ");
    ExcCode = EXCEPTIONS_cexception;
    break;
  } /* end switch (eCode) */

  if (_P3_EXC_CODE_INOUTERROR == eCode) {
    /* ignore input msg s, construct something from errno/_P3_err */
    /* sprintf(mess, "I/O Error. "); */
    char *sysMsg, *m2 = mess;
    char verb[32];
    char noun[16];
    int k = 0;

    if (_P3_err.verb) {
      k = sprintf (m2, "I/O error on %s of %s", getVerb(_P3_err.verb, verb),
                   getNoun(_P3_err.verb, noun));
    }
    else {
      k = sprintf (m2, "I/O error on file");
    }
    m2 += k;

    if (_P3_err.nam[0]) {       /* non-empty file name */
      k = sprintf (m2, " = '%s'", _P3_err.nam + 1);
      m2 += k;
    }

    if (errno) {
      sysMsg = strerror (errno);
      if (sysMsg)
        k = sprintf (m2, ": %s", sysMsg);
      else
        k = sprintf (m2, ": errno = %d, message not available", errno);
    }
    else if (_P3_err.notOpened) {
      k = sprintf (m2, ": file not open");
    }
    else if ((_P3_err.n > 0) && (5 != _P3_err.n)
             && (NULL != (sysMsg = strerror (_P3_err.n))) ) {
      k = sprintf (m2, ": %s", sysMsg);
    }
    else {
      k = sprintf (m2, ": IOResult = %d", _P3_err.n);
    }
    m2 += k;

    /* I/O-error caught - Reset _P3_err and errno: */
    p3ErrZap (&_P3_err);
    errno = 0;
    /* printf("RESETTING _P3_err AND errno\n"); */
  }
  else {
    if (s)
      strcat(mess, s);
    if (_P3_err.n) {
      /* should never happen!! _P3_err.n only gets set for I/O errors */
      sprintf(mess+strlen(mess), "   IoResult = %d", _P3_err.n);
      p3ErrZap (&_P3_err);
      errno = 0;
    }
  }

  /* We are ready to create and raise an exception now.
   * How we do that depends on the exception model */

#if _P3_EXC_MODEL_ == 1
  _P3_RAISE(_P3_Create_Exception2(ExcCode, mess));
#endif  /* if _P3_EXC_MODEL_ == 1 */

}  /* _P3_Exception */


int
_P3VariableCastError(const char *file, int line, size_t sizeof_a, size_t sizeof_b)
{
  char buf[1024];

  sprintf(buf, "Invalid variable typecast at (%s:%d):"
          " from size %u to size %u",
          file, line, (unsigned int) sizeof_b, (unsigned int) sizeof_a);
  _P3_Exception(_P3_EXC_CODE_INVALIDCAST, buf);
  return 0; /* yeah, right */
}

int _P3_Finalizing = 0; /* Set to 1 during Finalization sequence */

void _P3_halt(int rc)
{
  if (_P3_Finalizing) {
    /* halt called from finalization code. We have choices what to do:
    *  _P3_ABORT: Terminate, allow debugger to pinpoint error point,
    *  exit(rc):  Terminate "normally", but any other modules' finalization
    *             parts are skipped */

#if 0
    exit(rc);
#else
    {
      char tbuf[256];
      sprintf (tbuf, "_P3_halt(%d) called", rc);
      exit2R(tbuf);
    }
#endif

    /* or: _P3_ABORT(); */
    /* or: _P3_Exception("Halt called from finalization part - rc = ", 3, rc); */

  } else {
    _P3_Finalizing = 1; /* To guarantee against infinite loops here */

    /* _P3_DLL_UNWIND calls the exitproc chain, which has
     * _P3_Finalizing (in the main C file, not in p3io.?)
     *  as the last entry. That one calls finalization routines in modules */
    _P3_DLL_UNWIND();

#if 0
    exit(rc);
#else
    {
      char tbuf[256];
      sprintf (tbuf, "_P3_halt(%d) called", rc);
      exit2R(tbuf);
    }
#endif
  }
}


void _P3error_check(/*char* File, int Line*/)
  /* Called after I/O operations in $I+ state. */
{
  if (_P3_err.n) {
    _P3_Exception(_P3_EXC_CODE_INOUTERROR, "");
  }
}

long _P3rangechk(SYSTEM_longint i, SYSTEM_longint lb, SYSTEM_longint ub)
{
  char str[256];

  if (i < lb || i > ub) {
    sprintf(str, FMT_rangechk, i, lb, ub);
    _P3_Exception(_P3_EXC_CODE_RANGEERROR, str);
    return -1; /* Never makes it to here */
  }
  else
    return i;
}

/************** BEGINNING OF I/O ROUTINES *****************************/

/* INTEGER */

/* write integer with two modifiers as corresponding real */

/* Write unsigned integer with no modifiers */
void _P3write_u (_P3file *fil, SYSTEM_longword i)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  IOWRP(fprintf(fil->f,FMT_write_u0,i));
} /* _P3write_u */

/* Write integer with no modifiers */
void _P3write_i (_P3file *fil, SYSTEM_longint i)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  IOWRP(fprintf(fil->f,FMT_write_i0,i));
} /* _P3write_i */

/* Write nativeuint with no modifiers */
void _P3write_n (_P3file *fil, SYSTEM_nativeuint i)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  IOWRP(fprintf(fil->f,FMT_write_n0,i));
} /* _P3write_n */

/* Write uint64 with no modifiers */
void _P3write_y (_P3file *fil, SYSTEM_uint64 i)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  IOWRP(fprintf(fil->f,FMT_write_y0,i));
} /* _P3write_y */

/* Write int64 with no modifiers */
void _P3write_z (_P3file *fil, SYSTEM_int64 i)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  IOWRP(fprintf(fil->f,FMT_write_z0,i));
} /* _P3write_z */

/* Write unsigned integer in field of length len */
void _P3write_ux (_P3file *fil, SYSTEM_longword i, SYSTEM_longint len)
{
  char s[50];

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  sprintf(s,FMT_write_ux, len);
  IOWRP(fprintf(fil->f,s,i));
} /* _P3write_ux */

/* Write integer in field of length len */
void _P3write_ix (_P3file *fil, SYSTEM_longint i, SYSTEM_longint len)
{
  char s[50];

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  sprintf(s,FMT_write_ix, len);
  IOWRP(fprintf(fil->f,s,i));
} /* _P3write_ix */

/* Write nativeuint in field of length len */
void _P3write_nx (_P3file *fil, SYSTEM_nativeuint i, SYSTEM_longint len)
{
  char s[50];

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  sprintf(s,FMT_write_nx, len);
  /* printf("_nx format string: >%s<\n", s); */
  IOWRP(fprintf(fil->f,s,i));
} /* _P3write_nx */

/* Write uint64 in field of length len */
void _P3write_yx (_P3file *fil, SYSTEM_uint64 i, SYSTEM_longint len)
{
  char s[50];

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  sprintf(s,FMT_write_yx, len);
  /* printf("_yx format string: >%s<\n", s); */
  IOWRP(fprintf(fil->f,s,i));
} /* _P3write_yx */

/* Write int64 in field of length len */
void _P3write_zx (_P3file *fil, SYSTEM_int64 i, SYSTEM_longint len)
{
  char s[50];

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  sprintf(s,FMT_write_zx, len);
  /* printf("_zx format string: >%s<\n", s); */
  IOWRP(fprintf(fil->f,s,i));
} /* _P3write_zx */


/* Booleans (referenced as "<exp>? _P3true : _P3false") */
const SYSTEM_byte _P3true [] = "\4TRUE";
const SYSTEM_byte _P3false[] = "\5FALSE";

/* REALS  (double) */

/* write real with no modifiers, using default format */
void _P3write_r (SYSTEM_text *fil, SYSTEM_double d)
{
  char buf[30];
  char *s = buf;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
#if 0
  IOWRP(fprintf(fil->f,"%23.14E",d));
#else
  sprintf (buf,"%22.14E",d);
  if ('\0' == buf[22]) {
    /* assume libraries write in two ways: */
    if ('E' == buf[17]) { /* x.yyyE+eee */
      buf[23] = '\0';
      buf[22] = buf[21];
      buf[21] = buf[20];
      buf[20] = buf[19];
      buf[19] = '0';
    }
    else if (('E' == buf[18]) && ' ' == buf[0]) { /* x.yyyE+ee */
      buf[24] = '\0';
      buf[23] = buf[21];
      buf[22] = buf[20];
      buf[21] = '0';
      buf[20] = '0';
      s++;                      /* skip initial blank */
    }
    /* else we do not recognize this, just leave it alone */
  }
  IOWRP(fprintf(fil->f, "%s", s));
#endif
} /* _P3write_r */

/* write real with one modifier given, using default for second */
void _P3write_rx(SYSTEM_text *fil, SYSTEM_double   d, SYSTEM_longint len1)
{
  SYSTEM_byte s[256];

  _P3_Str_d1(d, len1, s, 255); /* Better not be > 255 at return */
  _P3_writefs0(fil, s);
} /* _P3write_rx */

/* Real with two modifiers */
void _P3write_ry(_P3file *fil, SYSTEM_double d,
                 SYSTEM_longint len1, SYSTEM_longint len2)
{
  SYSTEM_byte s[256];

  _P3_Str_d2(d, len1, len2, s, 255); /* Better not be > 255 at return */
  _P3_writefs0(fil, s);
} /* _P3write_ry */


/* CHAR, STRINGS */

/* Write a char */
void _P3write_c (_P3file *fil, SYSTEM_char c)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  IOWRP(putc(c,fil->f));
} /* _P3write_c */

/* Write a char in field of length len */
void _P3write_cx (_P3file *fil, SYSTEM_char c, SYSTEM_longint len)
{
  char s[50];

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  sprintf (s, FMT_write_cx, len);
  IOWRP(fprintf(fil->f,s,c));
} /* _P3write_cx */

#define CarriageReturn ((SYSTEM_byte)13)
#define LineFeed       ((SYSTEM_byte)10)
#define Tab            ((SYSTEM_byte)10)
#define CtrlZ          ((SYSTEM_byte)26)

void _P3_Readfs0 (_P3file *fil, SYSTEM_byte *s, SYSTEM_byte max)
{
  register int i, ch;
  FILE *f;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return;
  }
  f = fil->f;

  /* Read up to 'max' characters from a file into the string s.
   * If defined(P3UNIX) and the strings ends with CarriageReturn,
   * we'll return only max-1 chars
   */
  for (i=1; i<=max && (ch =getc(f)) != '\n' && (ch != EOF); i++)
    s[i] = ch;
  *s = i-1;
#if defined(P3UNIX) || defined(DJGPP)
  /* remove trailing CarriageReturn from PC files */
  if (*s > 0 && s[*s] == CarriageReturn)  --(*s);
#endif

  if (0 && ch != '\n' && ch != EOF)  /* Skip to, but not past, eoln */
    while ((ch = getc(f)) != '\n' && ch != EOF);

  if (ch == '\n')
    ch = ungetc(ch, f);

  if (ch == EOF) /* EOF or could be an error */
    if (ferror(f)) {
      _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_READ);
    }
} /* _P3_Readfs0 */


void _P3_writeln (void)
{
  if (0 > printf ("\n")) {
    _P3_IO_ERR(_P3_err, errno, &SYSTEM_output, _P3_ERR_VERB_WRITE);
  }
} /* _P3_writeln */


void _P3_writefn (SYSTEM_text *fil)
{
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  if (0 > fprintf (fil->f, "\n")) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_WRITE);
  }
} /* _P3_writefn */


void _P3_write_s0(const SYSTEM_byte *s)
{
#if 0
  int i, nChars;

  nChars = printf ("%.*s", (int)s[0], s+1);
  /* now handle special case where null byte stops printf early */
  for (i = 1+nChars;  i <= *s;  i++)
    putc(s[i],stdout);
  if (ferror(stdout)) {
    _P3_IO_ERR(_P3_err, errno, &SYSTEM_output, _P3_ERR_VERB_WRITE);
  }
#else
  char tbuf[256];
  memcpy(tbuf,s+1,*s);
  tbuf[*s] = '\0';
  msg2R (tbuf);
#endif
} /* _P3_write_s0 */

void _P3_writefs0 (SYSTEM_text *fil, const SYSTEM_byte *s)
{
  int i, nChars;
  FILE *f;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  f = fil->f;

  nChars = 0;
  /* now handle special case where null byte stops fprintf early,
   * or case where we don't use fprintf at all
   */
  for (i = 1+nChars;  i <= *s;  i++)
    putc(s[i],f);
  if (ferror(f)) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_WRITE);
  }
} /* _P3_writefs0 */


void _P3write_sx (_P3file *fil, const SYSTEM_byte *s, SYSTEM_longint len)
{
  int i, nBlanks, nChars;
  FILE *f;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_WRITE);
    return;
  }
  f = fil->f;
  for (i = *s, nBlanks = sizeof(blankBuf)-1;  i < len;  i += nBlanks) {
    if (len - i < nBlanks)
      nBlanks = len - i;
    fprintf (f, "%.*s", nBlanks, blankBuf);
  }
  nChars = fprintf (f, "%.*s", (int)*s, s+1);
  /* now handle special case where null byte stops fprintf early */
  for (i = 1+nChars;  i <= *s;  i++)
    putc(s[i],f);
  if (ferror(f)) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_WRITE);
  }
} /* _P3write_sx */


/* eoln, eof */

/* Must look ahead. */
SYSTEM_boolean _P3_eof(int Iplus, _P3file *fil, const char *File, int Line)
{
  int ch;
  FILE *f = fil->f;
  SYSTEM_boolean res;

  res = SYSTEM_false;
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_EOFCHK);
  }
  else {
    if (feof(f)) {
      res = 1;
    }
    else {
      ch = getc(f);
      if (ferror(f)) {
        _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_EOFCHK);
      }
      if (ch == EOF) {
        res = 1;
      }
      else {
        ungetc(ch,f);
        if (ferror(f)) {
          _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_EOFCHK);
        }
        res = 0;
      }
    }
  }

  if (Iplus) {
    _P3error_check(/*File, Line*/);
  }

  return res;
} /* _P3_eof */

void _P3_Seek (_P3file *fil, SYSTEM_longint offset, SYSTEM_longint origin)
{
  int rc, code;

  /* printf("Seek(1) ... offset %d, origin %d\n", offset, origin);  */
  /* Translate user's origin code to C library's, in case not same... */
  switch (origin) {
  case 0:
    code = SEEK_SET;            /* from beginning  */
    break;
  case 1:
    code = SEEK_CUR;            /* from current pos */
    break;
  case 2:
    code = SEEK_END;            /* from end of file */
    break;
  default:
    _P3_IO_ERR(_P3_err, EINVAL, fil, _P3_ERR_VERB_SEEK);
    return;
  } /* switch */

  /* printf("Seek(2) ... code %d, fil->f %d\n", code, fil->f); */

  if (_P3_ISOPEN(fil)) {        /* open file */
    rc = fseek(fil->f, offset * fil->block_size, code);

    if (rc == -1) {
      _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_SEEK);
    }
  }
  else {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_SEEK);
  }
} /* _P3_Seek */


SYSTEM_boolean _P3_seekeof(int Iplus, _P3file *fil, const char *File, int Line)
{
  int ch;
  FILE *f = fil->f;
  SYSTEM_boolean res;

  res = SYSTEM_false;
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_SEEKEOF);
  }
  else {
    if (feof(f))
      res = 1;
    else {
      while ((ch = getc(f)) == ' ' || ch == '\t' || ch == '\n') {
        if (ferror(f))
          break;
      }
      if (ferror(f)) {
        _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_SEEKEOF);
      }
      if (ch != EOF) {
        ungetc(ch,f);
        if (ferror(f)) {
          _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_SEEKEOF);
        }
        res = 0;
      }
      else
        res = 1;
    }
  }

  if (Iplus) {
    _P3error_check(/*File, Line*/);
  }
  return res;
} /* _P3_seekeof */

SYSTEM_boolean _P3_eoln(int Iplus, _P3file *fil, const char *File, int Line)
{
  int ch;
  FILE *f = fil->f;
  SYSTEM_boolean res;

  res = SYSTEM_false;
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_EOLNCHK);
  }
  else {
    if (feof(f)) {
      res = 1;
    }
    else {
      ch = getc(f);
      if (ferror(f)) {
        _P3_IO_NOTOPEN(_P3_err, errno, fil, _P3_ERR_VERB_EOLNCHK);
      }
      if (ch != EOF) {
        ungetc(ch,f);
        if (ferror(f)) {
          _P3_IO_NOTOPEN(_P3_err, errno, fil, _P3_ERR_VERB_EOLNCHK);
        }
      }
      res = (ch == '\n' || ch == EOF);
    }
  }

  if (Iplus) {
    _P3error_check(/*File, Line*/);
  }
  return res;
} /* _P3_eoln */

SYSTEM_boolean _P3_seekeoln (int Iplus, _P3file *fil, const char *File, int Line)
{
  int ch;
  FILE *f = fil->f;
  SYSTEM_boolean res;

  res = SYSTEM_false;
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_SEEKEOLN);
  }
  else {
    if (feof(f))
      res = 1;
    else {
      while ((ch = getc(f)) == ' ' || ch == '\t') {
        if (ferror(f))
          break;
      }
      if (ferror(f)) {
        _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_SEEKEOLN);
      }
      if (ch != EOF) {
        ungetc(ch,f);
        if (ferror(f)) {
          _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_SEEKEOLN);
        }
        res = (ch == '\n');
      }
      else
        res = 1;
    }
  }

  if (Iplus) {
    _P3error_check(/*File, Line*/);
  }
  return res;
} /* _P3_seekeoln */

#define DUMP_STUFF 0   /* Output stuff for Assign and _P3fileopn */
void _P3block_read_write (_P3file        *fil,
                          void           *buf,
                          size_t          count,
                          SYSTEM_longint *numread,
                          SYSTEM_boolean  wr)
{
  size_t num_transferred;
  FILE *f = fil->f;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, ((wr) ? _P3_ERR_VERB_WRITE : _P3_ERR_VERB_READ));
    return;
  }
#if 0
  printf("BlockRead/Write clld... count %d, f blksz %d, name = %s, sts = %d\n",
         count, fil->block_size, fil->nam+1, fil->status);
#endif

  if (wr)
    num_transferred = fwrite(buf, fil->block_size, count, f);
  else
    num_transferred = fread( buf, fil->block_size, count, f);

  if (NULL == numread) { /* 4th parm not given: Check for and report errors */
    if (ferror(f)) {
      _P3_IO_ERR(_P3_err, errno, fil, ((wr) ? _P3_ERR_VERB_WRITE : _P3_ERR_VERB_READ));
    }
    else {
      if (count != num_transferred) { /* Also an error indication to report */
        int eee;
        eee = feof(f) ? READ_AT_END_OF_FILE : EIO;
        _P3_IO_ERR(_P3_err, eee, fil, ((wr) ? _P3_ERR_VERB_WRITE : _P3_ERR_VERB_READ));
      }
    }
  }
  else /* numread is non-NULL */
    *numread = (SYSTEM_longint) num_transferred; /* always non-negative */

  return;
} /* _P3block_read_write */

void _P3rw_typed (_P3file *fil, void *buf, int wr)
{
  size_t count;
  FILE *f = fil->f;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, ((wr) ? _P3_ERR_VERB_WRITE : _P3_ERR_VERB_READ));
    return;
  }
  if (wr)
    count = fwrite(buf, fil->block_size, 1, f);
  else
    count = fread( buf, fil->block_size, 1, f);
  if (ferror(f)) {
    _P3_IO_ERR(_P3_err, errno, fil, ((wr) ? _P3_ERR_VERB_WRITE : _P3_ERR_VERB_READ));
  }
} /* _P3rw_typed */

void _P3_Assign (_P3file *fil, const SYSTEM_char *s)
{
  /* Initialize fil */
  fil->f = NULL;
  fil->status = _P3CLOSED;
  memmove(fil->nam, s, s[0]+1);
  fil->nam[s[0]+1] = '\0';

  if (DUMP_STUFF)
    printf ("Assigning name >%s< to FILE* at %p\n", fil->nam+1, fil);
} /* _P3_Assign */


/* **************************************************************************
 * _P3fileopn semantics:
 *    When called on text files from (given in s)           fopen str:
 *      Reset:     Open for Read-only.                         "r"
 *      Rewrite:   Open for Write-only.                        "w"
 *      Append:    Open for Write-only, at end                 "a"
 *    When called on typed or untyped files (i.e., non-text):
 *      Reset:     filemode = 0: Open for Read-only,           "r"
 *                 filemode = 1: Open for Write-only, at end   "w"
 *                 filemode = 2: Open for Update (ReadWrite).  "r+"
 *      Rewrite:   Open for Write-only.                        "w"
 *      Append:    Open for Write-only, at end                 "a"
 *
 * Parameters:
 *       s:
 *         _P3RESET,
 *         _P3REWRITE,
 *         _P3APPEND
 *         _P3UPDATE)
 *       t:
 *         0 for text file
 *         1 typed
 *         2 untyped with block_size parameter
 *       block_size:
 *         1 for text,
 *         sizeof(type) for typed
 *         128/user specified for untyped
 *       filemode: as follows:
 * ****************************************************************** */
#define _FM_RO 0 /* read-only  */
#define _FM_WO 1 /* write-only */
#define _FM_RW 2 /* read-write */
void _P3fileopn(_P3file *fil, SYSTEM_longint s, SYSTEM_longint t,
                SYSTEM_longint block_size)
{
  char  str[50];
  unsigned char verb = _P3_ERR_VERB_NONE; /* for exception handling */
  SYSTEM_byte action; /* 0: append, 4: reset, 8: rewrite, 12: update */

  if (DUMP_STUFF)
    printf ("Opening file for %s  >%s<    "
            "  at %p; s: %d, t: %d, blkz %d filemode %d\n",
            (s==_P3RESET?"Read ":
             s==_P3REWRITE?"Write":
             s==_P3APPEND?"Append" :
             s==_P3UPDATE? "Update" : " ??? "),
            fil->nam+1, fil, s, t, block_size, SYSTEM_filemode);

  fil->f = NULL;
  fil->status = _P3CLOSED;
  fil->block_size = block_size;

  if (SYSTEM_filemode != _FM_RO && SYSTEM_filemode != _FM_WO)
    SYSTEM_filemode = _FM_RW;

  /* Set action to perform on file */
  /* Weird case: s = _P3RESET and filemode = 1 (_FM_RO):
     Do append, not write which erases the file */
  switch (s) {
    case _P3APPEND:
      action = _P3APPEND ;
      verb = _P3_ERR_VERB_APPEND;
      break;
    case _P3REWRITE:
      action = _P3REWRITE;
      verb = _P3_ERR_VERB_REWRITE;
      break;
    case _P3RESET:
      verb = _P3_ERR_VERB_RESET;
      if (t == 0)
        action = _P3RESET;
      else {
        switch (SYSTEM_filemode) {
        case _FM_RO:
          action = _P3RESET;
          break;
        case _FM_WO:
          action = _P3APPEND;
          break;
        case _FM_RW:
          action = _P3UPDATE;
          break;
        } /* switch(SYSTEM_filemode) */
      }
      break;
  } /* switch(s) */

  if (DUMP_STUFF)
    printf("       action = %d\n", action);

  /* Done setting action */
  if (fil->nam[0] > 0) {  /* Non-empty file name */
    errno = 0;

    /* BEGIN SETTING fopen_string str */
#ifdef P3DOS
    if (t == 0) /* text file */
      switch (action) {
        case _P3APPEND : strcpy(str,"a"); break;
        case _P3RESET  : strcpy(str,"r"); break;
        case _P3REWRITE: strcpy(str,"w"); break;
/* NOTE: r+ below works for Paul's objtest.dpr, but w+ makes
   it write to stdout!!! FIGURE OUT WHY.                     */
        case _P3UPDATE : strcpy(str,"r+"); break;
      }
    else /* DOS, binary file */
      switch (action) {
        case _P3APPEND : strcpy(str,"ab"); break;
        case _P3RESET  : strcpy(str,"rb"); break;
        case _P3REWRITE: strcpy(str,"wb"); break;
/* NOTE: r+ below works for Paul's objtest.dpr, but w+ makes
   it write to stdout!!! FIGURE OUT WHY.                     */
        case _P3UPDATE : strcpy(str,"rb+"); break;
      }
#else
    /* Unix */
    switch (action) {
      case _P3APPEND : strcpy(str,"a"); break;
      case _P3RESET  : strcpy(str,"r"); break;
      case _P3REWRITE: strcpy(str,"w"); break;
/* NOTE: r+ below works for Paul's objtest.dpr, but w+ makes
   it write to stdout!!! FIGURE OUT WHY.                     */
      case _P3UPDATE : strcpy(str,"r+"); break;
    }
#endif

    /* END SETTING fopen_string str */

    if (DUMP_STUFF)
      printf("Fileopn %s s = %d, t = %d , str >%s< \n", fil->nam+1, s, t, str);

    fil->f = fopen((const char*)fil->nam+1, str);
    _P3_SETMODE(fil,s);         /* changed below if errors */

    if (fil->f == NULL) {
      _P3_IO_ERR(_P3_err, errno, fil, verb);
      fil->status = _P3CLOSED;
    }
    else {

#if defined(P3UNIX) || defined(DJGPP)
      int rCode;
      struct stat statBuf;

      rCode = fstat(fileno(fil->f), &statBuf);
      if (rCode != 0) {
        (void) fclose (fil->f);
        fil->f = NULL;
        fil->status = _P3CLOSED;
        _P3_IO_ERR(_P3_err, errno, fil, verb);
      }
      else if ((S_IFMT&statBuf.st_mode) == S_IFDIR) {
        (void) fclose (fil->f);
        fil->f = NULL;
        fil->status = _P3CLOSED;
        /* is a directory! */
        _P3_IO_ERR(_P3_err, EISDIR, fil, verb);
      }
#endif
    }

    /* if needed: if (*f != NULL) setvbuf(*f, NULL, _IOFBF, 4*1024); */

  }
  else { /* Empty file name. Assign to stdin/stdout, even if non-text file */
    exit2R("Assign to empty file name not allowed: should be silent instead");
  }
} /* _P3fileopn */

void _P3_Close(_P3file *fil)
{
  /* it is an error to close a non-opened file: includes closed or unassigned */
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_CLOSE);
    return;
  }
  errno = 0;
  if (fil->f) {
    if (fclose(fil->f)) {
      _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_CLOSE);
    }
  }
  fil->f = NULL;
  fil->status = _P3CLOSED;
} /* _P3_Close */

void _P3_Flush (_P3file *fil)
{
  /* it is an error to flush a non-opened file: includes closed or unassigned */
  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_FLUSH);
    return;
  }
  errno = 0;
  if (fflush(fil->f)) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_FLUSH);
  }
  return;
} /* _P3_Flush */

void _P3_Erase(_P3file *fil)
{
  /* in Kylix it works to erase an open file but not an unassigned one */
  if (_P3_NOTASSIGNED(fil)) {
    _P3_IO_ERR(_P3_err, ENOENT, fil, _P3_ERR_VERB_ERASE);
    return;
  }
  p3ErrZap (&_P3_err);
  /* NOTE: remove gives trouble on old cc's. Use unlink then. */
  if (remove((char *)(fil->nam+1))) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_ERASE);
  }
} /* _P3_Erase */

void SYSTEM_rmdir(const SYSTEM_char *s)
{
  char buf[256];
  memmove(buf, s+1, *s);
  buf[*s] = '\0';

  if (rmdir(buf))
    p3ErrCSet (&_P3_err, errno, s, _P3_ERR_VERB_RMDIR);
}

void SYSTEM_mkdir(const SYSTEM_char *s)
{
  char buf[256];
  memmove(buf, s+1, *s);
  buf[*s] = '\0';

#if defined(P3DOS)
  if (mkdir(buf))
    p3ErrCSet (&_P3_err, errno, s, _P3_ERR_VERB_MKDIR);
#else
  /* Check this is correct UNIX, i.e. 0777 */
  if (mkdir(buf,0777))
    p3ErrCSet (&_P3_err, errno, s, _P3_ERR_VERB_MKDIR);
#endif
}

void SYSTEM_chdir(const SYSTEM_char *s)
{
  char buf[256];
  memmove(buf, s+1, *s);
  buf[*s] = '\0';

  if (chdir(buf))
    p3ErrCSet (&_P3_err, errno, s, _P3_ERR_VERB_CHDIR);
}

/* SPD: P3_IO_ tested to here */

int SYSTEM_ioresult(void)
{ /* Returns _P3_err.n which is then cleared */
  int i = _P3_err.n;

  p3ErrZap (&_P3_err);
  return i;
}

/* A bunch of small functions that used to be macros in p3io.h
   but now have error checking: */

SYSTEM_longword _P3read_u (_P3file *fil)
{
  FILE *f;
  SYSTEM_longword res;
  int rc;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return 0;
  }
  f = fil->f;
  rc = fscanf(f, FMT_read_u, &res);
  if (rc == 0)   _P3_err.n = READ_CONVERSION_ERROR;
  if (rc == EOF) _P3_err.n = READ_AT_END_OF_FILE;
  if (ferror(f)) _P3_err.n = errno;
  if (_P3_err.n)
    _P3_IO_ERR(_P3_err, _P3_err.n, fil, _P3_ERR_VERB_READ);
  return res;
} /* _P3read_u */

SYSTEM_longint _P3read_i (_P3file *fil)
{
  FILE* f = fil->f;
  SYSTEM_longint res;
  int rc;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return 0;
  }
  rc = fscanf(f, FMT_read_i, &res);
  if (rc == 0)   _P3_err.n = READ_CONVERSION_ERROR;
  if (rc == EOF) _P3_err.n = READ_AT_END_OF_FILE;
  if (ferror(f)) _P3_err.n = errno;
  if (_P3_err.n)
    _P3_IO_ERR(_P3_err, _P3_err.n, fil, _P3_ERR_VERB_READ);
  return res;
} /* _P3read_i */

SYSTEM_nativeuint _P3read_n (_P3file *fil)
{
  FILE* f = fil->f;
  SYSTEM_nativeuint res;
  int rc;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return 0;
  }
  rc = fscanf(f, FMT_read_n, &res);
  if (rc == 0)   _P3_err.n = READ_CONVERSION_ERROR;
  if (rc == EOF) _P3_err.n = READ_AT_END_OF_FILE;
  if (ferror(f)) _P3_err.n = errno;
  if (_P3_err.n)
    _P3_IO_ERR(_P3_err, _P3_err.n, fil, _P3_ERR_VERB_READ);
  return res;
} /* _P3read_n */

SYSTEM_uint64 _P3read_y (_P3file *fil)
{
  FILE* f = fil->f;
  SYSTEM_uint64 res;
  int rc;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return 0;
  }
  rc = fscanf(f, FMT_read_y, &res);
  if (rc == 0)   _P3_err.n = READ_CONVERSION_ERROR;
  if (rc == EOF) _P3_err.n = READ_AT_END_OF_FILE;
  if (ferror(f)) _P3_err.n = errno;
  if (_P3_err.n)
    _P3_IO_ERR(_P3_err, _P3_err.n, fil, _P3_ERR_VERB_READ);
  return res;
} /* _P3read_y */

SYSTEM_int64 _P3read_z (_P3file *fil)
{
  FILE* f = fil->f;
  SYSTEM_int64 res;
  int rc;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return 0;
  }
  rc = fscanf(f, FMT_read_z, &res);
  if (rc == 0)   _P3_err.n = READ_CONVERSION_ERROR;
  if (rc == EOF) _P3_err.n = READ_AT_END_OF_FILE;
  if (ferror(f)) _P3_err.n = errno;
  if (_P3_err.n)
    _P3_IO_ERR(_P3_err, _P3_err.n, fil, _P3_ERR_VERB_READ);
  return res;
} /* _P3read_z */

SYSTEM_char _P3read_c (_P3file *fil)
{
  FILE* f = fil->f;
  int res;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return (SYSTEM_char) 0;
  }
  if ((res = getc(f)) < 0) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_READ);
    return (SYSTEM_char) 0;  /* or 255, maybe? */
  }
  return (SYSTEM_char) res;
} /* _P3read_c */


SYSTEM_double _P3read_d (_P3file *fil)
{
  FILE* f = fil->f;
  SYSTEM_double res;
  int rc;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return 0;
  }
  rc = fscanf(f, "%lf", &res);
  if (rc == 0)   _P3_err.n = READ_CONVERSION_ERROR;
  if (rc == EOF) _P3_err.n = READ_AT_END_OF_FILE;
  if (ferror(f)) _P3_err.n = errno;
  if (_P3_err.n)
    _P3_IO_ERR(_P3_err, _P3_err.n, fil, _P3_ERR_VERB_READ);
  return res;
} /* _P3read_d */


void _P3read_ln (_P3file *fil)
{
  FILE *f = fil->f;
  int ch;

  if (! _P3_ISOPEN(fil)) {
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_READ);
    return;
  }
  while ((ch = getc(f)) != '\n' && ch != EOF)
    ;                           /* do nothing */
  if (ferror(f)) {
    _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_READ);
  }
} /* _P3read_ln */

/************** END OF I/O ROUTINES *****************************/

/* implement Delphi setlength(s, newLen)
 * N.B.: Delphi just takes the low-order byte of newLen and assigns that,
 * regardless of the size of s
 */
void _P3setlength(SYSTEM_byte *s, SYSTEM_integer newLen,
                  SYSTEM_integer siz)
{
  /* if (newLen > siz)  newLen = siz;   */
  s[0] = (SYSTEM_byte) newLen;
} /* _P3setlength */

/* delphi: res := copy(s, i, cnt)
 * N.B.: res and s could be the same string
 */
SYSTEM_char *SYSTEM_copy(SYSTEM_byte *res, SYSTEM_byte max,
                         const SYSTEM_byte *s,
                         SYSTEM_integer i, SYSTEM_integer cnt)
{
  if ((i < 1) || (i > *s)) {
    *res = '\0';
    return res;
  }
  if (cnt > *s - i + 1)
    cnt = *s - i + 1;
  if (cnt > max)
    cnt = max;                  /* truncate resulting string */
  memmove(res+1, s+i, cnt);
  *res = (SYSTEM_byte) cnt;
  return res;
} /* SYSTEM_copy */

/* Delphi's system.delete */
void SYSTEM_delete(SYSTEM_byte *s, SYSTEM_integer i, SYSTEM_integer cnt)
{
  SYSTEM_integer j, k;

  if ((i <= 0) || (i > *s))
    return;
  if (cnt <= 0)
    return;
  if ((cnt > *s) || (i + cnt > *s)) {
    *s = i - 1;
    return;
  }

  j = *s - cnt + 1; /* last chars to move to */
  for (k = i;  k < j;  k++)
    s[k] = s[k+cnt];
  *s -= cnt;
} /* SYSTEM_delete */

SYSTEM_char SYSTEM_upcase(SYSTEM_char c)
{
  /* Don't try this on EBCDIC machines */
  if (c >= 'a' && c <= 'z')
    return c-('a'-'A');
  else
    return c;
}

SYSTEM_char SYSTEM_locase(SYSTEM_char c)
{
  /* Don't try this on EBCDIC machines */
  if (c >= 'A' && c <= 'Z')
    return c+('a'-'A');
  else
    return c;
}


/* implementation of Delphi's system.insert(sub, dest, i)
 */
void _P3_insert(const SYSTEM_byte *sub, SYSTEM_byte *dest,
                SYSTEM_byte maxDest, SYSTEM_integer i)
{
  int lenSub, m, n, k, stop;

  if (0 == (lenSub = *sub))     /* empty substring, nothing to do */
    return;
  m = *dest;
  if (i < 1)
    i = 1;                      /* like Delphi and Turbo */
  if (i > m)
    i = m+1;                    /* like Delphi and Turbo */

  if ((i-1 + lenSub) >= maxDest) {
    /* just copy in sub, ignore the tail */
    stop = maxDest - (i-1);
    for (k = 1;  k <= stop;  k++)
      dest[i+k-1] = sub[k];
    *dest = maxDest;
    return;
  }
  /* copy tail of dest to make space for sub */
  if ((n=m+lenSub) > maxDest)
    n = maxDest;
  *dest = n;
  for (k = n, stop = i + lenSub;  k >= stop;  --k)
    dest[k] = dest[k-lenSub];

  /* copy sub into the hole in dest */
  for (k = 1;  k <= lenSub;  k++)
    dest[i-1+k] = sub[k];
} /* _P3_insert */

/* return location of sub in s
 * on success: 0 < k <= length(s)
 * on failure: 0
 */
SYSTEM_integer SYSTEM_pos(const SYSTEM_byte *sub, const SYSTEM_byte *s)
{
  SYSTEM_longint j, k, stop;

  if (*sub == 1) { /* Common case, faster this way: (?) */
    for (k = 1;  k <= *s;  k++)
      if (sub[1] == s[k])
        return k;
    return 0;
  }

  stop = (SYSTEM_longint)*s - *sub + 1;
  for (k = 1;  k <= stop;  k++) {
    for (j = 1;  j <= *sub;  j++)
      if (sub[j] != s[k+j-1])
        goto next;

    return k;
  next:
    ;
  } /* k loop */

  return 0;
} /* SYSTEM_pos */


/* SET OPERATIONS */
/* The implementation requires >= 8 bits per byte */

_P3Tset_ptr _P3set_copy(SYSTEM_longint len, _P3Tset_ptr to,
                        const _P3set_elem *fr)
{
  SYSTEM_longint i;

#if defined(CATCH_EMPTY_SET)
  if (fr == _P3empty_set) {
    for (i = 0;  i < len;  i++)
      to[i] = (_P3set_elem)0;
    return to;
  }
#endif
  if (fr != to)
    for (i = 0;  i < len;  i++)
      to[i] = fr[i];
  return to;
} /* _P3set_copy */

/* set addition implemented as bitwise OR */
_P3Tset_ptr _P3set_p(SYSTEM_longint len, _P3Tset_ptr ret,
                     const _P3set_elem *s1, const _P3set_elem *s2)
{
  SYSTEM_longint i;

#if defined(CATCH_EMPTY_SET)
  if (s1 == _P3empty_set)
    return _P3set_copy(len, ret, s2);
  else if (s2 == _P3empty_set)
    return _P3set_copy(len, ret, s1);
#endif

  for (i = 0;  i < len;  i++)
    ret[i] = s1[i] | s2[i];

  return ret;
} /* _P3set_p */

/* set difference s1 - s2 implemented as bitwise s1 AND NOT s2  */
_P3Tset_ptr _P3set_m(SYSTEM_longint len, _P3Tset_ptr ret,
                     const _P3set_elem *s1, const _P3set_elem *s2)
{
  SYSTEM_longint i;

#if defined(CATCH_EMPTY_SET)
  if (s1 == _P3empty_set || s2 == _P3empty_set)
    return _P3set_copy(len, ret, s1);
#endif

  for (i = 0;  i < len;  i++)
    ret[i] = s1[i] & ~s2[i];

  return ret;
} /* _P3set_m */

/* set multiplication implemented as bitwise AND */
_P3Tset_ptr _P3set_t(SYSTEM_longint len, _P3Tset_ptr ret,
                     const _P3set_elem *s1, const _P3set_elem *s2)
{
  SYSTEM_longint i;

#if defined(CATCH_EMPTY_SET)
  if (s1 == _P3empty_set || s2 == _P3empty_set) {
    for (i = 0;  i < len;  i++)
      ret[i] = (_P3set_elem)0;
    return ret;
  }
#endif
  for (i = 0;  i < len;  i++)
    ret[i] = s1[i] & s2[i];

  return ret;
} /* _P3set_t */

_P3Tset_ptr _P3set_expand(SYSTEM_longint toLen, _P3Tset_ptr to,
                          SYSTEM_longint frLen, const _P3set_elem *fr)
{
  SYSTEM_longint i;

  /* printf("in _P3set_expand, toLen=%d  frLen=%d\n", toLen, frLen); */
  if (toLen < frLen) {
    for (i = 0;  i < toLen;  i++)
      to[i] = fr[i];
    return to;                  /* could this ever happen?? */
  }

  /* from should never be the empty set:
   * it's already full size so why would we expand it? */

  for (i = 0;  i < frLen;  i++)
    to[i] = fr[i];
  for (i = frLen;  i < toLen;  i++)
    to[i] = (_P3set_elem)0;
  return to;
} /* _P3set_expand */


/* check if  (i in s), where s is [0..mx]
 Has an inline variant, _P3SET_ic                       */
SYSTEM_boolean _P3set_i(SYSTEM_longint mx, SYSTEM_longint i,
                        const _P3set_elem *s)
{
#if defined(CATCH_EMPTY_SET)
  if (s == _P3empty_set)
    return 0;
#endif

  return (i < 0 || i > mx) ? 0 :
    ( (s)[i/_P3bits_per_elem]>>(i%_P3bits_per_elem) ) & 1;
} /* _P3set_i */

/* add the element i to the set s[0..mx] */
_P3Tset_ptr _P3set_add_elem(SYSTEM_longint mx, _P3Tset_ptr s,
                            _P3set_elem i)
{
#if defined(CATCH_EMPTY_SET)
  if (s == _P3empty_set)
    _P3_Exception(_P3_EXC_CODE_ASSERTIONFAILED,
                  "Constant empty set passed in to _P3set_add_elem");
#endif
  if (i <= mx) {
    s[i/_P3bits_per_elem] |= (1L << (i%_P3bits_per_elem));
  }
  return s;
} /* _P3set_add_elem */

/**********
void dump_set(int len, _P3Tset_ptr s)
{ int i;
  for (i=0; i<=len/_P3bits_per_elem; i++) printf(" %d", s[i]);
  printf("\n");
}
***********/

/* add the range [lo..up] to the set s[0..mx] */
_P3Tset_ptr _P3set_add_range(SYSTEM_longint mx, _P3Tset_ptr s,
                             _P3set_elem lo, _P3set_elem up)
{
  int i;

#if defined(CATCH_EMPTY_SET)
  if (s == _P3empty_set)
    _P3_Exception(_P3_EXC_CODE_ADDRANGE, "P_P3set_add_range");
#endif
  if (up > mx)
    up = mx;
  if (lo > up)
    return s;
  for (i = lo;  i <= up;  i++)
    s[i/_P3bits_per_elem] |= (1L << (i%_P3bits_per_elem));
  return s;
} /* _P3set_add_range */

/*  MATH ETC. */

double SYSTEM_int(SYSTEM_double x) /* Truncate towards zero */
{
  return x >= 0? floor(x) : ceil(x);
}

double SYSTEM_frac(SYSTEM_double x) /* Signed fraction; frac(x) = x - int(x) */
{
  return x - (x >= 0? floor(x) : ceil(x));
}

SYSTEM_int64 SYSTEM_round(double x)
{
  if (x >= 0)
    return (SYSTEM_int64) (x+0.5);
  else
    return (SYSTEM_int64) (x-0.5);
}

SYSTEM_int64 SYSTEM_sqr_i(SYSTEM_int64 i)
{
   return i*i;
}

SYSTEM_double SYSTEM_sqr_r(SYSTEM_double r)
{
  return r*r;
}

SYSTEM_integer SYSTEM_allocmemcount = 0, SYSTEM_allocmemsize   = 0;
SYSTEM_int64                             SYSTEM_allocmemsize64 = 0;

void _P3_new(void **p, SYSTEM_longint s)
{
  /* printf("_P3_new: %d\n", s);  */
  if (s > 0) {
    *p = (void*) malloc(s);
    if (!*p)
      _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
    SYSTEM_allocmemcount++;
    SYSTEM_allocmemsize64 += s;
    SYSTEM_allocmemsize = (SYSTEM_integer) SYSTEM_allocmemsize64;
  }
  else
    *p = NULL;
} /* _P3_new */

/* _P3_new64 only gets called on 64-bit machines!! */
void _P3_new64(void **p, SYSTEM_int64 s)
{
  /* printf("_P3_new64: " FMT_write_z0 "\n", s);  */
  if (s > 0) {
    *p = (void*) malloc((size_t)s);
    if (!*p)
      _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
    SYSTEM_allocmemcount++;
    SYSTEM_allocmemsize64 += s;
    SYSTEM_allocmemsize = (SYSTEM_integer) SYSTEM_allocmemsize64;
  }
  else
    *p = NULL;
} /* _P3_new64 */

void _P3_free(void *p, SYSTEM_longint s)
{
  if (p) {
    SYSTEM_allocmemcount--;
    SYSTEM_allocmemsize64 -= s; /* Note: freemem calls with s = 0 */
    SYSTEM_allocmemsize = (SYSTEM_integer) SYSTEM_allocmemsize64;
    free(p); /* Can't nullify p, only free the block */
  }
} /* _P3_free */

void _P3_free64(void *p, SYSTEM_int64 s)
{
  if (p) {
    SYSTEM_allocmemcount--;
    SYSTEM_allocmemsize64 -= s;
    SYSTEM_allocmemsize = (SYSTEM_integer) SYSTEM_allocmemsize64;
    free(p); /* Can't nullify p, only free the block */
  }
} /* _P3_free64 */

SYSTEM_tobject _P3_alloc_object(const SYSTEM_classdescriptor_t *CD)
{
  SYSTEM_tobject p;
/* printf("_P3_alloc_object: allocating %s of %ld bytes\n",
       CD->name+1, CD->CLS_size);   */

  _P3getmem(p, CD->CLS_size);
  memset(p, '\0', CD->CLS_size);
  p->CD = CD;
  return p;
} /* _P3_alloc_object */

void _P3_dealloc_object(void *obj)
{ /* SYSTEM_classreference_t clsref = ((SYSTEM_tobject)obj)->CD;
     printf("_P3_dealloc_object: de_allocating %s of %ld bytes\n",
       clsref->name+1, clsref->CLS_size); fflush(stdout); */
  if (obj)
    free(obj);
}

void SYSTEM_reallocmem(void **p, SYSTEM_longint s)
{
  /* printf("Call of ReAllocMem %d\n", s);  */
  /* How to keep track of sizes???   */

  if (s <= 0) {
    if (*p != NULL) {
      _P3_free(*p, 0);
      *p = NULL; /* to be sure - ReAllocMem must do this (as Delphi) */
    }
  }
  else {
    if (*p != NULL) {
      void *pNew;
      pNew = realloc(*p, s);
      if (!pNew)
        _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
      *p = pNew;
    }
    else {
      _P3_new(p, s);
    }
  }
} /* SYSTEM_reallocmem */

void SYSTEM_reallocmem64 (void **p, SYSTEM_int64 s)
{
  /* printf("Call of ReAllocMem64" FMT_write_z0 "\n", s);  */
  if (s <= 0) {
    if (*p != NULL) {
      _P3_free(*p, 0);
      *p = NULL; /* to be sure - ReAllocMem must do this (as Delphi) */
    }
  }
  else {
    if (*p != NULL) {
      void *pNew;
      pNew = realloc(*p, (size_t) s);
      if (!pNew)
        _P3_Exception(_P3_EXC_CODE_OUTOFMEMORY, "");
      *p = pNew;
    }
    else {
      _P3_new64(p, s);
    }
  }
} /* SYSTEM_reallocmem64 */

/* mess is a combo-string: leading length byte AND trailing nil byte */
void _P3assert(const SYSTEM_char *mess, const char *File, int Line)
{
  char buf[1024];

  if (0 == *mess) {
    sprintf(buf, "Assertion failure (%s:%d)", File, Line);
  }
  else {
    sprintf(buf, "%s (%s:%d)", mess+1, File, Line);
  }
  _P3_Exception(_P3_EXC_CODE_ASSERTIONFAILED, buf);
} /* _P3assert */



/* STRCPY and MEMCPY things */

SYSTEM_char *_P3_pcharn2str (SYSTEM_char *dst, SYSTEM_char dstSiz,
                             const SYSTEM_char *src, int srcLen)
{
  *dst = (srcLen > dstSiz) ? dstSiz : (SYSTEM_char) srcLen;
  memmove (dst+1, src, *dst);
  return dst;
} /* _P3_pcharn2str */

SYSTEM_char *_P3_pchar2str (SYSTEM_char *dst, SYSTEM_char dstSiz,
                            const SYSTEM_char *src)
{
  size_t srcLen;

  srcLen = strlen ((char *)src);
  *dst = (srcLen > dstSiz) ? dstSiz : (SYSTEM_char) srcLen;
  memmove (dst+1, src, *dst);
  return dst;
} /* _P3_pchar2str */

SYSTEM_char *_P3pa2str(SYSTEM_char *s, SYSTEM_char m,
                       const SYSTEM_char *p, SYSTEM_char n)
{
  memmove(s+1, p, *s = (n>m)?m:n);
  return s;
} /* _P3pa2str */

SYSTEM_char *_P3_ch2str(SYSTEM_char *st,SYSTEM_byte max,
                        SYSTEM_char ch)
{
  *st = 1;
  st[1] = ch;
  return st;
} /* _P3_ch2str */

SYSTEM_byte* _P3_strcpy(SYSTEM_byte *d, SYSTEM_integer max,
                        const SYSTEM_byte *s)
{
  SYSTEM_integer i;

  if (max > *s)
    max = *s;      /* Truncate if source too long */

  if (d > s) {                  /* Non-destructive copy ... */
    for (i = max; i > 0; i--)
      d[i] = s[i];
  }
  else if (d < s) {
    for (i = 1; i <= max; i++)
      d[i] = s[i];
  }
  /* else  -- who knows, it just might happen that d == s */
  *d = max;
  return d;
} /* _P3_strcpy */

/* Assumption: r==p1 or r and p1 are non-overlapping. Same with p2 */
SYSTEM_char *_P3_strcat(SYSTEM_char *r, SYSTEM_char max,
                        const SYSTEM_char *p1, const SYSTEM_char *p2)
{
  int len1 = *p1, len2 = *p2, i;

  if (len1 > max) {
    /* special case: implies r != p1  */
    for (i = 1;  i <= max;  i++)
      r[i] = p1[i];
    *r = max;
    return r;
  } /* special case */

  if (len1 + len2 > max)
    len2 = max - len1;          /* Truncate if too long */

  if (r == p1) {                /* Append to r */
    for (i = 1;  i <= len2;  i++)
      r[len1 + i] = p2[i];      /* OK even if r == p2, too */
  }
  else {
    /* General case, safe to assume no overlaps except r == p2 */
    for (i = len2;  i >= 1;  i--)
      r[len1 + i] = p2[i];      /* OK if r == p2 */
    for (i = 1;  i <= len1;  i++)
      r[i] = p1[i];
  }
  *r = len1 + len2;
  return r;
} /* _P3_strcat */

/* Block compares. Used for arrays and records.         */
/* Careful: Must return 0 or 1, not just 0 vs non-zero. */

/* returns 1 if strings are equal, 0 o/w */
SYSTEM_boolean _P3streq(const SYSTEM_char *p1, const SYSTEM_char *p2)
{
  SYSTEM_longword i, len1 = *p1;

  if (p1 == p2)        /* easy - of course they may both be garbage */
    return 1;

  for (i = 0;  i <= len1;  i++) /* start with length byte */
    if (p1[i] != p2[i])
      return 0;
  return 1;

  /* or simply: return (memcmp(p1,p2,*p1+1) == 0); but that's twice as slow! */
} /* _P3streq */

/* returns 1 if strings are equal, 0 o/w, ignoring case */
SYSTEM_boolean _P3streq_ic(const SYSTEM_char *p1, const SYSTEM_char *p2)
{
  SYSTEM_longword i, len1 = *p1;

  if (p1 == p2)        /* easy - of course they may both be garbage */
    return 1;

  if (p1[0] != p2[0])
    return 0;
  for (i = 1;  i <= len1;  i++) {
    char
      ch1 = p1[i],
      ch2 = p2[i];
    if (ch1 >= 'A' && ch1 <= 'Z')
      ch1 += 'a'-'A';
    if (ch2 >= 'A' && ch2 <= 'Z')
      ch2 += 'a'-'A';
    if (ch1 != ch2)
      return 0;
  }
  return 1;                     /* equal */
} /* _P3streq_ic */

/* Returns 0 if strings are equal, >0 if p1 > p2, <0 if p1 < p2. */
int _P3strcmp(const SYSTEM_char *p1, const SYSTEM_char *p2)
{
  int len, i;

  if (p1 == p2)        /* easy - of course they may both be garbage */
    return 0;

  /* ALTERNATIVE: Raw coding to avoid calling memcmp: 25% faster on Linux. */
  len = (*p1 > *p2) ? *p2 : *p1;
  for (i = 1;  i <= len;  i++) {
    if (p1[i] != p2[i]) {
      return p1[i] - p2[i];
    }
  }
  return *p1 - *p2;
} /* _P3strcmp */

/* Compare string st to packed array pa */
int _P3stpcmp(const SYSTEM_char *st, const SYSTEM_char *pa, int lgt)
{
  int len = (*st > lgt) ? lgt : *st;
  register int i;
  int stLen;

  stLen = *st++;
  for (i = 0;  i < len;  i++)
    if (*(st++) != *(pa++))
      return *(--st) - *(--pa);

  if (stLen > len)
    return 1;

  return 0;                     /* same len, same content */
} /* _P3stpcmp */

/* compare string to char */
int _P3stccmp(const SYSTEM_char *st, SYSTEM_char ch)
{
  if (*st == 0)
    return -1;
  if (*st == 1)
    return st[1]-ch;
  if (st[1] == ch)
    return 1;
  
  return st[1]-ch;
} /* _P3stccmp */


/* Str and Val stuff due to EK, Nov 99 */

/* buf should be pretty long, long enough so that the
 * buf[max] setting is safe */
static void _P3_Cstr2Pstr(SYSTEM_byte *s, char *buf, SYSTEM_byte sMax)
{
  int i;
  buf[sMax] = '\0';
  for (i = 0;  buf[i];  i++)
    s[i+1] = buf[i];
  s[0] = i;
} /* _P3_Cstr2Pstr */

void _P3_Str_i0(SYSTEM_integer i, SYSTEM_byte *s, SYSTEM_byte sMax)
{
  char buffer[1024];

  sprintf(buffer,FMT_Str_i0,i);
 _P3_Cstr2Pstr(s, buffer, sMax);
} /* _P3_Str_i0 */

void _P3_Str_i1(SYSTEM_integer i, SYSTEM_integer width, SYSTEM_byte *s,
                SYSTEM_byte sMax)
{
  char buffer[1024];
  char fmt[1024];

  sprintf(fmt,FMT_Str_i1,width);
  sprintf(buffer,fmt,i);
  _P3_Cstr2Pstr(s, buffer, sMax);
} /* _P3_Str_i1 */

void _P3_Str_d0(SYSTEM_double x, SYSTEM_byte *s, SYSTEM_byte sMax)
{
  char buffer[1024];

  sprintf(buffer,"%23.14E",x);
  _P3_Cstr2Pstr(s, buffer, sMax);
} /* _P3_Str_d0 */

void _P3_Str_d1(SYSTEM_double x, SYSTEM_integer width, SYSTEM_byte *s,
                SYSTEM_byte sMax)
{
  char fmt[1024];
  char buffer[1024];
  SYSTEM_longint decimals;

  if (width < 10)
    width = 10;

  decimals = width-8;    /* The 8 are non-decimals in: " 3.33E-001" */
  if (decimals > 18)
    decimals = 18;    /* Delphi outputs max 17, we have room for 18 */
  sprintf(fmt,FMT_Str_d1, width, decimals);

  sprintf(buffer,fmt,x);
  _P3_Cstr2Pstr(s, buffer, sMax);
} /* _P3_Str_d1 */

void _P3_Str_d2 (SYSTEM_double x, SYSTEM_integer width,
                 SYSTEM_integer decimals, SYSTEM_byte *s, SYSTEM_byte sMax)
{
  char buffer[1024];
  char fmt[1024];

  sprintf(fmt,FMT_Str_d2,width,decimals);
  if (fabs(x) > 1.0e37)
    sprintf(fmt,"%%%d.%dE",width,decimals);

  sprintf(buffer,fmt,x);
  _P3_Cstr2Pstr(s, buffer, sMax);
} /* _P3_Str_d2 */


/* valid strings look like: [+|-][0x|$]d+,
 * where d is a decimal digit unless preceded by the 0x or $,
 * in which case it is a hex digit
 */
void _P3_Val_i(const SYSTEM_byte *s, SYSTEM_integer *i, SYSTEM_integer *code)
{
  char buffer[256];
  SYSTEM_byte *end, *s2, *sd;
  long int li;
  int len;
  int sign = 1;

  len = s[0];
  strncpy (buffer, (char *)(s+1), len);
  buffer[len] = '\0';

  /* skip over blanks
   * - Kylix 3 does not treat any other chars as whitespace
   */
  for (s2 = (SYSTEM_byte *)buffer;  ' ' == *s2;  s2++);
  if ('+' == *s2) {
    sd = s2+1;
  }
  else if ('-' == *s2) {
    sign = -1;
    sd = s2+1;
  }
  else
    sd = s2;

  /* first check for the usual case - decimal digits */
  if (((*sd > '0') && (*sd <= '9')) ||
      (('0' == *sd) && ('\0' == sd[1] || isdigit(sd[1]))) ) {
    li = strtol((char *)s2, (char **)&end, 10);
    *i = li;
    if ('\0' == *end) {         /* reached the end, things went OK */
      *code = 0;
    }
    else
      *code = (SYSTEM_integer)(end - (SYSTEM_byte *)buffer + 1);
    return;
  }

  /* if not a decimal string,
   * must be either $ffff or 0xffff or an error */
  if ('$' == *sd) {
    if (isxdigit(sd[1])) {
      if (-1 == sign)
        *sd = '-';
      else
        sd++;
      li = strtol((char *)sd, (char **)&end, 16);
      *i = li;
      if ('\0' == *end) {               /* reached the end, things went OK */
        *code = 0;
      }
      else
        *code = (SYSTEM_integer)(end - (SYSTEM_byte *)buffer + 1);
    }
    else {
      *i = 0;
      sd++;
      *code = (SYSTEM_integer)(sd - (SYSTEM_byte *)buffer + 1);
    }
    return;
  }
  else if (('0' == *sd) &&
           (('x' == sd[1]) || ('X' == sd[1]))
           ) {
    li = strtol((char *)s2, (char **)&end, 16);
    *i = li;
    if ('\0' == *end) {         /* reached the end, things went OK */
      *code = 0;
    }
    else {
      /* we alread read the 0x, that is not an error for val */
      if (end < sd+2)
        end = sd+2;
      *code = (SYSTEM_integer)(end - (SYSTEM_byte *)buffer + 1);
    }
    return;
  }
  else {
    *i = 0;
    *code = (SYSTEM_integer)(sd - (SYSTEM_byte *)buffer + 1);
  }

  return;
} /* _P3_Val_i */

/* valid strings look like: [+|-][0x|$]d+,
 * where d is a decimal digit unless preceded by the 0x or $,
 * in which case it is a hex digit
 */
SYSTEM_integer _P3_Val_SPD(const SYSTEM_byte *s, SYSTEM_integer *code)
{
  char buffer[256];
  SYSTEM_byte *end, *s2, *sd;
  long int li;
  int len;
  int sign = 1;

  len = s[0];
  strncpy (buffer, (char *)(s+1), len);
  buffer[len] = '\0';

  /* skip over blanks
   * - Kylix 3 does not treat any other chars as whitespace
   */
  for (s2 = (SYSTEM_byte *)buffer;  ' ' == *s2;  s2++);
  if ('+' == *s2) {
    sd = s2+1;
  }
  else if ('-' == *s2) {
    sign = -1;
    sd = s2+1;
  }
  else
    sd = s2;

  /* first check for the usual case - decimal digits */
  if (((*sd > '0') && (*sd <= '9')) ||
      (('0' == *sd) && ('\0' == sd[1] || isdigit(sd[1]))) ) {
    li = strtol((char *)s2, (char **)&end, 10);
    if ('\0' == *end) {         /* reached the end, things went OK */
      *code = 0;
    }
    else
      *code = (SYSTEM_integer)(end - (SYSTEM_byte *)buffer + 1);
    return li;
  }

  /* if not a decimal string,
   * must be either $ffff or 0xffff or an error */
  if ('$' == *sd) {
    if (isxdigit(sd[1])) {
      if (-1 == sign)
        *sd = '-';
      else
        sd++;
      li = strtol((char *)sd, (char **)&end, 16);
      if ('\0' == *end) {               /* reached the end, things went OK */
        *code = 0;
      }
      else
        *code = (SYSTEM_integer)(end - (SYSTEM_byte *)buffer + 1);
    }
    else {
      li = 0;
      sd++;
      *code = (SYSTEM_integer)(sd - (SYSTEM_byte *)buffer + 1);
    }
    return li;
  }
  else if (('0' == *sd) &&
           (('x' == sd[1]) || ('X' == sd[1]))
           ) {
    li = strtol((char *)s2, (char **)&end, 16);
    if ('\0' == *end) {         /* reached the end, things went OK */
      *code = 0;
    }
    else {
      /* we already read the 0x, that is not an error for val */
      if (end < sd+2)
        end = sd+2;
      *code = (SYSTEM_integer)(end - (SYSTEM_byte *)buffer + 1);
    }
    return li;
  }
  else {
    li = 0;
    *code = (SYSTEM_integer)(sd - (SYSTEM_byte *)buffer + 1);
  }

  return li;
} /* _P3_Val_SPD */

void _P3_Val_d(const SYSTEM_byte *s, SYSTEM_double *d, SYSTEM_integer *code)
{
  SYSTEM_byte buffer[256];
  SYSTEM_byte *end;
  const SYSTEM_byte *stopper;
  const SYSTEM_byte *s0;        /* src ptr */
  SYSTEM_byte *s2, *sd;
  int len;
  int sign = 1;

  len = s[0];
#if 0
  strncpy ((char *)buffer, (char *)(s+1), len);
#else
  for (stopper = s+len, s0 = s+1, sd = buffer;
       s0 <= stopper;  s0++, sd++) {
    if (('D' == *s0) || ('d' == *s0))
      /* something that won't look the the exponent to MS strtod */
      *sd = 'Z';
    else
      *sd = *s0;
  }
#endif
  buffer[len] = '\0';

  /* skip over blanks
   * - Kylix 3 does not treat any other chars as whitespace
   */
  for (s2 = buffer;  ' ' == *s2;  s2++)
    ;
  if ('+' == *s2) {
    sd = s2+1;
  }
  else if ('-' == *s2) {
    sign = -1;
    sd = s2+1;
  }
  else
    sd = s2;

  /* guard against some special cases where strtod
   * doesn't do the right thing for val():
   *   the decimal string in front of the decimal exponent
   *      must be nonempty for strtod, not so for val
   *   hex input, starts with 0x or 0X
   *   nan or inf (case insensitive)
   */
  if (isdigit(*sd)) {
    if ('x' == tolower(sd[1])) {
      end = sd+1;
      *code = (SYSTEM_integer)(end - buffer + 1);
      *d = 0;
      return;
    }
    *d = strtod ((char *)s2, (char **)&end);
    if ('\0' == *end) {         /* reached the end, things went OK */
      *code = 0;
    }
    else {
      *code = (SYSTEM_integer)(end - buffer + 1);
    }
    return;
  } /* if digit after space and sign char */
  else if ('.' == *sd) {
    if ('\0' == sd[1]) {
      /* corner case of valid input not handled by strtod */
      *code = 0;
      *d = 0;
    }
    else {
      if ('e' == tolower(sd[1])) {
        *sd = '0';
      }
      *d = strtod ((char *)sd, (char **)&end);
      *d *= sign;
      if ('\0' == *end) {               /* reached the end, things went OK */
        *code = 0;
      }
      else {
        if (end <= sd)
          end = sd+1;
        *code = (SYSTEM_integer)(end - buffer + 1);
      }
    }
    return;
  }
  else {
    *d = 0;
    *code = (SYSTEM_integer)(sd - buffer + 1);
  }
  return;
} /* _P3_Val_d */


/****  Program initialization and command line arguments *****/

/* Save program parameters */
int _P3_argc;
char** _P3_argv;

/* Time and date procedures */
#include <sys/types.h>
#include <time.h>

#if defined(P3DOS)

#elif defined (P3UNIX) || defined(DJGPP)

#include <sys/param.h>
#include <sys/times.h>

#else

#endif

SYSTEM_longint _P3Filesize(int Iplus, SYSTEM_file *fil, const char *File, int Line)
{
  struct stat buf;
  int result;
  SYSTEM_longint res;

  if (_P3_ISOPEN(fil)) {        /* open file */
    result = fstat(fileno(fil->f), &buf);

    if (result == 0) {
      if (0 == fil->block_size)
        res = 1;
      else
        res = fil->block_size;
      res = buf.st_size / res;
    }
    else {
      res = -1;
      _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_FSIZE);
    }
  }
  else {                                /* not open */
    res = -1;
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_FSIZE);
  }

  if (Iplus) {
    _P3error_check(/*File, Line*/);
  }
  return res;
} /* _P3Filesize */


SYSTEM_longint _P3Filepos(int Iplus, SYSTEM_file *fil, const char *File, int Line)
{
  SYSTEM_longint res;

  if (! (_P3_ISOPEN(fil)) || fil->f == NULL || fil->block_size == 0) {
    res = -1;
    _P3_IO_NOTOPEN(_P3_err, EIO, fil, _P3_ERR_VERB_FPOS);
  }
  else {
    res = ftell(fil->f) / fil->block_size;
    if (-1 == res) {
      _P3_IO_ERR(_P3_err, errno, fil, _P3_ERR_VERB_FPOS);
    }
  }

  if (Iplus) {
    _P3error_check(/*File, Line*/);
  }
  return res;
} /* _P3Filepos */

/* Global variables: */

 /* NULL initializers changed in init, GCC objects otherwise */
SYSTEM_text SYSTEM_input     = {NULL, _P3OPEN | _P3RESET,   0, "\0"};
SYSTEM_text SYSTEM_output    = {NULL, _P3OPEN | _P3REWRITE, 0, "\0"};
SYSTEM_text SYSTEM_erroutput = {NULL, _P3OPEN | _P3REWRITE, 0, "\0"};

#if defined(__cplusplus)
const _P3set_elem _P3empty_set[ _P3set_max ] = {};
#else
const _P3set_elem _P3empty_set[ _P3set_max ] = {
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
#endif

/* Class/Interfaces stuff */

/* See p3io.h for declaration of SYSTEM_tobject and SYSTEM_classdescriptor_t */

static void * const SYSTEM_tobject_VT[] = {(void*)&SYSTEM_tobject_DOT_destroy};

const SYSTEM_classdescriptor_t SYSTEM_tobject_CD = {
  _P3str1("\007TObject"),
  NULL, NULL, 0,
  sizeof(SYSTEM_tobject),
  SYSTEM_tobject_VT,
  NULL
};

SYSTEM_boolean _P3is(const SYSTEM_tobject objref,
                     const SYSTEM_classdescriptor_t *clsref)
{
  if (!objref || !clsref)
    return SYSTEM_false;
  else {
    const struct SYSTEM_classdescriptor *cd;
    cd = objref->CD;

    while (cd && cd != clsref)
      cd = cd->ancestor;

    return (cd == clsref);
  }
} /* _P3is */

SYSTEM_tobject _P3as(SYSTEM_tobject objref,
                     const SYSTEM_classdescriptor_t *clsref,
                     const char *file, int line)
{
  char buf[1000];

  if (!objref)
    return NULL;
  else {
    if (clsref) {
      const struct SYSTEM_classdescriptor *cd;
      cd = objref->CD;

      while (cd && cd != clsref)
        cd = cd->ancestor;

      if (cd)
        return objref;
    }
  }

  /* Error return. Runtime error, or at least return nil */
  sprintf(buf, "Invalid class typecast at (%s:%d):"
          " '%s' not an instance of '%s'", file, line,
          (objref != NULL) && (objref->CD != NULL) ?
          (SYSTEM_byte*)((objref->CD->name)+1) :
          (SYSTEM_byte*)("nil"),
          (clsref != NULL) ?
          (SYSTEM_byte*)((clsref->name)+1) :
          (SYSTEM_byte*)("nil") );

  _P3_Exception(_P3_EXC_CODE_INVALIDCAST, buf);
  return NULL; /* Never reaches this */
} /* _P3as */

/* Some standard procs/funcs from Tobject: */
/*SYSTEM_tobject SYSTEM_tobject_DOT_destroy(SYSTEM_tobject obj)*/
Destructor(SYSTEM_tobject)  SYSTEM_tobject_DOT_destroy   (SYSTEM_tobject obj)
{
  EXCTST(printf("SYSTEM_tobject_DOT_destroy called\n");)
  return obj;
}

void _P3_abstract_call1(void *ref, ...)
{
  SYSTEM_tobject o = (SYSTEM_tobject) ref;
  char buf[512];

  sprintf(buf, "Call of abstract method from class '%s'",
               _P3get_CD(o)->name+1);

  _P3_Exception(_P3_EXC_CODE_ABSTRACTERROR, buf);
} /* _P3_abstract_call1 */

void
_P3_abstract_call3(SYSTEM_ansichar *result, SYSTEM_uint8 _len_ret,
                   SYSTEM_tobject o, ...)
{
  char buf[512];

  sprintf(buf, "Call of abstract method from class '%s'",
               _P3get_CD(o)->name+1);

  _P3_Exception(_P3_EXC_CODE_ABSTRACTERROR, buf);
} /* _P3_abstract_call3 */

void _P3_abstr_cl_call1(void *ref, ...)
{
  SYSTEM_classreference_t r = (SYSTEM_classreference_t) ref;
  char buf[512];

  sprintf(buf, "Call of abstract method from class '%s'", r->name+1);
  _P3_Exception(_P3_EXC_CODE_ABSTRACTERROR, buf);
} /* _P3_abstr_cl_call1 */

/* like _P3_abstr_cl_call,
 * but has two initial args for the shortstring return */
void
_P3_abstr_cl_call3(SYSTEM_ansichar *result, SYSTEM_uint8 _len_ret,
                   SYSTEM_classreference_t ref, ...)
{
  char buf[512];

  sprintf(buf, "Call of abstract method from class '%s'", ref->name+1);
  _P3_Exception(_P3_EXC_CODE_ABSTRACTERROR, buf);
} /* _P3_abstr_cl_call3 */

/* SYSTEM: TObject's free */
void SYSTEM_tobject_DOT_free(SYSTEM_tobject self)
{
  EXCTST(printf("P3PC TObject.free called\n"); fflush(stdout);)
  if (self != NULL) {
    _P3dealloc_object(VirtMethodCall((SYSTEM_tobject)self,
      SYSTEM_tobject_DOT_destroy_T, 0, ((SYSTEM_tobject)self)));
  }
}  /* free */

/* End class/Interfaces stuff */

void P3_PGM_init(char **_Argv, int _Argc, _P3void_procT PGM_Final)
{ int i;
  _P3_argc = _Argc-1;  _P3_argv = _Argv;

  /* Empty the empty set: */
  /* for (i=0; i<_P3set_max; i++) _P3empty_set[i] = (_P3set_elem)0; */

  /* default filemode: ReadWrite */
  SYSTEM_filemode = 2;

  /* Init std io - gcc doesn't like it in initializer above */
#if 0
  SYSTEM_input.f = stdin;
  SYSTEM_output.f = stdout;
  SYSTEM_erroutput.f = stderr;
#else
  SYSTEM_input.f = NULL;
  SYSTEM_output.f = NULL;
  SYSTEM_erroutput.f = NULL;
#endif

  _P3islibrary = 0; /* Note: Different in DLLs and application */

  /* Important: Set ExitProc to the finalization routine.
     This is the only way it'll be called!                     */
  SYSTEM_exitproc = PGM_Final;
#if (__linux__) || defined(__APPLE__)
  (void) atexit(SYSTEM_exitproc);
#else
  /* lets just see if any systems do NOT like atexit */
  (void) atexit(SYSTEM_exitproc);
#endif
} /* P3_PGM_init */

void P3_DLL_init(void)
{ int i;
  _P3_argc = -1;  _P3_argv = NULL; /* No arguments inside a DLL */

  /* Empty the empty set: */
  /* for (i=0; i<_P3set_max; i++) _P3empty_set[i] = (_P3set_elem)0; */

  /* default filemode: ReadWrite */
  SYSTEM_filemode = 2;

  /* Init std io - gcc doesn't like it in initializer above */
#if 0
  SYSTEM_input.f = stdin;
  SYSTEM_output.f = stdout;
  SYSTEM_erroutput.f = stderr;
#else
  SYSTEM_input.f = NULL;
  SYSTEM_output.f = NULL;
  SYSTEM_erroutput.f = NULL;
#endif

  _P3islibrary = 1; /* Note: Different in DLLs and application */

  /* Set ExitProc to NULL in DLLs */
  SYSTEM_exitproc = NULL;
} /* P3_DLL_init */


#if defined(_WIN32)
# include <process.h>
#else
# include <sys/wait.h>
#endif

#if ! defined(_WIN32)

/* Steve's routines from 23 Aug 2001 (3nd version)  */

/* followSymLink ():
 * replace input pathname (potentially a symlink) with actual file
 * returns:
 *   0      on success, or if input not a symlink
 *   errno  on failure
 */
static int
followSymLink (char absPathName[PATH_MAX])
{
  char newPathName[PATH_MAX];
  int rc;
  char *p;

  rc = readlink (absPathName, newPathName, PATH_MAX);
  if (-1 == rc) {
    if (EINVAL == errno) {
      return 0;                 /* not a symlink */
    }
    else {
      return errno;
    }
  }
  else {                        /* a symlink, successfully followed */
    newPathName[rc] = '\0';
    if ('/' == newPathName[0]) {
      strcpy (absPathName, newPathName);
    }
    else {
      p = strrchr(absPathName, '/');
      strcpy (p+1, newPathName);
    }
  }

  return 0;
} /* followSymLink */


/* getAbsPath ():
 * given a path to a file, returns the absolute path
 * it assumes the file has a slash in it, i.e. it's position
 * relative to the current dir is specified explicitly
 * return:
 *    absPath   on success
 *    NULL      on failure
 */
static char *
getAbsPath (char absPath[PATH_MAX], const char *fname)
{
  if ('/' == fname[0]) { /* easiest case: path name absolute already */
    strcpy (absPath, fname);
  }
  else {                        /* argument is relative to cwd */
    if (!getcwd(absPath, PATH_MAX)) {
      absPath[0]='\0';
      return NULL;
    }
    if ('.' == fname[0] && '/' == fname[1]) {
      strcat (absPath, fname+1);
    }
    else {
      strcat (absPath, "/");
      strcat (absPath, fname);
    }
  } /* relative path */

  return absPath;
} /* getAbsPath */

/* unixGetModuleFileName (prog, buf, bufLen)
 * get the full path of the executable being run
 * in the process, follow symlinks and make relative paths absolute
 * returns:
 *    0 on failure
 *    strlen of the returned buf otherwise
 */
static int
unixGetModuleFileName (const char *prog, char *buf, int bufLen)
{
  char pathName[PATH_MAX];      /* may be relative or not */
  char absPathName[PATH_MAX];
  char *path = NULL;
  char *savePtr;
  char *p, *dir;
  int len;
  int rc;
  struct stat statBuf;

  /* easy case: absolute path name */
  if ('/' == prog[0]) {
    strcpy (pathName, prog);
    rc = followSymLink (pathName);
    if (rc) return 0;
    if ('/' == *pathName) {     /* path already absolute */
      p = pathName;
    }
    else {
      p = getAbsPath (absPathName, pathName);
      if (NULL == p) return 0;
    }
    len = strlen (p);
    if (len+1 > bufLen) {
      len = bufLen - 1;
    }
    strncpy (buf, p, len+1);
    return len;
  }

  /* relative path name */
  if (strchr(prog, '/')) {
    if (!getcwd (pathName, PATH_MAX)) {
      buf[0]='\0';
      return 0;
    }
    if ('.' == prog[0] && '/' == prog[1]) {
      strcat (pathName, prog+1);
    }
    else {
      strcat (pathName, "/");
      strcat (pathName, prog);
    }
    rc = followSymLink (pathName);
    if (rc) return 0;
    if ('/' == *pathName) {     /* path already absolute */
      p = pathName;
    }
    else {
      p = getAbsPath (absPathName, pathName);
      if (NULL == p) return 0;
    }
    len = strlen (p);
    if (len+1 > bufLen) {
      len = bufLen - 1;
    }
    strncpy (buf, p, len+1);
    return len;
  } /* relative path */

  /* Third Step: Scan the path */
  p = getenv ("PATH");
  if (p) {
    len = strlen (p);
    path = (char *) malloc (len+1);
    if (NULL == path) {
      buf[0]='\0';
      errno = ENOMEM;
      return 0;
    }
    strcpy (path, p);
    errno = 0;
    for (dir = strtok_r (path, ":", &savePtr);
         dir;
         dir = strtok_r (NULL, ":", &savePtr)) {
      if ('.' == dir[0] && '\0' == dir[1]) {
        if (! getcwd (pathName, PATH_MAX)) {
          buf[0]='\0';
          goto wefail;
        }
      }
      else {
        strcpy (pathName, dir);
      }
      if (pathName[ strlen(pathName)-1 ] != '/') /* SSN added; this helps */
        strcat (pathName, "/");
      strcat (pathName, prog);
      if (!access (pathName, X_OK)
          && (0 == stat(pathName, &statBuf))
          && ! S_ISDIR(statBuf.st_mode)
          ) {
        rc = followSymLink (pathName);
        if (rc) goto wefail;
        if ('/' == *pathName) { /* path already absolute */
          p = pathName;
        }
        else {
          p = getAbsPath (absPathName, pathName);
          if (NULL == p) goto wefail;
        }
        len = strlen (p);
        if (len+1 > bufLen) {
          len = bufLen - 1;
        }
        strncpy (buf, p, len+1);
        free (path);
        return len;
      } /* if found in path */
    } /* end loop over dirs in path */
    free (path);
    errno = ENOENT;
    return 0;
  } /* if (path is set) */

 wefail:
  if (NULL != path) {
    free (path);
  }
  errno = ENOENT;
  return 0;
} /* unixGetModuleFileName */

#define Test_GMFN 0

static SYSTEM_char *
wrapUnixGMFN (SYSTEM_char *res, SYSTEM_byte max)
{
  char appDir[PATH_MAX+1];
  int rc, len, i;

  if (max > PATH_MAX)
    max = (SYSTEM_byte)PATH_MAX; /* max is max 255 */

  if (Test_GMFN)
    printf ("argv[0] is '%s'\n", _P3_argv[0]);
  rc = unixGetModuleFileName (_P3_argv[0], appDir, PATH_MAX);
  if (0 == rc) {
    if (Test_GMFN)
      printf ("wrapUnixGMFN(): Unable to get module file dir of '%s'\n",
              _P3_argv[0]);
  }
  else if (rc < max) {
    if (Test_GMFN)
      printf ("Module file dir is '%s'\n", appDir);
  }
  else {
    appDir[PATH_MAX-1] = '\0';
    if (Test_GMFN) {
      printf ("*** Warning: Module file dir was truncated:"
              " %d chars not enough!\n", PATH_MAX);
      printf ("Truncated module file dir is '%s'\n", appDir);
    }
    rc = 0;  /* just kill it, can't use it anyway */
  }

  /* Copy string (of length rc, sic) back to Pascal */
  len = (rc > max)? max : rc;
  for (i=1; i <= len; i++)
    res[i] = appDir[i-1];
  res[0] = (SYSTEM_char)len;

  return res;
} /* wrapUnixGMFN */

#endif  /* if ! WIN32 */

#ifdef _P3_LIBRARY
#define DLL_TEST_STUBS 0

/***************************************************************
        DLL / Shared Objects Stubs etc.
***************************************************************/


Extern_C int _P3_DllInit();
extern int  _P3_DllInit();  /* Redundant if Extern_C defined, otherwise necessary */
Extern_C void _P3_DllFini();
extern void _P3_DllFini();  /* Redundant if Extern_C defined, otherwise necessary */

#if defined(_WIN32)
BOOL WINAPI
DllMain(HINSTANCE hInst, DWORD reason, LPVOID lpReserved)
{
  switch (reason) {
  case DLL_PROCESS_ATTACH:
    if (DLL_TEST_STUBS)
      printf("SSNDLL: DllMain DLL_PROCESS_ATTACH, calling _init() now\n");
    (void) _P3_DllInit();
    break;
  case DLL_PROCESS_DETACH:
    if (DLL_TEST_STUBS)
      printf("SSNDLL: DllMain DLL_PROCESS_DETACH, calling _fini() now\n");
    _P3_DllFini();
    break;
  case DLL_THREAD_ATTACH:
  case DLL_THREAD_DETACH:
    /* ignored */
    if (DLL_TEST_STUBS)
      printf("SSNDLL: DllMain called, reason: %d\n", reason);
    break;
  } /* switch */

  return TRUE;
} /* DllMain */

#elif defined(DAR) || defined(DEG) || defined(DEI) || defined(DIG) || defined(DII) || defined(LEG) || defined(LNX)
Extern_C __attribute__((constructor))
void init(void)
{                               /* Called automatically at SO load */
  /* printf("--------------- DAR/DEG/DEI/DIG/DII init called ----------- \n"); */
  (void) _P3_DllInit();
}

Extern_C __attribute__((destructor))
void fini(void)
{                               /* Called automatically at SO unload */
  /* printf("----------------- DAR/DEG/DEI/DIG/DII fini  called ----------- \n"); */
  _P3_DllFini();
}

#elif  defined(AIX) || defined(BGP) || defined(LEI)
Extern_C void _init(void)
{                               /* Called automatically at SO load */
  /*printf("------------------------  _init  called ----------- \n");*/
  (void) _P3_DllInit();
}

Extern_C void _fini(void)
{                              /* Called automatically at SO unload */
  /*printf("------------------------  _fini  called ----------- \n");*/
  _P3_DllFini();
}

#elif defined(SIG) || defined(SOL) || defined(SOX)

Extern_C void P3_sol_initFunc (void)
{
  /*printf("------------------------ P3_sol_initFunc called ----------- \n");*/
  (void) _P3_DllInit();
}

Extern_C void P3_sol_finiFunc (void)
{
  /*printf("------------------------ P3_sol_finiFunc called ----------- \n");*/
  _P3_DllFini();
}

#else
#error "No shared library init/fini defined for this build"
  !! should die with an error right here;
#endif /* if (WIN32) .. elif .. else error endif */

#endif  /*  _P3_LIBRARY  */

/*************** NEW EXCEPTIONS - Object-based like Delphi/Kylix.         */
/***************                  Declarations for SYSTEM_exception class */

void * const SYSTEM_exception_VT[] = {(void*)&SYSTEM_tobject_DOT_destroy};

/* Descriptor, 'exception' */
const SYSTEM_classdescriptor_t SYSTEM_exception_CD = {
  _P3str1("\011exception"),
  &SYSTEM_tobject_CD, NULL, 0,
  sizeof(SYSTEM_exception_OD), SYSTEM_exception_VT, NULL
};

/*************** Constructor code: constructor create(const Mst: string): */

  Constructor(SYSTEM_exception ) SYSTEM_exception_DOT_create(
    SYSTEM_exception self,
    const SYSTEM_char *msg)
  {
    _P3strcpy(self->SYSTEM_exception_DOT_message,255,msg);
    return self;
  }  /* create */

/********   How to declare an exception object in C:
  static SYSTEM_exception exc;

  --- and how to create one:
  exc = ValueCast(SYSTEM_exception,SYSTEM_exception_DOT_create(ValueCast(
    SYSTEM_exception,_P3alloc_object(&SYSTEM_exception_CD)),
    _P3str1("\027A Test Exception Object")));            ********************/


/*************** END NEW EXCEPTIONS ***************************************/



/*************** PARAMSTRING STUFF ***************************************/

#if defined(_WIN32)
  /* Steve's routines for handling paramcount/paramstr in Windows. */
typedef SYSTEM_char *P3_pchar;

static void
pushChar (SYSTEM_char c, P3_pchar *r, SYSTEM_integer *len)
{
  if (*len < 255) {
    **r = c;
    ++*len;
    ++*r;
  }
}  /* pushChar */

static P3_pchar
getParamShortStr (P3_pchar p, SYSTEM_char *param)
{
  SYSTEM_integer len;
  P3_pchar s;
  P3_pchar r;

  for ( ; ; ) {
    while (*p != '\0' && *p <= ' ')
      ++p;
    s = p;
    ++s;
    if (*p == '\"' && *s == '\"')
      p += 2;
    else
      break;
  }

  len = 0;
  r = param+1;
  while (*p > ' ') {
    if (*p == _P3char('\"')) {
      p++;
      while (*p != '\0' && *p != '\"') {
        pushChar(*p,&r,&len);
        ++p;
      }
      if (*p != '\0')
        ++p;
    }
    else {
      pushChar(*p,&r,&len);
      ++p;
    }
  }
  _P3setlength(param,len,255);
  return p;
} /* getParamShortStr */

#else

#include "p3Custom2.h"

static SYSTEM_char *
UNIX_paramstr (SYSTEM_char *res, int k)
{
  /* Return parameter string k in the res string: assume sizeof(res)=255 */
  /* k == 0 is ok, returns name of program. */
  /* Note: paramstr(0) must return fully qualified name on all platforms */
  const int resMax = 255;

  if (k < 0 || k > _P3_argc)
    *res = 0;
  else

  if (0 == k) {
    int rc;
    SYSTEM_shortstring msg;

#if defined(__HOS_AIX__)
    rc = 999;
#else
    rc = xGetExecName (res, msg);
#endif
    switch (rc) {
    case 0:                     /* success */
      break;
    case 1:                     /* overflows short string */
      *res = '\0';              /* return empty string */
      break;
    default:                    /* other failure */
      (void) wrapUnixGMFN (res, resMax);
    }
  }
  else {
    int i;
    char *p = _P3_argv[k];
    i = strlen(p);
    if (i > resMax)
      i = resMax;
    memmove(res+1, p, i);
    *res = i;
  }
  return res;
} /* UNIX_paramstr */
#endif  /* if defined(_WIN32) .. else .. */

#if defined(_WIN32)
SYSTEM_integer
_P3_paramcount (void)
{
  SYSTEM_integer result;
  char *cmdLine;
  SYSTEM_char *p;
  SYSTEM_shortstring s;

  result = 0;
  cmdLine = GetCommandLine(); /* we work in the system memory: don't write */
  p = getParamShortStr((SYSTEM_char *)cmdLine, s);
  for ( ; ; ) {
    p = getParamShortStr(p, s);
    if (0 == SYSTEM_length(s))
      break;
    result++;
  }
  return result;
}  /* paramcount */
#endif  /* if defined(_WIN32) */

SYSTEM_char
*SYSTEM_paramstr (SYSTEM_char *result, SYSTEM_uint8 _len_ret,
                  SYSTEM_integer index)
{
  SYSTEM_shortstring s;

#if defined(_WIN32)

  SYSTEM_char *p;
  char buffer[261];   /* length take from Delphi's System.pas */
  DWORD fnLen;

  if (0 == index) {
    fnLen = GetModuleFileName (NULL, buffer, sizeof(buffer));
    (void) _P3pa2str(result, _len_ret, (SYSTEM_char *)buffer, fnLen);
    return result;
  }

  _P3strcpy(result,_len_ret,_P3str1("\000"));
  p = (SYSTEM_char *) GetCommandLine(); /* we work in the system memory: do not write */
  for ( ; ; ) {
    p = getParamShortStr(p, s);
    if (0 == index || 0 == SYSTEM_length(s))
      break;
    index--;
  }
  _P3strcpy(result,_len_ret, s);
#else
  _P3strcpy(result,_len_ret, UNIX_paramstr(s,index));
#endif /* if defined(_WIN32) */
  return result;
}  /* SYSTEM_paramstr */

/*************** END PARAMSTRING STUFF ***********************************/
