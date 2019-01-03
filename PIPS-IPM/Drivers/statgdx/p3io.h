
/*******************************************************/
/*         P3    RUNTIME SUPPORT ROUTINES              */
/*               (c) 1991   Soren Nielsen              */
/*         .h :  Include file.                         */
/*******************************************************/


#ifndef _P3IO_H
#define _P3IO_H

/************************** history ***********************

For recent history see the CHANGELOG.txt file

03/03/18  Remove any ifdef on EXCEPTION_MODEL, always use setjmp/longjmp.
03/03/14  Add declarations for SYSTEM_exception and all its stuff.
03/03/14  Add on-clause-macros (_P3_ON_CLAUSE, _P3_END_ON_CLAUSE,
          _P3_ON_ELSE, _P3_END_ON_ELSE).
02/09/06  Change _P3setlength back to procedure instead of inline.
02/08/30  Change P3_PGM_EXIT to do "return SYSTEM_exitcode".
02/08/23  Add macro for fillchar -> memset (from system.pas)
02/08/23  Alternative def of PointerCast for LNX to allow optimization
02/06/25  Bug correct: Remember to unwind jmp_buf list if no exception
02/06/22  _P3_Exception, cleanup, signal handlers, etc. etc.
02/06/15: Implement exception models _P3_EXCEPTION_MODEL...
02/04/16: SYSTEM_delete: Change byte parameters to integers.
01/11/27: Implement int64s on VIS as __64; impl. 64-bit I/O on all.
01/08/27: Implement enums as bytes; BYTEENUM macro
01/08/16: New integer types: (u)int8/16/32.
          Also: __alpha is back. 'integer' = 'longint' are 32 bits again.
01/08/14: dispose/freemem: Allow freeing non-lvalue, don't pass **.
01/07/13: Added runtime check of equal sizes in Variable Casts.
01/07/12: Added void to empty formal parm lists of proc/funcs.
00/10/26: Added definitions for _P3_DllExport
01/02/01: Added support for DJGPP: gpp for DOS.
00/12/17: Added SYSTEM_reallocmem.
00/11/04: Added backward compatibility for all standard proc/funcs
          declared in system.pas (e.g., define SYSTEM_ioresult _P3ioresult)
          The names in system.pas all changed from _P3... to SYSTEM_...
00/10/18: Removed special code for __alpha. USING 64-bit INTS NOW.
00/09/27: Various changes in integer output formats (%d to %ld or %lu):
          Add FMT_(routine) defines dependent upon integer size -> p3io.h
00/09/16: Add Prototype #define
00/09/09: Add #defines for calling virtual methods, VirtMethodCall etc.
00/08/22: Change c_gen to treat MININT specially (output SYSTEM_minint)
               Defined below as Steve suggested (always has been...)
00/05/18: Code for writeln(x:a:b) changed back to %f from %g.
00/05/16: Bug in readln (read from output...)
00/03/27: Printing of unformattet double chg from 17:10 to %23.14e.
00/03/27: Inserted #ifdef __alpha to make longs 32 bits, not 64.
               Works with both the Alpha's cc and gcc/g++.
99/12/15: Experiment using ZEROES to initialize anything. Works.
99/11/16: Correct _P3_insert (EK bug).
99/11/15: Merge EK's str and val code
99/09/26: Define standard string and set types when _P3_STRING_DEFS
99/09/26: Make sqr function in-line (_P3sqr_i and _P3sqr_r)
99/09/17: Implement _P3SET_ic for constant i in "i in set"
99/08/13: Change all occurrences of memcpy to memmove.
99/08/10: Add Filepos function.
99/08/05: Change SYSTEM_boolean to signed char for Watcom bug
              It compiled (!SYSTEM_false) to false otherwise!
99/07/30: Make _P3str2pa macro; compiler checks lengths OK
99/07/25: Set OS type to P3DOS or P3UNIX if not specified,
              instead of issuing error message.
99/06/27: Add SYSTEM_break, SYSTEM_continue.
          Add block_size parameter to untyped reset/rewrite.
          Add seek function.
99/06/26: Introduce IO_CHK flag; update I/O to set _P3ioresult
99/06/25: Clean-up _P3ioresult declarations
99/06/12: Merge Erwin's changes from compiling as C++
BRA
99/06/10: hi/lo functions
99/03/07: Added _P3_PGM_EXIT() macro -> exit(0)

SSN 11/14/98. Change to use _P3_Close instead of fclose
SSN 11/13/98. Code in _P3Uappend corrected. Also bug in _P3setcmpN
SSN 10/24/98. str2pa made macro; re-introduced in p3io
SSN 07/12/98. Bug fixes (stdout -> &_P3stdout).
SSN 09/09/97. Major rewrite to new file representation
                 (as structures) etc. Also, most type
                 names changed from _P3Txxx to SYSTEM_xxx
                 for consistency (they're declared in the
                 Pascal System module).
  nov 96 EK: added untyped files and blockread/write
  nov 96 EK: added seek
  dec 96 EK: incorporated changes by SD
  dec 96 EK: moved GAMS specific stuff into
             gamsext.h

****************************** end history ***************************/

#if ! defined(__cplusplus)
#  error "Detected a C compiler but we now require C++ compilation"
#endif

#if defined(_AIX)
   /* this includes _XOPEN_SOURCE, _POSIX_SOURCE, and other stuff too */
#  include <standards.h>
#endif
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <assert.h>
#include <errno.h>
#include <math.h>

#if defined(_WIN32)
#  define WIN32_LEAN_AND_MEAN   /* google it */
#  include <windows.h>
#  include <winbase.h>
#  include <direct.h>
#else
#  include <unistd.h>
#  include <pthread.h>
#  include <semaphore.h>
#endif

/* _P3_EXC_MODEL_ controls how P3 implements Delphi exceptions
 *  0 : setjmp/longjmp scheme, works identically with C and C++
 *  1 : use C++ try/catch
 * Also try a verbosity flag _P3_EXC_VERBOSE_ with values:
 *  0 : no debug output
 *  1 : some debug output
 *  2 : more debug output
 */
#define _P3_EXC_MODEL_  1
#if ! defined(_P3_EXC_VERBOSE_)
# define _P3_EXC_VERBOSE_ 0
#endif /* ! defined(_P3_EXC_VERBOSE_) */

#if defined(__cplusplus)
# include <iostream>
# include <cstddef>
# include <exception>
#else
/* if no C++, EXC_MODEL=1 is broken.  But zero is always broken! */
# error "NO valid EXC_MODEL without C++"
#endif

#undef P3OSDEFINED

/* Make sure DJGPP excludes P3DOS and P3UNIX */
#if defined(DJGPP)
#undef P3DOS
#undef P3UNIX
#define P3OSDEFINED
#endif

#if defined(P3DOS)
#if defined(P3OSDEFINED)
#error "P3<OSTYPE> redefined"
p3ostype_redefined__ERROR;
#endif
#define P3OSDEFINED
#endif

#if defined(P3UNIX)
#if defined(P3OSDEFINED)
#error "P3<OSTYPE> redefined"
p3ostype_redefined__ERROR;
#endif
#define P3OSDEFINED
#endif

#if defined(P3CMS) || defined(P3TSO) || defined(P3VMS) || defined(P3VMS)
#error "P3OS<TYPE> defined for unimplemented OS type"
P3OS_TYPE_defined_for_unimplemented_OS_type__ERROR;
#endif

#if ! defined(P3OSDEFINED)
/* Set to P3DOS or P3UNIX by default. _WIN32 is MS Visual C/C++ specific? */
#if defined(_WIN32)
#define P3OSDEFINED
#define P3DOS
#else
#define P3OSDEFINED
#define P3UNIX
/*   #error "P3<OSTYPE> not defined for any recognized OS type"  */
/*   P3OS_TYPE_not_defined_for_recognized_OS_type__ERROR;        */
#endif
#endif



/* Calling conventions: */
/* Adjusted by SPD: 26 Sep 02 to be set only if previously unset */
#if ! defined(REGISTER)
# define REGISTER   /* Default; leave empty here */
#endif
#if ! defined(CDECL)
# if defined(P3DOS)
#  define CDECL   __cdecl
# elif defined(P3UNIX) || defined(DJGPP)
#  define CDECL
# endif
#endif
#if ! defined(STDCALL)
# if defined(P3DOS)
#  define STDCALL   __stdcall
# elif defined(P3UNIX) || defined(DJGPP)
#  define STDCALL
# endif
#endif

#if defined(__cplusplus)
#define P3_INLINE inline
#else
#define P3_INLINE /* nothing */
#endif

#if ! defined(Extern_C)
# if defined(__cplusplus)
#   define Extern_C extern "C"
# else
#   define Extern_C      /* Nada */
# endif
#endif

#if ! defined(C_LINKAGE)
# if defined(__cplusplus)
#   define C_LINKAGE(FUNCDECL) extern "C" { FUNCDECL }
# else
#   define C_LINKAGE(FUNCDECL) FUNCDECL
# endif
#endif

#if ! defined(_P3_DllExport)
# if defined(_P3_LIBRARY) && defined(P3DOS) /* We're compiling a DOS DLL */
#   define _P3_DllExport __declspec( dllexport )
# else
#   define _P3_DllExport
# endif
#endif

#ifdef COULD_DO_THIS
#define int long  /* so as to NEVER use an int (16/32?)  */
#pragma warning(1:158)
#pragma warning(1:177)
#pragma warning(1:178)
#endif


/* This removes an irritating warning message from MS Visual C/C++: */
/* What's a better flag than _WIN32?                                */
#if defined(_WIN32)
#pragma warning(2:4761)
#endif

/* P3 TYPES: Definitions of Pascal types in C
 * Define the following integer types so that they hold EXACTLY:
 *  SYSTEM_shortint       -128 .. 127
 *  SYSTEM_smallint     -32768 .. 32767
 *  SYSTEM_longint -2147483648 .. 2147483647
 *  SYSTEM_byte              0 .. 255
 *  SYSTEM_word              0 .. 65535
 *  SYSTEM_longword          0 .. 4294967295
 *
 *  The following are generic types which may change in Delphi/Kylix,
 *  but are currently have have been for a long time
 *  SYSTEM_integer      SYSTEM_longint
 *  SYSTEM_cardinal     SYSTEM_longword

 *  Other system-defined scalar types:
 *  SYSTEM_boolean          0 .. 1
 *  SYSTEM_char (ascii)     chr(0) .. chr(255)
 */


/* As far as possible, directly declare Pascal standard types in C */
/* On boolean: If declared as an enum, has size 4 on RS/6000:
 * typedef enum{SYSTEM_false, SYSTEM_true} SYSTEM_boolean;  -- instead next 3:
 */
/* SSN: Add signed for Watcom: */
typedef  signed char SYSTEM_boolean;
#define SYSTEM_false ((SYSTEM_boolean)0)
#define SYSTEM_true  ((SYSTEM_boolean)1)
/* Can't do the following two in C - aren't considered to be constants...
  const SYSTEM_boolean SYSTEM_false = 0;
  const SYSTEM_boolean SYSTEM_true  = 1;
 * hence the casted constants above
 */

/* Declare the fundamental types: SYSTEM_[u]int{8|16|32|64},
 * along with FMT strings to go with them
 * apart from SYSTEM_[u]int64, these call all be the same,
 * regardless of the platform
 */
typedef  signed        char    SYSTEM_int8;
typedef  signed   short int    SYSTEM_int16;
typedef  signed         int    SYSTEM_int32;
typedef  unsigned      char    SYSTEM_uint8;
typedef  unsigned short int    SYSTEM_uint16;
typedef  unsigned       int    SYSTEM_uint32;

#define FMT_write_cx "%%%dc"
#define FMT_write_u0 "%u"
#define FMT_write_i0 "%d"
/*      FMT_write_y0        platform-dependent */
/*      FMT_write_z0        platform-dependent */
#define FMT_write_ux "%%%du"
#define FMT_write_ix "%%%dd"
/*      FMT_write_yx        platform-dependent */
/*      FMT_write_zx        platform-dependent */
#define FMT_read_u   "%u"
#define FMT_read_i   "%d"
/*      FMT_read_y          platform-dependent */
/*      FMT_read_z          platform-dependent */
#define FMT_Str_i0   "%d"
#define FMT_Str_i1   "%%%dd"
#define FMT_Str_d1   "%%%d.%dE"
#define FMT_Str_d2   "%%%d.%df"
#define FMT_rangechk "value=%d, range=[%d..%d]"

#if defined(P3DOS)
typedef          __int64 SYSTEM_int64;
typedef unsigned __int64 SYSTEM_uint64;
#  define FMT_write_y0 "%I64u"
#  define FMT_write_z0 "%I64d"
#  define FMT_write_yx "%%%dI64u"
#  define FMT_write_zx "%%%dI64d"
#  define FMT_read_y   "%I64u"
#  define FMT_read_z   "%I64d"

#elif defined(LEG) || defined(LEI)
typedef   signed long int SYSTEM_int64;
typedef unsigned long int SYSTEM_uint64;
#  define FMT_write_y0 "%lu"
#  define FMT_write_z0 "%ld"
#  define FMT_write_yx "%%%dlu"
#  define FMT_write_zx "%%%dld"
#  define FMT_read_y   "%lu"
#  define FMT_read_z   "%ld"

#else
typedef   signed long long int SYSTEM_int64;
typedef unsigned long long int SYSTEM_uint64;
#  define FMT_write_y0 "%llu"
#  define FMT_write_z0 "%lld"
#  define FMT_write_yx "%%%dllu"
#  define FMT_write_zx "%%%dlld"
#  define FMT_read_y   "%llu"
#  define FMT_read_z   "%lld"
#endif

#if defined(_WIN32)
# if defined(_WIN64)
typedef unsigned __int64  SYSTEM_nativeuint;
#  define FMT_write_n0 "%I64u"
#  define FMT_write_nx "%%%dI64u"
#  define FMT_read_n   "%I64u"
# else
typedef unsigned __int32  SYSTEM_nativeuint;
#  define FMT_write_n0 "%I32u"
#  define FMT_write_nx "%%%dI32u"
#  define FMT_read_n   "%I32u"
# endif
#else
typedef unsigned long int SYSTEM_nativeuint;
# define FMT_write_n0 "%lu"
# define FMT_write_nx "%%%dlu"
# define FMT_read_n   "%lu"
#endif

/* The following implement the integer type declarations in system.pas,
 * based on the "fundamental" types: SYSTEM_[u_]int{8|16|32}
 */
typedef SYSTEM_int32    SYSTEM_longint;
typedef SYSTEM_int16    SYSTEM_smallint;
typedef SYSTEM_int8     SYSTEM_shortint;
typedef SYSTEM_uint32   SYSTEM_longword;
typedef SYSTEM_uint16   SYSTEM_word;
typedef SYSTEM_uint8    SYSTEM_byte;

typedef SYSTEM_int32    SYSTEM_integer;
typedef SYSTEM_uint32   SYSTEM_cardinal;

typedef unsigned char   SYSTEM_ansichar; /* Keep unsigned so ord works!! */
typedef float           SYSTEM_single;
typedef double          SYSTEM_double;
typedef SYSTEM_double   SYSTEM_real;       /* Double/real are synonyms */

/* in P3, char and ansichar must always be synonyms */
#define SYSTEM_char SYSTEM_ansichar

typedef struct _P3file { /* Use this for all files */
  FILE            *f;
  SYSTEM_shortint status;
  SYSTEM_longint  block_size;
  SYSTEM_byte     nam[257];
} _P3file;
/* Status: 0 closed/error, 1 app, 2 read, 3 wrt */
typedef  _P3file        SYSTEM_text;     /* Pascal text file */
typedef  SYSTEM_text    SYSTEM_textfile; /* Pascal text file */
typedef  _P3file        SYSTEM_typed_file; /* typed  file */
typedef  _P3file        SYSTEM_file; /* Which one actually used? */
typedef  _P3file        SYSTEM_untypedfile; /* Untyped file */
typedef  _P3file        *_P3file_ptr;

typedef void            *SYSTEM_pointer; /* Can point to anything */
typedef void            SYSTEM_untyped;  /* Untyped parameters    */

typedef SYSTEM_char     SYSTEM_shortstring[256]; /* Used in system.pas */

#define ENUM1(x,i)  /* (SYSTEM_byte)   */(i)  /* just output the ordinal val */
#define ENUM2(x,i)  /* (SYSTEM_uint8)  */(i)  /* just output the ordinal val */
#define ENUM4(x,i)  /* (SYSTEM_uint16) */(i)  /* just output the ordinal val */
#define ENUM1a(x,i) /* (SYSTEM_byte)   */(x)  /* just use the enum literal */
#define ENUM2a(x,i) /* (SYSTEM_uint8)  */(x)  /* just use the enum literal */
#define ENUM4a(x,i) /* (SYSTEM_uint8)  */(x)  /* just use the enum literal */

/* SETS: Sizes and dimensions:
 *       The typedef for _P3set_elem sets everything. It should
 *       be SYSTEM_char for portability. Anything larger works on most
 *       machines but not SUNs (little-endians).
 */

#define _P3bits_per_byte 8 /* Don't run on machines with smaller bytes */
#define _P3setsize(i) ((i)/_P3bits_per_elem+1) /* bytes to hold 0..i bits */
#define _P3bits_per_elem  (_P3bits_per_byte*sizeof(_P3set_elem))
#define _P3set_max   (256/_P3bits_per_elem) /* Length of set in elem's */

typedef SYSTEM_char     _P3set_elem;
typedef _P3set_elem     *_P3Tset_ptr;
typedef _P3set_elem     _P3set255[_P3setsize(255)];

/* SETS: End of definitions */


/* To make proc/func declarations look a little nicer (?): */
#define Procedure void
#define Function(t) t
#define Constructor(t) t
#define Destructor(t)  t
#define Prototype typedef
#define FORWARD_REF_PTR(typ)  SYSTEM_pointer  /* Later replace by void*  */

#define cnstdef enum
#define CNST    /* empty. Could define as 'const' then change p3io.c/h */

#define SYSTEM_break(L)    goto L
#define SYSTEM_continue(L) goto L

/* codes for exceptions generated within P3 */
#define _P3_EXC_CODE_NONE         0
#define _P3_EXC_CODE_INVALIDCAST  1
#define _P3_EXC_CODE_INOUTERROR   2
#define _P3_EXC_CODE_RANGEERROR   3
#define _P3_EXC_CODE_ASSERTIONFAILED 4
#define _P3_EXC_CODE_ADDRANGE     5
#define _P3_EXC_CODE_OUTOFMEMORY  6
#define _P3_EXC_CODE_ABSTRACTERROR 7

/* TO DO: StrToInt64 in sysutils_p3.pas: exception on error? */

#define _P3_ERR_VERB_NONE     0
#define _P3_ERR_VERB_READ     1
#define _P3_ERR_VERB_WRITE    2
#define _P3_ERR_VERB_FLUSH    3
#define _P3_ERR_VERB_SEEK     4
#define _P3_ERR_VERB_EOFCHK   5
#define _P3_ERR_VERB_SEEKEOF  6
#define _P3_ERR_VERB_EOLNCHK  7
#define _P3_ERR_VERB_SEEKEOLN 8
#define _P3_ERR_VERB_FPOS     9
#define _P3_ERR_VERB_FSIZE    10
#define _P3_ERR_VERB_CLOSE    11
#define _P3_ERR_VERB_APPEND   12
#define _P3_ERR_VERB_REWRITE  13
#define _P3_ERR_VERB_RESET    14
#define _P3_ERR_VERB_ERASE    15
#define _P3_ERR_VERB_RMDIR    16
#define _P3_ERR_VERB_MKDIR    17
#define _P3_ERR_VERB_CHDIR    18
#define _P3_ERR_VERB_UNSPEC   19

typedef struct _P3err {
  int n;                        /* old _P3_errno */
  unsigned char verb;
  unsigned char notOpened;      /* handle this special case */
  unsigned char nam[257];
} _P3err_t;

extern _P3err_t _P3_err;
extern void _P3error_check(/*char*, int*/);

#define _Iplus_bgn()  /* I+ state, before I/O */
#define _Iplus_end()  /*           and after  */  _P3error_check(/*__FILE__, __LINE__*/);
#define _Iminus_bgn() /* I- state, before I/O */  if (!(_P3_err.n)) {
#define _Iminus_end() /*           and after  */  }
#define IOWRP(io) if ((io) < 0) _P3_err.n = errno;


#define _P3str1(s)  ((SYSTEM_char*)(s))
#define _P3str2(s)     s
#define _P3char(c)  ((SYSTEM_char)(c)) /* Otherwise '\377' is negative */


/*********  SIGNAL HANDLING  ***************************/

#include <signal.h>
#include <float.h>

#define _P3_ABORT() { signal( SIGABRT, SIG_DFL); abort(); }

/*********  END SIGNAL HANDLING  ***********************/

/* CAST Operators:

   ValueCast:      Used for scalar/pointer Rvalues
   PointerCast:    Used for pointer Lvalues (and for casts from "untyped")
   VariableCast:   Used for non-scalar/pointer casts; runtime size-check.  */


#define ValueCast(t, e)     (t)(e)
#define PointerCast(t, p)  *(t*)(p)  /* t pointers, p address */

/* Define VariableCast so (1) performs runtime size check; (2) is a lvalue */
/* Note that the runtime size check is a compile-time constant check; free */
/* Note: v below already has its address taken (or is an array name)       */
int _P3VariableCastError (const char *File, int Line,
                          size_t sizeof_a, size_t sizeof_b);

#define VariableCast(t, v, vtyp) \
  *((sizeof(t) != sizeof(vtyp)?  \
      _P3VariableCastError(__FILE__,__LINE__,sizeof(t),sizeof(vtyp)):0), \
      ((t*)(v)))  /* OK, so if it looks funny don't look at it. */

/* END CAST Operators */


/* ANSI prototypes of all routines defined in P3io.c  */

extern const SYSTEM_byte _P3true[];
extern const SYSTEM_byte _P3false[];
extern _P3Tset_ptr _P3set_p(SYSTEM_longint len, _P3Tset_ptr ret,
                            const _P3set_elem *s1, const _P3set_elem *s2);
extern _P3Tset_ptr _P3set_m(SYSTEM_longint len, _P3Tset_ptr ret,
                            const _P3set_elem *s1, const _P3set_elem *s2);
extern _P3Tset_ptr _P3set_t(SYSTEM_longint len, _P3Tset_ptr ret,
                            const _P3set_elem *s1, const _P3set_elem *s2);
extern _P3Tset_ptr _P3set_expand(SYSTEM_longint toLen, _P3Tset_ptr to,
                                 SYSTEM_longint frLen, const _P3set_elem *fr);
extern _P3Tset_ptr _P3set_copy(SYSTEM_longint len, _P3Tset_ptr to,
                               const _P3set_elem *fr);
#define _P3SET_p(d,l,s1,s2) _P3set_p(_P3setsize(l),d,s1,s2)
#define _P3SET_m(d,l,s1,s2) _P3set_m(_P3setsize(l),d,s1,s2)
#define _P3SET_t(d,l,s1,s2) _P3set_t(_P3setsize(l),d,s1,s2)
#define _P3SET_expand(d,L1,s,L2) _P3set_expand(_P3setsize(L1),d,_P3setsize(L2),s)
#define _P3SET_copy(s1,l,s2) _P3set_copy(_P3setsize(l),s1,s2)
#define _P3set1(s) ((_P3set_elem*)(s))  /* Cast literal string to set pointer */
#define _P3set2(s) s /*remove{}*/             /* Literal set in declaration */

extern SYSTEM_boolean _P3set_i(SYSTEM_longint len,
                               SYSTEM_longint i, const _P3set_elem *s);
extern _P3Tset_ptr _P3set_add_elem(SYSTEM_longint len, _P3Tset_ptr s,
                                   _P3set_elem i);
extern _P3Tset_ptr _P3set_add_range(SYSTEM_longint mx, _P3Tset_ptr s,
                                   _P3set_elem lo, _P3set_elem up);

extern double SYSTEM_int(SYSTEM_double x);
extern double SYSTEM_frac(SYSTEM_double x);
extern SYSTEM_int64 SYSTEM_round(SYSTEM_double x);

#define SYSTEM_fillchar(p, sz, val) memset((void*)(p),(int)(val),(size_t)(sz))
#define SYSTEM_abs_i(x) labs((SYSTEM_int64)(x))
#define SYSTEM_abs_r(x) fabs((double)(x))
extern SYSTEM_int64   SYSTEM_sqr_i(SYSTEM_int64 i);
extern SYSTEM_double  SYSTEM_sqr_r(SYSTEM_double x);


/* What's _P3hi supposed to do on a big-endian? */
#define SYSTEM_lo(x)     ((SYSTEM_byte)(x))
#define SYSTEM_hi(x)     ((SYSTEM_byte)((x)>>_P3bits_per_byte))

extern void _P3setlength(SYSTEM_byte *s, SYSTEM_integer newLen,
                         SYSTEM_integer siz);

extern SYSTEM_char SYSTEM_upcase(SYSTEM_char);
extern SYSTEM_char SYSTEM_locase(SYSTEM_char);
extern void _P3_new(void **p, SYSTEM_longint s);
extern void _P3_new64(void **p, SYSTEM_int64 s);
extern void _P3_free(void *p, SYSTEM_longint s);
extern void _P3_free64(void *p, SYSTEM_int64 s);
extern void SYSTEM_realloc(void** p, SYSTEM_longint lgt);

extern void SYSTEM_reallocmem(void **p, SYSTEM_longint s);
extern void SYSTEM_reallocmem64(void **p, SYSTEM_int64 s);
extern SYSTEM_char *_P3_pcharn2str (SYSTEM_char *dst, SYSTEM_char dstSiz,
                                    const SYSTEM_char *src, int srcLen);
extern SYSTEM_char *_P3_pchar2str (SYSTEM_char *dst, SYSTEM_char dstSiz,
                                   const SYSTEM_char *src);
/* only use this macro with a constant string src!! */
#define _P3conp2str(dst, dstSiz, src) _P3_pcharn2str(dst, dstSiz, (SYSTEM_char *)(src), sizeof(src)-1)
extern SYSTEM_char *_P3pa2str(SYSTEM_char *s, SYSTEM_char m,
                              const SYSTEM_char *p, SYSTEM_char n);
extern SYSTEM_char *_P3_ch2str(SYSTEM_char *st, SYSTEM_byte max,
                               SYSTEM_char ch);
/*extern SYSTEM_char *_P3str2pa( SYSTEM_char *p,  SYSTEM_longint m,
                               SYSTEM_char *s); */
#define _P3str2pa(p,m,s)  (SYSTEM_char*)_P3memcpy(p, m, s+1)
extern SYSTEM_char *_P3_strcat(SYSTEM_char *r,  SYSTEM_char max,
                               const SYSTEM_char *p1, const SYSTEM_char *p2);
extern SYSTEM_boolean _P3streq(const SYSTEM_char *p1, const SYSTEM_char *p2);
extern SYSTEM_boolean _P3streq_ic(const SYSTEM_char *p1, const SYSTEM_char *p2);
extern int _P3strcmp(const SYSTEM_char *p1, const SYSTEM_char *p2);
extern int _P3stpcmp(const SYSTEM_char *st, const SYSTEM_char *pa, int lgt);
extern int _P3stccmp(const SYSTEM_char *st, SYSTEM_char ch);
extern int _P3_argc;
extern char** _P3_argv;
extern SYSTEM_integer SYSTEM_allocmemcount, SYSTEM_allocmemsize;
extern SYSTEM_int64                         SYSTEM_allocmemsize64;
extern SYSTEM_char* SYSTEM_paramstr(SYSTEM_char *res,
                                    SYSTEM_byte max, int index);
extern _P3file SYSTEM_input, SYSTEM_output, SYSTEM_erroutput;
extern SYSTEM_char *SYSTEM_copy(SYSTEM_byte *res, SYSTEM_byte max,
                                const SYSTEM_byte *s,
                                SYSTEM_integer i, SYSTEM_integer cnt);
extern SYSTEM_integer SYSTEM_pos(const SYSTEM_byte *sub, const SYSTEM_byte *s);
extern void SYSTEM_delete(SYSTEM_byte *s, SYSTEM_integer i, SYSTEM_integer cnt);
extern void _P3_insert(const SYSTEM_byte *sub, SYSTEM_byte *dest,
                       SYSTEM_byte maxDest, SYSTEM_integer i);
extern SYSTEM_byte* _P3_strcpy(SYSTEM_byte *d, SYSTEM_integer max,
                               const SYSTEM_byte *s);
extern long _P3rangechk(SYSTEM_longint, SYSTEM_longint, SYSTEM_longint);


#define _P3strcat(r,max,b,c) _P3_strcat(r,max,b,c)
#define _P3ch2str(a,b,c) _P3_ch2str(a,b,c)
#define _P3insert(sub,dest,maxDest,i)  _P3_insert(sub,dest,maxDest,i)
/* #define _P3delete(s,i,c) SYSTEM_delete(s,i,c) */

#define _P3shl(a,b)  ((a) << (b))
#define _P3shr(a,b)  ((a) >> (b))
#define _P3inc0(x)       ++(x)
#define _P3inc1(x,i)   (x)+=(i)
#define _P3dec0(x)       --(x)
#define _P3dec1(x,i)   (x)-=(i)


extern int SYSTEM_ioresult(void); /* Test and reset _P3_err.n */

extern void _P3write_u (_P3file *fil, SYSTEM_longword i);
extern void _P3write_i (_P3file *fil, SYSTEM_longint  i);
extern void _P3write_n (_P3file *fil, SYSTEM_nativeuint i);
extern void _P3write_y (_P3file *fil, SYSTEM_uint64   i);
extern void _P3write_z (_P3file *fil, SYSTEM_int64    i);
extern void _P3write_ux(_P3file *fil, SYSTEM_longword i, SYSTEM_longint len);
extern void _P3write_ix(_P3file *fil, SYSTEM_longint  i, SYSTEM_longint len);
extern void _P3write_nx(_P3file *fil, SYSTEM_nativeuint i, SYSTEM_longint len);
extern void _P3write_yx(_P3file *fil, SYSTEM_uint64   i, SYSTEM_longint len);
extern void _P3write_zx(_P3file *fil, SYSTEM_int64    i, SYSTEM_longint len);
extern void _P3write_r (_P3file *fil, SYSTEM_double d);
extern void _P3write_rx(_P3file *fil, SYSTEM_double d, SYSTEM_longint len1);
extern void _P3write_ry(_P3file *fil, SYSTEM_double d,
                        SYSTEM_longint len1, SYSTEM_longint len2);
extern void _P3write_c (_P3file *fil, SYSTEM_char c);
extern void _P3write_cx(_P3file *fil, SYSTEM_char c, SYSTEM_longint len);
extern void _P3_writeln (void);
extern void _P3_writefn (SYSTEM_text *fil);
extern void _P3_write_s0(const SYSTEM_byte *s);
extern void _P3_writefs0(_P3file *fil, const SYSTEM_byte *s);
extern void _P3_Readfs0(SYSTEM_text *fil, SYSTEM_byte *s, SYSTEM_byte max);
extern void _P3write_sx(SYSTEM_text *fil,
                        const SYSTEM_byte *s, SYSTEM_longint len);

extern SYSTEM_longword _P3read_u(_P3file *fil);
extern SYSTEM_longint  _P3read_i(_P3file *fil);
extern SYSTEM_nativeuint _P3read_n(_P3file *fil);
extern SYSTEM_uint64   _P3read_y(_P3file *fil);
extern SYSTEM_int64    _P3read_z(_P3file *fil);
extern SYSTEM_char     _P3read_c(_P3file *fil);
extern SYSTEM_double   _P3read_d(_P3file *fil);

extern SYSTEM_boolean _P3_eof(int Iplus, _P3file *fil, const char *File, int Line);
extern SYSTEM_boolean _P3_eoln(int Iplus, _P3file *fil, const char *File, int Line);
extern SYSTEM_boolean _P3_seekeof(int Iplus, _P3file *fil, const char *File, int Line);
extern SYSTEM_boolean _P3_seekeoln(int Iplus, _P3file *fil, const char *File, int Line);
extern void _P3block_read_write(_P3file *fil, void *buf,
                                size_t count, SYSTEM_longint *numread,
                                SYSTEM_boolean wr);
extern void _P3rw_typed(_P3file *fil, void *buf, int wr);
extern void _P3_Assign(_P3file *f, const SYSTEM_char *s);
extern void _P3fileopn(_P3file *f, SYSTEM_longint s, SYSTEM_longint t,
                       SYSTEM_longint block_size);
extern void _P3_Close(_P3file *f);
extern void _P3_Flush(_P3file *f);

extern void SYSTEM_rmdir(const SYSTEM_char *s);
extern void SYSTEM_mkdir(const SYSTEM_char *s);
extern void SYSTEM_chdir(const SYSTEM_char *s);

/* INTEGER I/O: WRITE */

/* unsigned integer, 32 bit */
#define _P3write_u0(i)      _P3write_u (&SYSTEM_output,i)
#define _P3writefu0(i)      _P3write_u (_file_temp,i)
#define _P3write_u1(i, len) _P3write_ux(&SYSTEM_output,i,len)
#define _P3writefu1(i, len) _P3write_ux(_file_temp,i,len)
#define _P3write_u2(i,m,n)  _P3write_r2((double)(i),m,n)
#define _P3writefu2(i,m,n)  _P3writefr2((double)(i),m,n)

/* signed integer, 32 bit */
#define _P3write_i0(i)      _P3write_i (&SYSTEM_output,i)
#define _P3writefi0(i)      _P3write_i (_file_temp,i)
#define _P3write_i1(i, len) _P3write_ix(&SYSTEM_output,(SYSTEM_longint)(i),len)
#define _P3writefi1(i, len) _P3write_ix(_file_temp,(SYSTEM_longint)(i),len)
#define _P3write_i2(i,m,n)  _P3write_r2((double)(i),m,n)
#define _P3writefi2(i,m,n)  _P3writefr2((double)(i),m,n)

/* unsigned integer, native */
#define _P3write_n0(i)      _P3write_n (&SYSTEM_output,i)
#define _P3writefn0(i)      _P3write_n (_file_temp,i)
#define _P3write_n1(i,len)  _P3write_nx(&SYSTEM_output,i,len)
#define _P3writefn1(i,len)  _P3write_nx(_file_temp,i,len)
#define _P3write_n2(i,m,n)  _P3write_r2((double)(i),m,n)
#define _P3writefn2(i,m,n)  _P3writefr2((double)(i),m,n)

/* unsigned integer, 64 bit */
#define _P3write_y0(i)      _P3write_y (&SYSTEM_output,i)
#define _P3writefy0(i)      _P3write_y (_file_temp,i)
#define _P3write_y1(i,len)  _P3write_yx(&SYSTEM_output,i,len)
#define _P3writefy1(i,len)  _P3write_yx(_file_temp,i,len)
#define _P3write_y2(i,m,n)  _P3write_r2((double)(i),m,n)
#define _P3writefy2(i,m,n)  _P3writefr2((double)(i),m,n)

/* signed integer, 64 bit */
#define _P3write_z0(i)      _P3write_z (&SYSTEM_output,i)
#define _P3writefz0(i)      _P3write_z (_file_temp,i)
#define _P3write_z1(i,len)  _P3write_zx(&SYSTEM_output,i,len)
#define _P3writefz1(i,len)  _P3write_zx(_file_temp,i,len)
#define _P3write_z2(i,m,n)  _P3write_r2((double)(i),m,n)
#define _P3writefz2(i,m,n)  _P3writefr2((double)(i),m,n)

/* INTEGER I/O: READ */

#define _P3read__u0(i)    i = _P3read_u(&SYSTEM_input)
#define _P3read_fu0(i)    i = _P3read_u(_file_temp)

#define _P3read__i0(i)    i = _P3read_i(&SYSTEM_input)
#define _P3read_fi0(i)    i = _P3read_i(_file_temp)

#define _P3read__n0(i)    i = _P3read_n(&SYSTEM_input)
#define _P3read_fn0(i)    i = _P3read_n(_file_temp)

#define _P3read__y0(i)    i = _P3read_y(&SYSTEM_input)
#define _P3read_fy0(i)    i = _P3read_y(_file_temp)

#define _P3read__z0(i)    i = _P3read_z(&SYSTEM_input)
#define _P3read_fz0(i)    i = _P3read_z(_file_temp)

/* REALS  (double) */

/* write reals with zero, one, or two modifiers given */
#define _P3write_r0(d)     _P3write_r (&SYSTEM_output,(double)(d))
#define _P3writefr0(d)     _P3write_r (_file_temp,(double)d)
#define _P3write_d0(d)     _P3write_r (&SYSTEM_output,d)
#define _P3writefd0(d)     _P3write_r (_file_temp,d)
#define _P3write_r1(d,m)   _P3write_rx(&SYSTEM_output,(double)(d),m)
#define _P3writefr1(d,m)   _P3write_rx(_file_temp,(double)(d),m)
#define _P3write_d1(d,m)   _P3write_rx(&SYSTEM_output,(double)(d),m)
#define _P3writefd1(d,m)   _P3write_rx(_file_temp,(double)(d),m)
#define _P3write_r2(d,m,n) _P3write_ry(&SYSTEM_output,(double)d,m,n)
#define _P3writefr2(d,m,n) _P3write_ry(_file_temp,(double)d,m,n)
#define _P3write_d2(d,m,n) _P3write_ry(&SYSTEM_output,d,m,n)
#define _P3writefd2(d,m,n) _P3write_ry(_file_temp,d,m,n)

/* Read real */
#define _P3read__r0(r)    r = _P3read_d(&SYSTEM_input)
#define _P3read_fr0(r)    r = _P3read_d(_file_temp)
#define _P3read__d0(d)    d = _P3read_d(&SYSTEM_input)
#define _P3read_fd0(d)    d = _P3read_d(_file_temp)

/* CHAR */

#define _P3write_c0(ch)    _P3write_c (&SYSTEM_output,ch)
#define _P3writefc0(ch)    _P3write_c (_file_temp,ch)

#define _P3write_c1(ch,m)  _P3write_cx(&SYSTEM_output,ch,m)
#define _P3writefc1(ch,m)  _P3write_cx(_file_temp,ch,m)

/* Define print char with two modifiers: Ignore second modifier */
#define _P3write_c2(c,m,n) _P3write_cx(&SYSTEM_output,c,m)
#define _P3writefc2(c,m,n) _P3write_cx(_file_temp,c,m)

#define _P3read__c0(ch)    ch = _P3read_c(&SYSTEM_input)
#define _P3read_fc0(ch)    ch = _P3read_c(_file_temp)


/* STRINGS */

typedef SYSTEM_char _P3STR_3[4];
typedef SYSTEM_char _P3STR_7[8];
typedef SYSTEM_char _P3STR_15[16];
typedef SYSTEM_char _P3STR_31[32];
typedef SYSTEM_char _P3STR_63[64];
typedef SYSTEM_char _P3STR_95[96];
typedef SYSTEM_char _P3STR_127[128];
typedef SYSTEM_char _P3STR_159[160];
typedef SYSTEM_char _P3STR_191[192];
typedef SYSTEM_char _P3STR_223[224];
typedef SYSTEM_char _P3STR_255[256];
typedef SYSTEM_char _P3SET_3[_P3setsize(3)];
typedef SYSTEM_char _P3SET_7[_P3setsize(7)];
typedef SYSTEM_char _P3SET_15[_P3setsize(15)];
typedef SYSTEM_char _P3SET_31[_P3setsize(31)];
typedef SYSTEM_char _P3SET_63[_P3setsize(63)];
typedef SYSTEM_char _P3SET_95[_P3setsize(95)];
typedef SYSTEM_char _P3SET_127[_P3setsize(127)];
typedef SYSTEM_char _P3SET_159[_P3setsize(159)];
typedef SYSTEM_char _P3SET_191[_P3setsize(191)];
typedef SYSTEM_char _P3SET_223[_P3setsize(223)];
typedef SYSTEM_char _P3SET_255[_P3setsize(255)];

#define _P3read__s0(s,max)   _P3_Readfs0(&SYSTEM_input, s,max)
#define _P3read_fs0(s,max)   _P3_Readfs0(_file_temp, s,max)

#define _P3write_s0(st)  _P3_write_s0(st)
#define _P3writefs0(st)  _P3_writefs0(_file_temp,st)
#define _P3write_s1(st,m)  _P3write_sx(&SYSTEM_output, (st), m)
#define _P3writefs1(st,m)  _P3write_sx(_file_temp, (st), m)
#define _P3_copy(res,max,s,i,l)    _P3copy(res,max,s,i,l)

/* Booleans - inline test */
#define _P3write_b0(b)   _P3write_s0((b)? _P3true : _P3false)
#define _P3write_b1(b,m) _P3write_s1((b)? _P3true : _P3false, m)
#define _P3writefb0(b)   _P3writefs0((b)? _P3true : _P3false)
#define _P3writefb1(b,m) _P3writefs1((b)? _P3true : _P3false, m)


/* READLN / WRITELN  */

#if 0
#define _P3writeln()  IOWRP(printf("\n"))
#define _P3writefn()  IOWRP(fprintf(_file_temp->f,"\n"))
#else
#define _P3writeln()  _P3_writeln()
#define _P3writefn()  _P3_writefn(_file_temp)
#endif

extern void _P3read_ln(_P3file *fil);
#define _P3readlf()  _P3read_ln(_file_temp)
#define _P3readln()  _P3read_ln(&SYSTEM_input)

/* eoln, eof */

#define _P3eof_(I)  _P3_eof(I, &SYSTEM_input,__FILE__,__LINE__)
#define _P3eoff(I, f) _P3_eof(I, &(f),__FILE__,__LINE__)

#define _P3eoln(I)  _P3_eoln(I, &SYSTEM_input,__FILE__,__LINE__)
#define _P3eolf(I, f) _P3_eoln(I, &(f),__FILE__,__LINE__)

#define _P3seekeoln(I)   _P3_seekeoln(I,&SYSTEM_input,__FILE__,__LINE__)
#define _P3seefeoln(I,f) _P3_seekeoln(I,&(f),__FILE__,__LINE__)

#define _P3seekeof(I)    _P3_seekeof(I,&SYSTEM_input,__FILE__,__LINE__)
#define _P3seefeof(I,f)  _P3_seekeof(I,&(f),__FILE__,__LINE__)

/*  Read and Write, only for typed files. */
#define _P3fread(buf)   _P3rw_typed(_file_temp, buf, 0)
#define _P3frite(buf)   _P3rw_typed(_file_temp, buf, 1)


/* Blockread/write, only for untyped files: */
#define _P3blockR3(f,b,c)     _P3block_read_write(&(f),b,c,NULL,0);
#define _P3blockW3(f,b,c)     _P3block_read_write(&(f),b,c,NULL,1);
#define _P3blockR4(f,b,c,d)   _P3block_read_write(&(f),b,c,&(d),0);
#define _P3blockW4(f,b,c,d)   _P3block_read_write(&(f),b,c,&(d),1);

extern SYSTEM_longint _P3Filesize(int Iplus, _P3file *f, const char *File, int Line);
extern SYSTEM_longint _P3Filepos (int Iplus, _P3file *f, const char *File, int Line);
#define _P3filesize(I,f)  _P3Filesize(I,&(f),__FILE__,__LINE__)
#define _P3filepos(I,f)   _P3Filepos(I,&(f),__FILE__,__LINE__)


/* RESET, REWRITE, APPEND */

/* In name, and third parameter to fileopn:
(P3T.., 0): text file;  (P3F..., 1): typed file; (P3U..., 2): untyped file */
/* Also used for status field in P3file struct */

#define _P3_STATEMASK    3
/* we plan to add initialization code to set files unassigned */
#define _P3UNASSIGNED    0
#define _P3CLOSED        1
#define _P3OPEN          2
#define _P3_MODEMASK    12
#define _P3APPEND        0
#define _P3RESET         4
#define _P3REWRITE       8
#define _P3UPDATE       12

/* #define _P3_ISOPEN(f) (((f)->status & _P3_STATEMASK) == _P3OPEN) */
#define _P3_ISOPEN(f)      ((f)->status & _P3OPEN)
#define _P3_ISCLOSED(f)    ((f)->status & _P3CLOSED)
#define _P3_NOTASSIGNED(f) (((f)->status & _P3_STATEMASK) == _P3UNASSIGNED)
/* _P3_SETMODE has a dual function: set STATE to open and MODE = mode */
#define _P3_SETMODE(f,mode) (f)->status = _P3OPEN | (mode & _P3_MODEMASK)

#define _P3Treset(f)                _P3fileopn(&(f), _P3RESET  , 0, 1)
#define _P3Trewrite(f)              _P3fileopn(&(f), _P3REWRITE, 0, 1)
#define _P3Tappend(f)               _P3fileopn(&(f), _P3APPEND,  0, 1)
#define _P3Tupdate(f)               _P3fileopn(&(f), _P3UPDATE,  0, 1)
#define _P3Tclose(fil)              _P3_Close(&(fil))
#define _P3Tflush(fil)              _P3_Flush(&(fil))

#define _P3Freset(f,   elem_size)   _P3fileopn(&(f), _P3RESET  , 1, elem_size)
#define _P3Frewrite(f, elem_size)   _P3fileopn(&(f), _P3REWRITE, 1, elem_size)
#define _P3Fappend(f,  elem_size)   _P3fileopn(&(f), _P3APPEND , 1, elem_size)
#define _P3Fupdate(f,  elem_size)   _P3fileopn(&(f), _P3UPDATE,  1, elem_size)
#define _P3Fclose(fil)              _P3_Close(&(fil))

#define _P3Ureset(f,   block_size)  _P3fileopn(&(f), _P3RESET  , 2, block_size)
#define _P3Urewrite(f, block_size)  _P3fileopn(&(f), _P3REWRITE, 2, block_size)
#define _P3Uappend(f,  block_size)  _P3fileopn(&(f), _P3APPEND , 2, block_size)
#define _P3Uupdate(f,  block_size)  _P3fileopn(&(f), _P3UPDATE,  2, block_size)
#define _P3Uclose(fil)              _P3_Close(&(fil))


extern void _P3_Erase(SYSTEM_text *f);
#define _P3erase(f)    _P3_Erase(&(f))

extern void _P3_Seek(_P3file *fil, SYSTEM_longint pos, SYSTEM_longint from);
#define _P3Fseek(f,offset,origin) _P3_Seek(&(f),offset,(SYSTEM_longint)origin);

/* Communication between assign and reset/rewrite/append */

#define _P3PATHLEN 256 /* Max path length */

#define _P3assign(f,s)  _P3_Assign(&(f), s)

extern void _P3assert(const SYSTEM_char* mess, const char *File, int Line);
#define SYSTEM_assert(b,w)    {if (!(b)) _P3assert(w, __FILE__, __LINE__);}
#define SYSTEM_noassert(b,w)  /* ignore assertion */

/****************** END OF I/O DECLARATIONS AND MACROS *************/

/* SET OPERATIONS. Note: Bit-access operations (in-test,
   add element, add range) view sets as arrays of bytes,
   to ensure consistency with constant sets as strings */
#define _P3SET_i(mx,i,s)  _P3set_i(mx,i,s)
#define _P3SET_ic(mx,i,s)  /* mx and i are constant */    \
  (((i)<0||(i)>mx) ? 0 :                                  \
   (s)[(i)/_P3bits_per_elem]>>((i)%_P3bits_per_elem) & 1)
#define _P3SET_add_elem(mx,s,i) _P3set_add_elem(mx,s,i)
#define _P3SET_add_range(mx,s,lo,up) _P3set_add_range(mx,s,lo,up)

/*  MATH ETC. */

/* Keep the casts even though including math.h (Aug 99, DEC Alpha problem) */

#define SYSTEM_sin(x)         sin((double)(x))
#define SYSTEM_cos(x)         cos((double)(x))
#define SYSTEM_arctan(x)      atan((double)(x))
#define SYSTEM_sqrt(x)        sqrt((double)(x))
#define SYSTEM_exp(x)         exp((double)(x))
#define SYSTEM_ln(x)          log((double)(x))

#define SYSTEM_chr(x)         ((SYSTEM_byte)(x))
extern void _P3_halt(int);
extern int  _P3_Finalizing;
#define _P3halt1(x)       _P3_halt(x)
#define _P3halt0()        _P3_halt(0)
#define SYSTEM_length(x)      ((SYSTEM_byte)((x)[0]))  /* Cast necessary! */
/* extern void _P3setlength(SYSTEM_byte*, SYSTEM_integer, SYSTEM_integer);*/
#define SYSTEM_odd(x)         ((x)&1)

#if defined(_WIN32)
extern SYSTEM_integer _P3_paramcount(void);
#define SYSTEM_paramcount() _P3_paramcount()  /* Steve's routine in p3io.c */
#else
#define SYSTEM_paramcount()   _P3_argc
#endif

#define SYSTEM_trunc(x)    (SYSTEM_int64)(x)
#define SYSTEM_ord(x)          (x)
#define SYSTEM_pi()           3.1415926535897932385
#define SYSTEM_succ(i)         ((i)+1) /* Do not widen to 64 bits, ... */
#define SYSTEM_pred(i)         ((i)-1) /* ... neither does Kylix.      */

/* HEAP */
#define _P3new(x)         _P3_new((void**)&(x), sizeof(*(x)))
#define _P3getmem(x,i)    _P3_new((void**)&(x), i)
#define _P3freemem(x)     _P3_free((void*)(x), sizeof(*(x)))
#define _P3freemem0(x)    _P3_free((void*)(x), 0) /* When don't know size */
#define _P3freemem2(x, y) _P3_free((void*)(x), y)
#define _P3alloc_object(clsref) _P3_alloc_object(clsref)
#define _P3dealloc_object(p) _P3_dealloc_object(p)

/* STRCPY and MEMCPY things */
#define _P3strcpy(d,max,s)   _P3_strcpy(d,max,s)
#define _P3strclr(s1)        *(s1) = '\0'
#define _P3memcpy(s1,lgt,s2) memmove(s1,s2,lgt)  /* s1: dest, s2: src, lgt...*/
#define  SYSTEM_move(s,d,c)  memmove(d, s, c  )  /* dest, source, count      */


/* EK: STR variants */
void _P3_Str_i0(SYSTEM_integer i, SYSTEM_byte *s, SYSTEM_byte sMax);
void _P3_Str_i1(SYSTEM_integer i, SYSTEM_integer width, SYSTEM_byte *s,
                SYSTEM_byte sMax);
void _P3_Str_d0(SYSTEM_double x, SYSTEM_byte *s, SYSTEM_byte sMax);
void _P3_Str_d1(SYSTEM_double x, SYSTEM_integer width, SYSTEM_byte *s,
                SYSTEM_byte sMax);
void _P3_Str_d2(SYSTEM_double x, SYSTEM_integer width, SYSTEM_integer
                decimals, SYSTEM_byte *s, SYSTEM_byte sMax);
#define _P3str_i0(i,s,L)        _P3_Str_i0(i,s,L)
#define _P3str_i1(i,w,s,L)      _P3_Str_i1(i,w,s,L)
#define _P3str_i2(i,w,d,s,L)    _P3_Str_i1(i,w,s,L)
#define _P3str_d0(x,s,L)        _P3_Str_d0(x,s,L)
#define _P3str_d1(x,w,s,L)      _P3_Str_d1(x,w,s,L)
#define _P3str_d2(x,w,d,s,L)    _P3_Str_d2(x,w,d,s,L)


/* EK: VAL variants */
/* incomplete!!! Result can also be byte, single precision floating number
etc */
/* better to make a function ???? */
SYSTEM_integer _P3_Val_SPD(const SYSTEM_byte *s, SYSTEM_integer *code);
void _P3_Val_i(const SYSTEM_byte *s, SYSTEM_integer *i, SYSTEM_integer *code);
void _P3_Val_d(const SYSTEM_byte *s, SYSTEM_double *d, SYSTEM_integer *code);
#if 0
#define _P3val_i(s,i,c) _P3_Val_i(s,i,c)
#define _P3val_d(s,d,c) _P3_Val_d(s,d,c)
#else
#define _P3val_i(s,i,c) i=_P3_Val_SPD(s,c)
#define _P3val_d(s,d,c) _P3_Val_d(s,&(d),c)
#endif

/* Block compares. Used for arrays, sets and records.   */
/* Careful: Must return 0 or 1, not just 0 vs non-zero. */

#define _P3blkcmpE(p1,p2,l)   (memcmp(p1,p2,l)==0)
#define _P3blkcmpN(p1,p2,l)   (memcmp(p1,p2,l)!=0)
#define _P3blkcmpG(p1,p2,l)   (memcmp(p1,p2,l)>0)
#define _P3blkcmpL(p1,p2,l)   (memcmp(p1,p2,l)<0)
#define _P3blkcmpGE(p1,p2,l)  (memcmp(p1,p2,l)>=0)
#define _P3blkcmpLE(p1,p2,l)  (memcmp(p1,p2,l)<=0)

/* Next 4 used to inline smart 'in' tests: */
#define _P3SET_equal(i,j)     ((i) == (j))
#define _P3SET_in_1(i,j,k)    (((i) == (j)) || (k))
#define _P3SET_in_2(i,j,k,l)  (((i) >= (j) && (i) <= (k)) || (l))
#define _P3SET_in_3(i,j,k)    ((i) >= (j) && (i) <= (k))

#define _P3setcmpE(l,p1,p2)    (memcmp(p1,p2,l/_P3bits_per_byte+1)==0)
#define _P3setcmpN(l,p1,p2)    (memcmp(p1,p2,l/_P3bits_per_byte+1)!=0)

#define _P3strcmpE(p1,p2)   (_P3streq(p1,p2))
#define _P3strcmpN(p1,p2)   (!_P3streq(p1,p2))
#define _P3strcmpG(p1,p2)   (_P3strcmp(p1,p2)>0)
#define _P3strcmpL(p1,p2)   (_P3strcmp(p1,p2)<0)
#define _P3strcmpGE(p1,p2)  (_P3strcmp(p1,p2)>=0)
#define _P3strcmpLE(p1,p2)  (_P3strcmp(p1,p2)<=0)

#define _P3stpcmpE(p1,p2,l)   (_P3stpcmp(p1,p2,l)==0)
#define _P3stpcmpN(p1,p2,l)   (_P3stpcmp(p1,p2,l)!=0)
#define _P3stpcmpG(p1,p2,l)   (_P3stpcmp(p1,p2,l)>0)
#define _P3stpcmpL(p1,p2,l)   (_P3stpcmp(p1,p2,l)<0)
#define _P3stpcmpGE(p1,p2,l)  (_P3stpcmp(p1,p2,l)>=0)
#define _P3stpcmpLE(p1,p2,l)  (_P3stpcmp(p1,p2,l)<=0)

#define _P3stccmpE(p1,p2)   (_P3stccmp(p1,p2)==0)
#define _P3stccmpN(p1,p2)   (_P3stccmp(p1,p2)!=0)
#define _P3stccmpG(p1,p2)   (_P3stccmp(p1,p2)>0)
#define _P3stccmpL(p1,p2)   (_P3stccmp(p1,p2)<0)
#define _P3stccmpGE(p1,p2)  (_P3stccmp(p1,p2)>=0)
#define _P3stccmpLE(p1,p2)  (_P3stccmp(p1,p2)<=0)

extern SYSTEM_byte SYSTEM_filemode;

#if defined(__cplusplus)
extern "C" {
  typedef void (* _P3void_procT)(void);
  void _P3_Finalization(void);
}
#else
  typedef void (* _P3void_procT)(void);
  void _P3_Finalization(void);
#endif
extern SYSTEM_integer SYSTEM_exitcode;  /* system.ExitCode */
extern _P3void_procT  SYSTEM_exitproc;  /* system.ExitProc */

/* DLL stuff: */
extern SYSTEM_boolean _P3islibrary;     /* system.IsLibrary (Note: a function)*/
#define SYSTEM_islibrary()( _P3islibrary )
extern SYSTEM_integer SYSTEM_dll_refcount;

/* the Assigned procedure to test if (procedural) pointers are nil: */
#define _P3assigned(p)  ((p) != NULL)

/****  Program initialization and command line arguments *****/

/* User program initialization stuff: Save program parameters */

extern const _P3set_elem _P3empty_set[];

extern void P3_PGM_init(char **, int, _P3void_procT);
#define _P3_PGM_INIT()   P3_PGM_init(_Argv, _Argc, &_P3_Finalization);
/* #define _P3_PGM_EXIT()   {_P3_halt(0); ** Never gets here: ** return 0; }*/
#define _P3_PGM_EXIT()   { return SYSTEM_exitcode; }

extern void P3_DLL_init(void);
#define _P3_DLL_INIT()   P3_DLL_init(/*not used:  &_fini */);
#define _P3_DLL_UNWIND() {                     \
  while (SYSTEM_exitproc)    {                 \
    (*SYSTEM_exitproc)();                      \
  }                                            \
}

#define SYSTEM_maxint      2147483647  /* 2**31-1 */
#define SYSTEM_maxlongint  2147483647  /* 2**31-1 */
#define SYSTEM_maxlongword 4294967295U /* 2**32-1 */
#define SYSTEM_minint (-SYSTEM_maxint-1)
#if defined(P3DOS)
#  define SYSTEM_maxint64   9223372036854775807
#  define SYSTEM_minint64  -9223372036854775808 /* Visual doesn't like LL */
#  define SYSTEM_maxuint64  18446744073709551615
#else
#  define SYSTEM_maxint64    9223372036854775807LL
#  define SYSTEM_minint64  (-9223372036854775807LL-1) /* w/o LL gives warnings */
#  define SYSTEM_maxuint64  18446744073709551615ULL
#endif
#if defined(_WIN64)
# define SYSTEM_maxnativeuint 18446744073709551615U
#elif defined(_WIN32)
#  define SYSTEM_maxnativeuint 4294967295U /* 2**32-1 */
#elif defined(__WORDSIZE)
#  if 64 == __WORDSIZE
#    define SYSTEM_maxnativeuint 18446744073709551615UL
#  else
#    define SYSTEM_maxnativeuint 4294967295UL
#  endif
#elif defined(__SIZEOF_POINTER__)
#  if 8 == __SIZEOF_POINTER__
#    define SYSTEM_maxnativeuint 18446744073709551615UL
#  else
#    define SYSTEM_maxnativeuint 4294967295UL
#  endif
#elif (defined(__sparcv9) || defined(__HOS_AIX__))
/* assume AIX is 64-bit */
#  define SYSTEM_maxnativeuint 18446744073709551615UL
#elif defined(__sparc)
/* check __sparc after __sparcv9, both are defined for 64-bit */
#  define SYSTEM_maxnativeuint 4294967295UL
#else
#  error "no SYSTEM_maxnativeuint defined!!"
#endif

/*** Classes and Interfaces ***/

/* Each class type has an associated const of type SYSTEM_classdescriptor_t. */
/* A type identifier, or a class reference, is a SYSTEM_classreference_t.    */

struct SYSTEM_classdescriptor;
struct _P3InterfaceDescriptor;

typedef const struct SYSTEM_classdescriptor *SYSTEM_classreference_t;
typedef struct _P3InterfaceDescriptor *SYSTEM_interfacereference_t;

typedef struct SYSTEM_classdescriptor {
  const SYSTEM_char *name;         /* class' name */
  const struct SYSTEM_classdescriptor *ancestor; /* null for TObject */
  const SYSTEM_interfacereference_t *IF; /* Array of interfaces */
  SYSTEM_integer IF_count;         /* Number of interfaces */
  SYSTEM_longint CLS_size;         /* Size of an object */
  void * const *VT;                /* Virtual table */
  void **IT;                       /* Interface table */
} SYSTEM_classdescriptor_t;


typedef struct _P3InterfaceDescriptor {
  int i; /* nothing yet */
} _P3InterfaceDescriptor;

extern const SYSTEM_classdescriptor_t SYSTEM_tobject_CD;

typedef struct SYSTEM_tobject_OD_S *SYSTEM_tobject;
typedef struct SYSTEM_tobject_OD_S {
  const SYSTEM_classdescriptor_t *CD; /* = &SYSTEM_tobject_CD */
} SYSTEM_tobject_OD;

#define _P3get_CD(ref) ((ref)->CD) /* Get ref's ClassDescriptor */

/* A meta-class is just a reference to a class descriptor:
typedef SYSTEM_classreference_t SYSTEM_metaclass;  */


extern SYSTEM_boolean _P3is(const SYSTEM_tobject objref,
                            const SYSTEM_classdescriptor_t *clsref);
extern SYSTEM_tobject _P3as(SYSTEM_tobject objref,
                            const SYSTEM_classdescriptor_t *clsref,
                            const char *file, int line);
#define _P3_is(a,b)   _P3is((SYSTEM_tobject)(a), b)
#define _P3_as(a,b)   _P3as((SYSTEM_tobject)(a), b, __FILE__, __LINE__)


/* Implementations of methods from TObject:  */

#define SYSTEM_tobject_DOT_create(Self)       (Self)
#define SYSTEM_tobject_DOT_classname(a,b,c)   _P3_strcpy(a,b,(c)->name)
#define SYSTEM_tobject_DOT_instancesize(c)    ((c)->CLS_size)
#define SYSTEM_tobject_DOT_classnameis(a,b)   _P3streq_ic((a)->name, b)
#define SYSTEM_tobject_DOT_classtype(c)       _P3get_CD(c)
#define SYSTEM_tobject_DOT_classparent(c)     ((c)->ancestor)

/* Allocation/deallocation of objects:                                     */
/* constructors call _P3_alloc_object; destructors call _P3_dealloc_object */
extern SYSTEM_tobject _P3_alloc_object(const SYSTEM_classdescriptor_t *CD);
extern void _P3_dealloc_object(void* p);
/* Do not use the syntax in the next lines: kill Emacs auto-format
 * typedef Destructor(SYSTEM_tobject)(*SYSTEM_tobject_DOT_destroy_T)(SYSTEM_tobject);
 * extern  Destructor(SYSTEM_tobject)  SYSTEM_tobject_DOT_destroy   (SYSTEM_tobject);
 */
typedef SYSTEM_tobject (* SYSTEM_tobject_DOT_destroy_T)(SYSTEM_tobject);
extern  SYSTEM_tobject    SYSTEM_tobject_DOT_destroy   (SYSTEM_tobject);
extern void SYSTEM_tobject_DOT_free(SYSTEM_tobject);

/* NOTE: If you get segmentation errors/illegal references in the
         Virt...MethodCalls with a destructor it may be because the object
         is already released (_P3get_CD(ref) = NULL).                   */
#define VirtMethodCall(ref, typ, ix, prm) ((*(typ)(_P3get_CD(ref)->VT[ix])) prm)
#define VirtClassMethodCall(ref, typ, ix, prm) ((*(typ)((ref)->VT[ix])) prm)

extern void _P3_abstract_call1(void *ref, ...);
extern void _P3_abstract_call3(SYSTEM_ansichar *result, SYSTEM_uint8 _len_ret,
                               SYSTEM_tobject o, ...);
extern void _P3_abstr_cl_call1(void *ref, ...);
extern void _P3_abstr_cl_call3(SYSTEM_ansichar *result, SYSTEM_uint8 _len_ret,
                               SYSTEM_classreference_t ref, ...);
/*** End Classes and Interfaces ***/


/* Shared Objects/DLLs material: */
#ifdef DOS
  typedef HINSTANCE SYSTEM_handle;
#else
  typedef void* SYSTEM_handle;
#endif

/*************** NEW EXCEPTIONS - Object-based like Delphi/Kylix.         */
/***************                  Declarations for SYSTEM_exception class */
/*************** The following is generated by P3 so it'd better work!    */
/*************** I added comments so I can understand it myself...        */

/* A variable of type "exception" is a reference to an object descriptor: */
typedef struct SYSTEM_exception_OD_S* SYSTEM_exception; /* sy_class */

/* The object descriptor is a struct that contains the object's stuff: */
typedef struct SYSTEM_exception_OD_S { /* Objects of 'exception' */
  const SYSTEM_classdescriptor_t *CD; /* = &SYSTEM_exception_CD*/
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} SYSTEM_exception_OD;

/* Prototype for "exception"'s constructor create(const Msg: string): */
Constructor(SYSTEM_exception)
SYSTEM_exception_DOT_create(SYSTEM_exception self, const SYSTEM_char *msg);

/* Declaration of the Virtual Table, defined in p3io.c */
extern void * const SYSTEM_exception_VT[];

/* Class Descriptor for an exception. To identify it - contains name etc. */
extern const SYSTEM_classdescriptor_t SYSTEM_exception_CD;

/* Put message as the output of next one to get test output */
#define EXCTST(message) /* message */

extern void _P3_Exception (int n, const char *s);

#define _P3_exceptobject _P3_ExceptObject  /* local or global, depending on _P3_EXC_MODEL_ */
#define SYSTEM_exceptobject() _P3_ExceptObject

#if _P3_EXC_MODEL_ == 1
/* C++ try/catch */

class exWrap: public std::exception
{
 private:
  const char *mmm;
  SYSTEM_tobject payload;
 public:
  exWrap () : mmm(NULL), payload(NULL) {  }
  exWrap (const char *msg) : mmm(msg), payload(NULL) { }
  exWrap (SYSTEM_tobject o) : mmm(NULL), payload(o) { }
  exWrap (const char *msg, SYSTEM_tobject o) : mmm(msg), payload(o) { }
  virtual const char *what() const throw()
  {
    if (mmm)
      return mmm;
    else
      return "exWrap exception";
  } // what
  inline SYSTEM_tobject getPayload () { return payload; }
  inline void setPayload (SYSTEM_tobject o) { payload = o; }
}; // exWrap

class exFinal: public std::exception
{
 private:
  const char *mmm;
 public:
  exFinal () : mmm(NULL) {  }
  exFinal (const char *msg) : mmm(msg) { }
  virtual const char *what() const throw()
  {
    if (mmm)
      return mmm;
    else
      return "exFinal exception";
  } // what
}; // exFinal

#define _P3_TRY try {

#define _P3_RERAISE()    throw
#define _P3_RAISE(obj)   throw exWrap("_P3_RAISE",(SYSTEM_tobject)obj)

#if _P3_EXC_VERBOSE_ >= 2
# define _P3_EXCEPT } catch (exWrap& EW) { \
                      std::cout << "A p3-wrapper exception caught: " << EW.what() << std::endl; \
                      SYSTEM_tobject _P3_exceptobject = EW.getPayload();
#else
# define _P3_EXCEPT } catch (exWrap& EW) { \
                      SYSTEM_tobject _P3_exceptobject = EW.getPayload();
#endif /* if  _P3_EXC_VERBOSE_ >= 2 .. else .. */

#define _P3_END_TRY_EXCEPT  SYSTEM_tobject_DOT_free(_P3_exceptobject); }

#define _P3_ON_CLAUSE(class_CD)  if (_P3_is(_P3_exceptobject, class_CD))
#define _P3_END_ON_CLAUSE        else
#define _P3_ON_ELSE
#define _P3_END_ON_ELSE
#define _P3_EMPTY_ELSE_CLAUSE   _P3_RERAISE()

#if _P3_EXC_VERBOSE_ >= 2
# define _P3_FINALLY(LAB)  throw exFinal("_P3_FINALLY");  \
                         } catch (std::exception& EE) {   \
             std::cout << "A std::exception caught: " << EE.what() << std::endl;

# define _P3_END_TRY_FINALLY  \
    if (0 == dynamic_cast<exFinal*>(&EE)) { \
      std::cout << "really an exception: rethrow" << std::endl; \
      throw; \
     } \
   }

# define _P3_PGM_OUTER_CATCH  catch (exWrap& EW) { \
     std::cout << "A p3-wrapper exception caught in PGM_OUTER: " << EW.what() << std::endl; \
     SYSTEM_tobject _P3_exceptobject = EW.getPayload(); \
     std::cout << "HERE WE GO" << std::endl; \
     _P3_Std_Exception_Handler (_P3_exceptobject); \
   }
#else
# define _P3_FINALLY(LAB)  throw exFinal("_P3_FINALLY");  \
                         } catch (std::exception& EE) {

# define _P3_END_TRY_FINALLY  \
    if (0 == dynamic_cast<exFinal*>(&EE)) { \
      throw; \
     } \
   }

# define _P3_PGM_OUTER_CATCH  catch (exWrap& EW) { \
     SYSTEM_tobject _P3_exceptobject = EW.getPayload(); \
     _P3_Std_Exception_Handler (_P3_exceptobject); \
   }
#endif /* if _P3_EXC_VERBOSE_ >= 2 .. else .. */

#define _P3_PGM_OUTER_TRY    try

#define _P3_DLL_OUTER_TRY    _P3_DLL_OUTER_TRY    is_not_yet_implemented
#define _P3_DLL_OUTER_CATCH  _P3_DLL_OUTER_CATCH  is_not_yet_implemented

#endif  /* if _P3_EXC_MODEL_ == 1 */

/*************** END NEW EXCEPTIONS      ************************************/

#include "exceptions.h"


#endif /* ifndef _P3IO_H */
/* end of p3io.h */
