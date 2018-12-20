#include "p3io.h"
#include "p3platform.h"
#include "exceptions.h"
#include "math_p3.h"
#include "system_p3.h"
#include "p3ieeefp.h"

_arr_0P3IEEEFP P3IEEEFP_fpclasstext = {{13,'s','i','g','n','a','l','i','n','g',' ','N','a','N'}, {9,'q','u','i','e','t',' ','N','a','N'}, {17,'n','e','g','a','t','i','v','e',' ','i','n','f','i','n','i','t','y'}, {17,'p','o','s','i','t','i','v','e',' ','i','n','f','i','n','i','t','y'}, {30,'n','e','g','a','t','i','v','e',' ','d','e','n','o','r','m','a','l','i','z','e','d',' ','n','o','n','-','z','e','r','o'}, {30,'p','o','s','i','t','i','v','e',' ','d','e','n','o','r','m','a','l','i','z','e','d',' ','n','o','n','-','z','e','r','o'}, {4,'-','0','.','0'}, {4,'+','0','.','0'}, {28,'n','e','g','a','t','i','v','e',' ','n','o','r','m','a','l','i','z','e','d',' ','n','o','n','-','z','e','r','o'}, {28,'p','o','s','i','t','i','v','e',' ','n','o','r','m','a','l','i','z','e','d',' ','n','o','n','-','z','e','r','o'}};
SYSTEM_double P3IEEEFP_nanquiet;
SYSTEM_double P3IEEEFP_nansignaling;
SYSTEM_double P3IEEEFP_infpositive;
SYSTEM_double P3IEEEFP_infnegative;
/**** C code included from p3ieeefp.pas(99:1): 22 lines ****/
#if   defined(AIX)
# include <float.h>
# include <fpxcp.h>
# include <fptrap.h>
#elif defined(__APPLE__)
# include <fenv.h>
#elif defined(BGP) || defined(__linux__)
# include <fpu_control.h>
# include <fenv.h>
#elif defined(SOL)
# include <ieeefp.h>
#endif

#define ADD2MASK(X) _P3SET_p(result,_len_ret,result,X)
#define ISINFLAGS(X) (_P3SET_ic(5,X,flags))
#define EX_INVALIDOP   _P3set1("\001")
#define EX_DENORMAL    _P3set1("\002")
#define EX_ZERODIVIDE  _P3set1("\004")
#define EX_OVERFLOW    _P3set1("\010")
#define EX_UNDERFLOW   _P3set1("\020")
#define EX_PRECISION   _P3set1("\040")

typedef struct P3IEEEFP_ti64rec_S {
  union{
    struct{
      SYSTEM_double x;
    } _c1;
    struct{
      SYSTEM_int64 i64;
    } _c2;
  } _u;
} P3IEEEFP_ti64rec;

static SYSTEM_int64 P3IEEEFP_signmask;
static SYSTEM_int64 P3IEEEFP_expomask;
static SYSTEM_int64 P3IEEEFP_mantmask;
static SYSTEM_int64 P3IEEEFP_qnanmask;

static Procedure P3IEEEFP_double2i64(
  SYSTEM_double x,
  SYSTEM_int64 *i)
{
  P3IEEEFP_ti64rec i64rec;

  i64rec._u._c1.x = x;
  *i = i64rec._u._c2.i64;
}  /* double2i64 */

Function(P3IEEEFP_tfpclass ) P3IEEEFP_fpclass(
  SYSTEM_double x)
{
  P3IEEEFP_tfpclass result;
  SYSTEM_boolean isneg, isquiet;
  SYSTEM_int64 exponent;
  SYSTEM_int64 mantissa, ix;

  P3IEEEFP_double2i64(x,&ix);
  isneg = (ix & P3IEEEFP_signmask) == P3IEEEFP_signmask;
  exponent = ValueCast(SYSTEM_uint64,ix & P3IEEEFP_expomask) >> 52;
  isquiet = (ix & P3IEEEFP_qnanmask) == P3IEEEFP_qnanmask;
  mantissa = ix & P3IEEEFP_mantmask;
  if (exponent == 0) { 
    if (mantissa != 0) { 
      if (isneg) { 
        result = P3IEEEFP_fp_ndenorm;
      } else 
        result = P3IEEEFP_fp_pdenorm;
    } else 
      if (isneg) { 
        result = P3IEEEFP_fp_nzero;
      } else 
        result = P3IEEEFP_fp_pzero;
  } else 
    if (exponent == 2047) { 
      if (mantissa == 0) { 
        if (isneg) { 
          result = P3IEEEFP_fp_ninf;
        } else 
          result = P3IEEEFP_fp_pinf;
      } else 
        if (isquiet) { 
          result = P3IEEEFP_fp_qnan;
        } else 
          result = P3IEEEFP_fp_snan;
    } else 
      if (isneg) { 
        result = P3IEEEFP_fp_nnorm;
      } else 
        result = P3IEEEFP_fp_pnorm;
  return result;
}  /* fpclass */
typedef struct P3IEEEFP_fpenvrec_S {
  SYSTEM_word controlword;
  SYSTEM_word fill1;
  SYSTEM_word statusword;
  SYSTEM_word fill2;
  SYSTEM_word tagword;
  SYSTEM_word fill3;
  SYSTEM_longword instrpointer;
  SYSTEM_word codesegment;
  SYSTEM_word fill5;
  SYSTEM_longword operandaddress;
  SYSTEM_word fill7;
} P3IEEEFP_fpenvrec;


Function(SYSTEM_boolean ) P3IEEEFP_p3isfinite(
  SYSTEM_double x)
{
  SYSTEM_boolean result;
  P3IEEEFP_ti64rec i64rec;

  i64rec._u._c1.x = x;
  result = ValueCast(SYSTEM_uint64,i64rec._u._c2.i64 & 
    P3IEEEFP_expomask) >> 52 != 2047;
  return result;
}  /* p3isfinite */
static P3IEEEFP_ti64rec P3IEEEFP_i64rec;

/* unit p3ieeefp */
void _Init_Module_p3ieeefp(void)
{
  P3IEEEFP_signmask =  SYSTEM_minint;
  P3IEEEFP_signmask = ValueCast(SYSTEM_uint64,P3IEEEFP_signmask) << 32;
  P3IEEEFP_expomask = 2146435072;
  P3IEEEFP_expomask = ValueCast(SYSTEM_uint64,P3IEEEFP_expomask) << 32;
  P3IEEEFP_qnanmask = 524288;
  P3IEEEFP_qnanmask = ValueCast(SYSTEM_uint64,P3IEEEFP_qnanmask) << 32;
  P3IEEEFP_mantmask = ~(P3IEEEFP_signmask | P3IEEEFP_expomask);
  P3IEEEFP_i64rec._u._c2.i64 =  -1;
  P3IEEEFP_i64rec._u._c2.i64 = ValueCast(SYSTEM_uint64,P3IEEEFP_i64rec.
    _u._c2.i64) << 32;
  P3IEEEFP_nanquiet = P3IEEEFP_i64rec._u._c1.x;
  P3IEEEFP_i64rec._u._c2.i64 =  -524289;
  P3IEEEFP_i64rec._u._c2.i64 = ValueCast(SYSTEM_uint64,P3IEEEFP_i64rec.
    _u._c2.i64) << 32;
  P3IEEEFP_nansignaling = P3IEEEFP_i64rec._u._c1.x;
  P3IEEEFP_i64rec._u._c2.i64 = P3IEEEFP_expomask;
  P3IEEEFP_infpositive = P3IEEEFP_i64rec._u._c1.x;
  P3IEEEFP_i64rec._u._c2.i64 = P3IEEEFP_expomask | P3IEEEFP_signmask;
  P3IEEEFP_infnegative = P3IEEEFP_i64rec._u._c1.x;
  P3IEEEFP_i64rec._u._c2.i64 = 0;
} /* _Init_Module_p3ieeefp */

void _Final_Module_p3ieeefp(void)
{
} /* _Final_Module_p3ieeefp */

