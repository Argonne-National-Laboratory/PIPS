#include "p3io.h"
#include "math_p3.h"

SYSTEM_double MATH_P3_minsingle = 1.5e-45;
SYSTEM_double MATH_P3_maxsingle = 3.4e38;
SYSTEM_double MATH_P3_mindouble = 5.0e-324;
SYSTEM_double MATH_P3_maxdouble = 1.7e308;
typedef struct MATH_P3_ti64rec_S {
  union{
    struct{
      SYSTEM_double x;
    } _c1;
    struct{
      SYSTEM_int64 i64;
    } _c2;
  } _u;
} MATH_P3_ti64rec;

static SYSTEM_int64 MATH_P3_signmask;
static SYSTEM_int64 MATH_P3_expomask;
static SYSTEM_int64 MATH_P3_mantmask;

Function(SYSTEM_double ) MATH_P3_arccos(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(138:1): 1 lines ****/
  result = acos(x);
  return result;
}  /* arccos */

Function(SYSTEM_double ) MATH_P3_arcsin(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(149:1): 1 lines ****/
  result = asin(x);
  return result;
}  /* arcsin */

Function(SYSTEM_double ) MATH_P3_arctan2(
  SYSTEM_double y,
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(164:1): 1 lines ****/
  result = atan2(y,x);
  return result;
}  /* arctan2 */

Function(SYSTEM_double ) MATH_P3_tan(
  SYSTEM_double x)
{
  SYSTEM_double result;

  /**** C code included from math_p3.pas(180:1): 1 lines ****/
  result = tan(x);
  return result;
}  /* tan */

Function(SYSTEM_double ) MATH_P3_intpower(
  SYSTEM_double x,
  SYSTEM_integer i)
{
  SYSTEM_double result;
  SYSTEM_integer y;

  y = SYSTEM_abs_i(i);
  result = 1.0;
  while (y > 0) {
    while (!SYSTEM_odd(y)) {
      y = ValueCast(SYSTEM_uint32,y) >> 1;
      x = x * x;
    
}
    _P3dec0(y);
    result = result * x;
  
}
  if (i < 0) 
    result = 1.0 /  result;
  return result;
}  /* intpower */

Function(SYSTEM_double ) MATH_P3_power(
  SYSTEM_double x,
  SYSTEM_double y)
{
  SYSTEM_double result;

  if (y == 0.0) { 
    result = 1.0;
  } else 
    if (x == 0.0 && y > 0.0) { 
      result = 0.0;
    } else 
      if (SYSTEM_frac(y) == 0.0 && SYSTEM_abs_r(y) <= SYSTEM_maxint) { 
        result = MATH_P3_intpower(x,SYSTEM_trunc(y));
      } else 
        result = SYSTEM_exp(y * SYSTEM_ln(x));
  return result;
}  /* power */

Function(SYSTEM_double ) MATH_P3_roundto(
  SYSTEM_double x,
  SYSTEM_integer i)
{
  SYSTEM_double result;
  SYSTEM_double z;

  if (i == 0) { 
    if (x > 0) { 
      result = SYSTEM_int(x + 0.5);
    } else 
      result = SYSTEM_int(x - 0.5);
  } else 
    if (i > 0) {
      z = MATH_P3_intpower(10,i);
      if (x > 0) { 
        result = SYSTEM_int(x * z + 0.5) /  z;
      } else 
        result = SYSTEM_int(x * z - 0.5) /  z;
    } else {
      z = MATH_P3_intpower(10,-i);
      if (x > 0) { 
        result = SYSTEM_int(x /  z + 0.5) * z;
      } else 
        result = SYSTEM_int(x /  z - 0.5) * z;
    } 
  return result;
}  /* roundto */

Function(SYSTEM_boolean ) MATH_P3_isnan(
  SYSTEM_double avalue)
{
  SYSTEM_boolean result;
  MATH_P3_ti64rec i64rec;
  SYSTEM_int64 mantissa;

  result = SYSTEM_false;
  i64rec._u._c1.x = avalue;
  if (ValueCast(SYSTEM_uint64,i64rec._u._c2.i64 & MATH_P3_expomask) >> 52 == 2047) {
    mantissa = i64rec._u._c2.i64 & MATH_P3_mantmask;
    if (mantissa != 0) 
      result = SYSTEM_true;
  } 
  return result;
}  /* isnan */

Function(SYSTEM_boolean ) MATH_P3_isinfinite(
  SYSTEM_double avalue)
{
  SYSTEM_boolean result;
  MATH_P3_ti64rec i64rec;
  SYSTEM_int64 mantissa;

  result = SYSTEM_false;
  i64rec._u._c1.x = avalue;
  if (ValueCast(SYSTEM_uint64,i64rec._u._c2.i64 & MATH_P3_expomask) >> 52 == 2047) {
    mantissa = i64rec._u._c2.i64 & MATH_P3_mantmask;
    if (mantissa == 0) 
      result = SYSTEM_true;
  } 
  return result;
}  /* isinfinite */

/* unit math_p3 */
void _Init_Module_math_p3(void)
{
  MATH_P3_signmask =  SYSTEM_minint;
  MATH_P3_signmask = ValueCast(SYSTEM_uint64,MATH_P3_signmask) << 32;
  MATH_P3_expomask = 2146435072;
  MATH_P3_expomask = ValueCast(SYSTEM_uint64,MATH_P3_expomask) << 32;
  MATH_P3_mantmask = ~(MATH_P3_signmask | MATH_P3_expomask);
} /* _Init_Module_math_p3 */

void _Final_Module_math_p3(void)
{
} /* _Final_Module_math_p3 */

