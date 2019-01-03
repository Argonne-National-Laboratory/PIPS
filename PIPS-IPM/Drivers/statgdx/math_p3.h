#ifndef _P3___math_p3___H
#define _P3___math_p3___H

extern SYSTEM_double MATH_P3_minsingle;
extern SYSTEM_double MATH_P3_maxsingle;
extern SYSTEM_double MATH_P3_mindouble;
extern SYSTEM_double MATH_P3_maxdouble;

Function(SYSTEM_double ) MATH_P3_arccos(
  SYSTEM_double x);

Function(SYSTEM_double ) MATH_P3_arcsin(
  SYSTEM_double x);

Function(SYSTEM_double ) MATH_P3_arctan2(
  SYSTEM_double y,
  SYSTEM_double x);

Function(SYSTEM_double ) MATH_P3_tan(
  SYSTEM_double x);

Function(SYSTEM_double ) MATH_P3_intpower(
  SYSTEM_double x,
  SYSTEM_integer i);

Function(SYSTEM_double ) MATH_P3_power(
  SYSTEM_double x,
  SYSTEM_double y);

Function(SYSTEM_double ) MATH_P3_roundto(
  SYSTEM_double x,
  SYSTEM_integer i);

Function(SYSTEM_boolean ) MATH_P3_isnan(
  SYSTEM_double avalue);

Function(SYSTEM_boolean ) MATH_P3_isinfinite(
  SYSTEM_double avalue);
typedef SYSTEM_byte MATH_P3_tfpuexception; /* Anonymous */ enum{MATH_P3_exinvalidop,MATH_P3_exdenormalized,MATH_P3_exzerodivide,
  MATH_P3_exoverflow,MATH_P3_exunderflow,MATH_P3_exprecision};
typedef _P3SET_7 MATH_P3_tfpuexceptionmask;

Function(_P3set_elem *) MATH_P3_getexceptionmask(
  _P3set_elem *result,
  SYSTEM_uint8 _len_ret);

Function(_P3set_elem *) MATH_P3_setexceptionmask(
  _P3set_elem *result,
  SYSTEM_uint8 _len_ret,
  const _P3set_elem *mask);

Procedure MATH_P3_setexceptionmask2p3(void);

Procedure MATH_P3_clearexceptions(void);

extern void _Init_Module_math_p3(void);
extern void _Final_Module_math_p3(void);

#endif /* ! defined _P3___math_p3___H */
