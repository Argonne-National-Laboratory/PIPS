#ifndef _P3___p3ieeefp___H
#define _P3___p3ieeefp___H

typedef SYSTEM_byte P3IEEEFP_tfpclass; /* Anonymous */ enum{P3IEEEFP_fp_snan,P3IEEEFP_fp_qnan,P3IEEEFP_fp_ninf,
  P3IEEEFP_fp_pinf,P3IEEEFP_fp_ndenorm,P3IEEEFP_fp_pdenorm,
  P3IEEEFP_fp_nzero,P3IEEEFP_fp_pzero,P3IEEEFP_fp_nnorm,
  P3IEEEFP_fp_pnorm};
typedef SYSTEM_shortstring _arr_0P3IEEEFP[10];
extern _arr_0P3IEEEFP P3IEEEFP_fpclasstext;
extern SYSTEM_double P3IEEEFP_nanquiet;
extern SYSTEM_double P3IEEEFP_nansignaling;
extern SYSTEM_double P3IEEEFP_infpositive;
extern SYSTEM_double P3IEEEFP_infnegative;
typedef _P3SET_7 P3IEEEFP_tfpuexceptionflags;

Function(P3IEEEFP_tfpclass ) P3IEEEFP_fpclass(
  SYSTEM_double x);

Function(_P3set_elem *) P3IEEEFP_getexceptionflags(
  _P3set_elem *result,
  SYSTEM_uint8 _len_ret);

Procedure P3IEEEFP_clearexceptionflags(
  const _P3set_elem *flags);

Function(SYSTEM_boolean ) P3IEEEFP_p3isfinite(
  SYSTEM_double x);

extern void _Init_Module_p3ieeefp(void);
extern void _Final_Module_p3ieeefp(void);

#endif /* ! defined _P3___p3ieeefp___H */
