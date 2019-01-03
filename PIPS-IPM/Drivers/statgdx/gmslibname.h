#ifndef _P3___gmslibname___H
#define _P3___gmslibname___H


Function(SYSTEM_ansichar *) GMSLIBNAME_gamslibnameold(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) GMSLIBNAME_gamslibnamep3(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

extern void _Init_Module_gmslibname(void);
extern void _Final_Module_gmslibname(void);

#endif /* ! defined _P3___gmslibname___H */
