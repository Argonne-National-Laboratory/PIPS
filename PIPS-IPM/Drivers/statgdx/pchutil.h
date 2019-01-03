#ifndef _P3___pchutil___H
#define _P3___pchutil___H

typedef SYSTEM_uint8 _sub_0PCHUTIL;
typedef SYSTEM_ansichar PCHUTIL_shortstrbuf[256];

Function(SYSTEM_P3_pansichar ) PCHUTIL_strtopchar(
  const SYSTEM_ansichar *s);

Function(SYSTEM_P3_pansichar ) PCHUTIL_emptytopchar(void);

Function(SYSTEM_ansichar *) PCHUTIL_pchartostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pansichar p);

Procedure PCHUTIL_convertpchar(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *s);

Function(SYSTEM_integer ) PCHUTIL_pcharlen(
  SYSTEM_P3_pansichar p);

Procedure PCHUTIL_pcharconcatpchar(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  SYSTEM_P3_pansichar psrc);

Procedure PCHUTIL_pcharconcatstr(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  const SYSTEM_ansichar *src);

Procedure PCHUTIL_strpcopyn(
  SYSTEM_P3_pansichar pdest,
  const SYSTEM_ansichar *src,
  SYSTEM_integer n);

Function(SYSTEM_P3_pansichar ) PCHUTIL_strtostrbuf(
  const SYSTEM_ansichar *src,
  SYSTEM_ansichar *dest);

extern void _Init_Module_pchutil(void);
extern void _Final_Module_pchutil(void);

#endif /* ! defined _P3___pchutil___H */
