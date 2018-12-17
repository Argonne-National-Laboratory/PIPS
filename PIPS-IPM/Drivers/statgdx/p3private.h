#ifndef _P3___p3private___H
#define _P3___p3private___H

typedef SYSTEM_uint8 _sub_0P3PRIVATE;
typedef SYSTEM_ansichar P3PRIVATE_shortstrbuf[256];

Procedure P3PRIVATE_pcharconcatpchar(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  SYSTEM_P3_pansichar psrc);

Procedure P3PRIVATE_pcharconcatstr(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  const SYSTEM_ansichar *src);

Function(SYSTEM_ansichar *) P3PRIVATE_pchartostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pansichar psrc);

Function(SYSTEM_P3_pansichar ) P3PRIVATE_strtopchar(
  const SYSTEM_ansichar *src);

Function(SYSTEM_P3_pansichar ) P3PRIVATE_strtostrbuf(
  const SYSTEM_ansichar *src,
  SYSTEM_ansichar *dest);

Function(SYSTEM_ansichar *) P3PRIVATE_strbuftostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *src);

extern void _Init_Module_p3private(void);
extern void _Final_Module_p3private(void);

#endif /* ! defined _P3___p3private___H */
