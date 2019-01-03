#ifndef _P3___system_p3___H
#define _P3___system_p3___H

typedef SYSTEM_double SYSTEM_P3_tdatetime;
typedef SYSTEM_shortstring *SYSTEM_P3_pshortstring;
typedef SYSTEM_byte *SYSTEM_P3_pbyte;
typedef SYSTEM_ansichar *SYSTEM_P3_pansichar;
typedef SYSTEM_P3_pansichar SYSTEM_P3_pchar;
typedef SYSTEM_word *SYSTEM_P3_pword;
typedef SYSTEM_integer *SYSTEM_P3_pinteger;
typedef SYSTEM_double *SYSTEM_P3_pdouble;
typedef SYSTEM_longint *SYSTEM_P3_plongint;
typedef SYSTEM_pointer *SYSTEM_P3_ppointer;
typedef SYSTEM_longword SYSTEM_P3_thandle;
typedef SYSTEM_uint32 _sub_0SYSTEM_P3;
typedef SYSTEM_integer SYSTEM_P3_integerarray[251658240];
typedef SYSTEM_P3_integerarray *SYSTEM_P3_pintegerarray;
typedef SYSTEM_uint32 _sub_1SYSTEM_P3;
typedef SYSTEM_pointer SYSTEM_P3_pointerarray[268435455];
typedef SYSTEM_P3_pointerarray *SYSTEM_P3_ppointerarray;

Procedure SYSTEM_P3_getdir(
  SYSTEM_byte d,
  SYSTEM_ansichar *s);

Procedure SYSTEM_P3_fillchar(
  SYSTEM_untyped *p,
  SYSTEM_integer len,
  SYSTEM_byte v);

Procedure SYSTEM_P3_settextbuf(
  SYSTEM_text *f,
  SYSTEM_untyped *p);

Function(SYSTEM_integer ) SYSTEM_P3_paramcount(void);

Function(SYSTEM_ansichar *) SYSTEM_P3_paramstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer index);

extern void _Init_Module_system_p3(void);
extern void _Final_Module_system_p3(void);

#endif /* ! defined _P3___system_p3___H */
