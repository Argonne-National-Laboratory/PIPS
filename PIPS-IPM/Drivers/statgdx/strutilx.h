#ifndef _P3___strutilx___H
#define _P3___strutilx___H


Function(SYSTEM_P3_pshortstring ) STRUTILX_newstring(
  const SYSTEM_ansichar *s);

Procedure STRUTILX_disposestring(
  SYSTEM_P3_pshortstring ps);

Function(SYSTEM_ansichar *) STRUTILX_getstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pshortstring ps);

Procedure STRUTILX_strassign(
  SYSTEM_P3_pshortstring *p,
  const SYSTEM_ansichar *s);

Function(SYSTEM_P3_pshortstring ) STRUTILX_newstringm(
  const SYSTEM_ansichar *s,
  SYSTEM_int64 *m);

Procedure STRUTILX_disposestringm(
  SYSTEM_P3_pshortstring ps,
  SYSTEM_int64 *m);

Procedure STRUTILX_strassignm(
  SYSTEM_P3_pshortstring *p,
  const SYSTEM_ansichar *s,
  SYSTEM_int64 *m);

Function(SYSTEM_integer ) STRUTILX_strcmp(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_integer ) STRUTILX_pstrcmp(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2);

Function(SYSTEM_boolean ) STRUTILX_strequal(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_boolean ) STRUTILX_pstrequal(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2);

Function(SYSTEM_integer ) STRUTILX_strucmp(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_integer ) STRUTILX_strucmpnum(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_integer ) STRUTILX_pstrucmp(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2);

Function(SYSTEM_boolean ) STRUTILX_struequal(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_boolean ) STRUTILX_pstruequal(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2);

Function(SYSTEM_integer ) STRUTILX_lchpossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp);

Function(SYSTEM_integer ) STRUTILX_lchpos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_rchpossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp);

Function(SYSTEM_integer ) STRUTILX_rchpos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_lchsetpos(
  const _P3set_elem *cs,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_rchsetpos(
  const _P3set_elem *cs,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_lstrpossp(
  const SYSTEM_ansichar *pat,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp);

Function(SYSTEM_integer ) STRUTILX_lstrpos(
  const SYSTEM_ansichar *pat,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_lchupossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp);

Function(SYSTEM_integer ) STRUTILX_lchupos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_rchupossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp);

Function(SYSTEM_integer ) STRUTILX_rchupos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) STRUTILX_padleft(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer w);

Function(SYSTEM_ansichar *) STRUTILX_padright(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer w);

Function(SYSTEM_ansichar *) STRUTILX_padrightmod(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer m);

Function(SYSTEM_integer ) STRUTILX_padmodlength(
  const SYSTEM_ansichar *s,
  SYSTEM_integer m);

Function(SYSTEM_ansichar *) STRUTILX_inttostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) STRUTILX_inttostrex(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) STRUTILX_inttostrw(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n,
  SYSTEM_integer width);

Function(SYSTEM_ansichar *) STRUTILX_inttonicestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) STRUTILX_int64tonicestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n);

Function(SYSTEM_ansichar *) STRUTILX_inttonicestrw(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n,
  SYSTEM_integer width);

Function(SYSTEM_ansichar *) STRUTILX_int64tonicestrw(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n,
  SYSTEM_integer width);

Function(SYSTEM_ansichar *) STRUTILX_dbltostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v);

Function(SYSTEM_ansichar *) STRUTILX_dbltostrex(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v);

Function(SYSTEM_ansichar *) STRUTILX_gmsvaltostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v);

Function(SYSTEM_ansichar *) STRUTILX_mem64tonicestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n,
  SYSTEM_integer w);

Function(SYSTEM_ansichar *) STRUTILX_replacechar(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const _P3set_elem *chset,
  SYSTEM_ansichar _new,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) STRUTILX_deletechar(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const _P3set_elem *chset,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) STRUTILX_replacestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *old,
  const SYSTEM_ansichar *_new,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) STRUTILX_uppercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar ) STRUTILX_lowcase(
  SYSTEM_ansichar ch);

Function(SYSTEM_ansichar *) STRUTILX_lowercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRUTILX_integerwidth(
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) STRUTILX_fillstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_ansichar ch,
  SYSTEM_integer len);

Function(SYSTEM_ansichar *) STRUTILX_blankstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer len);

Function(SYSTEM_ansichar *) STRUTILX_extracttoken(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer *p);

Function(SYSTEM_integer ) STRUTILX_strasint(
  const SYSTEM_ansichar *s);

Function(SYSTEM_boolean ) STRUTILX_strasintex(
  const SYSTEM_ansichar *s,
  SYSTEM_integer *v);

Function(SYSTEM_boolean ) STRUTILX_strasintex2(
  const SYSTEM_ansichar *s,
  SYSTEM_integer *v);

Function(SYSTEM_boolean ) STRUTILX_strasdoubleex(
  const SYSTEM_ansichar *s,
  SYSTEM_double *v);

Function(SYSTEM_ansichar *) STRUTILX_excelcolstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer c);

Function(SYSTEM_integer ) STRUTILX_strexcelcol(
  const SYSTEM_ansichar *s);
extern SYSTEM_shortstring STRUTILX_blanks255;

extern void _Init_Module_strutilx(void);
extern void _Final_Module_strutilx(void);

#endif /* ! defined _P3___strutilx___H */
