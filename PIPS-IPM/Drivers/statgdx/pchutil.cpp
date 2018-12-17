#include "p3io.h"
#include "p3platform.h"
#include "p3utils.h"
#include "system_p3.h"
#include "p3process.h"
#include "p3library.h"
#include "math_p3.h"
#include "p3ieeefp.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "p3threads.h"
#include "idglobal_p3.h"
#include "gmsspecs.h"
#include "gmsgen.h"
#include "pchutil.h"


Function(SYSTEM_P3_pansichar ) PCHUTIL_strtopchar(
  const SYSTEM_ansichar *s)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer k;
  SYSTEM_integer l;

  l = SYSTEM_length(s);
  if (l == 0) { 
    result = NULL;
  } else {
    _P3getmem(result,l + 1);
    { register SYSTEM_int32 _stop = l;
      if ((k = 1) <=  _stop) do {
        (*ValueCast(GMSGEN_pansichararray,result))[k - 1] = s[k];
      } while (k++ !=  _stop);

    }
    (*ValueCast(GMSGEN_pansichararray,result))[l] = _P3char('\000');
  } 
  return result;
}  /* strtopchar */

Function(SYSTEM_P3_pansichar ) PCHUTIL_emptytopchar(void)
{
  SYSTEM_P3_pansichar result;

  _P3getmem(result,1);
  *result = _P3char('\000');
  return result;
}  /* emptytopchar */

Function(SYSTEM_ansichar *) PCHUTIL_pchartostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pansichar p)
{
  SYSTEM_integer k;

  if (p == NULL) { 
    _P3strclr(result);
  } else {
    k = 0;
    do {
      if ((*ValueCast(GMSGEN_pansichararray,p))[k] == _P3char('\000')) {
        _P3setlength(result,k,255);
        SYSTEM_break(BRK_1);
      } 
      result[k + 1] = (*ValueCast(GMSGEN_pansichararray,p))[k];
      k = k + 1;
      if (k >= 255) {
        _P3strcpy(result,_len_ret,_P3str1("\023PCharToStr Overflow"));
        SYSTEM_break(BRK_1);
      } 
    CNT_1:;
    } while (SYSTEM_true);
BRK_1:;
  } 
  return result;
}  /* pchartostr */

Procedure PCHUTIL_convertpchar(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *s)
{
  SYSTEM_integer k;

  if (p == NULL) { 
    _P3strclr(s);
  } else {
    k = 0;
    do {
      if ((*ValueCast(GMSGEN_pansichararray,p))[k] == _P3char('\000')) {
        _P3setlength(s,k,255);
        SYSTEM_break(BRK_2);
      } 
      s[k + 1] = (*ValueCast(GMSGEN_pansichararray,p))[k];
      k = k + 1;
      if (k >= 255) {
        _P3strcpy(s,255,_P3str1("\025ConvertPChar Overflow"));
        SYSTEM_break(BRK_2);
      } 
    CNT_2:;
    } while (SYSTEM_true);
BRK_2:;
  } 
}  /* convertpchar */

Function(SYSTEM_integer ) PCHUTIL_pcharlen(
  SYSTEM_P3_pansichar p)
{
  SYSTEM_integer result;

  result = 0;
  if (p != NULL) 
    while ((*ValueCast(GMSGEN_pansichararray,p))[result] != _P3char('\000')) {

      _P3inc0(result);
}
  return result;
}  /* pcharlen */

Procedure PCHUTIL_pcharconcatpchar(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  SYSTEM_P3_pansichar psrc)
{
  SYSTEM_integer k;

  if (psrc != NULL) {
    k = 0;
    while ((*ValueCast(GMSGEN_pansichararray,psrc))[k] != _P3char('\000')) {
      (*ValueCast(GMSGEN_pansichararray,pdest))[*w] = (*ValueCast(
        GMSGEN_pansichararray,psrc))[k];
      _P3inc0(*w);
      _P3inc0(k);
    
}
    (*ValueCast(GMSGEN_pansichararray,pdest))[*w] = _P3char('\000');
  } 
}  /* pcharconcatpchar */

Procedure PCHUTIL_pcharconcatstr(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  const SYSTEM_ansichar *src)
{
  SYSTEM_integer k;

  if (SYSTEM_length(src) > 0) {
    { register SYSTEM_int32 _stop = SYSTEM_length(src);
      if ((k = 1) <=  _stop) do {
        (*ValueCast(GMSGEN_pansichararray,pdest))[*w] = src[k];
        _P3inc0(*w);
      
      } while (k++ !=  _stop);

    }
    (*ValueCast(GMSGEN_pansichararray,pdest))[*w] = _P3char('\000');
  } 
}  /* pcharconcatstr */

Procedure PCHUTIL_strpcopyn(
  SYSTEM_P3_pansichar pdest,
  const SYSTEM_ansichar *src,
  SYSTEM_integer n)
{
  SYSTEM_integer k;
  SYSTEM_integer w;
  SYSTEM_integer m;

  w = 0;
  if (n < ValueCast(SYSTEM_int32,SYSTEM_length(src))) { 
    m = n;
  } else 
    m = SYSTEM_length(src);
  { register SYSTEM_int32 _stop = m;
    if ((k = 1) <=  _stop) do {
      (*ValueCast(GMSGEN_pansichararray,pdest))[w] = src[k];
      _P3inc0(w);
    
    } while (k++ !=  _stop);

  }
  (*ValueCast(GMSGEN_pansichararray,pdest))[w] = _P3char('\000');
}  /* strpcopyn */

Function(SYSTEM_P3_pansichar ) PCHUTIL_strtostrbuf(
  const SYSTEM_ansichar *src,
  SYSTEM_ansichar *dest)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer len;

  len = SYSTEM_length(src);
  SYSTEM_move(&src[1],dest,len);
  dest[len] = _P3char('\000');
  result = ValueCast(SYSTEM_P3_pansichar,&dest[0]);
  return result;
}  /* strtostrbuf */

/* unit pchutil */
void _Init_Module_pchutil(void)
{
} /* _Init_Module_pchutil */

void _Final_Module_pchutil(void)
{
} /* _Final_Module_pchutil */

