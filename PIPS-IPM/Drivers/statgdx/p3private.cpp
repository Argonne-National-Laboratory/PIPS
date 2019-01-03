#include "p3io.h"
#include "system_p3.h"
#include "p3private.h"


Procedure P3PRIVATE_pcharconcatpchar(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  SYSTEM_P3_pansichar psrc)
{
  SYSTEM_integer len;

  if (psrc != NULL) {
    /**** C code included from p3private.pas(43:1): 2 lines ****/
    len = strlen((char *) psrc);
    (void) strcat((char *) pdest+(*w), (char *) psrc);
    *w = *w + len;
  } 
}  /* pcharconcatpchar */

Procedure P3PRIVATE_pcharconcatstr(
  SYSTEM_P3_pansichar pdest,
  SYSTEM_integer *w,
  const SYSTEM_ansichar *src)
{
  SYSTEM_integer len;

  len = SYSTEM_length(src);
  if (len > 0) {
    /**** C code included from p3private.pas(65:1): 2 lines ****/
      memcpy((char *)pdest + *w, (char *)src+1, len);
      pdest[*w + len] = '\0';
    *w = *w + len;
  } 
}  /* pcharconcatstr */

Function(SYSTEM_P3_pansichar ) P3PRIVATE_strtopchar(
  const SYSTEM_ansichar *src)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer len;

  len = SYSTEM_length(src);
  _P3getmem(result,len + 1);
  /**** C code included from p3private.pas(85:1): 2 lines ****/
   memcpy((char *)result, (char *)src+1, len);
   result[len] = '\0';
  return result;
}  /* strtopchar */

Function(SYSTEM_P3_pansichar ) P3PRIVATE_strtostrbuf(
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

Function(SYSTEM_ansichar *) P3PRIVATE_strbuftostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *src)
{
  SYSTEM_integer len, i;

  len = 0;
  for (i = 0;i <= (SYSTEM_int32)254;++i) {
    if (src[i] != _P3char('\000')) {
      _P3inc0(len);
      result[len] = src[i];
    } else 
      SYSTEM_break(BRK_1);
  CNT_1:;
}
BRK_1:;
  _P3setlength(result,len,255);
  return result;
}  /* strbuftostr */

Function(SYSTEM_ansichar *) P3PRIVATE_pchartostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pansichar psrc)
{
  SYSTEM_integer len, i;

  len = 0;
  for (i = 0;i <= (SYSTEM_int32)254;++i) {
    if (*psrc != _P3char('\000')) {
      _P3inc0(len);
      result[len] = *psrc;
      _P3inc0(psrc);
    } else 
      SYSTEM_break(BRK_2);
  CNT_2:;
}
BRK_2:;
  _P3setlength(result,len,255);
  return result;
}  /* pchartostr */

/* unit p3private */
void _Init_Module_p3private(void)
{
} /* _Init_Module_p3private */

void _Final_Module_p3private(void)
{
} /* _Final_Module_p3private */

