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
#include "strutilx.h"

SYSTEM_shortstring STRUTILX_blanks255;
static GMSGEN_tcharset STRUTILX_whitespace = {255,255,255,255,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
static _P3STR_7 STRUTILX_maxint_s = {6,'m','a','x','i','n','t'};
static _P3STR_7 STRUTILX_minint_s = {6,'m','i','n','i','n','t'};
static SYSTEM_double STRUTILX_maxdouble_v = 1.0e299;
static _P3STR_15 STRUTILX_maxdouble_s = {9,'m','a','x','d','o','u','b','l','e'};
static SYSTEM_double STRUTILX_mindouble_v = -1.0e299;
static _P3STR_15 STRUTILX_mindouble_s = {9,'m','i','n','d','o','u','b','l','e'};
static SYSTEM_double STRUTILX_epsdouble_v = 1.0e-20;
static _P3STR_3 STRUTILX_epsdouble_s = {3,'e','p','s'};

Function(SYSTEM_P3_pshortstring ) STRUTILX_newstring(
  const SYSTEM_ansichar *s)
{
  SYSTEM_P3_pshortstring result;

  if (_P3strcmpE(s,_P3str1("\000"))) { 
    result = NULL;
  } else {
    _P3getmem(result,ValueCast(SYSTEM_int32,SYSTEM_length(s)) + 1);
    _P3strcpy(*result,255,s);
  } 
  return result;
}  /* newstring */

Procedure STRUTILX_disposestring(
  SYSTEM_P3_pshortstring ps)
{
  if (ps != NULL) 
    _P3freemem2(ps,ValueCast(SYSTEM_int32,SYSTEM_length(*ps)) + 1);
}  /* disposestring */

Function(SYSTEM_ansichar *) STRUTILX_getstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pshortstring ps)
{
  if (ps == NULL) { 
    _P3strclr(result);
  } else 
    _P3strcpy(result,_len_ret,*ps);
  return result;
}  /* getstring */

Function(SYSTEM_integer ) STRUTILX_strucmp(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_integer result;
  SYSTEM_byte k, l;

  l = SYSTEM_length(s1);
  if (l > SYSTEM_length(s2)) 
    l = SYSTEM_length(s2);
  { register SYSTEM_uint8 _stop = l;
    if ((k = 1) <=  _stop) do {
      result = SYSTEM_ord(SYSTEM_upcase(s1[k])) - SYSTEM_ord(
        SYSTEM_upcase(s2[k]));
      if (result != 0) 
        return result;
    
    } while (k++ !=  _stop);

  }
  result = ValueCast(SYSTEM_int32,SYSTEM_length(s1)) - SYSTEM_length(
    s2);
  return result;
}  /* strucmp */

Function(SYSTEM_integer ) STRUTILX_strucmpnum(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_integer result;
  SYSTEM_integer k, l, l1, l2;
  SYSTEM_shortstring ws1, ws2;

  l1 = SYSTEM_length(s1);
  while (l1 > 0 && _P3SET_in_3(s1[l1],_P3char('0'),_P3char('9'))) {

    l1 = l1 - 1;
}
  l2 = SYSTEM_length(s2);
  while (l2 > 0 && _P3SET_in_3(s2[l2],_P3char('0'),_P3char('9'))) {

    l2 = l2 - 1;
}
  {
    SYSTEM_shortstring _t1;
    SYSTEM_shortstring _t2;

    result = STRUTILX_strucmp(SYSTEM_copy(_t1,255,s1,1,l1),
      SYSTEM_copy(_t2,255,s2,1,l2));
  }
  if (result == 0) {
    SYSTEM_copy(ws1,255,s1,l1 + 1,255);
    SYSTEM_copy(ws2,255,s2,l2 + 1,255);
    l1 = SYSTEM_length(ws1);
    l2 = SYSTEM_length(ws2);
    if (l1 < l2) {
      l = l2;
      {
        SYSTEM_shortstring _t1;

        _P3strcat(ws1,255,SYSTEM_copy(_t1,255,STRUTILX_blanks255,1,
          l - l1),ws1);
      }
    } else 
      if (l1 == l2) { 
        l = l1;
      } else {
        l = l1;
        {
          SYSTEM_shortstring _t1;

          _P3strcat(ws2,255,SYSTEM_copy(_t1,255,
            STRUTILX_blanks255,1,l - l2),ws2);
        }
      } 
    { register SYSTEM_int32 _stop = l;
      if ((k = 1) <=  _stop) do {
        result = SYSTEM_ord(ws1[k]) - SYSTEM_ord(ws2[k]);
        if (result != 0) 
          SYSTEM_break(BRK_1);
      
CNT_1:;
      } while (k++ !=  _stop);
BRK_1:;

    }
  } 
  return result;
}  /* strucmpnum */

Function(SYSTEM_integer ) STRUTILX_pstrucmp(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2)
{
  SYSTEM_integer result;

  if (p1 != NULL && p2 != NULL) { 
    result = STRUTILX_strucmp(*p1,*p2);
  } else 
    result = SYSTEM_ord(p1 != NULL) - SYSTEM_ord(p2 != NULL);
  return result;
}  /* pstrucmp */

Procedure STRUTILX_strassign(
  SYSTEM_P3_pshortstring *p,
  const SYSTEM_ansichar *s)
{
  STRUTILX_disposestring(*p);
  *p = STRUTILX_newstring(s);
}  /* strassign */

Function(SYSTEM_integer ) STRUTILX_strcmp(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_integer result;
  SYSTEM_integer k, l;

  l = SYSTEM_length(s1);
  if (l > ValueCast(SYSTEM_int32,SYSTEM_length(s2))) 
    l = SYSTEM_length(s2);
  { register SYSTEM_int32 _stop = l;
    if ((k = 1) <=  _stop) do {
      result = SYSTEM_ord(s1[k]) - SYSTEM_ord(s2[k]);
      if (result != 0) 
        return result;
    
    } while (k++ !=  _stop);

  }
  result = ValueCast(SYSTEM_int32,SYSTEM_length(s1)) - SYSTEM_length(
    s2);
  return result;
}  /* strcmp */

Function(SYSTEM_integer ) STRUTILX_pstrcmp(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2)
{
  SYSTEM_integer result;

  if (p1 != NULL && p2 != NULL) { 
    result = STRUTILX_strcmp(*p1,*p2);
  } else 
    result = SYSTEM_ord(p1 != NULL) - SYSTEM_ord(p2 != NULL);
  return result;
}  /* pstrcmp */

Function(SYSTEM_boolean ) STRUTILX_struequal(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_boolean result;
  SYSTEM_integer k, l;

  result = SYSTEM_false;
  l = SYSTEM_length(s1);
  if (l != ValueCast(SYSTEM_int32,SYSTEM_length(s2))) 
    return result;
  for (k = l;k >= (SYSTEM_int32)1;--k) {
    if (SYSTEM_upcase(s1[k]) != SYSTEM_upcase(s2[k])) 
      return result;
  }
  result = SYSTEM_true;
  return result;
}  /* struequal */

Function(SYSTEM_boolean ) STRUTILX_strequal(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_boolean result;
  SYSTEM_integer k, l;

  result = SYSTEM_false;
  l = SYSTEM_length(s1);
  if (l != ValueCast(SYSTEM_int32,SYSTEM_length(s2))) 
    return result;
  for (k = l;k >= (SYSTEM_int32)1;--k) {
    if (s1[k] != s2[k]) 
      return result;
  }
  result = SYSTEM_true;
  return result;
}  /* strequal */

Function(SYSTEM_boolean ) STRUTILX_pstruequal(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2)
{
  SYSTEM_boolean result;
  SYSTEM_integer k, l;

  if (p1 == NULL || p2 == NULL) { 
    result = p1 == NULL && p2 == NULL;
  } else {
    result = SYSTEM_false;
    l = SYSTEM_ord((*p1)[0]);
    if (l != SYSTEM_ord((*p2)[0])) 
      return result;
    for (k = l;k >= (SYSTEM_int32)1;--k) {
      if (SYSTEM_upcase((*p1)[k]) != SYSTEM_upcase((*p2)[k])) 
        return result;
    }
    result = SYSTEM_true;
  } 
  return result;
}  /* pstruequal */

Function(SYSTEM_boolean ) STRUTILX_pstrequal(
  SYSTEM_P3_pshortstring p1,
  SYSTEM_P3_pshortstring p2)
{
  SYSTEM_boolean result;
  SYSTEM_integer k, l;

  if (p1 == NULL || p2 == NULL) { 
    result = p1 == NULL && p2 == NULL;
  } else {
    result = SYSTEM_false;
    l = SYSTEM_ord((*p1)[0]);
    if (l != SYSTEM_ord((*p2)[0])) 
      return result;
    for (k = l;k >= (SYSTEM_int32)1;--k) {
      if ((*p1)[k] != (*p2)[k]) 
        return result;
    }
    result = SYSTEM_true;
  } 
  return result;
}  /* pstrequal */

Function(SYSTEM_ansichar *) STRUTILX_inttostrw(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n,
  SYSTEM_integer width)
{
  SYSTEM_shortstring s;
  SYSTEM_boolean neg;
  SYSTEM_byte w, w2;
  SYSTEM_boolean dirty;

  neg = n < 0;
  if (neg) 
    n = -n;
  dirty = n < 0;
  if (dirty) 
    n = SYSTEM_maxint;
  w = 255;
  do {
    s[w] = ValueCast(SYSTEM_ansichar,n % 10 + 48);
    _P3dec0(w);
    n = n /  10;
  } while (!(n == 0));
  if (neg) {
    s[w] = _P3char('-');
    _P3dec0(w);
  } 
  if (dirty) 
    _P3inc0(s[255]);
  if (width > 255) 
    width = 255;
  n = width - (255 - w);
  w2 = 0;
  while (n > 0) {
    _P3inc0(w2);
    result[w2] = _P3char(' ');
    _P3dec0(n);
  
}
  while (w < 255) {
    _P3inc0(w2);
    _P3inc0(w);
    result[w2] = s[w];
  
}
  _P3setlength(result,w2,255);
  return result;
}  /* inttostrw */

Function(SYSTEM_ansichar *) STRUTILX_inttostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  STRUTILX_inttostrw(result,_len_ret,n,0);
  return result;
}  /* inttostr */

Function(SYSTEM_ansichar *) STRUTILX_inttostrex(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  if (n == SYSTEM_maxint) { 
    _P3strcpy(result,_len_ret,STRUTILX_maxint_s);
  } else 
    if (n ==  SYSTEM_minint) { 
      _P3strcpy(result,_len_ret,STRUTILX_minint_s);
    } else 
      STRUTILX_inttostrw(result,_len_ret,n,0);
  return result;
}  /* inttostrex */

Function(SYSTEM_ansichar *) STRUTILX_inttonicestrw(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n,
  SYSTEM_integer width)
{
  SYSTEM_shortstring s;
  SYSTEM_boolean neg;
  SYSTEM_byte w, w2;

  neg = n < 0;
  if (neg) 
    n = -n;
  w = 255;
  w2 = 0;
  do {
    s[w] = ValueCast(SYSTEM_ansichar,n % 10 + 48);
    _P3dec0(w);
    n = n /  10;
    _P3inc0(w2);
    if (w2 == 3) {
      if (n != 0) {
        s[w] = _P3char(',');
        _P3dec0(w);
      } 
      w2 = 0;
    } 
  } while (!(n == 0));
  if (neg) {
    s[w] = _P3char('-');
    _P3dec0(w);
  } 
  if (width > 255) 
    width = 255;
  n = width - (255 - w);
  w2 = 0;
  while (n > 0) {
    _P3inc0(w2);
    result[w2] = _P3char(' ');
    _P3dec0(n);
  
}
  while (w < 255) {
    _P3inc0(w2);
    _P3inc0(w);
    result[w2] = s[w];
  
}
  _P3setlength(result,w2,255);
  return result;
}  /* inttonicestrw */

Function(SYSTEM_ansichar *) STRUTILX_int64tonicestrw(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n,
  SYSTEM_integer width)
{
  SYSTEM_shortstring s;
  SYSTEM_boolean neg;
  SYSTEM_byte w, w2;

  neg = n < 0;
  if (neg) 
    n = -n;
  w = 255;
  w2 = 0;
  do {
    s[w] = ValueCast(SYSTEM_ansichar,n % 10 + 48);
    _P3dec0(w);
    n = n /  10;
    _P3inc0(w2);
    if (w2 == 3) {
      if (n != 0) {
        s[w] = _P3char(',');
        _P3dec0(w);
      } 
      w2 = 0;
    } 
  } while (!(n == 0));
  if (neg) {
    s[w] = _P3char('-');
    _P3dec0(w);
  } 
  if (width > 255) 
    width = 255;
  n = width - (255 - w);
  w2 = 0;
  while (n > 0) {
    _P3inc0(w2);
    result[w2] = _P3char(' ');
    _P3dec0(n);
  
}
  while (w < 255) {
    _P3inc0(w2);
    _P3inc0(w);
    result[w2] = s[w];
  
}
  _P3setlength(result,w2,255);
  return result;
}  /* int64tonicestrw */

Function(SYSTEM_ansichar *) STRUTILX_inttonicestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n)
{
  STRUTILX_inttonicestrw(result,_len_ret,n,0);
  return result;
}  /* inttonicestr */

Function(SYSTEM_ansichar *) STRUTILX_int64tonicestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n)
{
  STRUTILX_int64tonicestrw(result,_len_ret,n,0);
  return result;
}  /* int64tonicestr */

Function(SYSTEM_ansichar *) STRUTILX_dbltostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v)
{
  SYSTEM_shortstring s;
  SYSTEM_integer i, j, k, e;

  if (v == 0.0) {
    _P3strcpy(result,_len_ret,_P3str1("\0010"));
    return result;
  } 
  _P3str_d0(v,s,255);
  if (v < 0.0) 
    v = -v;
  k = STRUTILX_rchsetpos(_P3set1("\0\0\0\0\0\050\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"),
    s);
  j = STRUTILX_lchpos(_P3char('.'),s);
  if (v >= 1e-4 && v < 1e15) {
    {
      SYSTEM_shortstring _t1;

      _P3val_i(SYSTEM_copy(_t1,255,s,k,5),e,&i);
    }
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((i = k - 1) <=  _stop) do {
        s[i] = _P3char('0');
      } while (i++ !=  _stop);

    }
    if (e >= 0) {
      { register SYSTEM_int32 _stop = j + e;
        if ((i = j + 1) <=  _stop) do {
          s[i - 1] = s[i];
        } while (i++ !=  _stop);

      }
      s[j + e] = _P3char('.');
      { register SYSTEM_int32 _stop = j + e + 1;
        if ((i = SYSTEM_length(s)) >=  _stop) do {
          if (s[i] == _P3char('0')) {
            s[i] = _P3char(' ');
            if (i == j + e + 1) 
              s[j + e] = _P3char(' ');
          } else 
            SYSTEM_break(BRK_2);
CNT_2:;
        } while (i-- !=  _stop);
BRK_2:;

      }
    } else {
      s[j] = s[j - 1];
      s[j - 1] = _P3char('0');
      e = -e;
      { register SYSTEM_int32 _stop = j;
        if ((i = k - 2) >=  _stop) do {
          s[i + e] = s[i];
        } while (i-- !=  _stop);

      }
      { register SYSTEM_int32 _stop = j + e - 1;
        if ((i = j + 1) <=  _stop) do {
          s[i] = _P3char('0');
        } while (i++ !=  _stop);

      }
      s[j] = _P3char('.');
      _P3setlength(s,k + e - 2,255);
      { register SYSTEM_int32 _stop = j + e + 1;
        if ((i = SYSTEM_length(s)) >=  _stop) do {
          if (s[i] == _P3char('0')) { 
            s[i] = _P3char(' ');
          } else 
            SYSTEM_break(BRK_3);
CNT_3:;
        } while (i-- !=  _stop);
BRK_3:;

      }
    } 
  } else {
    if (s[k] == _P3char('+')) 
      s[k] = _P3char(' ');
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((i = k + 1) <=  _stop) do {
        if (s[i] == _P3char('0')) {
          s[i] = _P3char(' ');
          if (i == ValueCast(SYSTEM_int32,SYSTEM_length(s))) 
            s[k - 1] = _P3char(' ');
        } else 
          SYSTEM_break(BRK_4);
CNT_4:;
      } while (i++ !=  _stop);
BRK_4:;

    }
    { register SYSTEM_int32 _stop = j + 1;
      if ((i = k - 2) >=  _stop) do {
        if (s[i] == _P3char('0')) {
          s[i] = _P3char(' ');
          if (i == j + 1) 
            s[j] = _P3char(' ');
        } else 
          SYSTEM_break(BRK_5);
CNT_5:;
      } while (i-- !=  _stop);
BRK_5:;

    }
  } 
  j = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      if (s[i] != _P3char(' ')) {
        j = j + 1;
        result[j] = s[i];
      } 
    } while (i++ !=  _stop);

  }
  _P3setlength(result,j,255);
  return result;
}  /* dbltostr */

Function(SYSTEM_ansichar *) STRUTILX_dbltostrex(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v)
{
  if (v == STRUTILX_maxdouble_v) { 
    _P3strcpy(result,_len_ret,STRUTILX_maxdouble_s);
  } else 
    if (v == STRUTILX_mindouble_v) { 
      _P3strcpy(result,_len_ret,STRUTILX_mindouble_s);
    } else 
      STRUTILX_dbltostr(result,_len_ret,v);
  return result;
}  /* dbltostrex */

Function(SYSTEM_ansichar *) STRUTILX_fillstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_ansichar ch,
  SYSTEM_integer len)
{
  SYSTEM_integer k;

  if (len < 0) { 
    len = 0;
  } else 
    if (len > 255) 
      len = 255;
  _P3setlength(result,len,255);
  { register SYSTEM_int32 _stop = SYSTEM_length(result);
    if ((k = 1) <=  _stop) do {
      result[k] = ch;
    } while (k++ !=  _stop);

  }
  return result;
}  /* fillstr */

Function(SYSTEM_ansichar *) STRUTILX_padleft(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer w)
{
  SYSTEM_integer ww;

  ww = w - SYSTEM_length(s);
  if (ww <= 0) { 
    _P3strcpy(result,_len_ret,s);
  } else {
    if (ww + SYSTEM_length(s) > 255) 
      ww = 255 - SYSTEM_length(s);
    {
      SYSTEM_shortstring _t1;

      _P3strcat(result,_len_ret,SYSTEM_copy(_t1,255,
        STRUTILX_blanks255,1,ww),s);
    }
  } 
  return result;
}  /* padleft */

Function(SYSTEM_ansichar *) STRUTILX_padright(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer w)
{
  SYSTEM_integer ww;

  ww = w - SYSTEM_length(s);
  if (ww <= 0) { 
    _P3strcpy(result,_len_ret,s);
  } else 
    {
      SYSTEM_shortstring _t1;

      _P3strcat(result,_len_ret,s,SYSTEM_copy(_t1,255,
        STRUTILX_blanks255,1,ww));
    }
  return result;
}  /* padright */

Function(SYSTEM_integer ) STRUTILX_padmodlength(
  const SYSTEM_ansichar *s,
  SYSTEM_integer m)
{
  SYSTEM_integer result;

  result = SYSTEM_length(s);
  if (m > 0 && result % m != 0) 
    result = result + (m - result % m);
  return result;
}  /* padmodlength */

Function(SYSTEM_ansichar *) STRUTILX_padrightmod(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer m)
{
  {
    SYSTEM_shortstring _t1;

    _P3strcat(result,_len_ret,s,STRUTILX_blankstr(_t1,255,
      STRUTILX_padmodlength(s,m) - SYSTEM_length(s)));
  }
  return result;
}  /* padrightmod */

Function(SYSTEM_integer ) STRUTILX_lchpossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp)
{
  SYSTEM_integer result;
  SYSTEM_integer k;

  result = 0;
  if (sp <= 0) 
    sp = 1;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = sp) <=  _stop) do {
      if (s[k] == ch) {
        result = k;
        SYSTEM_break(BRK_6);
      } 
CNT_6:;
    } while (k++ !=  _stop);
BRK_6:;

  }
  return result;
}  /* lchpossp */

Function(SYSTEM_integer ) STRUTILX_lchpos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = STRUTILX_lchpossp(ch,s,1);
  return result;
}  /* lchpos */

Function(SYSTEM_integer ) STRUTILX_rchpossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp)
{
  SYSTEM_integer result;
  SYSTEM_integer l;
  SYSTEM_integer k;

  result = 0;
  l = SYSTEM_length(s);
  if (sp > l) 
    sp = l;
  for (k = sp;k >= (SYSTEM_int32)1;--k) {
    if (s[k] == ch) {
      result = k;
      SYSTEM_break(BRK_7);
    } 
  CNT_7:;
}
BRK_7:;
  return result;
}  /* rchpossp */

Function(SYSTEM_integer ) STRUTILX_rchpos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = STRUTILX_rchpossp(ch,s,SYSTEM_length(s));
  return result;
}  /* rchpos */

Function(SYSTEM_integer ) STRUTILX_lchupossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp)
{
  SYSTEM_integer result;
  SYSTEM_integer k;

  result = 0;
  ch = SYSTEM_upcase(ch);
  if (sp <= 0) 
    sp = 1;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = sp) <=  _stop) do {
      if (SYSTEM_upcase(s[k]) == ch) {
        result = k;
        SYSTEM_break(BRK_8);
      } 
CNT_8:;
    } while (k++ !=  _stop);
BRK_8:;

  }
  return result;
}  /* lchupossp */

Function(SYSTEM_integer ) STRUTILX_lchupos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = STRUTILX_lchupossp(ch,s,1);
  return result;
}  /* lchupos */

Function(SYSTEM_integer ) STRUTILX_rchupossp(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp)
{
  SYSTEM_integer result;
  SYSTEM_integer l;
  SYSTEM_integer k;

  result = 0;
  ch = SYSTEM_upcase(ch);
  l = SYSTEM_length(s);
  if (sp > l) 
    sp = l;
  for (k = sp;k >= (SYSTEM_int32)1;--k) {
    if (SYSTEM_upcase(s[k]) == ch) {
      result = k;
      SYSTEM_break(BRK_9);
    } 
  CNT_9:;
}
BRK_9:;
  return result;
}  /* rchupossp */

Function(SYSTEM_integer ) STRUTILX_rchupos(
  SYSTEM_ansichar ch,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = STRUTILX_rchupossp(ch,s,SYSTEM_length(s));
  return result;
}  /* rchupos */

Function(SYSTEM_integer ) STRUTILX_lchsetpos(
  const _P3set_elem *cs,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer k;

  result = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = 1) <=  _stop) do {
      if (_P3SET_i(255,s[k],cs)) {
        result = k;
        SYSTEM_break(BRK_10);
      } 
CNT_10:;
    } while (k++ !=  _stop);
BRK_10:;

  }
  return result;
}  /* lchsetpos */

Function(SYSTEM_integer ) STRUTILX_rchsetpos(
  const _P3set_elem *cs,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer k;

  result = 0;
  for (k = SYSTEM_length(s);k >= (SYSTEM_int32)1;--k) {
    if (_P3SET_i(255,s[k],cs)) {
      result = k;
      SYSTEM_break(BRK_11);
    } 
  CNT_11:;
}
BRK_11:;
  return result;
}  /* rchsetpos */

Function(SYSTEM_integer ) STRUTILX_lstrpossp(
  const SYSTEM_ansichar *pat,
  const SYSTEM_ansichar *s,
  SYSTEM_integer sp)
{
  SYSTEM_integer result;
  SYSTEM_integer p;
  SYSTEM_integer lp;
  SYSTEM_integer k;
  SYSTEM_integer ls;
  SYSTEM_ansichar pat1;

  result = 0;
  lp = SYSTEM_length(pat);
  ls = SYSTEM_length(s);
  if (lp == 0 || ls == 0 || sp <= 0 || sp + lp - 1 > 
    ls) 
    return result;
  pat1 = pat[1];
  if (lp == 1) { 
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((p = sp) <=  _stop) do {
        if (s[p] == pat1) {
          result = p;
          return result;
        } 
      } while (p++ !=  _stop);

    }
  } else 
    { register SYSTEM_int32 _stop = SYSTEM_length(s) - lp + 1;
      if ((p = sp) <=  _stop) do {
        if (s[p] != pat1) 
          SYSTEM_continue(CNT_12);
        result = p;
        { register SYSTEM_int32 _stop = lp;
          if ((k = 2) <=  _stop) do {
            if (pat[k] != s[p + k - 1]) {
              result = 0;
              SYSTEM_break(BRK_13);
            } 
CNT_13:;
          } while (k++ !=  _stop);
BRK_13:;

        }
        if (result != 0) 
          return result;
      
CNT_12:;
      } while (p++ !=  _stop);
BRK_12:;

    }
  return result;
}  /* lstrpossp */

Function(SYSTEM_integer ) STRUTILX_lstrpos(
  const SYSTEM_ansichar *pat,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = STRUTILX_lstrpossp(pat,s,1);
  return result;
}  /* lstrpos */

Function(SYSTEM_ansichar *) STRUTILX_replacechar(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const _P3set_elem *chset,
  SYSTEM_ansichar _new,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer k;

  _P3setlength(result,SYSTEM_length(s),255);
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = 1) <=  _stop) do {
      if (_P3SET_i(255,s[k],chset)) { 
        result[k] = _new;
      } else 
        result[k] = s[k];
    } while (k++ !=  _stop);

  }
  return result;
}  /* replacechar */

Function(SYSTEM_ansichar *) STRUTILX_deletechar(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const _P3set_elem *chset,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer k;
  SYSTEM_integer w;

  w = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = 1) <=  _stop) do {
      if (!_P3SET_i(255,s[k],chset)) {
        w = w + 1;
        result[w] = s[k];
      } 
    } while (k++ !=  _stop);

  }
  _P3setlength(result,w,255);
  return result;
}  /* deletechar */

Function(SYSTEM_ansichar *) STRUTILX_replacestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *old,
  const SYSTEM_ansichar *_new,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer i;
  SYSTEM_integer r;

  if (_P3strcmpE(old,_new) || _P3strcmpE(old,_P3str1("\000"))) { 
    _P3strcpy(result,_len_ret,s);
  } else {
    _P3strclr(result);
    r = 1;
    do {
      i = STRUTILX_lstrpossp(old,s,r);
      if (i <= 0) 
        SYSTEM_break(BRK_14);
      {
        SYSTEM_shortstring _t1;
        _P3STR_255 _t2;

        _P3strcat(result,_len_ret,_P3strcat(_t2,255,result,
          SYSTEM_copy(_t1,255,s,r,i - r)),_new);
      }
      r = i + SYSTEM_length(old);
    CNT_14:;
    } while (SYSTEM_true);
BRK_14:;
    {
      SYSTEM_shortstring _t1;

      _P3strcat(result,_len_ret,result,SYSTEM_copy(_t1,255,s,r,255));
    }
  } 
  return result;
}  /* replacestr */

Function(SYSTEM_ansichar *) STRUTILX_uppercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer k;

  _P3setlength(result,SYSTEM_length(s),255);
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = 1) <=  _stop) do {
      result[k] = SYSTEM_upcase(s[k]);
    } while (k++ !=  _stop);

  }
  return result;
}  /* uppercase */

Function(SYSTEM_ansichar ) STRUTILX_lowcase(
  SYSTEM_ansichar ch)
{
  SYSTEM_ansichar result;

  if (ch >= _P3char('A') && ch <= _P3char('Z')) { 
    result = ValueCast(SYSTEM_ansichar,SYSTEM_ord(ch) - 65 + 97);
  } else 
    result = ch;
  return result;
}  /* lowcase */

Function(SYSTEM_ansichar *) STRUTILX_lowercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer k;

  _P3setlength(result,SYSTEM_length(s),255);
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((k = 1) <=  _stop) do {
      result[k] = STRUTILX_lowcase(s[k]);
    } while (k++ !=  _stop);

  }
  return result;
}  /* lowercase */

Function(SYSTEM_integer ) STRUTILX_integerwidth(
  SYSTEM_integer n)
{
  SYSTEM_integer result;

  if (n >= 0) { 
    result = 0;
  } else {
    result = 1;
    n = -n;
  } 
  do {
    result = result + 1;
    n = n /  10;
  } while (!(n == 0));
  return result;
}  /* integerwidth */

Function(SYSTEM_ansichar *) STRUTILX_blankstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer len)
{
  SYSTEM_integer k;

  if (len <= 0) { 
    _P3strclr(result);
  } else {
    if (len > 255) 
      len = 255;
    _P3setlength(result,len,255);
    { register SYSTEM_int32 _stop = len;
      if ((k = 1) <=  _stop) do {
        result[k] = _P3char(' ');
      } while (k++ !=  _stop);

    }
  } 
  return result;
  SYSTEM_copy(result,_len_ret,STRUTILX_blanks255,1,len);
  return result;
}  /* blankstr */

Function(SYSTEM_ansichar *) STRUTILX_extracttoken(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s,
  SYSTEM_integer *p)
{
  SYSTEM_integer l;
  SYSTEM_ansichar stop;
  SYSTEM_integer rs;

  _P3strclr(result);
  if (*p <= 0) 
    return result;
  l = SYSTEM_length(s);
  while (*p <= l && s[*p] == _P3char(' ')) {

    *p = *p + 1;
}
  if (*p > l) 
    return result;
  if (!_P3SET_in_1(s[*p],_P3char('\''),_P3SET_equal(s[*p],_P3char('\"')))) { 
    stop = _P3char(' ');
  } else {
    stop = s[*p];
    *p = *p + 1;
  } 
  rs = *p;
  while (*p <= l && s[*p] != stop) {

    *p = *p + 1;
}
  SYSTEM_copy(result,_len_ret,s,rs,*p - rs);
  if (*p <= l && s[*p] == stop) 
    *p = *p + 1;
  return result;
}  /* extracttoken */

Function(SYSTEM_boolean ) STRUTILX_strasintex(
  const SYSTEM_ansichar *s,
  SYSTEM_integer *v)
{
  SYSTEM_boolean result;
  SYSTEM_integer k;

  result = SYSTEM_true;
  if (STRUTILX_struequal(s,STRUTILX_maxint_s)) { 
    *v = 2147483647;
  } else 
    if (STRUTILX_struequal(s,STRUTILX_minint_s)) { 
      *v =  SYSTEM_minint;
    } else 
      _P3_TRY {
        _P3val_i(s,*v,&k);
        result = k == 0;
      } _P3_EXCEPT {
{
          *v = 0;
          result = SYSTEM_false;
        } 
      } _P3_END_TRY_EXCEPT;
  return result;
}  /* strasintex */

Function(SYSTEM_boolean ) STRUTILX_strasintex2(
  const SYSTEM_ansichar *s,
  SYSTEM_integer *v)
{
  SYSTEM_boolean result;
  SYSTEM_double d;

  result = STRUTILX_strasintex(s,v);
  if (!result) {
    *v = 0;
    result = STRUTILX_strasdoubleex(s,&d);
    if (result) {
      result = d >=  SYSTEM_minint && d <= 2147483647 && SYSTEM_frac(d) == 0.0;
      if (result) 
        *v = SYSTEM_trunc(d);
    } 
  } 
  return result;
}  /* strasintex2 */

Function(SYSTEM_integer ) STRUTILX_strasint(
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer k;

  _P3val_i(s,result,&k);
  if (k != 0) 
    result = 0;
  return result;
}  /* strasint */

Function(SYSTEM_boolean ) STRUTILX_strasdoubleex(
  const SYSTEM_ansichar *s,
  SYSTEM_double *v)
{
  SYSTEM_boolean result;
  SYSTEM_integer k;
  SYSTEM_shortstring ws;

  result = SYSTEM_true;
  if (STRUTILX_struequal(s,STRUTILX_maxdouble_s)) { 
    *v = STRUTILX_maxdouble_v;
  } else 
    if (STRUTILX_struequal(s,STRUTILX_mindouble_s)) { 
      *v = STRUTILX_mindouble_v;
    } else 
      if (STRUTILX_struequal(s,STRUTILX_epsdouble_s)) { 
        *v = STRUTILX_epsdouble_v;
      } else {
        STRUTILX_replacechar(ws,255,_P3set1("\0\0\0\0\0\0\0\0\020\0\0\0\020\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0\0"),_P3char('E'),
          s);
        _P3_TRY {
          _P3val_d(ws,*v,&k);
          if (_P3SET_i(3,P3IEEEFP_fpclass(*v),_P3set1("\017"))) { 
            result = SYSTEM_false;
          } else 
            result = k == 0;
        } _P3_EXCEPT {

            result = SYSTEM_false;
        } _P3_END_TRY_EXCEPT;
} 
  return result;
}  /* strasdoubleex */

Function(SYSTEM_ansichar *) STRUTILX_gmsvaltostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v)
{
  if (SYSTEM_abs_r(v) <= GMSSPECS_valtiny) {
    _P3strcpy(result,_len_ret,_P3str1("\0010"));
    return result;
  } 
  if (v <= GMSSPECS_valbig && v >= GMSSPECS_valsmall) {
    SYSUTILS_P3_floattostr(result,_len_ret,v);
    return result;
  } 
  if (v < GMSSPECS_valsmall) { 
    _P3strcpy(result,_len_ret,_P3str1("\004-Inf"));
  } else 
    if (v == GMSSPECS_valpin) { 
      _P3strcpy(result,_len_ret,_P3str1("\004+Inf"));
    } else 
      if (v == GMSSPECS_valmin) { 
        _P3strcpy(result,_len_ret,_P3str1("\004-Inf"));
      } else 
        if (v == GMSSPECS_valeps) { 
          _P3strcpy(result,_len_ret,_P3str1("\003Eps"));
        } else 
          if (v == GMSSPECS_valna) { 
            _P3strcpy(result,_len_ret,_P3str1("\002Na"));
          } else 
            if (v == GMSSPECS_valund) { 
              _P3strcpy(result,_len_ret,_P3str1("\004Undf"));
            } else 
              _P3strcpy(result,_len_ret,_P3str1("\004+Inf"));
  return result;
}  /* gmsvaltostr */

Function(SYSTEM_ansichar *) STRUTILX_mem64tonicestr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n,
  SYSTEM_integer w)
{
  SYSTEM_integer d;
  SYSTEM_shortstring s;

  if (n < 16384) {
    d = 1;
    _P3strcpy(s,255,_P3str1("\002 b"));
  } else 
    if (n < 16777216) {
      d = 1024;
      _P3strcpy(s,255,_P3str1("\002Kb"));
    } else {
      d = 1048576;
      _P3strcpy(s,255,_P3str1("\002Mb"));
    } 
  {
    SYSTEM_shortstring _t1;
    SYSTEM_shortstring _t2;
    _P3STR_255 _t3;

    _P3strcat(result,_len_ret,_P3strcat(_t3,255,STRUTILX_padleft(
      _t1,255,STRUTILX_int64tonicestr(_t2,255,(n + d /  2) /  
      d),w - 3),_P3str1("\001 ")),s);
  }
  return result;
}  /* mem64tonicestr */

Function(SYSTEM_ansichar *) STRUTILX_excelcolstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer c)
{
  _P3strclr(result);
  if (c <= 0) 
    return result;
  do {
    c = c - 1;
    {
      _P3STR_3 _t1;

      _P3strcat(result,_len_ret,_P3ch2str(_t1,1,ValueCast(
        SYSTEM_ansichar,65 + c % 26)),result);
    }
    c = c /  26;
  } while (!(c == 0));
  return result;
}  /* excelcolstr */

Function(SYSTEM_integer ) STRUTILX_strexcelcol(
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer i, j;

  result = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      j = SYSTEM_ord(SYSTEM_upcase(s[i])) - 65;
      if (j < 0 || j > 25 || result >= 82595550) {
        result = 0;
        return result;
      } 
      result = result * 26 + j + 1;
    
    } while (i++ !=  _stop);

  }
  return result;
}  /* strexcelcol */

Function(SYSTEM_P3_pshortstring ) STRUTILX_newstringm(
  const SYSTEM_ansichar *s,
  SYSTEM_int64 *m)
{
  SYSTEM_P3_pshortstring result;
  SYSTEM_integer l;

  if (_P3strcmpE(s,_P3str1("\000"))) { 
    result = NULL;
  } else {
    l = ValueCast(SYSTEM_int32,SYSTEM_length(s)) + 1;
    _P3getmem(result,l);
    _P3strcpy(*result,255,s);
    *m = *m + l;
  } 
  return result;
}  /* newstringm */

Procedure STRUTILX_disposestringm(
  SYSTEM_P3_pshortstring ps,
  SYSTEM_int64 *m)
{
  SYSTEM_integer l;

  if (ps != NULL) {
    l = ValueCast(SYSTEM_int32,SYSTEM_length(*ps)) + 1;
    _P3freemem2(ps,l);
    *m = *m - l;
  } 
}  /* disposestringm */

Procedure STRUTILX_strassignm(
  SYSTEM_P3_pshortstring *p,
  const SYSTEM_ansichar *s,
  SYSTEM_int64 *m)
{
  STRUTILX_disposestringm(*p,m);
  *p = STRUTILX_newstringm(s,m);
}  /* strassignm */

/* unit strutilx */
void _Init_Module_strutilx(void)
{
  STRUTILX_fillstr(STRUTILX_blanks255,255,_P3char(' '),255);
} /* _Init_Module_strutilx */

void _Final_Module_strutilx(void)
{
} /* _Final_Module_strutilx */

