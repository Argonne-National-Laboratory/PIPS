#include "p3io.h"
#include "system_p3.h"


static Function(SYSTEM_P3_pansichar ) SYSTEM_P3_getparamstr(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *param)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer len;
  SYSTEM_P3_pansichar start, s;

  while (SYSTEM_true) {
    while (*p != _P3char('\000') && *p <= _P3char(' ')) {

      _P3inc0(p);
}
    start = p;
    _P3inc0(start);
    if (*p == _P3char('\"') && *start == _P3char('\"')) { 
      _P3inc1(p,2);
    } else 
      SYSTEM_break(BRK_1);
  
CNT_1:;
  }
BRK_1:;
  len = 0;
  start = p;
  while (*p > _P3char(' ')) {

    if (*p == _P3char('\"')) {
      _P3inc0(p);
      while (*p != _P3char('\000') && *p != _P3char('\"')) {
        _P3inc0(p);
        _P3inc0(len);
      
}
      if (*p != _P3char('\000')) 
        _P3inc0(p);
    } else {
      _P3inc0(p);
      _P3inc0(len);
    } 
}
  _P3setlength(param,len,255);
  p = start;
  s = ValueCast(SYSTEM_P3_pansichar,&param[1]);
  while (*p > _P3char(' ')) {

    if (*p == _P3char('\"')) {
      _P3inc0(p);
      while (*p != _P3char('\000') && *p != _P3char('\"')) {
        *s = *p;
        _P3inc0(s);
        _P3inc0(p);
      
}
      if (*p != _P3char('\000')) 
        _P3inc0(p);
    } else {
      *s = *p;
      _P3inc0(s);
      _P3inc0(p);
    } 
}
  result = p;
  return result;
}  /* getparamstr */

static Function(SYSTEM_P3_pansichar ) SYSTEM_P3_getparamshortstr(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *param);

static Procedure pushchar(
  SYSTEM_ansichar c,
  SYSTEM_P3_pansichar *_2r,
  SYSTEM_integer *_2len)
{
  if (*_2len < 255) {
    **_2r = c;
    _P3inc0(*_2len);
    _P3inc0(*_2r);
  } 
}  /* pushchar */

static Function(SYSTEM_P3_pansichar ) SYSTEM_P3_getparamshortstr(
  SYSTEM_P3_pansichar p,
  SYSTEM_ansichar *param)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer len;
  SYSTEM_P3_pansichar s;
  SYSTEM_P3_pansichar r;

  while (SYSTEM_true) {
    while (*p != _P3char('\000') && *p <= _P3char(' ')) {

      _P3inc0(p);
}
    s = p;
    _P3inc0(s);
    if (*p == _P3char('\"') && *s == _P3char('\"')) { 
      _P3inc1(p,2);
    } else 
      SYSTEM_break(BRK_2);
  
CNT_2:;
  }
BRK_2:;
  len = 0;
  r = ValueCast(SYSTEM_P3_pansichar,&param[1]);
  while (*p > _P3char(' ')) {

    if (*p == _P3char('\"')) {
      _P3inc0(p);
      while (*p != _P3char('\000') && *p != _P3char('\"')) {
        pushchar(*p,&r,&len);
        _P3inc0(p);
      
}
      if (*p != _P3char('\000')) 
        _P3inc0(p);
    } else {
      pushchar(*p,&r,&len);
      _P3inc0(p);
    } 
}
  _P3setlength(param,len,255);
  result = p;
  return result;
}  /* getparamshortstr */

Function(SYSTEM_integer ) SYSTEM_P3_paramcount(void)
{
  SYSTEM_integer result;

  /**** C code included from system_p3.pas(184:1): 19 lines ****/
#if defined(_WIN32)
{
  char *cmdLine;
  SYSTEM_char *p;
  SYSTEM_shortstring s;

  result = 0;
  cmdLine = GetCommandLine(); /* we work in the system memory: don't write */
  p = SYSTEM_P3_getparamshortstr((SYSTEM_char *)cmdLine, s);
  for ( ; ; ) {
    p = SYSTEM_P3_getparamshortstr(p, s);
    if (0 == SYSTEM_length(s))
      break;
    result++;
  }
}
#else
  result = SYSTEM_paramcount();
#endif /* if defined(_WIN32) */
  return result;
}  /* paramcount */

Function(SYSTEM_ansichar *) SYSTEM_P3_paramstr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer index)
{
  SYSTEM_shortstring s;

  /**** C code included from system_p3.pas(211:1): 26 lines ****/
/*  SYSTEM_shortstring s; */
#if defined(_WIN32)
{
  SYSTEM_char *p;
  char buffer[261];   /* length take from Delphi's System.pas */
  DWORD fnLen;

  if (0 == index) {
    fnLen = GetModuleFileName (NULL, buffer, sizeof(buffer));
    (void) _P3pa2str(result, _len_ret, (SYSTEM_char *)buffer, fnLen);
    return result;
  }

  _P3strcpy(result,_len_ret,_P3str1("\000"));
  p = (SYSTEM_char *)GetCommandLine(); /* we work in the system memory: do not write */
  for ( ; ; ) {
    p = SYSTEM_P3_getparamshortstr(p, s);
    if (0 == index || 0 == SYSTEM_length(s))
      break;
    index--;
  }
  _P3strcpy(result,_len_ret, s);
}
#else
  _P3strcpy(result,_len_ret, SYSTEM_paramstr(s,sizeof(s)-1,index));
#endif /* if defined(_WIN32) */
  return result;
}  /* paramstr */

Procedure SYSTEM_P3_getdir(
  SYSTEM_byte d,
  SYSTEM_ansichar *s)
{
  /**** C code included from system_p3.pas(243:1): 39 lines ****/
#if defined(_WIN32)
{
  char Drive[3];
  char DirBuf[MAX_PATH], SaveBuf[MAX_PATH];
  int len;

  if (d) {
    Drive[0] = 'A' + d - 1;
    Drive[1] = ':';
    Drive[2] = '\0';
    GetCurrentDirectory(sizeof(SaveBuf), SaveBuf);
    SetCurrentDirectory(Drive);
  }
  GetCurrentDirectory(sizeof(DirBuf), DirBuf);
  if (d) SetCurrentDirectory(SaveBuf);
  len = strlen(DirBuf);
  if (len > 255) len = 255;
  memcpy((char *)s+1, DirBuf, len);
  s[0] = (SYSTEM_char) len;
}
#else
{
  char DirBuf[512];
  char *p;
  int len;

  p = getcwd(DirBuf, sizeof(DirBuf));
  if (NULL == p) {
    /* this is unlikely, and not handled well here: exceptions? */
    s[0] = (SYSTEM_char) 0;
  }
  else {
    len = strlen(DirBuf);
    if (len > 255) len = 255;
    memcpy((char *)s+1, DirBuf, len);
    s[0] = (SYSTEM_char) len;
  }
}
#endif
}  /* getdir */

Procedure SYSTEM_P3_fillchar(
  SYSTEM_untyped *p,
  SYSTEM_integer len,
  SYSTEM_byte v)
{
  /**** C code included from system_p3.pas(288:1): 1 lines ****/
  (void) memset(p, (int)v, (size_t)len);
}  /* fillchar */

Procedure SYSTEM_P3_settextbuf(
  SYSTEM_text *f,
  SYSTEM_untyped *p)
{
}  /* settextbuf */

/* unit system_p3 */
void _Init_Module_system_p3(void)
{
} /* _Init_Module_system_p3 */

void _Final_Module_system_p3(void)
{
} /* _Final_Module_system_p3 */

