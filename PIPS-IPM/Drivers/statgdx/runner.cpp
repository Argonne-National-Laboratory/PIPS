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
#include "gmsobj.h"
#include "runner.h"


void * const RUNNER_tmsghandler_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'tmsghandler' */
const SYSTEM_classdescriptor_t RUNNER_tmsghandler_CD = {
  _P3str1("\013tmsghandler"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(RUNNER_tmsghandler_OD), RUNNER_tmsghandler_VT, NULL};


void * const RUNNER_trunner_VT[] = {(void*)&RUNNER_trunner_DOT_destroy};

/* Class descriptor for 'trunner' */
const SYSTEM_classdescriptor_t RUNNER_trunner_CD = {
  _P3str1("\007trunner"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(RUNNER_trunner_OD), RUNNER_trunner_VT, NULL};


Procedure RUNNER_trunner_DOT_paramsadd(
  RUNNER_trunner self,
  const SYSTEM_ansichar *v)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\011ParamsAdd"))) 
    return;
  GMSOBJ_txstrings_DOT_add(self->RUNNER_trunner_DOT_fparams,v);
  RUNNER_trunner_DOT_commandlinechanged(self);
}  /* paramsadd */

Procedure RUNNER_trunner_DOT_paramsclear(
  RUNNER_trunner self)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\013ParamsClear"))) 
    return;
  GMSOBJ_txlist_DOT_clear(ValueCast(GMSOBJ_txlist,self->
    RUNNER_trunner_DOT_fparams));
  RUNNER_trunner_DOT_commandlinechanged(self);
}  /* paramsclear */

static Procedure addstrtop(
  const SYSTEM_ansichar *s,
  SYSTEM_P3_pansichar *_2p,
  SYSTEM_integer *_2l,
  SYSTEM_integer *_2f,
  RUNNER_trunner *_2self)
{
  if (*_2f == 1) { 
    *_2l = *_2l + SYSTEM_length(s);
  } else {
    SYSTEM_move(&s[1],*_2p,SYSTEM_length(s));
    _P3inc1(*_2p,SYSTEM_length(s));
  } 
}  /* addstrtop */

Function(SYSTEM_P3_pansichar ) RUNNER_trunner_DOT_commandline(
  RUNNER_trunner self)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_integer l;
  SYSTEM_P3_pansichar p;
  SYSTEM_integer f;
  SYSTEM_ansichar q;
  SYSTEM_integer n;

  if (self->RUNNER_trunner_DOT_fcommandline != NULL) {
    result = self->RUNNER_trunner_DOT_fcommandline;
    return result;
  } 
  q = _P3char('\"');
  for (f = 1;f <= (SYSTEM_int32)2;++f) {
    if (f == 1) { 
      l = 0;
    } else {
      _P3getmem(p,l);
      result = p;
    } 
    if (SYSTEM_pos(_P3str1("\001 "),self->
      RUNNER_trunner_DOT_fexecutable) == 0) { 
      addstrtop(self->RUNNER_trunner_DOT_fexecutable,&p,&l,&f,&self);
    } else 
      {
        _P3STR_3 _t1;
        _P3STR_255 _t2;
        _P3STR_3 _t3;
        _P3STR_255 _t4;

        addstrtop(_P3strcat(_t4,255,_P3strcat(_t2,255,_P3ch2str(
          _t1,1,q),self->RUNNER_trunner_DOT_fexecutable),
          _P3ch2str(_t3,1,q)),&p,&l,&f,&self);
      }
    { register SYSTEM_int32 _stop = RUNNER_trunner_DOT_paramscount(
        self) - 1;
      if ((n = 0) <=  _stop) do {
        {
          SYSTEM_shortstring _t1;

          if (SYSTEM_pos(_P3str1("\001 "),GMSOBJ_txstrings_DOT_get(_t1,255,
            self->RUNNER_trunner_DOT_fparams,n)) == 0) { 
            {
              SYSTEM_shortstring _t1;
              _P3STR_255 _t2;

              addstrtop(_P3strcat(_t2,255,_P3str1("\001 "),
                GMSOBJ_txstrings_DOT_get(_t1,255,self->
                RUNNER_trunner_DOT_fparams,n)),&p,&l,&f,&self);
            }
          } else 
            {
              _P3STR_3 _t1;
              _P3STR_3 _t2;
              SYSTEM_shortstring _t3;
              _P3STR_255 _t4;
              _P3STR_3 _t5;
              _P3STR_255 _t6;

              addstrtop(_P3strcat(_t6,255,_P3strcat(_t4,255,
                _P3strcat(_t2,2,_P3str1("\001 "),_P3ch2str(_t1,1,
                q)),GMSOBJ_txstrings_DOT_get(_t3,255,self->
                RUNNER_trunner_DOT_fparams,n)),_P3ch2str(_t5,1,q)),&
                p,&l,&f,&self);
            }
        }
      } while (n++ !=  _stop);

    }
    addstrtop(_P3str1("\001\000"),&p,&l,&f,&self);
  
  }
  _P3dec1(p,l);
  SYSTEM_assert(result == p,_P3str1("\006Length"));
  self->RUNNER_trunner_DOT_fcommandline = result;
  return result;
}  /* commandline */

Constructor(RUNNER_trunner ) RUNNER_trunner_DOT_create(
  RUNNER_trunner self)
{
  ValueCast(RUNNER_trunner,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->RUNNER_trunner_DOT_fmsghandler = ValueCast(RUNNER_tmsghandler,
    RUNNER_tmsghandler_DOT_create(ValueCast(RUNNER_tmsghandler,
    _P3alloc_object(&RUNNER_tmsghandler_CD)),_P3str1("\006Runner")));
  self->RUNNER_trunner_DOT_fisrunning = SYSTEM_false;
  _P3strclr(self->RUNNER_trunner_DOT_fexecutable);
  self->RUNNER_trunner_DOT_fcommandline = NULL;
  self->RUNNER_trunner_DOT_fparams = ValueCast(GMSOBJ_txstrings,
    SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,_P3alloc_object(&
    GMSOBJ_txstrings_CD))));
  _P3strclr(self->RUNNER_trunner_DOT_fworkdir);
  return self;
}  /* create */

Destructor(RUNNER_trunner ) RUNNER_trunner_DOT_destroy(
  RUNNER_trunner self)
{
  if (self->RUNNER_trunner_DOT_fcommandline != NULL) 
    _P3freemem(self->RUNNER_trunner_DOT_fcommandline);
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    RUNNER_trunner_DOT_fparams));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    RUNNER_trunner_DOT_fmsghandler));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_boolean ) RUNNER_trunner_DOT_errorwhenrunning(
  RUNNER_trunner self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_boolean result;

  result = self->RUNNER_trunner_DOT_fisrunning;
  if (result) 
    {
      _P3STR_255 _t1;
      _P3STR_255 _t2;

      RUNNER_tmsghandler_DOT_errormessage(self->
        RUNNER_trunner_DOT_fmsghandler,RUNNER_ec_cannot_modify,
        _P3strcat(_t2,255,_P3strcat(_t1,255,_P3str1("\016Cannot modify "),
        s),_P3str1("\026when process is active")));
    }
  return result;
}  /* errorwhenrunning */

Function(SYSTEM_integer ) RUNNER_trunner_DOT_paramscount(
  RUNNER_trunner self)
{
  SYSTEM_integer result;

  result = self->RUNNER_trunner_DOT_fparams->GMSOBJ_txlist_DOT_fcount;
  return result;
}  /* paramscount */

Procedure RUNNER_trunner_DOT_setexecutable(
  RUNNER_trunner self,
  const SYSTEM_ansichar *v)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\012Executable"))) 
    return;
  _P3strcpy(self->RUNNER_trunner_DOT_fexecutable,255,v);
  RUNNER_trunner_DOT_commandlinechanged(self);
}  /* setexecutable */

Procedure RUNNER_trunner_DOT_setinherithandles(
  RUNNER_trunner self,
  SYSTEM_boolean v)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\016InheritHandles"))) 
    return;
  self->RUNNER_trunner_DOT_finherithandles = v;
}  /* setinherithandles */

Procedure RUNNER_trunner_DOT_setuseshell(
  RUNNER_trunner self,
  SYSTEM_boolean v)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\010UseShell"))) 
    return;
  self->RUNNER_trunner_DOT_fuseshell = v;
  RUNNER_trunner_DOT_commandlinechanged(self);
}  /* setuseshell */

Procedure RUNNER_trunner_DOT_setworkdir(
  RUNNER_trunner self,
  const SYSTEM_ansichar *v)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\007WorkDir"))) 
    return;
  _P3strcpy(self->RUNNER_trunner_DOT_fworkdir,255,v);
}  /* setworkdir */

Procedure RUNNER_trunner_DOT_commandlinechanged(
  RUNNER_trunner self)
{
  if (self->RUNNER_trunner_DOT_fcommandline != NULL) {
    _P3freemem(self->RUNNER_trunner_DOT_fcommandline);
    self->RUNNER_trunner_DOT_fcommandline = NULL;
  } 
}  /* commandlinechanged */

Function(SYSTEM_integer ) RUNNER_trunner_DOT_getverbose(
  RUNNER_trunner self)
{
  SYSTEM_integer result;

  result = self->RUNNER_trunner_DOT_fmsghandler->
    RUNNER_tmsghandler_DOT_fverbose;
  return result;
}  /* getverbose */

Procedure RUNNER_trunner_DOT_setverbose(
  RUNNER_trunner self,
  SYSTEM_integer v)
{
  self->RUNNER_trunner_DOT_fmsghandler->
    RUNNER_tmsghandler_DOT_fverbose = v;
}  /* setverbose */

Procedure RUNNER_trunner_DOT_setvisible(
  RUNNER_trunner self,
  RUNNER_tvisible v)
{
  if (RUNNER_trunner_DOT_errorwhenrunning(self,_P3str1("\007Visible"))) 
    return;
  self->RUNNER_trunner_DOT_fvisible = v;
}  /* setvisible */

static Function(SYSTEM_P3_pansichar ) RUNNER_strtopchar(
  const SYSTEM_ansichar *s)
{
  SYSTEM_P3_pansichar result;
  SYSTEM_P3_pansichar p;
  SYSTEM_integer l;

  l = SYSTEM_length(s);
  if (l == 0) { 
    result = NULL;
  } else {
    _P3getmem(p,l + 1);
    result = p;
    SYSTEM_move(&s[1],p,l);
    _P3inc1(p,l);
    *p = _P3char('\000');
  } 
  return result;
}  /* strtopchar */

static Function(SYSTEM_ansichar *) RUNNER_pchartostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_P3_pansichar p)
{
  SYSTEM_integer w;

  w = 0;
  while (w < 255 && *p != _P3char('\000')) {
    w = w + 1;
    result[w] = *p;
    _P3inc0(p);
  
}
  result[0] = ValueCast(SYSTEM_ansichar,w);
  return result;
}  /* pchartostr */

Function(SYSTEM_integer ) RUNNER_trunner_DOT_startandwait(
  RUNNER_trunner self)
{
  SYSTEM_integer result;

  if (self->RUNNER_trunner_DOT_fisrunning) {
    RUNNER_tmsghandler_DOT_errormessage(self->
      RUNNER_trunner_DOT_fmsghandler,RUNNER_ec_process_active,_P3str1("\036Cannot start an active process"));
    result = RUNNER_ec_process_active;
    return result;
  } 
  if (*RUNNER_trunner_DOT_commandline(self) == _P3char('\000')) {
    RUNNER_tmsghandler_DOT_errormessage(self->
      RUNNER_trunner_DOT_fmsghandler,RUNNER_ec_empty_cmd_line,_P3str1("\031No Command Line specified"));
    result = RUNNER_ec_empty_cmd_line;
    return result;
  } 
  if (self->RUNNER_trunner_DOT_fuseshell) {
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;

      RUNNER_tmsghandler_DOT_debugmessage(self->
        RUNNER_trunner_DOT_fmsghandler,_P3strcat(_t2,255,_P3str1("\013Use shell: "),
        RUNNER_pchartostr(_t1,255,RUNNER_trunner_DOT_commandline(
        self))));
    }
    result = P3PROCESS_p3systemp(RUNNER_trunner_DOT_commandline(self),&
      self->RUNNER_trunner_DOT_fprogrc);
  } else {
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;

      RUNNER_tmsghandler_DOT_debugmessage(self->
        RUNNER_trunner_DOT_fmsghandler,_P3strcat(_t2,255,_P3str1("\015Direct call: "),
        RUNNER_pchartostr(_t1,255,RUNNER_trunner_DOT_commandline(
        self))));
    }
    result = P3PROCESS_p3execp(RUNNER_trunner_DOT_commandline(self),&
      self->RUNNER_trunner_DOT_fprogrc);
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;
      SYSTEM_shortstring _t4;
      _P3STR_255 _t5;

      RUNNER_tmsghandler_DOT_debugmessage(self->
        RUNNER_trunner_DOT_fmsghandler,_P3strcat(_t5,255,_P3strcat(
        _t3,255,_P3strcat(_t2,255,_P3str1("\011Return = "),
        SYSUTILS_P3_inttostr(_t1,255,result)),_P3str1("\007  RC = ")),
        SYSUTILS_P3_inttostr(_t4,255,self->
        RUNNER_trunner_DOT_fprogrc)));
    }
  } 
  return result;
}  /* startandwait */

Constructor(RUNNER_tmsghandler ) RUNNER_tmsghandler_DOT_create(
  RUNNER_tmsghandler self,
  const SYSTEM_ansichar *msgpfx)
{
  ValueCast(RUNNER_tmsghandler,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  _P3strcpy(self->RUNNER_tmsghandler_DOT_fmsgpfx,255,msgpfx);
  self->RUNNER_tmsghandler_DOT_fverbose = 1;
  return self;
}  /* create */

Procedure RUNNER_tmsghandler_DOT_debugmessage(
  RUNNER_tmsghandler self,
  const SYSTEM_ansichar *s)
{
  if (self->RUNNER_tmsghandler_DOT_fverbose >= 2) 
    {
      _P3STR_255 _t1;
      _P3STR_255 _t2;

      RUNNER_tmsghandler_DOT_logmessage(self,_P3strcat(_t2,255,
        _P3strcat(_t1,255,self->RUNNER_tmsghandler_DOT_fmsgpfx,_P3str1("\002: ")),
        s));
    }
}  /* debugmessage */

Procedure RUNNER_tmsghandler_DOT_errormessage(
  RUNNER_tmsghandler self,
  SYSTEM_integer ec,
  const SYSTEM_ansichar *s)
{
  _Iplus_bgn();
  _P3write_s0(_P3str1("\013*** Error: "));
  _P3write_s0(s);
  _P3writeln();
  _Iplus_end();
}  /* errormessage */

Procedure RUNNER_tmsghandler_DOT_logmessage(
  RUNNER_tmsghandler self,
  const SYSTEM_ansichar *s)
{
  if (self->RUNNER_tmsghandler_DOT_fverbose >= 1) {
    _Iplus_bgn();
    _P3write_s0(s);
    _P3writeln();
    _Iplus_end();
    _Iplus_bgn();
    _P3Tflush(SYSTEM_output);
    _Iplus_end();
  } 
}  /* logmessage */

/* unit runner */
void _Init_Module_runner(void)
{
} /* _Init_Module_runner */

void _Final_Module_runner(void)
{
} /* _Final_Module_runner */

