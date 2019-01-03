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
#include "gmsglobx.h"
#include "gmsobj.h"
#include "paldoorg.h"

_P3STR_15 PALDOORG_cmexliccodes = {8,'0','1','0','2','0','3','0','4'};

void * const PALDOORG_tpalobject_VT[] = {(void*)&
  PALDOORG_tpalobject_DOT_destroy};

/* Class descriptor for 'tpalobject' */
const SYSTEM_classdescriptor_t PALDOORG_tpalobject_CD = {
  _P3str1("\012tpalobject"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(PALDOORG_tpalobject_OD), PALDOORG_tpalobject_VT, NULL};

static SYSTEM_integer PALDOORG_val1, PALDOORG_val2, PALDOORG_val3;
typedef struct PALDOORG_tusermem_S {
  SYSTEM_shortstring webpage;
  SYSTEM_integer nlnl;
} PALDOORG_tusermem;


static Procedure PALDOORG_blankfill(
  SYSTEM_ansichar *t,
  SYSTEM_integer maxt,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer i;

  _P3setlength(t,maxt,255);
  { register SYSTEM_int32 _stop = maxt;
    if ((i = 1) <=  _stop) do {
      if (i <= ValueCast(SYSTEM_int32,SYSTEM_length(s))) { 
        t[i] = s[i];
      } else 
        t[i] = _P3char(' ');
    } while (i++ !=  _stop);

  }
}  /* blankfill */

static Procedure PALDOORG_dlice(void)
{
  PALDOORG_val1 = 7;
  PALDOORG_val2 = 19;
  PALDOORG_val3 = 83;
}  /* dlice */

static Function(SYSTEM_boolean ) PALDOORG_platformsvals(
  SYSTEM_integer n,
  SYSTEM_integer *val1,
  SYSTEM_integer *val2,
  SYSTEM_integer *val3)
{
  SYSTEM_boolean result;

  result = SYSTEM_true;
  switch (n) {
    case 1: 
      *val1 = 3;
      *val2 = 97;
      *val3 = 13;
      break;
    case 2: 
      *val1 = 3;
      *val2 = 97;
      *val3 = 13;
      break;
    case 3: 
      *val1 = 7;
      *val2 = 19;
      *val3 = 83;
      break;
    case 4: 
      *val1 = 79;
      *val2 = 23;
      *val3 = 11;
      break;
    case 5: 
      *val1 = 11;
      *val2 = 79;
      *val3 = 19;
      break;
    case 6: 
      *val1 = 73;
      *val2 = 23;
      *val3 = 13;
      break;
    case 7: 
      *val1 = 7;
      *val2 = 83;
      *val3 = 19;
      break;
    case 8: 
      *val1 = 7;
      *val2 = 83;
      *val3 = 19;
      break;
    case 9: 
      *val1 = 7;
      *val2 = 19;
      *val3 = 83;
      break;
    case 10: 
      *val1 = 73;
      *val2 = 23;
      *val3 = 13;
      break;
    case 11: 
      *val1 = 79;
      *val2 = 23;
      *val3 = 11;
      break;
    case 12: 
      *val1 = 79;
      *val2 = 23;
      *val3 = 11;
      break;
    case 13: 
      *val1 = 5;
      *val2 = 89;
      *val3 = 19;
      break;
    default:
      result = SYSTEM_false;
  }
  return result;
}  /* platformsvals */

Procedure PALDOORG_tpalobject_DOT_gutsofcreate(
  PALDOORG_tpalobject self)
{
  self->PALDOORG_tpalobject_DOT_juliantoday = SYSTEM_trunc(
    SYSUTILS_P3_now() - 1);
  PALDOORG_blankfill(self->PALDOORG_tpalobject_DOT_gdlsysnam,16,_P3str1("\015Uninitialized"));
  _P3strclr(self->PALDOORG_tpalobject_DOT_gdlauditline);
  PALDOORG_blankfill(self->PALDOORG_tpalobject_DOT_gdlrelcpr,70,_P3str1("\075Copyright (C) 1987-2015 GAMS Development. All rights reserved"));
  PALDOORG_blankfill(self->PALDOORG_tpalobject_DOT_gdlreldat,21,_P3str1("\025ALFA Released 16Feb15"));
  PALDOORG_blankfill(self->PALDOORG_tpalobject_DOT_gdllicdat,12,_P3str1("\014Jan 26, 2015"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_gdlrelmaj,2,_P3str1("\00224"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_gdlrelmin,1,_P3str1("\0015"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_gdlrelgold,1,_P3str1("\0010"));
  self->PALDOORG_tpalobject_DOT_gdllicjul = 42029;
  _P3strcpy(self->PALDOORG_tpalobject_DOT_gdlrelplc,3,_P3str1("\003LEX"));
  PALDOORG_blankfill(self->PALDOORG_tpalobject_DOT_gdlrelplt,22,_P3str1("\017x86 64bit/Linux"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_gdlbldcod,3,_P3str1("\003LEG"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_gdlrevision,12,_P3str1("\006r50613"));
  _P3strcat(self->PALDOORG_tpalobject_DOT_gdlsysver,3,self->
    PALDOORG_tpalobject_DOT_gdlrelmaj,self->
    PALDOORG_tpalobject_DOT_gdlrelmin);
  if (PALDOORG_tpalobject_DOT_palisalfa(self)) 
    P3UTILS_p3nopopups();
  PALDOORG_tpalobject_DOT_makedemolicense(self);
  self->PALDOORG_tpalobject_DOT_rowcnt = 301;
  self->PALDOORG_tpalobject_DOT_colcnt = 301;
  self->PALDOORG_tpalobject_DOT_nzcnt = 2001;
  self->PALDOORG_tpalobject_DOT_nlnzcnt = 1001;
  self->PALDOORG_tpalobject_DOT_disccnt = 51;
}  /* gutsofcreate */

Constructor(PALDOORG_tpalobject ) PALDOORG_tpalobject_DOT_create(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *msg)
{
  ValueCast(PALDOORG_tpalobject,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  PALDOORG_tpalobject_DOT_gutsofcreate(self);
  _P3strclr(msg);
  return self;
}  /* create */

Constructor(PALDOORG_tpalobject ) PALDOORG_tpalobject_DOT_createx(
  PALDOORG_tpalobject self)
{
  ValueCast(PALDOORG_tpalobject,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  PALDOORG_tpalobject_DOT_gutsofcreate(self);
  return self;
}  /* createx */

Destructor(PALDOORG_tpalobject ) PALDOORG_tpalobject_DOT_destroy(
  PALDOORG_tpalobject self)
{
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure PALDOORG_tpalobject_DOT_palsetauditline(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *auditline)
{
  {
    SYSTEM_shortstring _t1;

    PALDOORG_blankfill(self->PALDOORG_tpalobject_DOT_gdlauditline,78,
      SYSTEM_copy(_t1,255,auditline,7,ValueCast(SYSTEM_int32,
      SYSTEM_length(auditline)) - 12));
  }
  SYSTEM_copy(self->PALDOORG_tpalobject_DOT_gdlsysnam,16,self->
    PALDOORG_tpalobject_DOT_gdlauditline,1,16);
}  /* palsetauditline */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palauditrun(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;
  SYSTEM_integer i;
  SYSTEM_shortstring audit;

  result = SYSTEM_false;
  if (SYSTEM_P3_paramcount() > 0) 
    {
      SYSTEM_shortstring _t1;

      if (SYSTEM_length(SYSTEM_P3_paramstr(_t1,255,1)) == 5) {
        SYSTEM_P3_paramstr(audit,255,1);
        for (i = 1;i <= (SYSTEM_int32)5;++i) {
          audit[i] = SYSTEM_upcase(audit[i]);
        }
        if (_P3strcmpE(audit,_P3str1("\005AUDIT"))) 
          result = SYSTEM_true;
      } 
    }
  return result;
}  /* palauditrun */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlauditline);
  return result;
}  /* palgetauditline */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetcpr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlrelcpr);
  return result;
}  /* palgetcpr */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgetver(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;
  SYSTEM_integer rc;

  _P3val_i(self->PALDOORG_tpalobject_DOT_gdlsysver,result,&rc);
  return result;
}  /* palgetver */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetrevision(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlrevision);
  return result;
}  /* palgetrevision */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetrel(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  {
    _P3STR_3 _t1;

    _P3strcat(result,_len_ret,_P3strcat(_t1,3,self->
      PALDOORG_tpalobject_DOT_gdlrelmaj,_P3str1("\001.")),self->
      PALDOORG_tpalobject_DOT_gdlrelmin);
  }
  return result;
}  /* palgetrel */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetgold(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlrelgold);
  return result;
}  /* palgetgold */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetcod(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlrelplc);
  return result;
}  /* palgetcod */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgethdr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlrelplt);
  return result;
}  /* palgethdr */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgetjul(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  result = self->PALDOORG_tpalobject_DOT_gdllicjul;
  return result;
}  /* palgetjul */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetreldat(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlreldat);
  return result;
}  /* palgetreldat */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetbldcod(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdlbldcod);
  return result;
}  /* palgetbldcod */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetlicdat(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_gdllicdat);
  return result;
}  /* palgetlicdat */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palisbeta(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  {
    SYSTEM_shortstring _t1;

    result = _P3strcmpE(_P3str1("\004BETA"),SYSTEM_copy(_t1,255,
      self->PALDOORG_tpalobject_DOT_gdlreldat,1,4));
  }
  return result;
}  /* palisbeta */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palisalfa(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  {
    SYSTEM_shortstring _t1;

    result = _P3strcmpE(_P3str1("\004ALFA"),SYSTEM_copy(_t1,255,
      self->PALDOORG_tpalobject_DOT_gdlreldat,1,4));
  }
  return result;
}  /* palisalfa */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgettoday(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  result = self->PALDOORG_tpalobject_DOT_juliantoday;
  return result;
}  /* palgettoday */

Procedure PALDOORG_tpalobject_DOT_palauditfields(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *auditline,
  SYSTEM_ansichar *v1,
  SYSTEM_ansichar *v2,
  SYSTEM_ansichar *v3)
{
  {
    SYSTEM_shortstring _t2;

    SYSUTILS_P3_trim(v1,255,SYSTEM_copy(_t2,255,auditline,1,16));
  }
  {
    SYSTEM_shortstring _t2;

    SYSUTILS_P3_trim(v2,255,SYSTEM_copy(_t2,255,auditline,18,12));
  }
  {
    SYSTEM_shortstring _t2;

    SYSUTILS_P3_trim(v3,255,SYSTEM_copy(_t2,255,auditline,31,255));
  }
}  /* palauditfields */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_palgetshortauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_shortstring v1, v2, v3;

  {
    SYSTEM_shortstring _t1;

    PALDOORG_tpalobject_DOT_palauditfields(self,
      PALDOORG_tpalobject_DOT_palgetauditline(_t1,255,self),v1,v2,
      v3);
  }
  {
    _P3STR_255 _t1;

    _P3strcat(result,_len_ret,_P3strcat(_t1,255,v1,_P3str1("\001 ")),
      v3);
  }
  return result;
}  /* palgetshortauditline */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensereadu(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *filename,
  SYSTEM_ansichar *msg,
  SYSTEM_integer *rc)
{
  SYSTEM_boolean result;
  cnstdef {bufclean = 325};
  cnstdef {bufmax = 490};
  SYSTEM_untypedfile f;
  SYSTEM_integer numread;
  typedef SYSTEM_uint16 _sub_1PALLICENSEREADU;
  typedef SYSTEM_ansichar _arr_0PALLICENSEREADU[490];
  _arr_0PALLICENSEREADU buf;
  SYSTEM_integer i, j, off, len6;

  result = SYSTEM_false;
  _Iplus_bgn();
  _P3assign(f,filename);
  _Iplus_end();
  SYSTEM_filemode = 0;
  _Iminus_bgn();
  _P3Ureset(f,1);
  _Iminus_end();
  *rc = SYSTEM_ioresult();
  if (*rc != 0) {
    _P3strcat(msg,255,_P3str1("\031Cannot open license file "),
      filename);
    return result;
  } 
  self->PALDOORG_tpalobject_DOT_licenselevel = 0;
  self->PALDOORG_tpalobject_DOT_licensestatus = 1;
  _Iplus_bgn();
  _P3blockR4(f,buf,sizeof(_arr_0PALLICENSEREADU),numread);
  _Iplus_end();
  j = 0;
  { register SYSTEM_int32 _stop = numread;
    if ((i = 1) <=  _stop) do {
      if (SYSTEM_ord(buf[i - 1]) > 32) {
        _P3inc0(j);
        buf[j - 1] = buf[i - 1];
      } 
    } while (i++ !=  _stop);

  }
  for (i = j + 1;i <= bufmax;++i) {
    buf[i - 1] = _P3char('_');
  }
  _P3setlength(self->PALDOORG_tpalobject_DOT_license1,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license2,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license3,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license4,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license5,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license6,0,65);
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    self->PALDOORG_tpalobject_DOT_license1[i] = buf[i - 1];
  }
  off = PALDOORG_maxlicense;
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    self->PALDOORG_tpalobject_DOT_license2[i] = buf[i + off - 1];
  }
  off = off + PALDOORG_maxlicense;
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    self->PALDOORG_tpalobject_DOT_license3[i] = buf[i + off - 1];
  }
  off = off + PALDOORG_maxlicense;
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    self->PALDOORG_tpalobject_DOT_license4[i] = buf[i + off - 1];
  }
  off = off + PALDOORG_maxlicense;
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    self->PALDOORG_tpalobject_DOT_license5[i] = buf[i + off - 1];
  }
  off = off + PALDOORG_maxlicense;
  len6 = j - bufclean;
  if (len6 > 0) {
    if (len6 > 65) 
      len6 = PALDOORG_maxlicense;
    { register SYSTEM_int32 _stop = len6;
      if ((i = 1) <=  _stop) do {
        self->PALDOORG_tpalobject_DOT_license6[i] = buf[i + off - 1];
      } while (i++ !=  _stop);

    }
    _P3setlength(self->PALDOORG_tpalobject_DOT_license6,len6,65);
  } 
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    if (self->PALDOORG_tpalobject_DOT_license1[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license1[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license2[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license2[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license3[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license3[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license4[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license4[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license5[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license5[i] = _P3char(' ');
  
  }
  {
    SYSTEM_shortstring _t1;

    self->PALDOORG_tpalobject_DOT_licenselevel = SYSUTILS_P3_strtoint(
      SYSTEM_copy(_t1,255,self->PALDOORG_tpalobject_DOT_license3,9,2));
  }
  self->PALDOORG_tpalobject_DOT_licenseversion = 
    PALDOORG_tpalobject_DOT_pallicensegetversion(self);
  self->PALDOORG_tpalobject_DOT_licenseactsub = 
    PALDOORG_tpalobject_DOT_licensegetmaxsubsys(self);
  self->PALDOORG_tpalobject_DOT_currentsearchhelper = self->
    PALDOORG_tpalobject_DOT_licenseactsub + 1;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
  _Iplus_bgn();
  _P3Uclose(f);
  _Iplus_end();
  result = SYSTEM_true;
  return result;
}  /* pallicensereadu */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_msgadd(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *msg)
{
  SYSTEM_integer result;

  if (self->PALDOORG_tpalobject_DOT_ml == NULL) 
    self->PALDOORG_tpalobject_DOT_ml = ValueCast(SYSTEM_pointer,
      SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,
      _P3alloc_object(&GMSOBJ_txstrings_CD))));
  result = GMSOBJ_txstrings_DOT_add(ValueCast(GMSOBJ_txstrings,self->
    PALDOORG_tpalobject_DOT_ml),msg);
  return result;
}  /* msgadd */

Procedure PALDOORG_tpalobject_DOT_makedemolicense(
  PALDOORG_tpalobject self)
{
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license1,65,_P3str1("\101GAMS Development Corporation, Washington, DC   G871201/0000CA-ANY"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license2,65,_P3str1("\101Free Demo,  202-342-0180,  sales@gams.com,  www.gams.com         "));
  {
    SYSTEM_shortstring _t1;

    _P3strcat(self->PALDOORG_tpalobject_DOT_license3,65,
      STRUTILX_inttostr(_t1,255,SYSTEM_round(10000000 + 12 * 
      SYSUTILS_P3_now())),_P3str1("\071000S                                                   00"));
  }
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license4,65,_P3str1("\10101234567000000                                                   "));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license5,65,_P3str1("\101DC0000       Ref: Generated by Base            A Demo            "));
  _P3strclr(self->PALDOORG_tpalobject_DOT_license6);
  self->PALDOORG_tpalobject_DOT_licensestatus = 1;
  self->PALDOORG_tpalobject_DOT_licenselevel = 0;
  self->PALDOORG_tpalobject_DOT_licenseversion = 2;
  self->PALDOORG_tpalobject_DOT_licenseactsub = PALDOORG_licensemaxsub;
  self->PALDOORG_tpalobject_DOT_currentsearchhelper = self->
    PALDOORG_tpalobject_DOT_licenseactsub + 1;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
}  /* makedemolicense */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetversion(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  if (self->PALDOORG_tpalobject_DOT_license1[55] == _P3char('/')) { 
    result = 2;
  } else 
    if (self->PALDOORG_tpalobject_DOT_license1[55] == _P3char(':')) { 
      result = 1;
    } else 
      result = 0;
  return result;
}  /* pallicensegetversion */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_licensegetmaxsubsys(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  if (self->PALDOORG_tpalobject_DOT_licenseversion > 1) { 
    result = PALDOORG_licensemaxsub;
  } else 
    result = 19;
  return result;
}  /* licensegetmaxsubsys */

Procedure PALDOORG_tpalobject_DOT_pallicenseregistergams(
  PALDOORG_tpalobject self,
  SYSTEM_integer linenr,
  const SYSTEM_ansichar *liceline)
{
  switch (linenr) {
    case 1: 
      _P3strcpy(self->PALDOORG_tpalobject_DOT_license1,65,liceline);
      break;
    case 2: 
      _P3strcpy(self->PALDOORG_tpalobject_DOT_license2,65,liceline);
      break;
    case 3: 
      _P3strcpy(self->PALDOORG_tpalobject_DOT_license3,65,liceline);
      break;
    case 4: 
      _P3strcpy(self->PALDOORG_tpalobject_DOT_license4,65,liceline);
      break;
    case 5: 
      _P3strcpy(self->PALDOORG_tpalobject_DOT_license5,65,liceline);
      break;
    case 6: 
      _P3strcpy(self->PALDOORG_tpalobject_DOT_license6,65,liceline);
      break;
    default: break;
  }
}  /* pallicenseregistergams */

Procedure PALDOORG_tpalobject_DOT_pallicenseregistergamsdone(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer i;

  self->PALDOORG_tpalobject_DOT_licenselevel = 0;
  self->PALDOORG_tpalobject_DOT_licensestatus = 1;
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    if (ValueCast(SYSTEM_int32,SYSTEM_length(self->
      PALDOORG_tpalobject_DOT_license1)) < i) 
      self->PALDOORG_tpalobject_DOT_license1[i] = _P3char(' ');
    if (ValueCast(SYSTEM_int32,SYSTEM_length(self->
      PALDOORG_tpalobject_DOT_license2)) < i) 
      self->PALDOORG_tpalobject_DOT_license2[i] = _P3char(' ');
    if (ValueCast(SYSTEM_int32,SYSTEM_length(self->
      PALDOORG_tpalobject_DOT_license3)) < i) 
      self->PALDOORG_tpalobject_DOT_license3[i] = _P3char(' ');
    if (ValueCast(SYSTEM_int32,SYSTEM_length(self->
      PALDOORG_tpalobject_DOT_license4)) < i) 
      self->PALDOORG_tpalobject_DOT_license4[i] = _P3char(' ');
    if (ValueCast(SYSTEM_int32,SYSTEM_length(self->
      PALDOORG_tpalobject_DOT_license5)) < i) 
      self->PALDOORG_tpalobject_DOT_license5[i] = _P3char(' ');
  
  }
  _P3setlength(self->PALDOORG_tpalobject_DOT_license1,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license2,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license3,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license4,
    PALDOORG_maxlicense,65);
  _P3setlength(self->PALDOORG_tpalobject_DOT_license5,
    PALDOORG_maxlicense,65);
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    if (self->PALDOORG_tpalobject_DOT_license1[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license1[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license2[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license2[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license3[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license3[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license4[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license4[i] = _P3char(' ');
    if (self->PALDOORG_tpalobject_DOT_license5[i] == _P3char('_')) 
      self->PALDOORG_tpalobject_DOT_license5[i] = _P3char(' ');
  
  }
  {
    SYSTEM_shortstring _t1;

    self->PALDOORG_tpalobject_DOT_licenselevel = SYSUTILS_P3_strtoint(
      SYSTEM_copy(_t1,255,self->PALDOORG_tpalobject_DOT_license3,9,2));
  }
  self->PALDOORG_tpalobject_DOT_licenseversion = 
    PALDOORG_tpalobject_DOT_pallicensegetversion(self);
  self->PALDOORG_tpalobject_DOT_licenseactsub = 
    PALDOORG_tpalobject_DOT_licensegetmaxsubsys(self);
  self->PALDOORG_tpalobject_DOT_currentsearchhelper = self->
    PALDOORG_tpalobject_DOT_licenseactsub + 1;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
}  /* pallicenseregistergamsdone */

Procedure PALDOORG_tpalobject_DOT_pallicenseregistersystem(
  PALDOORG_tpalobject self,
  SYSTEM_integer numcodes,
  const SYSTEM_ansichar *codes,
  SYSTEM_integer magicnum,
  SYSTEM_integer globalflag)
{
  self->PALDOORG_tpalobject_DOT_subsyssecondary = numcodes;
  _P3strcpy(self->PALDOORG_tpalobject_DOT_subsyscode,20,codes);
  self->PALDOORG_tpalobject_DOT_checksum = magicnum;
  self->PALDOORG_tpalobject_DOT_isglobal = globalflag;
}  /* pallicenseregistersystem */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_tampercheck(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;
  _P3STR_31 platformheader;
  SYSTEM_integer sysdatenumber, code1, i;

  sysdatenumber = self->PALDOORG_tpalobject_DOT_gdllicjul;
  _P3strcpy(platformheader,16,self->
    PALDOORG_tpalobject_DOT_gdlrelplt);
  {
    _P3STR_3 _t1;

    _P3val_i(_P3strcat(_t1,3,self->
      PALDOORG_tpalobject_DOT_gdlrelmaj,self->
      PALDOORG_tpalobject_DOT_gdlrelmin),code1,&i);
  }
  code1 = sysdatenumber /  PALDOORG_val1 + code1 * PALDOORG_val2;
  for (i = 1;i <= (SYSTEM_int32)16;++i) {
    code1 = code1 + SYSTEM_ord(platformheader[i]) * i * PALDOORG_val3;
  }
  { register SYSTEM_int32 _stop = self->
      PALDOORG_tpalobject_DOT_subsyssecondary - 1;
    if ((i = 0) <=  _stop) do {
      code1 = code1 + SYSTEM_ord(self->
        PALDOORG_tpalobject_DOT_subsyscode[i * 2 + 1]) * 97 + 
        SYSTEM_ord(self->PALDOORG_tpalobject_DOT_subsyscode[i * 2 + 2]) * 7;
    } while (i++ !=  _stop);

  }
  result = self->PALDOORG_tpalobject_DOT_checksum != code1;
  return result;
}  /* tampercheck */

static Function(SYSTEM_boolean ) PALDOORG_havedata(
  SYSTEM_untyped *buf,
  SYSTEM_integer len,
  SYSTEM_pointer usermem)
{
  SYSTEM_boolean result;
  SYSTEM_integer k;
  SYSTEM_P3_pansichar p;
  typedef PALDOORG_tusermem *_ptr_0HAVEDATA;
  _ptr_0HAVEDATA rp;
  GMSGEN_pansichararray q;

  result = SYSTEM_false;
  rp = ValueCast(_ptr_0HAVEDATA,usermem);
  q = ValueCast(GMSGEN_pansichararray,rp->webpage);
  p = ValueCast(SYSTEM_P3_pansichar,buf);
  if (rp->nlnl < 0) { 
    { register SYSTEM_int32 _stop = len;
      if ((k = 1) <=  _stop) do {
        _Iplus_bgn();
        _P3write_c0(*p);
        _Iplus_end();
        _P3inc0(p);
      
      } while (k++ !=  _stop);

    }
  } else 
    { register SYSTEM_int32 _stop = len;
      if ((k = 1) <=  _stop) do {
        if (rp->nlnl > 3) {
          if (ValueCast(SYSTEM_int32,(*q)[0]) == 255) 
            return result;
          _P3inc0((*q)[0]);
          (*q)[ValueCast(SYSTEM_int32,(*q)[0])] = *p;
        } else 
          if (*p == _P3char('\015') && (rp->nlnl == 0 || rp->nlnl == 2) || *
            p == _P3char('\012') && (rp->nlnl == 1 || rp->nlnl == 3)) { 
            _P3inc0(rp->nlnl);
          } else 
            rp->nlnl = 0;
        _P3inc0(p);
      
      } while (k++ !=  _stop);

    }
  result = SYSTEM_true;
  return result;
}  /* havedata */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_licensecheckv1to3(
  PALDOORG_tpalobject self,
  SYSTEM_integer v1,
  SYSTEM_integer v2,
  SYSTEM_integer v3)
{
  SYSTEM_boolean result;
  SYSTEM_longint code1, code2;
  SYSTEM_integer i;
  SYSTEM_shortstring msg, nodelocktype, hostid, acthostid;
  PALDOORG_tusermem u;

  if (self->PALDOORG_tpalobject_DOT_licenselevel == 4) { 
    code1 = 159;
  } else 
    code1 = 0;
  code2 = 0;
  for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
    code1 = code1 + SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license1[
      i]) * i * v1 + SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license2[
      i]) * i * v2 + SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license3[
      i]) * i * v3;
  }
  for (i = 1;i <= (SYSTEM_int32)8;++i) {
    code1 = code1 - SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license3[
      i]) * i * v3;
    code2 = code2 * 10 + SYSTEM_ord(self->
      PALDOORG_tpalobject_DOT_license3[i]) - 48;
  
  }
  if (code1 == code2 && _P3SET_in_1(self->
    PALDOORG_tpalobject_DOT_license1[55],_P3char(':'),_P3SET_equal(
    self->PALDOORG_tpalobject_DOT_license1[55],_P3char('/')))) {
    code1 = 0;
    code2 = 0;
    for (i = 1;i <= (SYSTEM_int32)PALDOORG_maxlicense;++i) {
      code1 = code1 + SYSTEM_ord(self->
        PALDOORG_tpalobject_DOT_license3[i]) * i * v1 + SYSTEM_ord(
        self->PALDOORG_tpalobject_DOT_license4[i]) * i * v2 + 
        SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license5[i]) * i * v3;
    }
    for (i = 1;i <= (SYSTEM_int32)8;++i) {
      code1 = code1 - SYSTEM_ord(self->
        PALDOORG_tpalobject_DOT_license4[i]) * i * v2;
      code2 = code2 * 10 + SYSTEM_ord(self->
        PALDOORG_tpalobject_DOT_license4[i]) - 48;
    
    }
  } 
  if (code1 == code2 && self->PALDOORG_tpalobject_DOT_licenselevel == 4) {
    result = SYSTEM_true;
    i = STRUTILX_lchpossp(_P3char('@'),self->
      PALDOORG_tpalobject_DOT_license2,1);
    if (i == 0) 
      return result;
    {
      SYSTEM_shortstring _t2;

      SYSUTILS_P3_trim(nodelocktype,255,SYSTEM_copy(_t2,255,self->
        PALDOORG_tpalobject_DOT_license2,1,i));
    }
    if (SYSUTILS_P3_sametext(_P3str1("\004aws@"),nodelocktype)) {
      i = STRUTILX_lchpossp(_P3char('@'),self->
        PALDOORG_tpalobject_DOT_license2,1);
      SYSTEM_copy(hostid,255,self->PALDOORG_tpalobject_DOT_license2,
        i + 1,STRUTILX_lchpossp(_P3char(' '),self->
        PALDOORG_tpalobject_DOT_license2,i + 1) - i - 1);
      u.webpage[0] = _P3char('\000');
      u.nlnl = 0;
      P3UTILS_p3getfromurl(_P3str1("\017169.254.169.254"),_P3str1("\037latest/meta-data/product-codes/"),80,ValueCast(
        P3UTILS_thavedatacb,&PALDOORG_havedata),ValueCast(
        SYSTEM_pointer,&u),msg);
      if (_P3strcmpN(msg,_P3str1("\000"))) {
        _Iplus_bgn();
        _P3writeln();
        _Iplus_end();
        _Iplus_bgn();
        _P3write_s0(_P3str1("\003***"));
        _P3writeln();
        _Iplus_end();
        _Iplus_bgn();
        _P3write_s0(_P3str1("\117*** Problems retrieving http://169.254.169.254/latest/meta-data/product-codes/."));
        _P3writeln();
        _Iplus_end();
        _Iplus_bgn();
        {
          _P3STR_255 _t1;

          _P3write_s0(_P3strcat(_t1,255,_P3str1("\004*** "),msg));
          _P3writeln();
        }
        _Iplus_end();
        _Iplus_bgn();
        _P3write_s0(_P3str1("\003***"));
        _P3writeln();
        _Iplus_end();
        _Iplus_bgn();
        _P3write_s0(_P3str1("\035*** Responds from web server:"));
        _P3writeln();
        _Iplus_end();
        _Iplus_bgn();
        _P3write_s0(_P3str1("\003***"));
        _P3writeln();
        _Iplus_end();
        u.nlnl =  -1;
        P3UTILS_p3getfromurl(_P3str1("\017169.254.169.254"),_P3str1("\037latest/meta-data/product-codes/"),80,ValueCast(
          P3UTILS_thavedatacb,&PALDOORG_havedata),ValueCast(
          SYSTEM_pointer,&u),msg);
        _Iplus_bgn();
        _P3write_s0(_P3str1("\003***"));
        _P3writeln();
        _Iplus_end();
        return result;
      } 
      SYSUTILS_P3_trim(acthostid,255,u.webpage);
      if (!SYSUTILS_P3_sametext(hostid,acthostid)) 
        return result;
    } else 
      return result;
  } 
  result = code1 != code2;
  return result;
}  /* licensecheckv1to3 */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetplatform(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_copy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license1,63,3);
  return result;
}  /* pallicensegetplatform */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensevalidation(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  result = PALDOORG_tpalobject_DOT_licensecheckinternal(self);
  return result;
}  /* pallicensevalidation */

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicensevalidateforplatform(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *pf)
{
  SYSTEM_boolean result;
  SYSTEM_integer v1, v2, v3;

  PALDOORG_platformsvals(GMSGLOBX_platformslookup(pf),&v1,&v2,&v3);
  result = PALDOORG_tpalobject_DOT_licensecheckv1to3(self,v1,v2,v3);
  return result;
}  /* pallicensevalidateforplatform */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_licensecheckinternal(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;
  SYSTEM_integer v1, v2, v3;

  result = PALDOORG_tpalobject_DOT_licensecheckv1to3(self,
    PALDOORG_val1,PALDOORG_val2,PALDOORG_val3);
  {
    SYSTEM_shortstring _t1;

    if (result && SYSUTILS_P3_sametext(_P3str1("\003GEN"),
      PALDOORG_tpalobject_DOT_pallicensegetplatform(_t1,255,self))) 
      if (PALDOORG_platformsvals(GMSGLOBX_platformslookup(_P3str1("\003GEN")),&
        v1,&v2,&v3)) 
        result = PALDOORG_tpalobject_DOT_licensecheckv1to3(self,v1,v2,
          v3);
  }
  {
    SYSTEM_shortstring _t1;
    SYSTEM_shortstring _t2;

    if (result && P3PLATFORM_osfiletype() == P3PLATFORM_osfilewin && (
      SYSUTILS_P3_sametext(_P3str1("\003LNX"),
      PALDOORG_tpalobject_DOT_pallicensegetplatform(_t1,255,self)) || 
      SYSUTILS_P3_sametext(_P3str1("\003DAR"),
      PALDOORG_tpalobject_DOT_pallicensegetplatform(_t2,255,self))) && 
      _P3SET_i(4,P3UTILS_p3getwindowsversion(),_P3set1("\030"))) {
      {
        SYSTEM_shortstring _t1;

        if (SYSUTILS_P3_sametext(_P3str1("\003LNX"),
          PALDOORG_tpalobject_DOT_pallicensegetplatform(_t1,255,
          self)) && PALDOORG_platformsvals(GMSGLOBX_platformslookup(_P3str1("\003LNX")),&
          v1,&v2,&v3)) 
          result = PALDOORG_tpalobject_DOT_licensecheckv1to3(self,v1,
            v2,v3);
      }
      {
        SYSTEM_shortstring _t1;

        if (SYSUTILS_P3_sametext(_P3str1("\003DAR"),
          PALDOORG_tpalobject_DOT_pallicensegetplatform(_t1,255,
          self)) && PALDOORG_platformsvals(GMSGLOBX_platformslookup(_P3str1("\003DAR")),&
          v1,&v2,&v3)) 
          result = PALDOORG_tpalobject_DOT_licensecheckv1to3(self,v1,
            v2,v3);
      }
    } 
  }
  return result;
}  /* licensecheckinternal */

static Function(SYSTEM_integer ) PALDOORG_lnumtoint(
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  if (SYSTEM_length(s) == 0 || SYSTEM_length(s) > 1) { 
    result = 0;
  } else 
    if (SYSTEM_ord(s[1]) >= 48 && SYSTEM_ord(s[1]) <= 57) { 
      result = SYSTEM_ord(s[1]) - 48;
    } else 
      if (SYSTEM_ord(s[1]) >= 65 && SYSTEM_ord(s[1]) <= 90) { 
        result = SYSTEM_ord(s[1]) - 65 + 10;
      } else 
        if (SYSTEM_ord(s[1]) >= 97 && SYSTEM_ord(s[1]) <= 122) { 
          result = SYSTEM_ord(s[1]) - 97 + 36;
        } else 
          result = 0;
  return result;
}  /* lnumtoint */

static Function(SYSTEM_double ) PALDOORG_gamsencodedate(
  SYSTEM_integer year,
  SYSTEM_integer month,
  SYSTEM_integer day)
{
  SYSTEM_double result;
  SYSTEM_double a, b, c, d;
  SYSTEM_word yr, mo;

  a = year;
  b = month;
  c = day;
  d = SYSTEM_int((b - 1) /  12);
  a = a + d;
  b = b - d * 12;
  if (b <= 0) {
    a = a - 1;
    b = 12 + b;
  } 
  if (a < 1 || a > 9999) {
    result = 0;
    return result;
  } 
  yr = SYSTEM_trunc(a);
  mo = SYSTEM_trunc(b);
  result = SYSUTILS_P3_encodedate(yr,mo,1) + c - 2;
  return result;
}  /* gamsencodedate */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgetjuliandays(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer yr, month, day;

  {
    SYSTEM_shortstring _t1;

    yr = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,s,1,2));
  }
  if (yr < 87) { 
    yr = 2000 + yr;
  } else 
    yr = 1900 + yr;
  {
    SYSTEM_shortstring _t1;

    month = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,s,3,2));
  }
  {
    SYSTEM_shortstring _t1;

    day = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,s,5,2));
  }
  result = SYSTEM_trunc(PALDOORG_gamsencodedate(yr,month,day));
  return result;
}  /* palgetjuliandays */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetsubeval(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  if (self->PALDOORG_tpalobject_DOT_licenseversion < 2) { 
    result = SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license3[66 - 
      self->PALDOORG_tpalobject_DOT_currentsearchhelper]) - 48;
  } else 
    {
      SYSTEM_shortstring _t1;

      result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
        PALDOORG_tpalobject_DOT_license4,7 + 2 * self->
        PALDOORG_tpalobject_DOT_currentsearchhelper + 1,1));
    }
  return result;
}  /* pallicensegetsubeval */

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetsubmaint(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  result = 0;
  switch (self->PALDOORG_tpalobject_DOT_licenseversion) {
    case 1: 
      {
        SYSTEM_shortstring _t1;

        result = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,self->
          PALDOORG_tpalobject_DOT_license4,7 + 2 * self->
          PALDOORG_tpalobject_DOT_currentsearchhelper,2));
      }
      break;
    case 2: 
      {
        SYSTEM_shortstring _t1;

        result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
          PALDOORG_tpalobject_DOT_license4,7 + 2 * self->
          PALDOORG_tpalobject_DOT_currentsearchhelper,1));
      }
      break;
    default: break;
  }
  return result;
}  /* pallicensegetsubmaint */

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_licensechecksubinternal(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *msg,
  SYSTEM_integer numcodes,
  const SYSTEM_ansichar *codes)
{
  SYSTEM_boolean result;
  SYSTEM_shortstring lcode;
  SYSTEM_integer i, eval, maint, over;

  result = SYSTEM_false;
  _P3strclr(msg);
  if (PALDOORG_tpalobject_DOT_pallicensegetversion(self) < 1) {
    _P3strcpy(msg,255,_P3str1("\067Old style license file - cannot run with current system"));
    result = SYSTEM_true;
    return result;
  } 
  {
    SYSTEM_shortstring _t1;

    if (SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,self->
      PALDOORG_tpalobject_DOT_license3,9,2)) == 0 && 
      _P3strcmpE(codes,_P3str1("\0020S"))) {
      _P3strcpy(msg,255,_P3str1("\021Demo license only"));
      result = SYSTEM_false;
      return result;
    } 
  }
  {
    SYSTEM_shortstring _t1;

    if (SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,self->
      PALDOORG_tpalobject_DOT_license3,9,2)) > 0 && 
      PALDOORG_tpalobject_DOT_licensecheckinternal(self)) {
      _P3strcpy(msg,255,_P3str1("\031License validation failed"));
      result = SYSTEM_true;
      return result;
    } 
  }
  { register SYSTEM_int32 _stop = numcodes;
    if ((i = 1) <=  _stop) do {
      SYSTEM_copy(lcode,255,codes,i * 2 - 1,2);
      if (_P3strcmpE(lcode,_P3str1("\002FR"))) {
        _P3strclr(msg);
        return result;
      } 
      self->PALDOORG_tpalobject_DOT_currentsearchhelper = 0;
      while (PALDOORG_tpalobject_DOT_pallicensegetsubnext(self)) {

        {
          SYSTEM_shortstring _t1;

          if (_P3strcmpE(lcode,SYSTEM_copy(_t1,255,self->
            PALDOORG_tpalobject_DOT_license3,7 + 2 * self->
            PALDOORG_tpalobject_DOT_currentsearchhelper,2))) {
            eval = PALDOORG_tpalobject_DOT_pallicensegetsubeval(self);
            maint = PALDOORG_tpalobject_DOT_pallicensegetsubmaint(self);
            if (eval == 0) {
              if (maint == 0) {
                _P3strclr(msg);
                return result;
              } 
              {
                SYSTEM_shortstring _t1;

                over = self->PALDOORG_tpalobject_DOT_gdllicjul - 
                  PALDOORG_tpalobject_DOT_palgetjuliandays(self,
                  SYSTEM_copy(_t1,255,self->
                  PALDOORG_tpalobject_DOT_license1,49,6)) - 30 * 
                  maint;
              }
              if (over < 30) {
                _P3strclr(msg);
                return result;
              } 
              if ((PALDOORG_tpalobject_DOT_palisbeta(self) || 
                PALDOORG_tpalobject_DOT_palisalfa(self)) && self->
                PALDOORG_tpalobject_DOT_juliantoday - self->
                PALDOORG_tpalobject_DOT_gdllicjul < 60) {
                _P3strcpy(msg,255,_P3str1("\054BETA system running with expired maintenance"));
                return result;
              } 
              {
                SYSTEM_shortstring _t1;
                _P3STR_255 _t2;

                _P3strcat(msg,255,_P3strcat(_t2,255,_P3str1("\007Module "),
                  STRUTILX_inttostr(_t1,255,over)),_P3str1("\054 days too young for this license - demo only"));
              }
            } else {
              {
                SYSTEM_shortstring _t1;

                over = self->PALDOORG_tpalobject_DOT_juliantoday - 
                  PALDOORG_tpalobject_DOT_palgetjuliandays(self,
                  SYSTEM_copy(_t1,255,self->
                  PALDOORG_tpalobject_DOT_license1,49,6)) - 30 * 
                  eval;
              }
              if (over < 0) {
                _P3strclr(msg);
                return result;
              } 
              {
                SYSTEM_shortstring _t1;
                _P3STR_255 _t2;

                _P3strcat(msg,255,_P3strcat(_t2,255,_P3str1("\027Evaluation has expired "),
                  STRUTILX_inttostr(_t1,255,over)),_P3str1("\011 days ago"));
              }
              if (over <= 30) 
                return result;
            } 
          } 
        }
}
    
    } while (i++ !=  _stop);

  }
  result = SYSTEM_true;
  if (_P3strcmpE(msg,_P3str1("\000"))) 
    _P3strcpy(msg,255,_P3str1("\020No license found"));
  return result;
}  /* licensechecksubinternal */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensecheck(
  PALDOORG_tpalobject self,
  SYSTEM_integer m,
  SYSTEM_integer n,
  SYSTEM_integer nz,
  SYSTEM_integer nlnz,
  SYSTEM_integer ndisc)
{
  SYSTEM_boolean result;
  SYSTEM_shortstring msg;

  self->PALDOORG_tpalobject_DOT_rowcnt = m;
  self->PALDOORG_tpalobject_DOT_colcnt = n;
  self->PALDOORG_tpalobject_DOT_nzcnt = nz;
  self->PALDOORG_tpalobject_DOT_nlnzcnt = nlnz;
  self->PALDOORG_tpalobject_DOT_disccnt = ndisc;
  if (self->PALDOORG_tpalobject_DOT_checksum != 0 && 
    PALDOORG_tpalobject_DOT_tampercheck(self)) {
    PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\035*** Module has been modified."));
    PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\052*** Continue to run in demonstration mode."));
    PALDOORG_tpalobject_DOT_makedemolicense(self);
  } 
  if (self->PALDOORG_tpalobject_DOT_licenselevel > 0 && 
    PALDOORG_tpalobject_DOT_licensecheckinternal(self)) {
    PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\043*** License file validation failed."));
    PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\052*** Continue to run in demonstration mode."));
    PALDOORG_tpalobject_DOT_makedemolicense(self);
  } 
  if (self->PALDOORG_tpalobject_DOT_isglobal == 1 && (self->
    PALDOORG_tpalobject_DOT_rowcnt > 10 || self->
    PALDOORG_tpalobject_DOT_colcnt > 10) || self->
    PALDOORG_tpalobject_DOT_isglobal == 0 && (self->
    PALDOORG_tpalobject_DOT_disccnt > 50 || self->
    PALDOORG_tpalobject_DOT_nlnzcnt > 1000 || self->
    PALDOORG_tpalobject_DOT_nzcnt > 2000 || self->
    PALDOORG_tpalobject_DOT_rowcnt > 300 || self->
    PALDOORG_tpalobject_DOT_colcnt > 300)) {
    {
      SYSTEM_shortstring _t1;

      if (_P3strcmpE(SYSTEM_copy(_t1,255,self->
        PALDOORG_tpalobject_DOT_subsyscode,1,2),_P3str1("\00200")) && 
        self->PALDOORG_tpalobject_DOT_licenselevel == 0) {
        result = SYSTEM_true;
        _P3strcpy(msg,255,_P3str1("\027Model exceeds demo size"));
      } else 
        result = PALDOORG_tpalobject_DOT_licensechecksubinternal(self,
          msg,self->PALDOORG_tpalobject_DOT_subsyssecondary,self->
          PALDOORG_tpalobject_DOT_subsyscode);
    }
    if (result) 
      if (_P3strcmpN(msg,_P3str1("\000"))) 
        {
          _P3STR_255 _t1;

          PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t1,255,_P3str1("\004*** "),
            msg));
        }
    if (result) {
      PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\074*** To update your license, please contact your distributor."));
      if (self->PALDOORG_tpalobject_DOT_isglobal == 1) {
        if (self->PALDOORG_tpalobject_DOT_rowcnt > 10) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\053*** Too many equations for global demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of equations          : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_rowcnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demoglobalmaxrow)));
          }
        } 
        if (self->PALDOORG_tpalobject_DOT_colcnt > 10) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\053*** Too many variables for global demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of variables          : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_colcnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demoglobalmaxcol)));
          }
        } 
      } else {
        if (self->PALDOORG_tpalobject_DOT_rowcnt > 300) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\044*** Too many equations for demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of equations          : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_rowcnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demomaxrow)));
          }
        } 
        if (self->PALDOORG_tpalobject_DOT_colcnt > 300) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\044*** Too many variables for demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of variables          : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_colcnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demomaxcol)));
          }
        } 
        if (self->PALDOORG_tpalobject_DOT_nzcnt > 2000) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\053*** Too many nonzero elements for demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of nonzeros           : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_nzcnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demomaxnz)));
          }
        } 
        if (self->PALDOORG_tpalobject_DOT_nlnzcnt > 1000) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\065*** Too many nonlinear nonzero elements for demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of nonlinear nonzeros : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_nlnzcnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demomaxnlnz)));
          }
        } 
        if (self->PALDOORG_tpalobject_DOT_disccnt > 50) {
          PALDOORG_tpalobject_DOT_msgadd(self,_P3str1("\055*** Too many discrete variables for demo mode"));
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Number of discrete variables : "),
              STRUTILX_inttostr(_t1,255,self->
              PALDOORG_tpalobject_DOT_disccnt)));
          }
          {
            SYSTEM_shortstring _t1;
            _P3STR_255 _t2;

            PALDOORG_tpalobject_DOT_msgadd(self,_P3strcat(_t2,255,_P3str1("\045***   Maximum allowed              : "),
              STRUTILX_inttostr(_t1,255,PALDOORG_demomaxdisc)));
          }
        } 
      } 
    } 
  } else 
    result = SYSTEM_false;
  return result;
}  /* pallicensecheck */

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicensechecksubsys(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *codes)
{
  SYSTEM_boolean result;
  SYSTEM_shortstring msg;

  result = PALDOORG_tpalobject_DOT_licensechecksubinternal(self,msg,
    SYSTEM_ord(codes[0]) /  2,codes);
  return result;
}  /* pallicensechecksubsys */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensegetmessage(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *msg)
{
  SYSTEM_boolean result;

  if (self->PALDOORG_tpalobject_DOT_ml == NULL) {
    _P3strclr(msg);
    result = SYSTEM_false;
    return result;
  } 
  GMSOBJ_txstrings_DOT_get(msg,255,ValueCast(GMSOBJ_txstrings,self->
    PALDOORG_tpalobject_DOT_ml),0);
  GMSOBJ_txlist_DOT_delete(ValueCast(GMSOBJ_txlist,self->
    PALDOORG_tpalobject_DOT_ml),0);
  if ((ValueCast(GMSOBJ_txstrings,self->PALDOORG_tpalobject_DOT_ml))->
    GMSOBJ_txlist_DOT_fcount == 0) 
    SYSUTILS_P3_freeandnil(&PointerCast(GMSOBJ_txstrings,&self->
      PALDOORG_tpalobject_DOT_ml));
  result = SYSTEM_true;
  return result;
}  /* pallicensegetmessage */

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicenseisdemocheckout(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  if (self->PALDOORG_tpalobject_DOT_isglobal == 1) { 
    result = !(self->PALDOORG_tpalobject_DOT_rowcnt > 10 || self->
      PALDOORG_tpalobject_DOT_colcnt > 10);
  } else 
    result = !(self->PALDOORG_tpalobject_DOT_disccnt > 50 || self->
      PALDOORG_tpalobject_DOT_nlnzcnt > 1000 || self->
      PALDOORG_tpalobject_DOT_nzcnt > 2000 || self->
      PALDOORG_tpalobject_DOT_rowcnt > 300 || self->
      PALDOORG_tpalobject_DOT_colcnt > 300);
  return result;
}  /* pallicenseisdemocheckout */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicenseisacademic(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  result = self->PALDOORG_tpalobject_DOT_license1[60] == _P3char('A');
  return result;
}  /* pallicenseisacademic */

Procedure PALDOORG_tpalobject_DOT_pallicensedemolimits(
  PALDOORG_tpalobject self,
  SYSTEM_integer *m,
  SYSTEM_integer *n,
  SYSTEM_integer *nz,
  SYSTEM_integer *nlnz,
  SYSTEM_integer *ndisc)
{
  *m = PALDOORG_demomaxrow;
  *n = PALDOORG_demomaxcol;
  *nz = PALDOORG_demomaxnz;
  *nlnz = PALDOORG_demomaxnlnz;
  *ndisc = PALDOORG_demomaxdisc;
}  /* pallicensedemolimits */

Procedure PALDOORG_tpalobject_DOT_pallicenseglobaldemolimits(
  PALDOORG_tpalobject self,
  SYSTEM_integer *m,
  SYSTEM_integer *n)
{
  *m = PALDOORG_demoglobalmaxrow;
  *n = PALDOORG_demoglobalmaxcol;
}  /* pallicenseglobaldemolimits */

Procedure PALDOORG_tpalobject_DOT_pallicensesetsubsearch(
  PALDOORG_tpalobject self)
{
  self->PALDOORG_tpalobject_DOT_currentsearchhelper = 0;
}  /* pallicensesetsubsearch */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensegetsubnext(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  if (self->PALDOORG_tpalobject_DOT_currentsearchhelper > self->
    PALDOORG_tpalobject_DOT_licenseactsub) 
    return result;
  _P3inc0(self->PALDOORG_tpalobject_DOT_currentsearchhelper);
  if (self->PALDOORG_tpalobject_DOT_license3[7 + 2 * self->
    PALDOORG_tpalobject_DOT_currentsearchhelper] == _P3char(' ')) { 
    self->PALDOORG_tpalobject_DOT_currentsearchhelper = self->
      PALDOORG_tpalobject_DOT_licenseactsub + 1;
  } else 
    result = SYSTEM_true;
  return result;
}  /* pallicensegetsubnext */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetjullice(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  {
    SYSTEM_shortstring _t1;

    result = PALDOORG_tpalobject_DOT_palgetjuliandays(self,SYSTEM_copy(
      _t1,255,self->PALDOORG_tpalobject_DOT_license1,49,6));
  }
  return result;
}  /* pallicensegetjullice */

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetevaldate(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  if (self->PALDOORG_tpalobject_DOT_licenseversion < 2) { 
    {
      SYSTEM_shortstring _t1;

      result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
        PALDOORG_tpalobject_DOT_license3,65,1));
    }
  } else 
    {
      SYSTEM_shortstring _t1;

      result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
        PALDOORG_tpalobject_DOT_license4,10,1));
    }
  if (result == 0) { 
    result = SYSTEM_maxlongint;
  } else 
    result = PALDOORG_tpalobject_DOT_pallicensegetjullice(self) + 30 * 
      result;
  return result;
}  /* pallicensegetevaldate */

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetmaintdate(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  switch (self->PALDOORG_tpalobject_DOT_licenseversion) {
    case 1: 
      {
        SYSTEM_shortstring _t1;

        result = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,self->
          PALDOORG_tpalobject_DOT_license4,9,2));
      }
      break;
    case 2: 
      {
        SYSTEM_shortstring _t1;

        result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
          PALDOORG_tpalobject_DOT_license4,9,1));
      }
      break;
    default:
      result = 0;
  }
  if (result == 0) { 
    result = SYSTEM_maxlongint;
  } else 
    result = PALDOORG_tpalobject_DOT_pallicensegetjulbase(self) + 30 * 
      result;
  return result;
}  /* pallicensegetmaintdate */

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetsubmaintdate(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  result = PALDOORG_tpalobject_DOT_pallicensegetjulbase(self) + 30 * 
    PALDOORG_tpalobject_DOT_pallicensegetsubmaint(self);
  return result;
}  /* pallicensegetsubmaintdate */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetjulbase(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  {
    SYSTEM_shortstring _t1;

    result = PALDOORG_tpalobject_DOT_palgetjuliandays(self,SYSTEM_copy(
      _t1,255,self->PALDOORG_tpalobject_DOT_license1,49,6));
  }
  switch (self->PALDOORG_tpalobject_DOT_licenseversion) {
    case 2: 
      {
        SYSTEM_shortstring _t1;

        result = result - 30 * (SYSUTILS_P3_strtoint(SYSTEM_copy(
          _t1,255,self->PALDOORG_tpalobject_DOT_license1,56,2)) - 1);
      }
      break;
    default: break;
  }
  return result;
}  /* pallicensegetjulbase */

Procedure PALDOORG_tpalobject_DOT_pallicenseclear(
  PALDOORG_tpalobject self)
{
  self->PALDOORG_tpalobject_DOT_licensestatus = 1;
  self->PALDOORG_tpalobject_DOT_licenselevel =  -1;
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license1,65,_P3str1("\020license1 not set"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license2,65,_P3str1("\020license2 not set"));
  _P3strclr(self->PALDOORG_tpalobject_DOT_license3);
  _P3strclr(self->PALDOORG_tpalobject_DOT_license4);
  _P3strclr(self->PALDOORG_tpalobject_DOT_license5);
  _P3strclr(self->PALDOORG_tpalobject_DOT_license6);
  self->PALDOORG_tpalobject_DOT_currentsearchhelper = 29;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
}  /* pallicenseclear */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensechecksize(
  PALDOORG_tpalobject self,
  SYSTEM_integer m,
  SYSTEM_integer n,
  SYSTEM_integer nz,
  SYSTEM_integer nlnz,
  SYSTEM_integer ndisc)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  if (self->PALDOORG_tpalobject_DOT_licenselevel == 0) 
    if (ndisc > 50 || nlnz > 1000 || nz > 2000 || m > 300 || 
      n > 300) 
      result = SYSTEM_true;
  return result;
}  /* pallicensechecksize */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensereadw(
  PALDOORG_tpalobject self,
  SYSTEM_text *f)
{
  SYSTEM_integer result;

  result = 0;
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3read_fs0(self->PALDOORG_tpalobject_DOT_license1,65);
    _P3readlf();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3read_fs0(self->PALDOORG_tpalobject_DOT_license2,65);
    _P3readlf();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3read_fs0(self->PALDOORG_tpalobject_DOT_license3,65);
    _P3readlf();
  }
  _Iplus_end();
  if (self->PALDOORG_tpalobject_DOT_license1[55] != _P3char('-')) {
    _Iplus_bgn();
    {
      _P3file_ptr _file_temp = f;

      _P3read_fs0(self->PALDOORG_tpalobject_DOT_license4,65);
      _P3readlf();
    }
    _Iplus_end();
    _Iplus_bgn();
    {
      _P3file_ptr _file_temp = f;

      _P3read_fs0(self->PALDOORG_tpalobject_DOT_license5,65);
      _P3readlf();
    }
    _Iplus_end();
  } 
  return result;
}  /* pallicensereadw */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensescipw(
  PALDOORG_tpalobject self,
  SYSTEM_text *f)
{
  SYSTEM_integer result;
  SYSTEM_shortstring temp;

  result = 0;
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3read_fs0(temp,255);
    _P3readlf();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3readlf();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3readlf();
  }
  _Iplus_end();
  if (temp[55] != _P3char('-')) {
    _Iplus_bgn();
    {
      _P3file_ptr _file_temp = f;

      _P3readlf();
    }
    _Iplus_end();
    _Iplus_bgn();
    {
      _P3file_ptr _file_temp = f;

      _P3readlf();
    }
    _Iplus_end();
  } 
  return result;
}  /* pallicensescipw */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensedisplay(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self,
  SYSTEM_integer i)
{
  switch (i) {
    case 1: 
      _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license1);
      break;
    case 2: 
      _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license2);
      break;
    case 3: 
      _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license3);
      break;
    case 4: 
      _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license4);
      break;
    case 5: 
      _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license5);
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\041**** Unknown license display line"));
  }
  return result;
}  /* pallicensedisplay */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensegetdc(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_integer i;

  i = SYSTEM_pos(_P3str1("\001 "),self->
    PALDOORG_tpalobject_DOT_license5);
  if (i > 0) { 
    SYSTEM_copy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license5,1,
      i - 1);
  } else 
    _P3strclr(result);
  return result;
}  /* pallicensegetdc */

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicensegetacademic(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  result = self->PALDOORG_tpalobject_DOT_license1[60] == _P3char('A');
  return result;
}  /* pallicensegetacademic */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensemessage(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  _P3strcpy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license6);
  return result;
}  /* pallicensemessage */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palolderlicense3(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  {
    SYSTEM_shortstring _t1;

    _P3strcat(result,_len_ret,STRUTILX_inttostr(_t1,255,
      SYSTEM_round(10000000 + (12 * SYSUTILS_P3_now() - 1))),_P3str1("\071000S                                                   00"));
  }
  return result;
}  /* palolderlicense3 */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicenseexist(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  result = self->PALDOORG_tpalobject_DOT_licenselevel >= 0;
  return result;
}  /* pallicenseexist */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensewritec1(
  PALDOORG_tpalobject self,
  SYSTEM_text *f)
{
  SYSTEM_integer result;

  result = 0;
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3writefs0(self->PALDOORG_tpalobject_DOT_license1);
    _P3writefc0(_P3char('0'));
    _P3writefn();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3writefs0(self->PALDOORG_tpalobject_DOT_license2);
    _P3writefc0(_P3char('0'));
    _P3writefn();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3writefs0(self->PALDOORG_tpalobject_DOT_license3);
    _P3writefc0(_P3char('0'));
    _P3writefn();
  }
  _Iplus_end();
  return result;
}  /* pallicensewritec1 */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensewritec2(
  PALDOORG_tpalobject self,
  SYSTEM_text *f)
{
  SYSTEM_integer result;

  result = 0;
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3writefs0(self->PALDOORG_tpalobject_DOT_license4);
    _P3writefc0(_P3char('0'));
    _P3writefn();
  }
  _Iplus_end();
  _Iplus_bgn();
  {
    _P3file_ptr _file_temp = f;

    _P3writefs0(self->PALDOORG_tpalobject_DOT_license5);
    _P3writefc0(_P3char('0'));
    _P3writefn();
  }
  _Iplus_end();
  return result;
}  /* pallicensewritec2 */

Procedure PALDOORG_tpalobject_DOT_pallicensedemo(
  PALDOORG_tpalobject self)
{
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license1,65,_P3str1("\101GAMS Development Corporation, Washington, DC   G871201/0000CA-ANY"));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license2,65,_P3str1("\101Free Demo,  202-342-0180,  sales@gams.com,  www.gams.com         "));
  {
    SYSTEM_shortstring _t1;

    _P3strcat(self->PALDOORG_tpalobject_DOT_license3,65,
      STRUTILX_inttostr(_t1,255,SYSTEM_round(10000000 + 12 * 
      SYSUTILS_P3_now())),_P3str1("\071000S                                                   00"));
  }
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license4,65,_P3str1("\10101234567000000                                                   "));
  _P3strcpy(self->PALDOORG_tpalobject_DOT_license5,65,_P3str1("\101DC0000       Ref: Generated by Base            A Demo            "));
  _P3strclr(self->PALDOORG_tpalobject_DOT_license6);
  self->PALDOORG_tpalobject_DOT_licensestatus = 1;
  self->PALDOORG_tpalobject_DOT_licenselevel = 0;
  self->PALDOORG_tpalobject_DOT_licenseversion = 2;
  self->PALDOORG_tpalobject_DOT_licenseactsub = PALDOORG_licensemaxsub;
  self->PALDOORG_tpalobject_DOT_currentsearchhelper = self->
    PALDOORG_tpalobject_DOT_licenseactsub + 1;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
  self->PALDOORG_tpalobject_DOT_optionsearchhelper = 0;
}  /* pallicensedemo */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_paltampercheck(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;

  result = PALDOORG_tpalobject_DOT_tampercheck(self);
  return result;
}  /* paltampercheck */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetlevel(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  {
    SYSTEM_shortstring _t1;

    result = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,self->
      PALDOORG_tpalobject_DOT_license3,9,2));
  }
  return result;
}  /* pallicensegetlevel */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetinstitution(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSUTILS_P3_trim(result,_len_ret,self->
    PALDOORG_tpalobject_DOT_license2);
  return result;
}  /* pallicensegetinstitution */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetinstdc(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_shortstring temp;
  static _P3STR_95 blank = {64,' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '};

  PALDOORG_tpalobject_DOT_pallicensegetdc(temp,255,self);
  {
    SYSTEM_shortstring _t1;
    SYSTEM_shortstring _t2;
    _P3STR_255 _t3;

    _P3strcat(result,_len_ret,SYSTEM_copy(_t1,255,_P3strcat(_t3,255,
      PALDOORG_tpalobject_DOT_pallicensegetinstitution(_t2,255,self),
      blank),1,65 - SYSTEM_length(temp)),temp);
  }
  return result;
}  /* pallicensegetinstdc */

Function(SYSTEM_longint ) PALDOORG_tpalobject_DOT_pallicensegetkey(
  PALDOORG_tpalobject self)
{
  SYSTEM_longint result;
  SYSTEM_integer i;

  {
    SYSTEM_shortstring _t1;

    _P3val_i(SYSTEM_copy(_t1,255,self->
      PALDOORG_tpalobject_DOT_license3,1,8),result,&i);
  }
  if (i != 0) 
    result = 0;
  return result;
}  /* pallicensegetkey */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palsecuritycheck(
  PALDOORG_tpalobject self)
{
  SYSTEM_boolean result;
  SYSTEM_shortstring msg;

  result = PALDOORG_tpalobject_DOT_licensechecksubinternal(self,msg,1,_P3str1("\0020S"));
  return result;
}  /* palsecuritycheck */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_licensegetsubstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_copy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license3,7 + 2 * 
    self->PALDOORG_tpalobject_DOT_currentsearchhelper,2);
  return result;
}  /* licensegetsubstring */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_licensegetsubeval(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  if (self->PALDOORG_tpalobject_DOT_licenseversion < 2) { 
    result = SYSTEM_ord(self->PALDOORG_tpalobject_DOT_license3[66 - 
      self->PALDOORG_tpalobject_DOT_currentsearchhelper]) - 48;
  } else 
    {
      SYSTEM_shortstring _t1;

      result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
        PALDOORG_tpalobject_DOT_license4,7 + 2 * self->
        PALDOORG_tpalobject_DOT_currentsearchhelper + 1,1));
    }
  return result;
}  /* licensegetsubeval */

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_licensegetsubmaint(
  PALDOORG_tpalobject self)
{
  SYSTEM_integer result;

  result = 0;
  switch (self->PALDOORG_tpalobject_DOT_licenseversion) {
    case 1: 
      {
        SYSTEM_shortstring _t1;

        result = SYSUTILS_P3_strtoint(SYSTEM_copy(_t1,255,self->
          PALDOORG_tpalobject_DOT_license4,7 + 2 * self->
          PALDOORG_tpalobject_DOT_currentsearchhelper,2));
      }
      break;
    case 2: 
      {
        SYSTEM_shortstring _t1;

        result = PALDOORG_lnumtoint(SYSTEM_copy(_t1,255,self->
          PALDOORG_tpalobject_DOT_license4,7 + 2 * self->
          PALDOORG_tpalobject_DOT_currentsearchhelper,1));
      }
      break;
    default: break;
  }
  return result;
}  /* licensegetsubmaint */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensegetdates(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *lcode,
  SYSTEM_longint *eval,
  SYSTEM_longint *maint)
{
  SYSTEM_boolean result;

  while (PALDOORG_tpalobject_DOT_pallicensegetsubnext(self)) {
    result = SYSTEM_true;
    PALDOORG_tpalobject_DOT_licensegetsubstring(lcode,255,self);
    *eval = PALDOORG_tpalobject_DOT_licensegetsubeval(self);
    *maint = PALDOORG_tpalobject_DOT_licensegetsubmaint(self);
    if (*eval == 0) {
      *eval = SYSTEM_maxlongint;
      if (*maint == 0) { 
        *maint = SYSTEM_maxlongint;
      } else 
        *maint = PALDOORG_tpalobject_DOT_pallicensegetjulbase(self) + 30 * *
          maint;
    } else {
      *eval = PALDOORG_tpalobject_DOT_pallicensegetjullice(self) + 30 * *
        eval;
      *maint = *eval;
    } 
    return result;
  
}
  result = SYSTEM_false;
  _P3strclr(lcode);
  *eval = 0;
  *maint = 0;
  return result;
}  /* pallicensegetdates */

Function(SYSTEM_longint ) 
  PALDOORG_tpalobject_DOT_pallicensegetclipcomponent(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *lcode)
{
  SYSTEM_longint result;
  SYSTEM_longint i, k;

  i = GMSGLOBX_clipcodeslookup(lcode);
  result = 0;
  for (k = 1;k <= (SYSTEM_int32)GMSGLOBX_maxcomponentsolvermap;++
    k) {
    if (i == GMSGLOBX_clipcomponentmapmap(k,1)) {
      result = GMSGLOBX_clipcomponentmapmap(k,2);
      return result;
    } 
  }
  return result;
}  /* pallicensegetclipcomponent */

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensegetid(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_copy(result,_len_ret,self->PALDOORG_tpalobject_DOT_license1,48,18);
  return result;
}  /* pallicensegetid */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetlicensee(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  {
    SYSTEM_shortstring _t2;

    SYSUTILS_P3_trim(result,_len_ret,SYSTEM_copy(_t2,255,self->
      PALDOORG_tpalobject_DOT_license1,1,47));
  }
  return result;
}  /* pallicensegetlicensee */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetleveltext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  GMSGLOBX_gamslicensestext(result,_len_ret,self->
    PALDOORG_tpalobject_DOT_licenselevel + 1);
  return result;
}  /* pallicensegetleveltext */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetvendor(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  {
    SYSTEM_shortstring _t2;

    GMSGLOBX_vendorstext(result,_len_ret,GMSGLOBX_vendorslookup(
      SYSTEM_copy(_t2,255,self->PALDOORG_tpalobject_DOT_license1,48,1)));
  }
  return result;
}  /* pallicensegetvendor */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetplatformtext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_shortstring xxx;

  PALDOORG_tpalobject_DOT_pallicensegetplatform(xxx,255,self);
  if (_P3strcmpE(xxx,_P3str1("\003ANY"))) { 
    _P3strcpy(result,_len_ret,_P3str1("\030Any platform demo system"));
  } else 
    GMSGLOBX_platformstext(result,_len_ret,GMSGLOBX_platformslookup(
      xxx));
  return result;
}  /* pallicensegetplatformtext */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetmudtext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  SYSTEM_integer j, k;

  j = SYSTEM_pos(_P3str1("\003 g "),self->
    PALDOORG_tpalobject_DOT_license5);
  if (j == 0) 
    j = SYSTEM_pos(_P3str1("\003_g_"),self->
      PALDOORG_tpalobject_DOT_license5);
  if (j == 0) { 
    _P3strclr(result);
  } else {
    {
      SYSTEM_shortstring _t1;

      k = GMSGLOBX_gamslicensetypeslookup(SYSTEM_copy(_t1,255,self->
        PALDOORG_tpalobject_DOT_license5,j + 3,1));
    }
    if (k == 0) { 
      _P3strclr(result);
    } else 
      GMSGLOBX_gamslicensetypestext(result,_len_ret,k);
  } 
  return result;
}  /* pallicensegetmudtext */

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensestatusmessage(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self)
{
  switch (self->PALDOORG_tpalobject_DOT_licensestatus) {
    case 0: 
      _P3strcpy(result,_len_ret,_P3str1("\035license validated succesfully"));
      break;
    case 1: 
      _P3strcpy(result,_len_ret,_P3str1("\026license status unknown"));
      break;
    case 2: 
      _P3strcpy(result,_len_ret,_P3str1("\045could not open specified license file"));
      break;
    case 3: 
      _P3strcpy(result,_len_ret,_P3str1("\031license validation failed"));
      break;
    case 4: 
      _P3strcpy(result,_len_ret,_P3str1("\023tamper check failed"));
      break;
    case 5: 
      _P3strcpy(result,_len_ret,_P3str1("\033license file format too old"));
      break;
    case 6: 
      _P3strcpy(result,_len_ret,_P3str1("\036evaluation license has expired"));
      break;
    case 7: 
      _P3strcpy(result,_len_ret,_P3str1("\027maintenance has expired"));
      break;
    default:
      _P3strcpy(result,_len_ret,_P3str1("\030*** no message found ***"));
  }
  return result;
}  /* pallicensestatusmessage */

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensechecksubx(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *sname,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_integer *daysleft)
{
  SYSTEM_shortstring codes;
  SYSTEM_boolean result;
  SYSTEM_shortstring msg;
  SYSTEM_integer eval;

  _P3strcpy(codes,255,_ftmp1);
  result = SYSTEM_true;
  if (ValueCast(SYSTEM_int32,SYSTEM_length(codes)) % 2 != 0 || 
    SYSTEM_length(codes) == 0) 
    return result;
  if (self->PALDOORG_tpalobject_DOT_licenselevel == 0) 
    return result;
  if ((SYSUTILS_P3_sametext(_P3str1("\006soplex"),sname) || 
    SYSUTILS_P3_sametext(_P3str1("\011osisoplex"),sname) || 
    SYSUTILS_P3_sametext(_P3str1("\004scip"),sname) || 
    SYSUTILS_P3_sametext(_P3str1("\010coinscip"),sname)) && 
    PALDOORG_tpalobject_DOT_pallicensegetacademic(self)) 
    _P3strcpy(codes,255,PALDOORG_cmexliccodes);
  if (!PALDOORG_tpalobject_DOT_licensechecksubinternal(self,msg,ValueCast(
    SYSTEM_int32,SYSTEM_length(codes)) /  2,codes)) {
    result = SYSTEM_false;
    eval = PALDOORG_tpalobject_DOT_pallicensegetsubeval(self);
    if (eval > 0) { 
      *daysleft = -self->PALDOORG_tpalobject_DOT_juliantoday + 
        PALDOORG_tpalobject_DOT_pallicensegetjullice(self) + 30 * 
        eval;
    } else 
      *daysleft = 0;
  } 
  return result;
}  /* pallicensechecksubx */

/* unit paldoorg */
void _Init_Module_paldoorg(void)
{
  PALDOORG_dlice();
} /* _Init_Module_paldoorg */

void _Final_Module_paldoorg(void)
{
} /* _Final_Module_paldoorg */

