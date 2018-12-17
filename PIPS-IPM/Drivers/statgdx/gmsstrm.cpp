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
#include "clibtypes.h"
#include "gmslibname.h"
#include "xcompress.h"
#include "gmsstrm.h"

SYSTEM_integer GMSSTRM_buffersize = 32768;
_arr_0GMSSTRM GMSSTRM_rwtypetext = {{4,'B','y','t','e'}, {4,'B','o','o','l'}, {4,'C','h','a','r'}, {4,'W','o','r','d'}, {7,'I','n','t','e','g','e','r'}, {5,'I','n','t','6','4'}, {6,'D','o','u','b','l','e'}, {6,'S','t','r','i','n','g'}, {5,'P','C','h','a','r'}, {7,'P','S','t','r','i','n','g'}};

void * const GMSSTRM_txstream_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy, (void*)&_P3_abstract_call1, (void*)&
  _P3_abstract_call1, (void*)&_P3_abstract_call1, (void*)&
  _P3_abstract_call1, (void*)&_P3_abstract_call1, (void*)&
  GMSSTRM_txstream_DOT_readdouble, (void*)&
  GMSSTRM_txstream_DOT_readinteger, (void*)&
  GMSSTRM_txstream_DOT_readword, (void*)&
  GMSSTRM_txstream_DOT_readint64};

/* Class descriptor for 'txstream' */
const SYSTEM_classdescriptor_t GMSSTRM_txstream_CD = {
  _P3str1("\010txstream"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSSTRM_txstream_OD), GMSSTRM_txstream_VT, NULL};


void * const GMSSTRM_txfilestream_VT[] = {(void*)&
  GMSSTRM_txfilestream_DOT_destroy, (void*)&
  GMSSTRM_txfilestream_DOT_getposition, (void*)&
  GMSSTRM_txfilestream_DOT_setposition, (void*)&
  GMSSTRM_txfilestream_DOT_getsize, (void*)&
  GMSSTRM_txfilestream_DOT_read, (void*)&
  GMSSTRM_txfilestream_DOT_write, (void*)&
  GMSSTRM_txstream_DOT_readdouble, (void*)&
  GMSSTRM_txstream_DOT_readinteger, (void*)&
  GMSSTRM_txstream_DOT_readword, (void*)&
  GMSSTRM_txstream_DOT_readint64};

/* Class descriptor for 'txfilestream' */
const SYSTEM_classdescriptor_t GMSSTRM_txfilestream_CD = {
  _P3str1("\014txfilestream"), 
  &GMSSTRM_txstream_CD, NULL, 0, 
  sizeof(GMSSTRM_txfilestream_OD), GMSSTRM_txfilestream_VT, NULL};


void * const GMSSTRM_tbufferedfilestream_VT[] = {(void*)&
  GMSSTRM_tbufferedfilestream_DOT_destroy, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_getposition, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_setposition, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_getsize, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_read, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_write, (void*)&
  GMSSTRM_txstream_DOT_readdouble, (void*)&
  GMSSTRM_txstream_DOT_readinteger, (void*)&
  GMSSTRM_txstream_DOT_readword, (void*)&
  GMSSTRM_txstream_DOT_readint64};

/* Class descriptor for 'tbufferedfilestream' */
const SYSTEM_classdescriptor_t GMSSTRM_tbufferedfilestream_CD = {
  _P3str1("\023tbufferedfilestream"), 
  &GMSSTRM_txfilestream_CD, NULL, 0, 
  sizeof(GMSSTRM_tbufferedfilestream_OD), 
    GMSSTRM_tbufferedfilestream_VT, NULL};


void * const GMSSTRM_tmibufferedstream_VT[] = {(void*)&
  GMSSTRM_tbufferedfilestream_DOT_destroy, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_getposition, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_setposition, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_getsize, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_read, (void*)&
  GMSSTRM_tbufferedfilestream_DOT_write, (void*)&
  GMSSTRM_tmibufferedstream_DOT_readdouble, (void*)&
  GMSSTRM_tmibufferedstream_DOT_readinteger, (void*)&
  GMSSTRM_tmibufferedstream_DOT_readword, (void*)&
  GMSSTRM_tmibufferedstream_DOT_readint64};

/* Class descriptor for 'tmibufferedstream' */
const SYSTEM_classdescriptor_t GMSSTRM_tmibufferedstream_CD = {
  _P3str1("\021tmibufferedstream"), 
  &GMSSTRM_tbufferedfilestream_CD, NULL, 0, 
  sizeof(GMSSTRM_tmibufferedstream_OD), GMSSTRM_tmibufferedstream_VT, NULL};


void * const GMSSTRM_tgzipinputstream_VT[] = {(void*)&
  GMSSTRM_tgzipinputstream_DOT_destroy};

/* Class descriptor for 'tgzipinputstream' */
const SYSTEM_classdescriptor_t GMSSTRM_tgzipinputstream_CD = {
  _P3str1("\020tgzipinputstream"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSSTRM_tgzipinputstream_OD), GMSSTRM_tgzipinputstream_VT, NULL};


void * const GMSSTRM_tbinarytextfileio_VT[] = {(void*)&
  GMSSTRM_tbinarytextfileio_DOT_destroy};

/* Class descriptor for 'tbinarytextfileio' */
const SYSTEM_classdescriptor_t GMSSTRM_tbinarytextfileio_CD = {
  _P3str1("\021tbinarytextfileio"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSSTRM_tbinarytextfileio_OD), GMSSTRM_tbinarytextfileio_VT, NULL};

cnstdef {GMSSTRM_paranoid = SYSTEM_false};
static SYSTEM_byte GMSSTRM_signature_header = 255;
static _P3STR_7 GMSSTRM_signature_gams = {6,'*','G','A','M','S','*'};
cnstdef {GMSSTRM_verify_offset = 100};
typedef SYSTEM_uint8 _sub_3GMSSTRM;
typedef SYSTEM_byte _arr_2GMSSTRM[8];
typedef struct GMSSTRM_tdoublevar_S {
  union{
    struct{
      SYSTEM_double v;
    } _c1;
    struct{
      _arr_2GMSSTRM va;
    } _c2;
  } _u;
} GMSSTRM_tdoublevar;


Procedure GMSSTRM_txstream_DOT_parwrite(
  GMSSTRM_txstream self,
  GMSSTRM_rwtype t)
{
  SYSTEM_byte b;

  b = SYSTEM_ord(t);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&b,1));
}  /* parwrite */

Procedure GMSSTRM_txstream_DOT_parcheck(
  GMSSTRM_txstream self,
  GMSSTRM_rwtype t)
{
  SYSTEM_byte b;
  SYSTEM_shortstring s;

  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&b,1));
  if (ValueCast(SYSTEM_int32,b) != SYSTEM_ord(t)) {
    {
      _P3STR_63 _t1;

      _P3strcat(s,255,_P3strcat(_t1,39,_P3str1("\040Stream check failed: Expected = "),
        GMSSTRM_rwtypetext[t]),_P3str1("\010 Read = "));
    }
    if (ValueCast(SYSTEM_int32,b) > 9) { 
      {
        _P3STR_255 _t1;
        SYSTEM_shortstring _t2;

        _P3strcat(s,255,_P3strcat(_t1,255,s,_P3str1("\003???")),
          SYSUTILS_P3_inttostr(_t2,255,b));
      }
    } else 
      _P3strcat(s,255,s,GMSSTRM_rwtypetext[ValueCast(GMSSTRM_rwtype,
        b)]);
    _P3_RAISE(ValueCast(SYSTEM_exception,SYSTEM_exception_DOT_create(ValueCast(
      SYSTEM_exception,_P3alloc_object(&SYSTEM_exception_CD)),s)));
  } 
}  /* parcheck */

Procedure GMSSTRM_txstream_DOT_writestring(
  GMSSTRM_txstream self,
  const SYSTEM_ansichar *s)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_string);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&s[0],ValueCast(
    SYSTEM_int32,SYSTEM_length(s)) + 1));
}  /* writestring */

Function(SYSTEM_ansichar *) GMSSTRM_txstream_DOT_readstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSSTRM_txstream self)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_string);
  if (VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&
    result[0],1)) == 0 || SYSTEM_ord(result[0]) == 0) { 
    _P3strclr(result);
  } else 
    VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result[1],
      SYSTEM_ord(result[0])));
  return result;
}  /* readstring */

Procedure GMSSTRM_txstream_DOT_writepstring(
  GMSSTRM_txstream self,
  SYSTEM_P3_pshortstring ps)
{
  SYSTEM_byte l;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_pstring);
  if (ps != NULL) { 
    VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&(*ps)[0],
      SYSTEM_ord((*ps)[0]) + 1));
  } else {
    l = 0;
    VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&l,1));
  } 
}  /* writepstring */

Procedure GMSSTRM_txstream_DOT_readpstring(
  GMSSTRM_txstream self,
  SYSTEM_P3_pshortstring *ps)
{
  SYSTEM_byte l;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_pstring);
  if (*ps != NULL) 
    _P3freemem(*ps);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&l,1));
  if (l == 0) { 
    *ps = NULL;
  } else {
    _P3getmem(*ps,ValueCast(SYSTEM_int32,l) + 1);
    VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&(**ps)[1],
      l));
    _P3setlength(**ps,l,255);
  } 
}  /* readpstring */

Procedure GMSSTRM_txstream_DOT_writedouble(
  GMSSTRM_txstream self,
  SYSTEM_double x)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_double);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&x,sizeof(
    SYSTEM_double)));
}  /* writedouble */

Function(SYSTEM_double ) GMSSTRM_txstream_DOT_readdouble(
  GMSSTRM_txstream self)
{
  SYSTEM_double result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_double);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_double)));
  return result;
}  /* readdouble */

Procedure GMSSTRM_txstream_DOT_writeinteger(
  GMSSTRM_txstream self,
  SYSTEM_integer n)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_integer);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&n,sizeof(
    SYSTEM_int32)));
}  /* writeinteger */

Function(SYSTEM_integer ) GMSSTRM_txstream_DOT_readinteger(
  GMSSTRM_txstream self)
{
  SYSTEM_integer result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_integer);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_int32)));
  return result;
}  /* readinteger */

Procedure GMSSTRM_txstream_DOT_writebyte(
  GMSSTRM_txstream self,
  SYSTEM_byte b)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_byte);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&b,sizeof(
    SYSTEM_uint8)));
}  /* writebyte */

Function(SYSTEM_byte ) GMSSTRM_txstream_DOT_readbyte(
  GMSSTRM_txstream self)
{
  SYSTEM_byte result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_byte);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_uint8)));
  return result;
}  /* readbyte */

Procedure GMSSTRM_txstream_DOT_writeword(
  GMSSTRM_txstream self,
  SYSTEM_word w)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_word);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&w,sizeof(
    SYSTEM_uint16)));
}  /* writeword */

Function(SYSTEM_word ) GMSSTRM_txstream_DOT_readword(
  GMSSTRM_txstream self)
{
  SYSTEM_word result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_word);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_uint16)));
  return result;
}  /* readword */

Procedure GMSSTRM_txstream_DOT_writeint64(
  GMSSTRM_txstream self,
  SYSTEM_int64 n)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_int64);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&n,sizeof(
    SYSTEM_int64)));
}  /* writeint64 */

Function(SYSTEM_int64 ) GMSSTRM_txstream_DOT_readint64(
  GMSSTRM_txstream self)
{
  SYSTEM_int64 result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_int64);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_int64)));
  return result;
}  /* readint64 */

Procedure GMSSTRM_txstream_DOT_writebool(
  GMSSTRM_txstream self,
  SYSTEM_boolean b)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_bool);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&b,sizeof(
    SYSTEM_boolean)));
}  /* writebool */

Procedure GMSSTRM_txstream_DOT_writepchar(
  GMSSTRM_txstream self,
  SYSTEM_P3_pansichar p,
  SYSTEM_integer l)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_pchar);
  GMSSTRM_txstream_DOT_writeinteger(self,l);
  if (l > 0) 
    VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,p,l));
}  /* writepchar */

Function(SYSTEM_boolean ) GMSSTRM_txstream_DOT_readbool(
  GMSSTRM_txstream self)
{
  SYSTEM_boolean result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_bool);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_boolean)));
  return result;
}  /* readbool */

Procedure GMSSTRM_txstream_DOT_writechar(
  GMSSTRM_txstream self,
  SYSTEM_ansichar c)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parwrite(self,GMSSTRM_rw_char);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_write_T, 5, (self,&c,sizeof(
    SYSTEM_ansichar)));
}  /* writechar */

Function(SYSTEM_ansichar ) GMSSTRM_txstream_DOT_readchar(
  GMSSTRM_txstream self)
{
  SYSTEM_ansichar result;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_char);
  VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,&result,sizeof(
    SYSTEM_ansichar)));
  return result;
}  /* readchar */

Procedure GMSSTRM_txstream_DOT_readpchar(
  GMSSTRM_txstream self,
  SYSTEM_P3_pansichar *p,
  SYSTEM_integer *l)
{
  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(self,GMSSTRM_rw_pchar);
  *l = VirtMethodCall(self, GMSSTRM_txstream_DOT_readinteger_T, 7, (
    self));
  if (*l <= 0) { 
    *p = NULL;
  } else {
    _P3getmem(*p,*l);
    VirtMethodCall(self, GMSSTRM_txstream_DOT_read_T, 4, (self,*p,*l));
  } 
}  /* readpchar */

Constructor(GMSSTRM_txfilestream ) GMSSTRM_txfilestream_DOT_create(
  GMSSTRM_txfilestream self,
  const SYSTEM_ansichar *afilename,
  SYSTEM_word amode)
{
  P3UTILS_tp3fileopenaction fmode;

  ValueCast(GMSSTRM_txfilestream,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  _P3strcpy(self->GMSSTRM_txfilestream_DOT_ffilename,255,afilename);
  GMSSTRM_txfilestream_DOT_setpassword(self,_P3str1("\000"));
  self->GMSSTRM_txfilestream_DOT_flastioresult = 0;
  switch (amode) {
    case GMSSTRM_fmcreate: 
      fmode = P3UTILS_p3openwrite;
      break;
    case GMSSTRM_fmopenwrite: 
      fmode = P3UTILS_p3openwrite;
      break;
    case GMSSTRM_fmopenread: 
      fmode = P3UTILS_p3openread;
      break;
    case GMSSTRM_fmopenreadwrite: 
      fmode = P3UTILS_p3openreadwrite;
      break;
    default:
      {
        SYSTEM_shortstring _t1;
        _P3STR_255 _t2;

        SYSTEM_assert(SYSTEM_false,_P3strcat(_t2,255,_P3str1("\026TXFileStream.Create = "),
          SYSUTILS_P3_inttostr(_t1,255,amode)));
      }
  }
  GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3fileopen(
    self->GMSSTRM_txfilestream_DOT_ffilename,fmode,&self->
    GMSSTRM_txfilestream_DOT_fs));
  self->GMSSTRM_txfilestream_DOT_fileisopen = self->
    GMSSTRM_txfilestream_DOT_flastioresult == 0;
  self->GMSSTRM_txfilestream_DOT_physposition = 0;
  return self;
}  /* create */

Destructor(GMSSTRM_txfilestream ) GMSSTRM_txfilestream_DOT_destroy(
  GMSSTRM_txfilestream self)
{
  if (self->GMSSTRM_txfilestream_DOT_fileisopen) 
    GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3fileclose(&
      self->GMSSTRM_txfilestream_DOT_fs));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSSTRM_txfilestream_DOT_applypassword(
  GMSSTRM_txfilestream self,
  GMSGEN_pbytedataarray pr,
  GMSGEN_pbytedataarray pw,
  SYSTEM_integer len,
  SYSTEM_int64 offs)
{
  SYSTEM_integer l, n;
  SYSTEM_integer fpwnxt;

  l = SYSTEM_length(self->GMSSTRM_txfilestream_DOT_fpassword);
  fpwnxt = offs % l;
  { register SYSTEM_int32 _stop = len - 1;
    if ((n = 0) <=  _stop) do {
      fpwnxt = fpwnxt + 1;
      if (fpwnxt > l) 
        fpwnxt = 1;
      (*pw)[n] = (*pr)[n] ^ SYSTEM_ord(self->
        GMSSTRM_txfilestream_DOT_fpassword[fpwnxt]);
    
    } while (n++ !=  _stop);

  }
}  /* applypassword */

Function(SYSTEM_longword ) GMSSTRM_txfilestream_DOT_read(
  GMSSTRM_txfilestream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count)
{
  SYSTEM_longword result;
  GMSGEN_pbytedataarray pw, pr;

  if (SYSTEM_length(self->GMSSTRM_txfilestream_DOT_fpassword) == 0) { 
    GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3fileread(
      self->GMSSTRM_txfilestream_DOT_fs,buffer,count,&result));
  } else {
    pw = ValueCast(GMSGEN_pbytedataarray,buffer);
    _P3getmem(pr,count);
    GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3fileread(
      self->GMSSTRM_txfilestream_DOT_fs,&(*pr)[0],count,&result));
    GMSSTRM_txfilestream_DOT_applypassword(self,pr,pw,count,self->
      GMSSTRM_txfilestream_DOT_physposition);
    _P3freemem(pr);
  } 
  _P3inc1(self->GMSSTRM_txfilestream_DOT_physposition,result);
  return result;
}  /* read */

Function(SYSTEM_longword ) GMSSTRM_txfilestream_DOT_write(
  GMSSTRM_txfilestream self,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword count)
{
  SYSTEM_longword result;
  GMSGEN_pbytedataarray pw, pr;

  if (SYSTEM_length(self->GMSSTRM_txfilestream_DOT_fpassword) == 0) { 
    GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3filewrite(
      self->GMSSTRM_txfilestream_DOT_fs,buffer,count,&result));
  } else {
    pr = ValueCast(GMSGEN_pbytedataarray,buffer);
    _P3getmem(pw,count);
    GMSSTRM_txfilestream_DOT_applypassword(self,pr,pw,count,self->
      GMSSTRM_txfilestream_DOT_physposition);
    GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3filewrite(
      self->GMSSTRM_txfilestream_DOT_fs,&(*pw)[0],count,&result));
    _P3freemem(pw);
  } 
  _P3inc1(self->GMSSTRM_txfilestream_DOT_physposition,result);
  return result;
}  /* write */

Procedure GMSSTRM_txfilestream_DOT_setpassword(
  GMSSTRM_txfilestream self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer k;
  SYSTEM_integer w;
  SYSTEM_byte b;
  SYSTEM_boolean bb;

  if (SYSTEM_length(s) == 0) { 
    _P3strclr(self->GMSSTRM_txfilestream_DOT_fpassword);
  } else {
    w = 0;
    bb = SYSTEM_false;
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((k = 1) <=  _stop) do {
        if (s[k] != _P3char(' ')) { 
          bb = SYSTEM_false;
        } else {
          if (bb) 
            SYSTEM_continue(CNT_1);
          bb = SYSTEM_true;
        } 
        b = SYSTEM_ord(s[k]);
        if ((ValueCast(SYSTEM_int32,b) & 1) == 0) { 
          b = b >> 1;
        } else 
          b = 128 + (b >> 1);
        w = w + 1;
        self->GMSSTRM_txfilestream_DOT_fpassword[w] = ValueCast(
          SYSTEM_ansichar,b);
      
CNT_1:;
      } while (k++ !=  _stop);
BRK_1:;

    }
    _P3setlength(self->GMSSTRM_txfilestream_DOT_fpassword,w,255);
  } 
}  /* setpassword */

Function(SYSTEM_boolean ) GMSSTRM_txfilestream_DOT_getusespassword(
  GMSSTRM_txfilestream self)
{
  SYSTEM_boolean result;

  result = SYSTEM_length(self->GMSSTRM_txfilestream_DOT_fpassword) > 0;
  return result;
}  /* getusespassword */

static Function(SYSTEM_ansichar ) randch(
  SYSTEM_integer *_2seed,
  GMSSTRM_txfilestream *_2self)
{
  SYSTEM_ansichar result;

  *_2seed = *_2seed * 12347 + 1023 & 134217727;
  result = ValueCast(SYSTEM_ansichar,*_2seed & 255);
  return result;
}  /* randch */

Function(SYSTEM_ansichar *) GMSSTRM_txfilestream_DOT_randstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSSTRM_txfilestream self,
  SYSTEM_integer l)
{
  SYSTEM_integer seed;
  SYSTEM_integer n;

  _P3setlength(result,l,255);
  seed = 1234 * l;
  { register SYSTEM_int32 _stop = l;
    if ((n = 1) <=  _stop) do {
      result[n] = randch(&seed,&self);
    } while (n++ !=  _stop);

  }
  return result;
}  /* randstring */

Function(SYSTEM_int64 ) GMSSTRM_txfilestream_DOT_getsize(
  GMSSTRM_txfilestream self)
{
  SYSTEM_int64 result;

  GMSSTRM_txfilestream_DOT_setlastioresult(self,P3UTILS_p3filegetsize(
    self->GMSSTRM_txfilestream_DOT_fs,&result));
  return result;
}  /* getsize */

Function(SYSTEM_int64 ) GMSSTRM_txfilestream_DOT_getposition(
  GMSSTRM_txfilestream self)
{
  SYSTEM_int64 result;

  result = self->GMSSTRM_txfilestream_DOT_physposition;
  return result;
}  /* getposition */

Procedure GMSSTRM_txfilestream_DOT_setposition(
  GMSSTRM_txfilestream self,
  SYSTEM_int64 p)
{
  SYSTEM_int64 np;

  self->GMSSTRM_txfilestream_DOT_physposition = p;
  GMSSTRM_txfilestream_DOT_setlastioresult(self,
    P3UTILS_p3filesetpointer(self->GMSSTRM_txfilestream_DOT_fs,p,&np,
    P3UTILS_p3_file_begin));
}  /* setposition */

Procedure GMSSTRM_txfilestream_DOT_setlastioresult(
  GMSSTRM_txfilestream self,
  SYSTEM_integer v)
{
  if (self->GMSSTRM_txfilestream_DOT_flastioresult == 0) 
    self->GMSSTRM_txfilestream_DOT_flastioresult = v;
}  /* setlastioresult */

Function(SYSTEM_integer ) GMSSTRM_txfilestream_DOT_getlastioresult(
  GMSSTRM_txfilestream self)
{
  SYSTEM_integer result;

  result = self->GMSSTRM_txfilestream_DOT_flastioresult;
  self->GMSSTRM_txfilestream_DOT_flastioresult = 0;
  return result;
}  /* getlastioresult */

Constructor(GMSSTRM_tbufferedfilestream ) 
  GMSSTRM_tbufferedfilestream_DOT_createwithpath(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode,
  const SYSTEM_ansichar *loadpath)
{
  SYSTEM_shortstring msg;

  ValueCast(GMSSTRM_tbufferedfilestream,
    GMSSTRM_txfilestream_DOT_create(ValueCast(GMSSTRM_txfilestream,
    self),filename,mode));
  GMSSTRM_tbufferedfilestream_DOT_setloadpath(self,loadpath);
  if (!XCOMPRESS_zlibdllloaded()) 
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;

      XCOMPRESS_loadzliblibrary(_P3strcat(_t2,255,
        GMSSTRM_tbufferedfilestream_DOT_getloadpath(_t1,255,self),_P3str1("\010gmszlib1")),
        msg);
    }
  self->GMSSTRM_tbufferedfilestream_DOT_fcancompress = 
    XCOMPRESS_zlibdllloaded();
  self->GMSSTRM_tbufferedfilestream_DOT_bufsize = GMSSTRM_buffersize;
  _P3getmem(self->GMSSTRM_tbufferedfilestream_DOT_bufptr,self->
    GMSSTRM_tbufferedfilestream_DOT_bufsize);
  self->GMSSTRM_tbufferedfilestream_DOT_cbufsize = SYSTEM_round(ValueCast(
    SYSTEM_double,self->GMSSTRM_tbufferedfilestream_DOT_bufsize * 12) /  10) + 20;
  _P3getmem(self->GMSSTRM_tbufferedfilestream_DOT_cbufptr,sizeof(
    GMSSTRM_tcompressheader) + self->
    GMSSTRM_tbufferedfilestream_DOT_cbufsize);
  self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_nrread = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_nrwritten = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_fcompress = SYSTEM_false;
  return self;
}  /* createwithpath */

Constructor(GMSSTRM_tbufferedfilestream ) 
  GMSSTRM_tbufferedfilestream_DOT_create(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode)
{
  ValueCast(GMSSTRM_tbufferedfilestream,
    GMSSTRM_tbufferedfilestream_DOT_createwithpath(self,filename,mode,_P3str1("\000")));
  return self;
}  /* create */

Destructor(GMSSTRM_tbufferedfilestream ) 
  GMSSTRM_tbufferedfilestream_DOT_destroy(
  GMSSTRM_tbufferedfilestream self)
{
  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten > 0) 
    GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self);
  _P3freemem2(self->GMSSTRM_tbufferedfilestream_DOT_bufptr,self->
    GMSSTRM_tbufferedfilestream_DOT_bufsize);
  _P3freemem2(self->GMSSTRM_tbufferedfilestream_DOT_cbufptr,self->
    GMSSTRM_tbufferedfilestream_DOT_cbufsize);
  GMSSTRM_txfilestream_DOT_destroy(ValueCast(GMSSTRM_txfilestream,self));
  return self;
}  /* destroy */

Function(SYSTEM_boolean ) GMSSTRM_tbufferedfilestream_DOT_fillbuffer(
  GMSSTRM_tbufferedfilestream self)
{
  SYSTEM_boolean result;
  CLIBTYPES_clib_ulong xlen;
  SYSTEM_word wlen;
  SYSTEM_word rlen;

  if (!self->GMSSTRM_tbufferedfilestream_DOT_fcompress) { 
    self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 
      GMSSTRM_txfilestream_DOT_read(ValueCast(GMSSTRM_txfilestream,
      self),*self->GMSSTRM_tbufferedfilestream_DOT_bufptr,self->
      GMSSTRM_tbufferedfilestream_DOT_bufsize);
  } else 
    if (!self->GMSSTRM_tbufferedfilestream_DOT_fcancompress) {
      self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
      self->GMSSTRM_txfilestream_DOT_flastioresult =  -100044;
    } else {
      rlen = GMSSTRM_txfilestream_DOT_read(ValueCast(
        GMSSTRM_txfilestream,self),&self->
        GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader,sizeof(
        GMSSTRM_tcompressheader));
      if (rlen < sizeof(GMSSTRM_tcompressheader)) { 
        self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
      } else {
        wlen = ValueCast(SYSTEM_int32,self->
          GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxb1 << 8) + 
          self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxb2;
        if (self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.
          cxtyp == 0) { 
          self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 
            GMSSTRM_txfilestream_DOT_read(ValueCast(
            GMSSTRM_txfilestream,self),*self->
            GMSSTRM_tbufferedfilestream_DOT_bufptr,wlen);
        } else {
          GMSSTRM_txfilestream_DOT_read(ValueCast(GMSSTRM_txfilestream,
            self),&self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->
            cxdata,wlen);
          xlen = self->GMSSTRM_tbufferedfilestream_DOT_bufsize;
          XCOMPRESS_uncompress(self->
            GMSSTRM_tbufferedfilestream_DOT_bufptr,&xlen,ValueCast(
            SYSTEM_pointer,&self->
            GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxdata),wlen);
          self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = xlen;
        } 
      } 
    } 
  self->GMSSTRM_tbufferedfilestream_DOT_nrread = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_nrwritten = 0;
  result = self->GMSSTRM_tbufferedfilestream_DOT_nrloaded > 0;
  return result;
}  /* fillbuffer */

Function(SYSTEM_boolean ) GMSSTRM_tbufferedfilestream_DOT_flushbuffer(
  GMSSTRM_tbufferedfilestream self)
{
  SYSTEM_boolean result;
  SYSTEM_longword actwritten;
  CLIBTYPES_clib_ulong len;

  result = SYSTEM_true;
  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten == 0) 
    return result;
  if (!(self->GMSSTRM_tbufferedfilestream_DOT_fcompress && self->
    GMSSTRM_tbufferedfilestream_DOT_fcancompress)) {
    actwritten = GMSSTRM_txfilestream_DOT_write(ValueCast(
      GMSSTRM_txfilestream,self),*self->
      GMSSTRM_tbufferedfilestream_DOT_bufptr,self->
      GMSSTRM_tbufferedfilestream_DOT_nrwritten);
    result = self->GMSSTRM_tbufferedfilestream_DOT_nrwritten == 
      actwritten;
  } else {
    len = self->GMSSTRM_tbufferedfilestream_DOT_cbufsize - sizeof(
      GMSSTRM_tcompressheader);
    XCOMPRESS_compress(ValueCast(SYSTEM_pointer,&self->
      GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxdata),&len,self->
      GMSSTRM_tbufferedfilestream_DOT_bufptr,self->
      GMSSTRM_tbufferedfilestream_DOT_nrwritten);
    if (len < ValueCast(SYSTEM_int64,self->
      GMSSTRM_tbufferedfilestream_DOT_nrwritten)) {
      self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxtyp = 1;
      self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxb1 = ValueCast(
        SYSTEM_uint64,len) >> 8;
      self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxb2 = 
        len & 255;
      len = len + sizeof(GMSSTRM_tcompressheader);
      actwritten = GMSSTRM_txfilestream_DOT_write(ValueCast(
        GMSSTRM_txfilestream,self),&self->
        GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxtyp,len);
      result = len == ValueCast(SYSTEM_int64,actwritten);
    } else {
      self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxtyp = 0;
      self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxb1 = 
        self->GMSSTRM_tbufferedfilestream_DOT_nrwritten >> 8;
      self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.cxb2 = 
        self->GMSSTRM_tbufferedfilestream_DOT_nrwritten & 255;
      GMSSTRM_txfilestream_DOT_write(ValueCast(GMSSTRM_txfilestream,
        self),&self->GMSSTRM_tbufferedfilestream_DOT_cbufptr->cxheader.
        cxtyp,sizeof(GMSSTRM_tcompressheader));
      actwritten = GMSSTRM_txfilestream_DOT_write(ValueCast(
        GMSSTRM_txfilestream,self),*self->
        GMSSTRM_tbufferedfilestream_DOT_bufptr,self->
        GMSSTRM_tbufferedfilestream_DOT_nrwritten);
      result = self->GMSSTRM_tbufferedfilestream_DOT_nrwritten == 
        actwritten;
    } 
  } 
  self->GMSSTRM_tbufferedfilestream_DOT_nrwritten = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_nrread = 0;
  return result;
}  /* flushbuffer */

Function(SYSTEM_int64 ) GMSSTRM_tbufferedfilestream_DOT_getposition(
  GMSSTRM_tbufferedfilestream self)
{
  SYSTEM_int64 result;

  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten == 0) { 
    result = self->GMSSTRM_txfilestream_DOT_physposition - self->
      GMSSTRM_tbufferedfilestream_DOT_nrloaded + self->
      GMSSTRM_tbufferedfilestream_DOT_nrread;
  } else {
    if (self->GMSSTRM_tbufferedfilestream_DOT_fcompress) 
      GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self);
    result = self->GMSSTRM_txfilestream_DOT_physposition + self->
      GMSSTRM_tbufferedfilestream_DOT_nrwritten;
  } 
  return result;
}  /* getposition */

Procedure GMSSTRM_tbufferedfilestream_DOT_setposition(
  GMSSTRM_tbufferedfilestream self,
  SYSTEM_int64 p)
{
  SYSTEM_int64 startofbuf;

  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten > 0) {
    if (p == self->GMSSTRM_txfilestream_DOT_physposition + self->
      GMSSTRM_tbufferedfilestream_DOT_nrwritten && !self->
      GMSSTRM_tbufferedfilestream_DOT_fcompress) 
      return;
    GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self);
  } 
  if (self->GMSSTRM_tbufferedfilestream_DOT_nrloaded > 0 && !self->
    GMSSTRM_tbufferedfilestream_DOT_fcompress) {
    startofbuf = self->GMSSTRM_txfilestream_DOT_physposition - self->
      GMSSTRM_tbufferedfilestream_DOT_nrloaded;
    if (p >= startofbuf && p < self->
      GMSSTRM_txfilestream_DOT_physposition) {
      self->GMSSTRM_tbufferedfilestream_DOT_nrread = p - startofbuf;
      return;
    } 
  } 
  GMSSTRM_txfilestream_DOT_setposition(ValueCast(GMSSTRM_txfilestream,
    self),p);
  self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
  self->GMSSTRM_tbufferedfilestream_DOT_nrread = 0;
}  /* setposition */

Function(SYSTEM_int64 ) GMSSTRM_tbufferedfilestream_DOT_getsize(
  GMSSTRM_tbufferedfilestream self)
{
  SYSTEM_int64 result;
  SYSTEM_integer siz;

  result = GMSSTRM_txfilestream_DOT_getsize(ValueCast(
    GMSSTRM_txfilestream,self));
  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten > 0) {
    siz = self->GMSSTRM_txfilestream_DOT_physposition + self->
      GMSSTRM_tbufferedfilestream_DOT_nrwritten;
    if (siz > result) 
      result = siz;
  } 
  return result;
}  /* getsize */

Function(SYSTEM_longword ) GMSSTRM_tbufferedfilestream_DOT_read(
  GMSSTRM_tbufferedfilestream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count)
{
  SYSTEM_longword result;
  GMSGEN_pbytedataarray usrptr;
  SYSTEM_integer usrreadcnt;
  SYSTEM_longword nrbytes;

  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten > 0) 
    GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self);
  if (count <= self->GMSSTRM_tbufferedfilestream_DOT_nrloaded - self->
    GMSSTRM_tbufferedfilestream_DOT_nrread) {
    if (count <= 32) { 
      GMSOBJ_cmove(&(*self->GMSSTRM_tbufferedfilestream_DOT_bufptr)[
        self->GMSSTRM_tbufferedfilestream_DOT_nrread],buffer,count);
    } else 
      SYSTEM_move(&(*self->GMSSTRM_tbufferedfilestream_DOT_bufptr)[
        self->GMSSTRM_tbufferedfilestream_DOT_nrread],buffer,count);
    _P3inc1(self->GMSSTRM_tbufferedfilestream_DOT_nrread,count);
    result = count;
  } else {
    result = 0;
    usrptr = ValueCast(GMSGEN_pbytedataarray,buffer);
    usrreadcnt = 0;
    while (count > 0) {
      if (self->GMSSTRM_tbufferedfilestream_DOT_nrread >= self->
        GMSSTRM_tbufferedfilestream_DOT_nrloaded && !
        GMSSTRM_tbufferedfilestream_DOT_fillbuffer(self)) 
        SYSTEM_break(BRK_2);
      nrbytes = self->GMSSTRM_tbufferedfilestream_DOT_nrloaded - self->
        GMSSTRM_tbufferedfilestream_DOT_nrread;
      if (nrbytes > count) 
        nrbytes = count;
      if (nrbytes <= 32) { 
        GMSOBJ_cmove(&(*self->GMSSTRM_tbufferedfilestream_DOT_bufptr)[
          self->GMSSTRM_tbufferedfilestream_DOT_nrread],&(*usrptr)[
          usrreadcnt],nrbytes);
      } else 
        SYSTEM_move(&(*self->GMSSTRM_tbufferedfilestream_DOT_bufptr)[
          self->GMSSTRM_tbufferedfilestream_DOT_nrread],&(*usrptr)[
          usrreadcnt],nrbytes);
      _P3inc1(self->GMSSTRM_tbufferedfilestream_DOT_nrread,nrbytes);
      _P3inc1(usrreadcnt,nrbytes);
      _P3dec1(count,nrbytes);
      _P3inc1(result,nrbytes);
    
CNT_2:;
    }
BRK_2:;
  } 
  return result;
}  /* read */

Function(SYSTEM_longword ) GMSSTRM_tbufferedfilestream_DOT_write(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword count)
{
  SYSTEM_longword result;
  GMSGEN_pbytedataarray usrptr;
  SYSTEM_integer usrwritecnt;
  SYSTEM_longword nrbytes;

  if (self->GMSSTRM_tbufferedfilestream_DOT_nrloaded > 0) {
    GMSSTRM_txfilestream_DOT_setposition(ValueCast(
      GMSSTRM_txfilestream,self),self->
      GMSSTRM_txfilestream_DOT_physposition - self->
      GMSSTRM_tbufferedfilestream_DOT_nrloaded + self->
      GMSSTRM_tbufferedfilestream_DOT_nrread);
    self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
    self->GMSSTRM_tbufferedfilestream_DOT_nrread = 0;
  } 
  if (count <= self->GMSSTRM_tbufferedfilestream_DOT_bufsize - self->
    GMSSTRM_tbufferedfilestream_DOT_nrwritten) {
    if (count <= 32) { 
      GMSOBJ_cmove(buffer,&(*self->
        GMSSTRM_tbufferedfilestream_DOT_bufptr)[self->
        GMSSTRM_tbufferedfilestream_DOT_nrwritten],count);
    } else 
      SYSTEM_move(buffer,&(*self->
        GMSSTRM_tbufferedfilestream_DOT_bufptr)[self->
        GMSSTRM_tbufferedfilestream_DOT_nrwritten],count);
    _P3inc1(self->GMSSTRM_tbufferedfilestream_DOT_nrwritten,count);
    result = count;
  } else {
    usrptr = ValueCast(GMSGEN_pbytedataarray,buffer);
    usrwritecnt = 0;
    result = 0;
    while (count > 0) {
      nrbytes = self->GMSSTRM_tbufferedfilestream_DOT_bufsize - self->
        GMSSTRM_tbufferedfilestream_DOT_nrwritten;
      if (nrbytes > count) 
        nrbytes = count;
      if (nrbytes <= 32) { 
        GMSOBJ_cmove(&(*usrptr)[usrwritecnt],&(*self->
          GMSSTRM_tbufferedfilestream_DOT_bufptr)[self->
          GMSSTRM_tbufferedfilestream_DOT_nrwritten],nrbytes);
      } else 
        SYSTEM_move(&(*usrptr)[usrwritecnt],&(*self->
          GMSSTRM_tbufferedfilestream_DOT_bufptr)[self->
          GMSSTRM_tbufferedfilestream_DOT_nrwritten],nrbytes);
      _P3inc1(self->GMSSTRM_tbufferedfilestream_DOT_nrwritten,nrbytes);
      _P3inc1(usrwritecnt,nrbytes);
      _P3dec1(count,nrbytes);
      _P3inc1(result,nrbytes);
      if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten >= self->
        GMSSTRM_tbufferedfilestream_DOT_bufsize && !
        GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self)) 
        SYSTEM_break(BRK_3);
    
CNT_3:;
    }
BRK_3:;
  } 
  return result;
}  /* write */

Function(SYSTEM_boolean ) GMSSTRM_tbufferedfilestream_DOT_iseof(
  GMSSTRM_tbufferedfilestream self)
{
  SYSTEM_boolean result;

  if (self->GMSSTRM_tbufferedfilestream_DOT_nrread < self->
    GMSSTRM_tbufferedfilestream_DOT_nrloaded) { 
    result = SYSTEM_false;
  } else 
    result = VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self))) >= VirtMethodCall(ValueCast(
      GMSSTRM_txstream,self), GMSSTRM_txstream_DOT_getsize_T, 3, (ValueCast(
      GMSSTRM_txstream,self)));
  return result;
}  /* iseof */
static SYSTEM_word GMSSTRM_pat_word = 4660;
static SYSTEM_integer GMSSTRM_pat_integer = 305419896;
static SYSTEM_double GMSSTRM_pat_double = 3.1415926535897932385;
cnstdef {GMSSTRM_pat_bad_order = 254};
cnstdef {GMSSTRM_pat_bad_size = 255};

Function(SYSTEM_ansichar *) 
  GMSSTRM_tbufferedfilestream_DOT_getloadpath(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSSTRM_tbufferedfilestream self)
{
  if (_P3strcmpE(self->GMSSTRM_tbufferedfilestream_DOT_floadpath,_P3str1("\000"))) { 
    _P3strclr(result);
  } else 
    SYSUTILS_P3_includetrailingpathdelimiter(result,_len_ret,self->
      GMSSTRM_tbufferedfilestream_DOT_floadpath);
  return result;
}  /* getloadpath */

Procedure GMSSTRM_tbufferedfilestream_DOT_setloadpath(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_ansichar *s)
{
  SYSUTILS_P3_excludetrailingpathdelimiter(self->
    GMSSTRM_tbufferedfilestream_DOT_floadpath,255,s);
}  /* setloadpath */

Constructor(GMSSTRM_tmibufferedstream ) 
  GMSSTRM_tmibufferedstream_DOT_createwithpath(
  GMSSTRM_tmibufferedstream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode,
  const SYSTEM_ansichar *loadpath)
{
  SYSTEM_byte b;
  SYSTEM_word w;
  SYSTEM_integer i;
  SYSTEM_double d;
  GMSSTRM_tdoublevar x;

  ValueCast(GMSSTRM_tmibufferedstream,
    GMSSTRM_tbufferedfilestream_DOT_createwithpath(ValueCast(
    GMSSTRM_tbufferedfilestream,self),filename,mode,loadpath));
  if (self->GMSSTRM_txfilestream_DOT_flastioresult != 0) 
    return self;
  if (mode != GMSSTRM_fmcreate) { 
    GMSSTRM_tmibufferedstream_DOT_determinebyteorder(self);
  } else {
    b = sizeof(SYSTEM_uint16);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      b,sizeof(SYSTEM_uint8)));
    w = GMSSTRM_pat_word;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      w,sizeof(SYSTEM_uint16)));
    b = sizeof(SYSTEM_longint);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      b,sizeof(SYSTEM_uint8)));
    i = GMSSTRM_pat_integer;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      i,sizeof(SYSTEM_int32)));
    b = sizeof(SYSTEM_double);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      b,sizeof(SYSTEM_uint8)));
    d = GMSSTRM_pat_double;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      d,sizeof(SYSTEM_double)));
  } 
  x._u._c1.v = 1.0;
  self->GMSSTRM_tmibufferedstream_DOT_normalorder = x._u._c2.va[0] == 0;
  return self;
}  /* createwithpath */

Constructor(GMSSTRM_tmibufferedstream ) 
  GMSSTRM_tmibufferedstream_DOT_create(
  GMSSTRM_tmibufferedstream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode)
{
  ValueCast(GMSSTRM_tmibufferedstream,
    GMSSTRM_tmibufferedstream_DOT_createwithpath(self,filename,mode,_P3str1("\000")));
  return self;
}  /* create */

Function(SYSTEM_word ) GMSSTRM_tmibufferedstream_DOT_readword(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_word result;
  SYSTEM_word w;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(ValueCast(GMSSTRM_txstream,self),
      GMSSTRM_rw_word);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_word == 0) { 
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      result,sizeof(SYSTEM_uint16)));
  } else {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      w,sizeof(SYSTEM_uint16)));
    GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
      SYSTEM_pointer,&w),ValueCast(SYSTEM_pointer,&result),sizeof(
      SYSTEM_uint16));
  } 
  return result;
}  /* readword */

Function(SYSTEM_integer ) GMSSTRM_tmibufferedstream_DOT_readinteger(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(ValueCast(GMSSTRM_txstream,self),
      GMSSTRM_rw_integer);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_integer == 0) { 
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      result,sizeof(SYSTEM_int32)));
  } else {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      n,sizeof(SYSTEM_int32)));
    GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
      SYSTEM_pointer,&n),ValueCast(SYSTEM_pointer,&result),sizeof(
      SYSTEM_int32));
  } 
  return result;
}  /* readinteger */

Function(SYSTEM_double ) GMSSTRM_tmibufferedstream_DOT_readdouble(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_double result;
  SYSTEM_double f;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(ValueCast(GMSSTRM_txstream,self),
      GMSSTRM_rw_double);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_double == 0) { 
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      result,sizeof(SYSTEM_double)));
  } else {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      f,sizeof(SYSTEM_double)));
    GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
      SYSTEM_pointer,&f),ValueCast(SYSTEM_pointer,&result),sizeof(
      SYSTEM_double));
  } 
  return result;
}  /* readdouble */

Function(SYSTEM_int64 ) GMSSTRM_tmibufferedstream_DOT_readint64(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_int64 result;
  SYSTEM_int64 n;

  if (GMSSTRM_paranoid) 
    GMSSTRM_txstream_DOT_parcheck(ValueCast(GMSSTRM_txstream,self),
      GMSSTRM_rw_int64);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_integer == 0) { 
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      result,sizeof(SYSTEM_int64)));
  } else {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      n,sizeof(SYSTEM_int64)));
    GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
      SYSTEM_pointer,&n),ValueCast(SYSTEM_pointer,&result),sizeof(
      SYSTEM_int64));
  } 
  return result;
}  /* readint64 */

Procedure GMSSTRM_tmibufferedstream_DOT_determinebyteorder(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_word w, w2;
  SYSTEM_integer n, n2;
  SYSTEM_double f, f2;

  VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
    GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
    self->GMSSTRM_tmibufferedstream_DOT_size_word,sizeof(SYSTEM_uint8)));
  if (self->GMSSTRM_tmibufferedstream_DOT_size_word != sizeof(
    SYSTEM_uint16)) {
    self->GMSSTRM_tmibufferedstream_DOT_order_word = 
      GMSSTRM_pat_bad_size;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self),VirtMethodCall(ValueCast(GMSSTRM_txstream,
      self), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self))) + self->
      GMSSTRM_tmibufferedstream_DOT_size_word));
  } else {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      w,sizeof(SYSTEM_uint16)));
    self->GMSSTRM_tmibufferedstream_DOT_order_word = 0;
    if (w != GMSSTRM_pat_word) {
      self->GMSSTRM_tmibufferedstream_DOT_order_word = 1;
      GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
        SYSTEM_pointer,&w),ValueCast(SYSTEM_pointer,&w2),sizeof(
        SYSTEM_uint16));
      if (w2 != GMSSTRM_pat_word) 
        self->GMSSTRM_tmibufferedstream_DOT_order_word = 
          GMSSTRM_pat_bad_order;
    } 
  } 
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
    GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
    self->GMSSTRM_tmibufferedstream_DOT_size_integer,sizeof(
    SYSTEM_uint8)));
  if (self->GMSSTRM_tmibufferedstream_DOT_size_integer != sizeof(
    SYSTEM_longint)) {
    self->GMSSTRM_tmibufferedstream_DOT_order_integer = 
      GMSSTRM_pat_bad_size;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self),VirtMethodCall(ValueCast(GMSSTRM_txstream,
      self), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self))) + self->
      GMSSTRM_tmibufferedstream_DOT_size_integer));
  } else {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      n,sizeof(SYSTEM_int32)));
    self->GMSSTRM_tmibufferedstream_DOT_order_integer = 0;
    if (n != GMSSTRM_pat_integer) {
      self->GMSSTRM_tmibufferedstream_DOT_order_integer = 1;
      GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
        SYSTEM_pointer,&n),ValueCast(SYSTEM_pointer,&n2),sizeof(
        SYSTEM_int32));
      if (n2 != GMSSTRM_pat_integer) 
        self->GMSSTRM_tmibufferedstream_DOT_order_integer = 
          GMSSTRM_pat_bad_order;
    } 
  } 
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
    GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
    self->GMSSTRM_tmibufferedstream_DOT_size_double,sizeof(
    SYSTEM_uint8)));
  if (self->GMSSTRM_tmibufferedstream_DOT_size_double != sizeof(
    SYSTEM_double)) {
    self->GMSSTRM_tmibufferedstream_DOT_order_double = 
      GMSSTRM_pat_bad_size;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self),VirtMethodCall(ValueCast(GMSSTRM_txstream,
      self), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self))) + self->
      GMSSTRM_tmibufferedstream_DOT_size_double));
  } else {
    self->GMSSTRM_tmibufferedstream_DOT_order_double = 0;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      f,sizeof(SYSTEM_double)));
    if (f != GMSSTRM_pat_double) {
      self->GMSSTRM_tmibufferedstream_DOT_order_double = 1;
      GMSSTRM_tmibufferedstream_DOT_reversebytes(self,ValueCast(
        SYSTEM_pointer,&f),ValueCast(SYSTEM_pointer,&f2),sizeof(
        SYSTEM_double));
      if (f2 != GMSSTRM_pat_double) 
        self->GMSSTRM_tmibufferedstream_DOT_order_double = 
          GMSSTRM_pat_bad_order;
    } 
  } 
}  /* determinebyteorder */

Function(SYSTEM_integer ) GMSSTRM_tmibufferedstream_DOT_goodbyteorder(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_integer result;

  result = 0;
  if (self->GMSSTRM_tmibufferedstream_DOT_order_word == 
    GMSSTRM_pat_bad_size) 
    _P3inc1(result,1);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_word == 
    GMSSTRM_pat_bad_order) 
    _P3inc1(result,2);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_integer == 
    GMSSTRM_pat_bad_size) 
    _P3inc1(result,4);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_integer == 
    GMSSTRM_pat_bad_order) 
    _P3inc1(result,8);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_double == 
    GMSSTRM_pat_bad_size) 
    _P3inc1(result,16);
  if (self->GMSSTRM_tmibufferedstream_DOT_order_double == 
    GMSSTRM_pat_bad_order) 
    _P3inc1(result,32);
  return result;
}  /* goodbyteorder */

Procedure GMSSTRM_tmibufferedstream_DOT_reversebytes(
  GMSSTRM_tmibufferedstream self,
  SYSTEM_pointer psrc,
  SYSTEM_pointer pdest,
  SYSTEM_integer sz)
{
  SYSTEM_integer k;

  _P3inc1(PointerCast(SYSTEM_P3_pbyte,&pdest),sz - 1);
  { register SYSTEM_int32 _stop = sz;
    if ((k = 1) <=  _stop) do {
      *ValueCast(SYSTEM_P3_pbyte,pdest) = *ValueCast(SYSTEM_P3_pbyte,
        psrc);
      _P3inc0(PointerCast(SYSTEM_P3_pbyte,&psrc));
      _P3dec0(PointerCast(SYSTEM_P3_pbyte,&pdest));
    
    } while (k++ !=  _stop);

  }
}  /* reversebytes */

Function(SYSTEM_boolean ) GMSSTRM_tmibufferedstream_DOT_wordsneedflip(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_boolean result;

  result = self->GMSSTRM_tmibufferedstream_DOT_order_word != 0;
  return result;
}  /* wordsneedflip */

Function(SYSTEM_boolean ) GMSSTRM_tmibufferedstream_DOT_intsneedflip(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_boolean result;

  result = self->GMSSTRM_tmibufferedstream_DOT_order_integer != 0;
  return result;
}  /* intsneedflip */

Procedure GMSSTRM_tbufferedfilestream_DOT_setcompression(
  GMSSTRM_tbufferedfilestream self,
  SYSTEM_boolean v)
{
  if ((self->GMSSTRM_tbufferedfilestream_DOT_fcompress || v) && self->
    GMSSTRM_tbufferedfilestream_DOT_nrwritten > 0) 
    GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self);
  if (self->GMSSTRM_tbufferedfilestream_DOT_fcompress != v) {
    self->GMSSTRM_tbufferedfilestream_DOT_nrloaded = 0;
    self->GMSSTRM_tbufferedfilestream_DOT_nrread = 0;
  } 
  self->GMSSTRM_tbufferedfilestream_DOT_fcompress = v;
}  /* setcompression */

Function(SYSTEM_double ) GMSSTRM_tmibufferedstream_DOT_readgmsdouble(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_double result;
  SYSTEM_integer i;
  GMSSTRM_tdoublevar z;
  SYSTEM_byte b;
  SYSTEM_integer c;

  b = GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,self));
  if ((ValueCast(SYSTEM_int32,b) & 128) == 0) { 
    switch (b) {
      case 1: 
        result = GMSSPECS_valund;
        break;
      case 2: 
        result = GMSSPECS_valna;
        break;
      case 3: 
        result = GMSSPECS_valpin;
        break;
      case 4: 
        result = GMSSPECS_valmin;
        break;
      case 5: 
        result = GMSSPECS_valeps;
        break;
      case 6: 
        result = GMSSPECS_valacr * 
          GMSSTRM_tmibufferedstream_DOT_readgmsinteger(self);
        break;
      case 7: 
        result = 0.0;
        break;
      case 8: 
        result = 1.0;
        break;
      case 9: 
        result = -1.0;
        break;
      default:
        result = 0.0;
    }
  } else {
    c = ValueCast(SYSTEM_int32,b) & 127;
    if (self->GMSSTRM_tmibufferedstream_DOT_normalorder) { 
      for (i = 1;i <= (SYSTEM_int32)8;++i) {
        if (c == 0) { 
          z._u._c2.va[i - 1] = GMSSTRM_txstream_DOT_readbyte(ValueCast(
            GMSSTRM_txstream,self));
        } else {
          z._u._c2.va[i - 1] = 0;
          _P3dec0(c);
        } 
      }
    } else 
      for (i = 8;i >= (SYSTEM_int32)1;--i) {
        if (c == 0) { 
          z._u._c2.va[i - 1] = GMSSTRM_txstream_DOT_readbyte(ValueCast(
            GMSSTRM_txstream,self));
        } else {
          z._u._c2.va[i - 1] = 0;
          _P3dec0(c);
        } 
      }
    result = z._u._c1.v;
  } 
  return result;
}  /* readgmsdouble */

Function(SYSTEM_integer ) GMSSTRM_tmibufferedstream_DOT_readgmsinteger(
  GMSSTRM_tmibufferedstream self)
{
  SYSTEM_integer result;
  typedef SYSTEM_uint8 _sub_1READGMSINTEGER;
  typedef SYSTEM_byte _arr_0READGMSINTEGER[5];
  _arr_0READGMSINTEGER w;
  SYSTEM_byte b;
  typedef SYSTEM_uint8 _sub_2READGMSINTEGER;
  _sub_2READGMSINTEGER c;
  SYSTEM_boolean neg;

  VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
    GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
    b,1));
  w[0] = ValueCast(SYSTEM_int32,b) & 15;
  neg = b >= 128;
  c = ValueCast(SYSTEM_int32,b >> 4) & 7;
  if (c > 0) 
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,self),&
      w[1],c));
  result = 0;
  while (c >= 1) {
    result = ValueCast(SYSTEM_uint32,result) << 8 | w[c];
    c = ValueCast(SYSTEM_int32,c) - 1;
  
}
  result = ValueCast(SYSTEM_uint32,result) << 4 | w[0];
  if (neg) 
    result = -result;
  return result;
}  /* readgmsinteger */

Procedure GMSSTRM_tmibufferedstream_DOT_writegmsdouble(
  GMSSTRM_tmibufferedstream self,
  SYSTEM_double d)
{
  SYSTEM_integer i;
  GMSSTRM_tdoublevar z;
  GMSSPECS_tgmsvalue gv;
  SYSTEM_byte b;
  SYSTEM_integer c;

  gv = GMSSPECS_mapval(d);
  b = SYSTEM_ord(gv);
  if (gv == GMSSPECS_xvreal) 
    if (d == 0.0) { 
      b = 7;
    } else 
      if (d == 1.0) { 
        b = 8;
      } else 
        if (d == -1.0) 
          b = 9;
  if (b != 0) {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
      GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
      b,1));
    if (gv == GMSSPECS_xvacr) 
      GMSSTRM_tmibufferedstream_DOT_writegmsinteger(self,SYSTEM_round(
        d /  GMSSPECS_valacr));
  } else {
    z._u._c1.v = d;
    if (self->GMSSTRM_tmibufferedstream_DOT_normalorder) {
      c = 0;
      for (i = 1;i <= (SYSTEM_int32)8;++i) {
        if (z._u._c2.va[i - 1] == 0) { 
          _P3inc0(c);
        } else 
          SYSTEM_break(BRK_4);
      CNT_4:;
}
BRK_4:;
      b = 128 | c;
      VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
        GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,
        self),&b,1));
      VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
        GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,
        self),&z._u._c2.va[c + 1 - 1],8 - c));
    } else {
      c = 0;
      for (i = 8;i >= (SYSTEM_int32)1;--i) {
        if (z._u._c2.va[i - 1] == 0) { 
          _P3inc0(c);
        } else 
          SYSTEM_break(BRK_5);
      CNT_5:;
}
BRK_5:;
      b = 128 | c;
      VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
        GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,
        self),&b,1));
      for (i = 8 - c;i >= (SYSTEM_int32)1;--i) {
        VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
          GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,
          self),&z._u._c2.va[i - 1],1));
      }
    } 
  } 
}  /* writegmsdouble */

Procedure GMSSTRM_tmibufferedstream_DOT_writegmsinteger(
  GMSSTRM_tmibufferedstream self,
  SYSTEM_integer n)
{
  typedef SYSTEM_uint8 _sub_1WRITEGMSINTEGER;
  typedef SYSTEM_byte _arr_0WRITEGMSINTEGER[5];
  _arr_0WRITEGMSINTEGER w;
  SYSTEM_integer c;
  SYSTEM_byte b;

  if (n >= 0) { 
    b = 0;
  } else {
    b = 128;
    n = -n;
  } 
  b = b | n & 15;
  n = ValueCast(SYSTEM_uint32,n) >> 4;
  c = 0;
  while (n != 0) {
    c = c + 1;
    w[c] = n & 255;
    n = ValueCast(SYSTEM_uint32,n) >> 8;
  
}
  w[0] = b | ValueCast(SYSTEM_uint32,c) << 4;
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self), 
    GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(GMSSTRM_txstream,self),&
    w[0],c + 1));
}  /* writegmsinteger */

Constructor(GMSSTRM_tbinarytextfileio ) 
  GMSSTRM_tbinarytextfileio_DOT_openforread(
  GMSSTRM_tbinarytextfileio self,
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *password,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg)
{
  SYSTEM_byte b1, b2;
  SYSTEM_ansichar ch;
  SYSTEM_boolean hascomp, haspswd;
  SYSTEM_shortstring src, targ;

  ValueCast(GMSSTRM_tbinarytextfileio,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSSTRM_tbinarytextfileio_DOT_fs = ValueCast(
    GMSSTRM_tbufferedfilestream,GMSSTRM_tbufferedfilestream_DOT_create(ValueCast(
    GMSSTRM_tbufferedfilestream,_P3alloc_object(&
    GMSSTRM_tbufferedfilestream_CD)),fn,GMSSTRM_fmopenread));
  *errnr = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
    GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
  if (*errnr != 0) {
    SYSUTILS_P3_syserrormessage(errmsg,255,*errnr);
    *errnr = GMSSTRM_strmerrorioresult;
    goto _Lerrorexit_67;
  } 
  b1 = GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs));
  b2 = GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs));
  if (b1 == 31 && b2 == 139) {
    self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature = 
      GMSSTRM_fsign_gzip;
    SYSUTILS_P3_freeandnil(&self->GMSSTRM_tbinarytextfileio_DOT_fs);
    self->GMSSTRM_tbinarytextfileio_DOT_gzfs = ValueCast(
      GMSSTRM_tgzipinputstream,GMSSTRM_tgzipinputstream_DOT_create(ValueCast(
      GMSSTRM_tgzipinputstream,_P3alloc_object(&
      GMSSTRM_tgzipinputstream_CD)),fn,errmsg));
    if (_P3strcmpN(errmsg,_P3str1("\000"))) 
      *errnr = 1;
    return self;
  } 
  _P3setlength(src,b2,255);
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs), GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(
    GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs),&src[1],
    b2));
  if (!(b1 == GMSSTRM_signature_header && _P3strcmpE(src,
    GMSSTRM_signature_gams))) {
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs), 
      GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs),0));
    self->GMSSTRM_tbinarytextfileio_DOT_frewindpoint = 0;
    self->GMSSTRM_tbinarytextfileio_DOT_fmajorversionread = 0;
    self->GMSSTRM_tbinarytextfileio_DOT_fminorversionread = 0;
    self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature = 
      GMSSTRM_fsign_text;
    _P3strclr(errmsg);
    return self;
  } 
  *errnr = GMSSTRM_strmerrorgamsheader;
  _P3strcpy(errmsg,255,_P3str1("\025GAMS header not found"));
  self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature = ValueCast(
    GMSSTRM_tfilesignature,GMSSTRM_txstream_DOT_readbyte(ValueCast(
    GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs)) - 65);
  {
    SYSTEM_shortstring _t2;

    GMSSTRM_txstream_DOT_readstring(_t2,255,ValueCast(
      GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
  }
  self->GMSSTRM_tbinarytextfileio_DOT_fmajorversionread = 
    GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs));
  self->GMSSTRM_tbinarytextfileio_DOT_fminorversionread = 
    GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs));
  ch = ValueCast(SYSTEM_ansichar,GMSSTRM_txstream_DOT_readbyte(ValueCast(
    GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs)));
  if (ch == _P3char('P')) { 
    haspswd = SYSTEM_true;
  } else 
    if (ch == _P3char('p')) { 
      haspswd = SYSTEM_false;
    } else 
      goto _Lerrorexit_67;
  ch = ValueCast(SYSTEM_ansichar,GMSSTRM_txstream_DOT_readbyte(ValueCast(
    GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs)));
  if (ch == _P3char('C')) { 
    hascomp = SYSTEM_true;
  } else 
    if (ch == _P3char('c')) { 
      hascomp = SYSTEM_false;
    } else 
      goto _Lerrorexit_67;
  if (haspswd && _P3strcmpE(password,_P3str1("\000"))) {
    *errnr = GMSSTRM_strmerrornopassword;
    _P3strcpy(errmsg,255,_P3str1("\026A Password is required"));
    goto _Lerrorexit_67;
  } 
  *errnr = GMSSTRM_strmerrorintegrity;
  _P3strcpy(errmsg,255,_P3str1("\026Integrity check failed"));
  if (haspswd) {
    GMSSTRM_txfilestream_DOT_setpassword(ValueCast(
      GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),
      password);
    GMSSTRM_txstream_DOT_readstring(src,255,ValueCast(
      GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
    GMSSTRM_txfilestream_DOT_applypassword(ValueCast(
      GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),ValueCast(
      GMSGEN_pbytedataarray,&src[1]),ValueCast(
      GMSGEN_pbytedataarray,&targ[1]),SYSTEM_length(src),
      GMSSTRM_verify_offset);
    _P3setlength(targ,SYSTEM_length(src),255);
    {
      SYSTEM_shortstring _t1;

      if (_P3strcmpN(targ,GMSSTRM_txfilestream_DOT_randstring(_t1,255,ValueCast(
        GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),
        SYSTEM_length(src)))) 
        goto _Lerrorexit_67;
    }
  } 
  self->GMSSTRM_tbinarytextfileio_DOT_frewindpoint = VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs), 
    GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(GMSSTRM_txstream,
    self->GMSSTRM_tbinarytextfileio_DOT_fs)));
  GMSSTRM_tbufferedfilestream_DOT_setcompression(self->
    GMSSTRM_tbinarytextfileio_DOT_fs,SYSTEM_true);
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs), 
    GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(GMSSTRM_txstream,
    self->GMSSTRM_tbinarytextfileio_DOT_fs),self->
    GMSSTRM_tbinarytextfileio_DOT_frewindpoint));
  if (!hascomp) 
    GMSSTRM_tbufferedfilestream_DOT_setcompression(self->
      GMSSTRM_tbinarytextfileio_DOT_fs,SYSTEM_false);
  {
    SYSTEM_shortstring _t1;

    if (_P3strcmpN(GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs)),
      GMSSTRM_signature_gams)) 
      goto _Lerrorexit_67;
  }
  *errnr = GMSSTRM_strmerrornoerror;
  _P3strclr(errmsg);
  return self;
  _Lerrorexit_67:;
  SYSUTILS_P3_freeandnil(&self->GMSSTRM_tbinarytextfileio_DOT_fs);
  return self;
}  /* openforread */

Constructor(GMSSTRM_tbinarytextfileio ) 
  GMSSTRM_tbinarytextfileio_DOT_openforwrite(
  GMSSTRM_tbinarytextfileio self,
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *producer,
  const SYSTEM_ansichar *password,
  GMSSTRM_tfilesignature signature,
  SYSTEM_boolean comp,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg)
{
  SYSTEM_shortstring src, targ;

  ValueCast(GMSSTRM_tbinarytextfileio,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature = signature;
  self->GMSSTRM_tbinarytextfileio_DOT_frw = GMSSTRM_fm_write;
  self->GMSSTRM_tbinarytextfileio_DOT_fs = ValueCast(
    GMSSTRM_tbufferedfilestream,GMSSTRM_tbufferedfilestream_DOT_create(ValueCast(
    GMSSTRM_tbufferedfilestream,_P3alloc_object(&
    GMSSTRM_tbufferedfilestream_CD)),fn,GMSSTRM_fmcreate));
  if (signature != GMSSTRM_fsign_text || _P3strcmpN(password,_P3str1("\000")) || 
    comp) {
    GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),GMSSTRM_signature_header);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),GMSSTRM_signature_gams);
    GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),SYSTEM_ord(signature) + 65);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),producer);
    GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),1);
    GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),1);
    if (_P3strcmpE(password,_P3str1("\000"))) { 
      GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
        GMSSTRM_tbinarytextfileio_DOT_fs),112);
    } else 
      GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
        GMSSTRM_tbinarytextfileio_DOT_fs),80);
    if (comp) { 
      GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
        GMSSTRM_tbinarytextfileio_DOT_fs),67);
    } else 
      GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
        GMSSTRM_tbinarytextfileio_DOT_fs),99);
    if (_P3strcmpN(password,_P3str1("\000"))) {
      GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self->
        GMSSTRM_tbinarytextfileio_DOT_fs);
      GMSSTRM_txfilestream_DOT_setpassword(ValueCast(
        GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),
        password);
      GMSSTRM_txfilestream_DOT_randstring(src,255,ValueCast(
        GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),
        SYSTEM_length(password));
      GMSSTRM_txfilestream_DOT_applypassword(ValueCast(
        GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),ValueCast(
        GMSGEN_pbytedataarray,&src[1]),ValueCast(
        GMSGEN_pbytedataarray,&targ[1]),SYSTEM_length(src),
        GMSSTRM_verify_offset);
      GMSSTRM_txfilestream_DOT_setpassword(ValueCast(
        GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),_P3str1("\000"));
      _P3setlength(targ,SYSTEM_length(src),255);
      GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
        GMSSTRM_tbinarytextfileio_DOT_fs),targ);
    } 
    if (comp) { 
      GMSSTRM_tbufferedfilestream_DOT_setcompression(self->
        GMSSTRM_tbinarytextfileio_DOT_fs,SYSTEM_true);
    } else 
      GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self->
        GMSSTRM_tbinarytextfileio_DOT_fs);
    GMSSTRM_txfilestream_DOT_setpassword(ValueCast(
      GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs),
      password);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs),GMSSTRM_signature_gams);
  } 
  *errnr = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
    GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
  if (*errnr == 0) {
    *errnr = GMSSTRM_strmerrornoerror;
    _P3strclr(errmsg);
  } else {
    SYSUTILS_P3_syserrormessage(errmsg,255,*errnr);
    SYSUTILS_P3_freeandnil(&self->GMSSTRM_tbinarytextfileio_DOT_fs);
  } 
  return self;
}  /* openforwrite */

Destructor(GMSSTRM_tbinarytextfileio ) 
  GMSSTRM_tbinarytextfileio_DOT_destroy(
  GMSSTRM_tbinarytextfileio self)
{
  if (self->GMSSTRM_tbinarytextfileio_DOT_fs != NULL) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
      GMSSTRM_tbinarytextfileio_DOT_fs));
  if (_P3assigned(self->GMSSTRM_tbinarytextfileio_DOT_gzfs)) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
      GMSSTRM_tbinarytextfileio_DOT_gzfs));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_integer ) GMSSTRM_tbinarytextfileio_DOT_read(
  GMSSTRM_tbinarytextfileio self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer count)
{
  SYSTEM_integer result;

  if (self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature == 
    GMSSTRM_fsign_gzip) { 
    result = GMSSTRM_tgzipinputstream_DOT_read(self->
      GMSSTRM_tbinarytextfileio_DOT_gzfs,buffer,count);
  } else 
    result = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs), GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(
      GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs),buffer,
      count));
  return result;
}  /* read */

Function(SYSTEM_integer ) GMSSTRM_tbinarytextfileio_DOT_write(
  GMSSTRM_tbinarytextfileio self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer count)
{
  SYSTEM_integer result;

  SYSTEM_assert(self->GMSSTRM_tbinarytextfileio_DOT_frw == 
    GMSSTRM_fm_write,_P3str1("\026TBinaryTextFileIO.Read"));
  if (self->GMSSTRM_tbinarytextfileio_DOT_fs == NULL) { 
    result =  -1;
  } else 
    result = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GMSSTRM_tbinarytextfileio_DOT_fs), GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(
      GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs),buffer,
      count));
  return result;
}  /* write */

Function(SYSTEM_boolean ) GMSSTRM_tbinarytextfileio_DOT_usespassword(
  GMSSTRM_tbinarytextfileio self)
{
  SYSTEM_boolean result;

  result = self->GMSSTRM_tbinarytextfileio_DOT_fs != NULL && 
    GMSSTRM_txfilestream_DOT_getusespassword(ValueCast(
    GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
  return result;
}  /* usespassword */

Procedure GMSSTRM_tbinarytextfileio_DOT_rewind(
  GMSSTRM_tbinarytextfileio self)
{
  SYSTEM_assert(self->GMSSTRM_tbinarytextfileio_DOT_frw == 
    GMSSTRM_fm_read,_P3str1("\031TBinaryTextFileIO.ReWind1"));
  SYSTEM_assert(self->GMSSTRM_tbinarytextfileio_DOT_fs != NULL,_P3str1("\031TBinaryTextFileIO.ReWind2"));
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GMSSTRM_tbinarytextfileio_DOT_fs), 
    GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(GMSSTRM_txstream,
    self->GMSSTRM_tbinarytextfileio_DOT_fs),self->
    GMSSTRM_tbinarytextfileio_DOT_frewindpoint));
  if (self->GMSSTRM_tbinarytextfileio_DOT_fs->
    GMSSTRM_tbufferedfilestream_DOT_fcompress) 
    {
      SYSTEM_shortstring _t2;

      GMSSTRM_txstream_DOT_readstring(_t2,255,ValueCast(
        GMSSTRM_txstream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
    }
}  /* rewind */

Function(SYSTEM_integer ) 
  GMSSTRM_tbinarytextfileio_DOT_getlastioresult(
  GMSSTRM_tbinarytextfileio self)
{
  SYSTEM_integer result;

  result = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
    GMSSTRM_txfilestream,self->GMSSTRM_tbinarytextfileio_DOT_fs));
  return result;
}  /* getlastioresult */

Procedure GMSSTRM_compresstextfile(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *fo,
  const SYSTEM_ansichar *password,
  SYSTEM_boolean comp,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg)
{
  cnstdef {bufsize = 4096};
  GMSSTRM_tbinarytextfileio fin;
  GMSSTRM_tbinarytextfileio fout;
  typedef SYSTEM_uint16 _sub_2COMPRESSTEXTFILE;
  typedef SYSTEM_byte _arr_1COMPRESSTEXTFILE[4096];
  _arr_1COMPRESSTEXTFILE buffer;
  SYSTEM_integer nrread;

  fout = NULL;
  fin = ValueCast(GMSSTRM_tbinarytextfileio,
    GMSSTRM_tbinarytextfileio_DOT_openforread(ValueCast(
    GMSSTRM_tbinarytextfileio,_P3alloc_object(&
    GMSSTRM_tbinarytextfileio_CD)),fn,_P3str1("\000"),errnr,errmsg));
  if (_P3strcmpN(errmsg,_P3str1("\000"))) 
    goto _Lalldone_75;
  fout = ValueCast(GMSSTRM_tbinarytextfileio,
    GMSSTRM_tbinarytextfileio_DOT_openforwrite(ValueCast(
    GMSSTRM_tbinarytextfileio,_P3alloc_object(&
    GMSSTRM_tbinarytextfileio_CD)),fo,_P3str1("\020CompressTextFile"),
    password,GMSSTRM_fsign_text,comp,errnr,errmsg));
  if (_P3strcmpN(errmsg,_P3str1("\000"))) 
    goto _Lalldone_75;
  do {
    nrread = GMSSTRM_tbinarytextfileio_DOT_read(fin,buffer,bufsize);
    if (nrread == 0) 
      SYSTEM_break(BRK_6);
    GMSSTRM_tbinarytextfileio_DOT_write(fout,buffer,nrread);
  CNT_6:;
  } while (!(nrread < 4096));
BRK_6:;
  _Lalldone_75:;
  if (fin != NULL) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,fin));
  if (fout != NULL) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,fout));
}  /* compresstextfile */

Procedure GMSSTRM_uncompresstextfile(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *fo,
  const SYSTEM_ansichar *password,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg)
{
  cnstdef {bufsize_2 = 4096};
  GMSSTRM_tbinarytextfileio fin;
  GMSSTRM_tbinarytextfileio fout;
  typedef SYSTEM_uint16 _sub_2UNCOMPRESSTEXTFILE;
  typedef SYSTEM_byte _arr_1UNCOMPRESSTEXTFILE[4096];
  _arr_1UNCOMPRESSTEXTFILE buffer;
  SYSTEM_integer nrread;

  fout = NULL;
  fin = ValueCast(GMSSTRM_tbinarytextfileio,
    GMSSTRM_tbinarytextfileio_DOT_openforread(ValueCast(
    GMSSTRM_tbinarytextfileio,_P3alloc_object(&
    GMSSTRM_tbinarytextfileio_CD)),fn,password,errnr,errmsg));
  if (_P3strcmpN(errmsg,_P3str1("\000"))) 
    goto _Lalldone_76;
  fout = ValueCast(GMSSTRM_tbinarytextfileio,
    GMSSTRM_tbinarytextfileio_DOT_openforwrite(ValueCast(
    GMSSTRM_tbinarytextfileio,_P3alloc_object(&
    GMSSTRM_tbinarytextfileio_CD)),fo,_P3str1("\000"),_P3str1("\000"),
    GMSSTRM_fsign_text,SYSTEM_false,errnr,errmsg));
  if (_P3strcmpN(errmsg,_P3str1("\000"))) 
    goto _Lalldone_76;
  do {
    nrread = GMSSTRM_tbinarytextfileio_DOT_read(fin,buffer,bufsize_2);
    if (nrread == 0) 
      SYSTEM_break(BRK_7);
    GMSSTRM_tbinarytextfileio_DOT_write(fout,buffer,nrread);
  CNT_7:;
  } while (!(nrread < 4096));
BRK_7:;
  _Lalldone_76:;
  if (fin != NULL) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,fin));
  if (fout != NULL) 
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,fout));
}  /* uncompresstextfile */

Procedure GMSSTRM_tbinarytextfileio_DOT_readline(
  GMSSTRM_tbinarytextfileio self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer *len,
  SYSTEM_integer maxinp,
  SYSTEM_ansichar *lastchar)
{
  typedef SYSTEM_ansichar *_ptr_0READLINE;
  _ptr_0READLINE pbuf;

  if (self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature == 
    GMSSTRM_fsign_gzip) { 
    GMSSTRM_tgzipinputstream_DOT_readline(self->
      GMSSTRM_tbinarytextfileio_DOT_gzfs,buffer,len,maxinp,lastchar);
  } else {
    *len = 0;
    pbuf = ValueCast(_ptr_0READLINE,buffer);
    while (!(*lastchar == _P3char('\015') || *lastchar == _P3char('\012') || *
      lastchar == _P3char('\032') || *len == maxinp)) {
      *len = *len + 1;
      *pbuf = *lastchar;
      _P3inc0(pbuf);
      { register GMSSTRM_tbufferedfilestream_OD *_W2=self->
        GMSSTRM_tbinarytextfileio_DOT_fs;
        if (1 <= _W2->GMSSTRM_tbufferedfilestream_DOT_nrloaded - 
          _W2->GMSSTRM_tbufferedfilestream_DOT_nrread) {
          *lastchar = ValueCast(SYSTEM_ansichar,(*_W2->
            GMSSTRM_tbufferedfilestream_DOT_bufptr)[_W2->
            GMSSTRM_tbufferedfilestream_DOT_nrread]);
          _P3inc0(_W2->GMSSTRM_tbufferedfilestream_DOT_nrread);
        } else 
          if (VirtMethodCall(ValueCast(GMSSTRM_txstream,_W2), 
            GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(GMSSTRM_txstream,
            _W2),lastchar,1)) <= 0) 
            *lastchar = _P3char('\032');

      }
    
}
  } 
}  /* readline */

Procedure GMSSTRM_compressfromstdin(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *password,
  SYSTEM_boolean comp,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg)
{
  cnstdef {bufsize_3 = 4096};
  GMSSTRM_tbinarytextfileio fout;
  typedef SYSTEM_uint16 _sub_2COMPRESSFROMSTDIN;
  typedef SYSTEM_byte _arr_1COMPRESSFROMSTDIN[4096];
  _arr_1COMPRESSFROMSTDIN buffer;
  typedef SYSTEM_uint16 _sub_4COMPRESSFROMSTDIN;
  typedef SYSTEM_byte _arr_3COMPRESSFROMSTDIN[4096];
  _arr_3COMPRESSFROMSTDIN inbuffer;
  SYSTEM_integer nrread;
  SYSTEM_ansichar ch;

  fout = ValueCast(GMSSTRM_tbinarytextfileio,
    GMSSTRM_tbinarytextfileio_DOT_openforwrite(ValueCast(
    GMSSTRM_tbinarytextfileio,_P3alloc_object(&
    GMSSTRM_tbinarytextfileio_CD)),fn,_P3str1("\015CompressStdIn"),
    password,GMSSTRM_fsign_text,comp,errnr,errmsg));
  if (*errnr != 0 || _P3strcmpN(errmsg,_P3str1("\000"))) 
    goto _Lalldone_78;
  SYSTEM_P3_settextbuf(&SYSTEM_input,inbuffer);
  do {
    nrread = 0;
    while (!_P3eoff(1,SYSTEM_input)) {
      _Iplus_bgn();
      {
        _P3file_ptr _file_temp = &SYSTEM_input;

        _P3read_fc0(ch);
      }
      _Iplus_end();
      nrread = nrread + 1;
      buffer[nrread - 1] = SYSTEM_ord(ch);
      if (nrread == 4096) 
        SYSTEM_break(BRK_8);
    
CNT_8:;
    }
BRK_8:;
    if (nrread == 0) 
      SYSTEM_break(BRK_9);
    GMSSTRM_tbinarytextfileio_DOT_write(fout,buffer,nrread);
  CNT_9:;
  } while (!(nrread < 4096));
BRK_9:;
  *errnr = GMSSTRM_tbinarytextfileio_DOT_getlastioresult(fout);
  if (*errnr != 0) 
    SYSUTILS_P3_syserrormessage(errmsg,255,*errnr);
  _Lalldone_78:;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,fout));
}  /* compressfromstdin */

Procedure GMSSTRM_uncompresstostdout(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *password,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg)
{
  cnstdef {bufsize_4 = 4096};
  GMSSTRM_tbinarytextfileio fin;
  typedef SYSTEM_uint16 _sub_2UNCOMPRESSTOSTDOUT;
  typedef SYSTEM_byte _arr_1UNCOMPRESSTOSTDOUT[4096];
  _arr_1UNCOMPRESSTOSTDOUT buffer;
  typedef SYSTEM_uint16 _sub_4UNCOMPRESSTOSTDOUT;
  typedef SYSTEM_byte _arr_3UNCOMPRESSTOSTDOUT[4096];
  _arr_3UNCOMPRESSTOSTDOUT outbuffer;
  SYSTEM_integer nrread;
  SYSTEM_integer n;

  fin = ValueCast(GMSSTRM_tbinarytextfileio,
    GMSSTRM_tbinarytextfileio_DOT_openforread(ValueCast(
    GMSSTRM_tbinarytextfileio,_P3alloc_object(&
    GMSSTRM_tbinarytextfileio_CD)),fn,password,errnr,errmsg));
  if (*errnr != 0 || _P3strcmpN(errmsg,_P3str1("\000"))) 
    goto _Lalldone_79;
  SYSTEM_P3_settextbuf(&SYSTEM_output,outbuffer);
  do {
    nrread = GMSSTRM_tbinarytextfileio_DOT_read(fin,buffer,bufsize_4);
    if (nrread == 0) 
      SYSTEM_break(BRK_10);
    { register SYSTEM_int32 _stop = nrread;
      if ((n = 1) <=  _stop) do {
        _Iplus_bgn();
        _P3write_c0(ValueCast(SYSTEM_ansichar,buffer[n - 1]));
        _Iplus_end();
      
      } while (n++ !=  _stop);

    }
  CNT_10:;
  } while (!(nrread < 4096));
BRK_10:;
  *errnr = GMSSTRM_tbinarytextfileio_DOT_getlastioresult(fin);
  if (*errnr != 0) 
    SYSUTILS_P3_syserrormessage(errmsg,255,*errnr);
  _Lalldone_79:;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,fin));
}  /* uncompresstostdout */

Function(SYSTEM_ansichar ) GMSSTRM_tbinarytextfileio_DOT_readcharacter(
  GMSSTRM_tbinarytextfileio self)
{
  SYSTEM_ansichar result;

  if (self->GMSSTRM_tbinarytextfileio_DOT_ffilesignature == 
    GMSSTRM_fsign_gzip) { 
    if (GMSSTRM_tgzipinputstream_DOT_read(self->
      GMSSTRM_tbinarytextfileio_DOT_gzfs,&result,1) == 0) 
      result = _P3char('\032');
  } else 
    result = GMSSTRM_tbufferedfilestream_DOT_readcharacter(self->
      GMSSTRM_tbinarytextfileio_DOT_fs);
  return result;
}  /* readcharacter */

Constructor(GMSSTRM_tgzipinputstream ) 
  GMSSTRM_tgzipinputstream_DOT_create(
  GMSSTRM_tgzipinputstream self,
  const SYSTEM_ansichar *fn,
  SYSTEM_ansichar *errmsg)
{
  ValueCast(GMSSTRM_tgzipinputstream,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSSTRM_tgzipinputstream_DOT_pgz = XCOMPRESS_gzreadopen(fn);
  if (self->GMSSTRM_tgzipinputstream_DOT_pgz == NULL) { 
    _P3strcpy(errmsg,255,_P3str1("\020Cannot open file"));
  } else {
    _P3strclr(errmsg);
    self->GMSSTRM_tgzipinputstream_DOT_bufsize = GMSSTRM_buffersize;
    _P3getmem(self->GMSSTRM_tgzipinputstream_DOT_bufptr,self->
      GMSSTRM_tgzipinputstream_DOT_bufsize);
    self->GMSSTRM_tgzipinputstream_DOT_nrread = 0;
    self->GMSSTRM_tgzipinputstream_DOT_nrloaded = 0;
  } 
  return self;
}  /* create */

Destructor(GMSSTRM_tgzipinputstream ) 
  GMSSTRM_tgzipinputstream_DOT_destroy(
  GMSSTRM_tgzipinputstream self)
{
  XCOMPRESS_gzreadclose(&self->GMSSTRM_tgzipinputstream_DOT_pgz);
  if (self->GMSSTRM_tgzipinputstream_DOT_bufptr != NULL) 
    _P3freemem(self->GMSSTRM_tgzipinputstream_DOT_bufptr);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

static Function(SYSTEM_boolean ) fillbuffer(
  GMSSTRM_tgzipinputstream *_2self)
{
  SYSTEM_boolean result;

  (*_2self)->GMSSTRM_tgzipinputstream_DOT_nrloaded = XCOMPRESS_gzread((*
    _2self)->GMSSTRM_tgzipinputstream_DOT_pgz,*(*_2self)->
    GMSSTRM_tgzipinputstream_DOT_bufptr,(*_2self)->
    GMSSTRM_tgzipinputstream_DOT_bufsize);
  (*_2self)->GMSSTRM_tgzipinputstream_DOT_nrread = 0;
  result = (*_2self)->GMSSTRM_tgzipinputstream_DOT_nrloaded > 0;
  return result;
}  /* fillbuffer */

Function(SYSTEM_longword ) GMSSTRM_tgzipinputstream_DOT_read(
  GMSSTRM_tgzipinputstream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count)
{
  SYSTEM_longword result;
  GMSGEN_pbytedataarray usrptr;
  SYSTEM_integer usrreadcnt;
  SYSTEM_longword nrbytes;

  if (count <= self->GMSSTRM_tgzipinputstream_DOT_nrloaded - self->
    GMSSTRM_tgzipinputstream_DOT_nrread) {
    if (count <= 32) { 
      GMSOBJ_cmove(&(*self->GMSSTRM_tgzipinputstream_DOT_bufptr)[self->
        GMSSTRM_tgzipinputstream_DOT_nrread],buffer,count);
    } else 
      SYSTEM_move(&(*self->GMSSTRM_tgzipinputstream_DOT_bufptr)[self->
        GMSSTRM_tgzipinputstream_DOT_nrread],buffer,count);
    _P3inc1(self->GMSSTRM_tgzipinputstream_DOT_nrread,count);
    result = count;
  } else {
    result = 0;
    usrptr = ValueCast(GMSGEN_pbytedataarray,buffer);
    usrreadcnt = 0;
    while (count > 0) {
      if (self->GMSSTRM_tgzipinputstream_DOT_nrread >= self->
        GMSSTRM_tgzipinputstream_DOT_nrloaded && !fillbuffer(&self)) 
        SYSTEM_break(BRK_11);
      nrbytes = self->GMSSTRM_tgzipinputstream_DOT_nrloaded - self->
        GMSSTRM_tgzipinputstream_DOT_nrread;
      if (nrbytes > count) 
        nrbytes = count;
      if (nrbytes <= 32) { 
        GMSOBJ_cmove(&(*self->GMSSTRM_tgzipinputstream_DOT_bufptr)[
          self->GMSSTRM_tgzipinputstream_DOT_nrread],&(*usrptr)[
          usrreadcnt],nrbytes);
      } else 
        SYSTEM_move(&(*self->GMSSTRM_tgzipinputstream_DOT_bufptr)[self->
          GMSSTRM_tgzipinputstream_DOT_nrread],&(*usrptr)[usrreadcnt],
          nrbytes);
      _P3inc1(self->GMSSTRM_tgzipinputstream_DOT_nrread,nrbytes);
      _P3inc1(usrreadcnt,nrbytes);
      _P3dec1(count,nrbytes);
      _P3inc1(result,nrbytes);
    
CNT_11:;
    }
BRK_11:;
  } 
  return result;
}  /* read */

Procedure GMSSTRM_tgzipinputstream_DOT_readline(
  GMSSTRM_tgzipinputstream self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer *len,
  SYSTEM_integer maxinp,
  SYSTEM_ansichar *lastchar)
{
  typedef SYSTEM_ansichar *_ptr_0READLINE;
  _ptr_0READLINE pbuf;

  *len = 0;
  pbuf = ValueCast(_ptr_0READLINE,buffer);
  while (!(*lastchar == _P3char('\015') || *lastchar == _P3char('\012') || *
    lastchar == _P3char('\032') || *len == maxinp)) {
    *len = *len + 1;
    *pbuf = *lastchar;
    _P3inc0(pbuf);
    if (self->GMSSTRM_tgzipinputstream_DOT_nrloaded - self->
      GMSSTRM_tgzipinputstream_DOT_nrread >= 1) {
      *lastchar = ValueCast(SYSTEM_ansichar,(*self->
        GMSSTRM_tgzipinputstream_DOT_bufptr)[self->
        GMSSTRM_tgzipinputstream_DOT_nrread]);
      _P3inc0(self->GMSSTRM_tgzipinputstream_DOT_nrread);
    } else 
      if (GMSSTRM_tgzipinputstream_DOT_read(self,lastchar,1) <= 0) 
        *lastchar = _P3char('\032');
  
}
}  /* readline */

Function(SYSTEM_ansichar ) 
  GMSSTRM_tbufferedfilestream_DOT_readcharacter(
  GMSSTRM_tbufferedfilestream self)
{
  SYSTEM_ansichar result;

  if (self->GMSSTRM_tbufferedfilestream_DOT_nrwritten > 0) 
    GMSSTRM_tbufferedfilestream_DOT_flushbuffer(self);
  if (self->GMSSTRM_tbufferedfilestream_DOT_nrread >= self->
    GMSSTRM_tbufferedfilestream_DOT_nrloaded && !
    GMSSTRM_tbufferedfilestream_DOT_fillbuffer(self)) { 
    result = _P3char('\032');
  } else {
    result = ValueCast(SYSTEM_ansichar,(*self->
      GMSSTRM_tbufferedfilestream_DOT_bufptr)[self->
      GMSSTRM_tbufferedfilestream_DOT_nrread]);
    _P3inc0(self->GMSSTRM_tbufferedfilestream_DOT_nrread);
  } 
  return result;
}  /* readcharacter */

/* unit gmsstrm */
void _Init_Module_gmsstrm(void)
{
} /* _Init_Module_gmsstrm */

void _Final_Module_gmsstrm(void)
{
} /* _Final_Module_gmsstrm */

