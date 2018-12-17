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
#include "gxdefs.h"
#include "gmsglob.h"
#include "gmsgen.h"
#include "strutilx.h"
#include "gmsobj.h"
#include "gmsdata.h"
#include "clibtypes.h"
#include "gmslibname.h"
#include "xcompress.h"
#include "gmsstrm.h"
#include "strhash.h"
#include "gmsheapnew.h"
#include "datastorage.h"
#include "gdlaudit.h"
#include "runner.h"
#include "gxfile.h"

_P3STR_7 GXFILE_baduel_prefix = {4,'?','L','_','_'};
_P3STR_7 GXFILE_badstr_prefix = {6,'?','S','t','r','_','_'};
_P3STR_15 GXFILE_strgdxcompress = {11,'G','D','X','C','O','M','P','R','E','S','S'};
_P3STR_15 GXFILE_strgdxconvert = {10,'G','D','X','C','O','N','V','E','R','T'};

void * const GXFILE_tdfilter_VT[] = {(void*)&
  GXFILE_tdfilter_DOT_destroy};

/* Class descriptor for 'tdfilter' */
const SYSTEM_classdescriptor_t GXFILE_tdfilter_CD = {
  _P3str1("\010tdfilter"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GXFILE_tdfilter_OD), GXFILE_tdfilter_VT, NULL};


void * const GXFILE_tfilterlist_VT[] = {(void*)&
  GXFILE_tfilterlist_DOT_destroy};

/* Class descriptor for 'tfilterlist' */
const SYSTEM_classdescriptor_t GXFILE_tfilterlist_CD = {
  _P3str1("\013tfilterlist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GXFILE_tfilterlist_OD), GXFILE_tfilterlist_VT, NULL};

GXFILE_tgxmodeset GXFILE_anywritemode = {252,1,0};
GXFILE_tgxmodeset GXFILE_anyreadmode = {2,240,0};

void * const GXFILE_tintegermapping_VT[] = {(void*)&
  GXFILE_tintegermapping_DOT_destroy};

/* Class descriptor for 'tintegermapping' */
const SYSTEM_classdescriptor_t GXFILE_tintegermapping_CD = {
  _P3str1("\017tintegermapping"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GXFILE_tintegermapping_OD), GXFILE_tintegermapping_VT, NULL};


void * const GXFILE_tueltable_VT[] = {(void*)&
  GXFILE_tueltable_DOT_destroy, (void*)&STRHASH_txstrhashlist_DOT_hash, (void*)&
  STRHASH_txstrhashlist_DOT_entryequal, (void*)&
  STRHASH_txstrhashlist_DOT_compare, (void*)&
  STRHASH_txstrhashlist_DOT_freeitem};

/* Class descriptor for 'tueltable' */
const SYSTEM_classdescriptor_t GXFILE_tueltable_CD = {
  _P3str1("\011tueltable"), 
  &STRHASH_txstrhashlist_CD, NULL, 0, 
  sizeof(GXFILE_tueltable_OD), GXFILE_tueltable_VT, NULL};


void * const GXFILE_tacronym_VT[] = {(void*)&
  GXFILE_tacronym_DOT_destroy};

/* Class descriptor for 'tacronym' */
const SYSTEM_classdescriptor_t GXFILE_tacronym_CD = {
  _P3str1("\010tacronym"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GXFILE_tacronym_OD), GXFILE_tacronym_VT, NULL};


void * const GXFILE_tacronymlist_VT[] = {(void*)&
  GXFILE_tacronymlist_DOT_destroy};

/* Class descriptor for 'tacronymlist' */
const SYSTEM_classdescriptor_t GXFILE_tacronymlist_CD = {
  _P3str1("\014tacronymlist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GXFILE_tacronymlist_OD), GXFILE_tacronymlist_VT, NULL};


void * const GXFILE_tgxfileobj_VT[] = {(void*)&
  GXFILE_tgxfileobj_DOT_destroy};

/* Class descriptor for 'tgxfileobj' */
const SYSTEM_classdescriptor_t GXFILE_tgxfileobj_CD = {
  _P3str1("\012tgxfileobj"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GXFILE_tgxfileobj_OD), GXFILE_tgxfileobj_VT, NULL};

SYSTEM_shortstring GXFILE_dllloadpath;
typedef struct GXFILE_uint64_S {
  union{
    struct{
      SYSTEM_int64 i;
    } _c1;
    struct{
      SYSTEM_pointer p;
    } _c2;
  } _u;
} GXFILE_uint64;


Function(SYSTEM_integer )  STDCALL 
  GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_fc(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *indx,
  const SYSTEM_double *vals,
  SYSTEM_pointer uptr)
{
  SYSTEM_integer result;
  GXDEFS_tdatastorefiltproc_f local_gdxdatareadrawfastfilt_dp;
  GXFILE_uint64 local_uptr;

  if (self->GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_callbyref) {
    PointerCast(SYSTEM_pointer,&local_gdxdatareadrawfastfilt_dp) = ValueCast(
      SYSTEM_pointer,self->
      GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp);
    local_uptr._u._c1.i = 0;
    local_uptr._u._c2.p = uptr;
    result = (*local_gdxdatareadrawfastfilt_dp)(indx,vals,&local_uptr.
      _u._c1.i);
  } else 
    result = (*self->GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp)(
      indx,vals,uptr);
  return result;
}  /* gdxdatareadrawfastfilt_dp_fc */

Procedure  STDCALL GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_fc(
  GXFILE_tgxfileobj self,
  SYSTEM_integer rawindex,
  SYSTEM_integer mappedindex,
  SYSTEM_pointer uptr)
{
  GXDEFS_tdomainindexproc_f local_gdxgetdomainelements_dp;
  GXFILE_uint64 local_uptr;

  if (self->GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_callbyref) {
    PointerCast(SYSTEM_pointer,&local_gdxgetdomainelements_dp) = ValueCast(
      SYSTEM_pointer,self->
      GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp);
    local_uptr._u._c1.i = 0;
    local_uptr._u._c2.p = uptr;
    (*local_gdxgetdomainelements_dp)(&rawindex,&mappedindex,&
      local_uptr._u._c1.i);
  } else 
    (*self->GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp)(rawindex,
      mappedindex,uptr);
}  /* gdxgetdomainelements_dp_fc */
cnstdef {GXFILE_version = 7};
cnstdef {GXFILE_gdxheadernr = 123};
static _P3STR_7 GXFILE_gdxheaderid = {7,'G','A','M','S','G','D','X'};
static SYSTEM_double GXFILE_bigfloat = 1e20;
cnstdef {GXFILE_mark_boi = 19510624};
static _P3STR_7 GXFILE_mark_uel = {5,'_','U','E','L','_'};
static _P3STR_7 GXFILE_mark_symb = {6,'_','S','Y','M','B','_'};
static _P3STR_7 GXFILE_mark_data = {6,'_','D','A','T','A','_'};
static _P3STR_7 GXFILE_mark_sett = {6,'_','S','E','T','T','_'};
static _P3STR_7 GXFILE_mark_acro = {6,'_','A','C','R','O','_'};
static _P3STR_7 GXFILE_mark_doms = {6,'_','D','O','M','S','_'};
cnstdef {GXFILE_mytrue = 1};
cnstdef {GXFILE_myfalse = 0};
cnstdef {GXFILE_index_initial =  -256};
typedef _P3STR_15 _arr_7GXFILE[18];
static _arr_7GXFILE GXFILE_fmode_str = {{11,'F','i','l','e','N','o','t','O','p','e','n'}, {11,'R','e','a','d','C','o','m','m','a','n','d'}, {12,'W','r','i','t','e','C','o','m','m','a','n','d'}, {13,'W','r','i','t','e','-','D','o','m','-','R','a','w'}, {13,'W','r','i','t','e','-','D','o','m','-','M','a','p'}, {13,'W','r','i','t','e','-','D','o','m','-','S','t','r'}, {9,'W','r','i','t','e','-','R','a','w'}, {9,'W','r','i','t','e','-','M','a','p'}, {9,'W','r','i','t','e','-','S','t','r'}, {9,'R','e','g','i','s','-','R','a','w'}, {9,'R','e','g','i','s','-','M','a','p'}, {9,'R','e','g','i','s','-','S','t','r'}, {8,'R','e','a','d','-','R','a','w'}, {8,'R','e','a','d','-','M','a','p'}, {9,'R','e','a','d','_','M','a','p','R'}, {8,'R','e','a','d','-','S','t','r'}, {12,'R','e','g','i','s','-','F','i','l','t','e','r'}, {10,'R','e','a','d','-','S','l','i','c','e'}};
cnstdef {GXFILE_err_noerror = 0};
cnstdef {GXFILE_err_nofile =  -100000};
cnstdef {GXFILE_err_fileerror =  -100001};
cnstdef {GXFILE_err_badmode =  -100002};
cnstdef {GXFILE_err_baddimension =  -100003};
cnstdef {GXFILE_err_badelementindex =  -100004};
cnstdef {GXFILE_err_badsymbolindex =  -100005};
cnstdef {GXFILE_err_elementsequence =  -100006};
cnstdef {GXFILE_err_duplicatesymbol =  -100007};
cnstdef {GXFILE_err_datanotsorted =  -100008};
cnstdef {GXFILE_err_dataduplicate =  -100009};
cnstdef {GXFILE_err_unknownfilter =  -100010};
cnstdef {GXFILE_err_badstringformat =  -100011};
cnstdef {GXFILE_err_badidentformat =  -100012};
cnstdef {GXFILE_err_uelconflict =  -100013};
cnstdef {GXFILE_err_duplicatespecval =  -100014};
cnstdef {GXFILE_err_baderrorrecord =  -100015};
cnstdef {GXFILE_err_duplicateuel =  -100016};
cnstdef {GXFILE_err_baduelstr =  -100017};
cnstdef {GXFILE_err_undefuel =  -100018};
cnstdef {GXFILE_err_uelsecondwrite =  -100019};
cnstdef {GXFILE_err_uelnotempty =  -100020};
cnstdef {GXFILE_err_bad_filter_nr =  -100021};
cnstdef {GXFILE_err_bad_filter_indx =  -100022};
cnstdef {GXFILE_err_filter_unmapped =  -100023};
cnstdef {GXFILE_err_obsolete_function =  -100024};
cnstdef {GXFILE_err_rawnotsorted =  -100025};
cnstdef {GXFILE_err_bad_alias_dim =  -100026};
cnstdef {GXFILE_err_baddatamarker_data =  -100029};
cnstdef {GXFILE_err_baddatamarker_dim =  -100030};
cnstdef {GXFILE_err_open_boi =  -100031};
cnstdef {GXFILE_err_open_fileheader =  -100032};
cnstdef {GXFILE_err_open_fileversion =  -100033};
cnstdef {GXFILE_err_open_filemarker =  -100034};
cnstdef {GXFILE_err_open_symbolmarker1 =  -100035};
cnstdef {GXFILE_err_open_symbolmarker2 =  -100036};
cnstdef {GXFILE_err_open_uelmarker1 =  -100037};
cnstdef {GXFILE_err_open_uelmarker2 =  -100038};
cnstdef {GXFILE_err_open_textmarker1 =  -100039};
cnstdef {GXFILE_err_open_textmarker2 =  -100040};
cnstdef {GXFILE_err_baddataformat =  -100041};
cnstdef {GXFILE_err_next_error =  -100042};
cnstdef {GXFILE_err_out_of_memory =  -100043};
cnstdef {GXFILE_err_zlib_not_found =  -100044};
cnstdef {GXFILE_err_open_acromarker1 =  -100045};
cnstdef {GXFILE_err_open_acromarker2 =  -100046};
cnstdef {GXFILE_err_badacroindex =  -100047};
cnstdef {GXFILE_err_badacronumber =  -100048};
cnstdef {GXFILE_err_badacroname =  -100049};
cnstdef {GXFILE_err_acrodupemap =  -100050};
cnstdef {GXFILE_err_acrobadaddition =  -100051};
cnstdef {GXFILE_err_unknowndomain =  -100052};
cnstdef {GXFILE_err_baddomain =  -100053};
cnstdef {GXFILE_err_nodomaindata =  -100054};
cnstdef {GXFILE_err_aliassetexpected =  -100055};
cnstdef {GXFILE_err_baddatatype =  -100056};
cnstdef {GXFILE_err_nosymbolforcomment =  -100057};
cnstdef {GXFILE_err_domainviolation =  -100058};
cnstdef {GXFILE_err_filealreadyopen =  -100059};
cnstdef {GXFILE_err_filetooldforappend =  -100060};
cnstdef {GXFILE_err_open_domsmarker1 =  -100061};
cnstdef {GXFILE_err_open_domsmarker2 =  -100062};
cnstdef {GXFILE_err_open_domsmarker3 =  -100063};
cnstdef {GXFILE_err_gdxcopy =  -100100};
cnstdef {GXFILE_err_parameter =  -100101};
cnstdef {GXFILE_err_dll_not_found =  -100102};
cnstdef {GXFILE_err_create_dir =  -100103};
cnstdef {GXFILE_err_file_open =  -100104};
cnstdef {GXFILE_err_file_write =  -100105};
cnstdef {GXFILE_err_uel_length =  -100106};
cnstdef {GXFILE_err_uel_register =  -100107};
cnstdef {GXFILE_err_expl_text =  -100108};
cnstdef {GXFILE_err_dimension =  -100109};
cnstdef {GXFILE_err_write_symbol =  -100110};
cnstdef {GXFILE_err_close_file =  -100111};
cnstdef {GXFILE_err_cannot_delete =  -100112};
cnstdef {GXFILE_err_cannot_rename =  -100113};

static Function(SYSTEM_integer ) GXFILE_getenvcompressflag(void)
{
  SYSTEM_integer result;
  SYSTEM_shortstring s;

  SYSUTILS_P3_getenvironmentvariable(s,255,GXFILE_strgdxcompress);
  {
    SYSTEM_shortstring _t2;

    STRUTILX_uppercase(s,255,SYSTEM_copy(_t2,255,s,1,1));
  }
  if (_P3strcmpE(s,_P3str1("\000")) || _P3stccmpE(s,_P3char('N')) || 
    _P3stccmpE(s,_P3char('0'))) { 
    result = 0;
  } else 
    result = 1;
  return result;
}  /* getenvcompressflag */

static Function(SYSTEM_integer ) GXFILE_convertgdxfile(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *mycomp)
{
  SYSTEM_integer result;
  SYSTEM_shortstring conv;
  SYSTEM_shortstring comp;
  RUNNER_trunner r;

  result = 0;
  {
    SYSTEM_shortstring _t2;
    SYSTEM_shortstring _t3;

    SYSUTILS_P3_trim(conv,255,STRUTILX_uppercase(_t2,255,
      SYSUTILS_P3_getenvironmentvariable(_t3,255,
      GXFILE_strgdxconvert)));
  }
  if (_P3strcmpE(conv,_P3str1("\000"))) 
    _P3strcpy(conv,255,_P3str1("\002V7"));
  if (_P3strcmpE(conv,_P3str1("\002V5"))) { 
    _P3strclr(comp);
  } else 
    if (GXFILE_getenvcompressflag() == 0) { 
      _P3strcpy(comp,255,_P3str1("\001U"));
    } else 
      _P3strcpy(comp,255,_P3str1("\001C"));
  {
    _P3STR_255 _t1;
    _P3STR_255 _t2;

    if (STRUTILX_struequal(_P3strcat(_t1,255,conv,comp),_P3strcat(
      _t2,255,_P3str1("\002V7"),mycomp))) 
      return result;
  }
  r = ValueCast(RUNNER_trunner,RUNNER_trunner_DOT_create(ValueCast(
    RUNNER_trunner,_P3alloc_object(&RUNNER_trunner_CD))));
  if (_P3strcmpE(GXFILE_dllloadpath,_P3str1("\000"))) { 
    RUNNER_trunner_DOT_setexecutable(r,_P3str1("\007gdxcopy"));
  } else 
    {
      _P3STR_3 _t1;
      _P3STR_255 _t2;
      _P3STR_255 _t3;

      RUNNER_trunner_DOT_setexecutable(r,_P3strcat(_t3,255,
        _P3strcat(_t2,255,GXFILE_dllloadpath,_P3ch2str(_t1,1,
        SYSUTILS_P3_pathdelim)),_P3str1("\007gdxcopy")));
    }
  {
    _P3STR_255 _t1;
    _P3STR_255 _t2;

    RUNNER_trunner_DOT_paramsadd(r,_P3strcat(_t2,255,_P3strcat(_t1,255,_P3str1("\001-"),
      conv),comp));
  }
  RUNNER_trunner_DOT_paramsadd(r,_P3str1("\010-Replace"));
  RUNNER_trunner_DOT_paramsadd(r,fn);
  result = RUNNER_trunner_DOT_startandwait(r);
  if (result == 0) 
    if (r->RUNNER_trunner_DOT_fprogrc != 0) 
      result = GXFILE_err_gdxcopy - r->RUNNER_trunner_DOT_fprogrc;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,r));
  return result;
}  /* convertgdxfile */

static Function(SYSTEM_ansichar *) GXFILE_makegoodexpltext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer i;
  SYSTEM_ansichar q;
  SYSTEM_ansichar ch;

  q = _P3char('\000');
  _P3setlength(result,SYSTEM_length(s),255);
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      ch = s[i];
      if (!_P3SET_in_1(ch,_P3char('\''),_P3SET_equal(ch,_P3char('\"')))) { 
        if (ch < _P3char(' ')) 
          ch = _P3char('?');
      } else {
        if (q == _P3char('\000')) 
          q = ch;
        ch = q;
      } 
      result[i] = ch;
    
    } while (i++ !=  _stop);

  }
  return result;
}  /* makegoodexpltext */

static Function(SYSTEM_boolean ) GXFILE_isgoodident(
  const SYSTEM_ansichar *s)
{
  SYSTEM_boolean result;
  SYSTEM_integer n;

  result = SYSTEM_length(s) > 0 && SYSTEM_length(s) <= 
    GMSSPECS_maxnamelen && _P3SET_in_2(s[1],_P3char('A'),_P3char('Z'),
    _P3SET_in_3(s[1],_P3char('a'),_P3char('z')));
  if (result) 
    { register SYSTEM_int32 _stop = SYSTEM_length(s);
      if ((n = 2) <=  _stop) do {
        if (!_P3SET_i(122,s[n],
          _P3set1("\0\0\0\0\0\0\377\003\376\377\377\207\376\377\377\007"))) {
          result = SYSTEM_false;
          SYSTEM_break(BRK_1);
        } 
CNT_1:;
      } while (n++ !=  _stop);
BRK_1:;

    }
  return result;
}  /* isgoodident */

static Function(SYSTEM_integer ) GXFILE_imax(
  SYSTEM_integer a,
  SYSTEM_integer b)
{
  SYSTEM_integer result;

  if (a >= b) { 
    result = a;
  } else 
    result = b;
  return result;
}  /* imax */

static Function(GXFILE_tgdxelemsize ) GXFILE_getintegersize(
  SYSTEM_integer n)
{
  GXFILE_tgdxelemsize result;

  if (n <= 0) { 
    result = GXFILE_sz_integer;
  } else 
    if (n <= 255) { 
      result = GXFILE_sz_byte;
    } else 
      if (n <= 65535) { 
        result = GXFILE_sz_word;
      } else 
        result = GXFILE_sz_integer;
  return result;
}  /* getintegersize */

Constructor(GXFILE_tintegermapping ) GXFILE_tintegermapping_DOT_create(
  GXFILE_tintegermapping self)
{
  ValueCast(GXFILE_tintegermapping,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GXFILE_tintegermapping_DOT_fcapacity = 0;
  self->GXFILE_tintegermapping_DOT_fhighestindex = 0;
  self->GXFILE_tintegermapping_DOT_pmap = NULL;
  return self;
}  /* create */

Destructor(GXFILE_tintegermapping ) GXFILE_tintegermapping_DOT_destroy(
  GXFILE_tintegermapping self)
{
  if (self->GXFILE_tintegermapping_DOT_pmap != NULL) 
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      GXFILE_tintegermapping_DOT_pmap),0);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GXFILE_tintegermapping_DOT_setmapping(
  GXFILE_tintegermapping self,
  SYSTEM_integer f,
  SYSTEM_integer t)
{
  SYSTEM_integer n;
  SYSTEM_integer delta;

  if (f >= self->GXFILE_tintegermapping_DOT_fcapacity) {
    delta = 0;
    do {
      if (self->GXFILE_tintegermapping_DOT_fcapacity == 0) { 
        _P3inc1(delta,1024);
      } else 
        if (self->GXFILE_tintegermapping_DOT_fcapacity <= 32768) { 
          _P3inc1(delta,self->GXFILE_tintegermapping_DOT_fcapacity);
        } else 
          _P3inc1(delta,self->GXFILE_tintegermapping_DOT_fcapacity /  4);
    } while (!(f < self->GXFILE_tintegermapping_DOT_fcapacity + delta));
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      GXFILE_tintegermapping_DOT_pmap),(self->
      GXFILE_tintegermapping_DOT_fcapacity + delta) * sizeof(
      SYSTEM_longint));
    { register SYSTEM_int32 _stop = self->
        GXFILE_tintegermapping_DOT_fcapacity + delta - 1;
      if ((n = self->GXFILE_tintegermapping_DOT_fcapacity) <=  _stop) do {
        (*self->GXFILE_tintegermapping_DOT_pmap)[n] =  -1;
      } while (n++ !=  _stop);

    }
    self->GXFILE_tintegermapping_DOT_fcapacity = self->
      GXFILE_tintegermapping_DOT_fcapacity + delta;
  } 
  (*self->GXFILE_tintegermapping_DOT_pmap)[f] = t;
  if (f > self->GXFILE_tintegermapping_DOT_fhighestindex) 
    self->GXFILE_tintegermapping_DOT_fhighestindex = f;
}  /* setmapping */

Function(SYSTEM_integer ) GXFILE_tintegermapping_DOT_getmapping(
  GXFILE_tintegermapping self,
  SYSTEM_integer f)
{
  SYSTEM_integer result;

  if (f >= 0 && f < self->GXFILE_tintegermapping_DOT_fcapacity) { 
    result = (*self->GXFILE_tintegermapping_DOT_pmap)[f];
  } else 
    result =  -1;
  return result;
}  /* getmapping */

Function(SYSTEM_int64 ) GXFILE_tintegermapping_DOT_memoryused(
  GXFILE_tintegermapping self)
{
  SYSTEM_int64 result;

  result = self->GXFILE_tintegermapping_DOT_fcapacity * sizeof(
    SYSTEM_longint);
  return result;
}  /* memoryused */

Constructor(GXFILE_tueltable ) GXFILE_tueltable_DOT_create(
  GXFILE_tueltable self)
{
  ValueCast(GXFILE_tueltable,STRHASH_txstrhashlist_DOT_create(ValueCast(
    STRHASH_txstrhashlist,self)));
  self->STRHASH_txstrhashlist_DOT_onebased = SYSTEM_true;
  self->GXFILE_tueltable_DOT_usruel2ent = ValueCast(
    GXFILE_tintegermapping,GXFILE_tintegermapping_DOT_create(ValueCast(
    GXFILE_tintegermapping,_P3alloc_object(&GXFILE_tintegermapping_CD))));
  GXFILE_tueltable_DOT_resetmaptouserstatus(self);
  return self;
}  /* create */

Destructor(GXFILE_tueltable ) GXFILE_tueltable_DOT_destroy(
  GXFILE_tueltable self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tueltable_DOT_usruel2ent));
  STRHASH_txstrhashlist_DOT_destroy(ValueCast(STRHASH_txstrhashlist,
    self));
  return self;
}  /* destroy */

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_addusrnew(
  GXFILE_tueltable self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer en;

  en = STRHASH_txstrhashlist_DOT_addobject(ValueCast(
    STRHASH_txstrhashlist,self),s,ValueCast(SYSTEM_tobject,
    GMSOBJ_copyint2ptr( -1)));
  result = GMSOBJ_copyptr2int(STRHASH_txstrhashlist_DOT_getobject(ValueCast(
    STRHASH_txstrhashlist,self),en));
  if (result < 0) {
    result = self->GXFILE_tueltable_DOT_usruel2ent->
      GXFILE_tintegermapping_DOT_fhighestindex + 1;
    STRHASH_txstrhashlist_DOT_setobject(ValueCast(
      STRHASH_txstrhashlist,self),en,ValueCast(SYSTEM_tobject,
      GMSOBJ_copyint2ptr(result)));
    GXFILE_tintegermapping_DOT_setmapping(self->
      GXFILE_tueltable_DOT_usruel2ent,result,en);
  } 
  GXFILE_tueltable_DOT_resetmaptouserstatus(self);
  return result;
}  /* addusrnew */

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_addusrindxnew(
  GXFILE_tueltable self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer uelnr)
{
  SYSTEM_integer result;
  SYSTEM_integer en;

  en = STRHASH_txstrhashlist_DOT_addobject(ValueCast(
    STRHASH_txstrhashlist,self),s,ValueCast(SYSTEM_tobject,
    GMSOBJ_copyint2ptr( -1)));
  result = GMSOBJ_copyptr2int(STRHASH_txstrhashlist_DOT_getobject(ValueCast(
    STRHASH_txstrhashlist,self),en));
  if (result < 0) {
    result = uelnr;
    STRHASH_txstrhashlist_DOT_setobject(ValueCast(
      STRHASH_txstrhashlist,self),en,ValueCast(SYSTEM_tobject,
      GMSOBJ_copyint2ptr(result)));
    GXFILE_tintegermapping_DOT_setmapping(self->
      GXFILE_tueltable_DOT_usruel2ent,result,en);
  } else 
    if (result != uelnr) 
      result =  -1;
  GXFILE_tueltable_DOT_resetmaptouserstatus(self);
  return result;
}  /* addusrindxnew */

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_getusermap(
  GXFILE_tueltable self,
  SYSTEM_integer n)
{
  SYSTEM_integer result;

  result = GMSOBJ_copyptr2int(STRHASH_txstrhashlist_DOT_getobject(ValueCast(
    STRHASH_txstrhashlist,self),n));
  return result;
}  /* getusermap */

Procedure GXFILE_tueltable_DOT_loadfromstream(
  GXFILE_tueltable self,
  GMSSTRM_txstream s)
{
  SYSTEM_integer n;

  STRHASH_txstrhashlist_DOT_loadfromstream(ValueCast(
    STRHASH_txstrhashlist,self),s);
  if (self->GXFILE_tueltable_DOT_usruel2ent != NULL) {
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
      GXFILE_tueltable_DOT_usruel2ent));
    self->GXFILE_tueltable_DOT_usruel2ent = ValueCast(
      GXFILE_tintegermapping,GXFILE_tintegermapping_DOT_create(ValueCast(
      GXFILE_tintegermapping,_P3alloc_object(&
      GXFILE_tintegermapping_CD))));
  } 
  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_fcount;
    if ((n = 1) <=  _stop) do {
      STRHASH_txstrhashlist_DOT_setobject(ValueCast(
        STRHASH_txstrhashlist,self),n,ValueCast(SYSTEM_tobject,
        GMSOBJ_copyint2ptr( -1)));
    } while (n++ !=  _stop);

  }
  GXFILE_tueltable_DOT_resetmaptouserstatus(self);
}  /* loadfromstream */

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_newusruel(
  GXFILE_tueltable self,
  SYSTEM_integer en)
{
  SYSTEM_integer result;

  result = GMSOBJ_copyptr2int(STRHASH_txstrhashlist_DOT_getobject(ValueCast(
    STRHASH_txstrhashlist,self),en));
  if (result < 0) {
    result = self->GXFILE_tueltable_DOT_usruel2ent->
      GXFILE_tintegermapping_DOT_fhighestindex + 1;
    STRHASH_txstrhashlist_DOT_setobject(ValueCast(
      STRHASH_txstrhashlist,self),en,ValueCast(SYSTEM_tobject,
      GMSOBJ_copyint2ptr(result)));
    GXFILE_tintegermapping_DOT_setmapping(self->
      GXFILE_tueltable_DOT_usruel2ent,result,en);
  } 
  GXFILE_tueltable_DOT_resetmaptouserstatus(self);
  return result;
}  /* newusruel */

Function(GXFILE_tuelusermapstatus ) 
  GXFILE_tueltable_DOT_getmaptouserstatus(
  GXFILE_tueltable self)
{
  GXFILE_tuelusermapstatus result;
  SYSTEM_integer n;
  SYSTEM_integer v;
  SYSTEM_integer lv;
  SYSTEM_boolean c;

  if (self->GXFILE_tueltable_DOT_fmaptouserstatus == 
    GXFILE_map_unknown) {
    lv =  -1;
    c = SYSTEM_true;
    self->GXFILE_tueltable_DOT_fmaptouserstatus = GXFILE_map_sortgrow;
    { register SYSTEM_int32 _stop = self->
        STRHASH_txstrhashlist_DOT_fcount;
      if ((n = 1) <=  _stop) do {
        v = GXFILE_tueltable_DOT_getusermap(self,n);
        if (v < 0) { 
          c = SYSTEM_false;
        } else 
          if (v > lv) {
            lv = v;
            if (!c) 
              self->GXFILE_tueltable_DOT_fmaptouserstatus = 
                GXFILE_map_sorted;
          } else {
            self->GXFILE_tueltable_DOT_fmaptouserstatus = 
              GXFILE_map_unsorted;
            SYSTEM_break(BRK_2);
          } 
      
CNT_2:;
      } while (n++ !=  _stop);
BRK_2:;

    }
    if (self->GXFILE_tueltable_DOT_fmaptouserstatus == 
      GXFILE_map_sortgrow && c) 
      self->GXFILE_tueltable_DOT_fmaptouserstatus = 
        GXFILE_map_sortfull;
  } 
  result = self->GXFILE_tueltable_DOT_fmaptouserstatus;
  return result;
}  /* getmaptouserstatus */

Procedure GXFILE_tueltable_DOT_resetmaptouserstatus(
  GXFILE_tueltable self)
{
  self->GXFILE_tueltable_DOT_fmaptouserstatus = GXFILE_map_unknown;
}  /* resetmaptouserstatus */

Procedure GXFILE_tgxfileobj_DOT_initerrors(
  GXFILE_tgxfileobj self)
{
  self->GXFILE_tgxfileobj_DOT_errcnt = 0;
  self->GXFILE_tgxfileobj_DOT_errcnttotal = 0;
  self->GXFILE_tgxfileobj_DOT_lasterror = GXFILE_err_noerror;
  self->GXFILE_tgxfileobj_DOT_lastreperror = GXFILE_err_noerror;
}  /* initerrors */

Procedure GXFILE_tgxfileobj_DOT_seterror(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n)
{
  if (n == 0) 
    return;
  if (self->GXFILE_tgxfileobj_DOT_lasterror == 0) 
    self->GXFILE_tgxfileobj_DOT_lasterror = n;
  _P3inc0(self->GXFILE_tgxfileobj_DOT_errcnt);
  _P3inc0(self->GXFILE_tgxfileobj_DOT_errcnttotal);
}  /* seterror */

Procedure GXFILE_tgxfileobj_DOT_reporterror(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n)
{
  SYSTEM_shortstring s;

  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_errors && n != 
    self->GXFILE_tgxfileobj_DOT_lastreperror) {
    if (_P3strcmpN(self->GXFILE_tgxfileobj_DOT_majcontext,_P3str1("\000"))) {
      _Iplus_bgn();
      _P3write_s0(_P3str1("\024Error after call to "));
      _P3write_s0(self->GXFILE_tgxfileobj_DOT_majcontext);
      _P3writeln();
      _Iplus_end();
    } 
    GXFILE_tgxfileobj_DOT_gdxerrorstr(self,n,s);
    _Iplus_bgn();
    _P3write_s0(_P3str1("\010Error = "));
    _P3write_i0(n);
    _P3write_s0(_P3str1("\003 : "));
    _P3write_s0(s);
    _P3writeln();
    _Iplus_end();
  } 
  GXFILE_tgxfileobj_DOT_seterror(self,n);
  self->GXFILE_tgxfileobj_DOT_lastreperror = n;
}  /* reporterror */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_errorcondition(
  GXFILE_tgxfileobj self,
  SYSTEM_boolean cnd,
  SYSTEM_integer n)
{
  SYSTEM_boolean result;

  result = !cnd;
  if (result) 
    GXFILE_tgxfileobj_DOT_reporterror(self,n);
  return result;
}  /* errorcondition */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenwriteex(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *producer,
  SYSTEM_integer compr,
  SYSTEM_integer *errnr)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_fmode != GXFILE_f_not_open) {
    *errnr = GXFILE_err_filealreadyopen;
    return result;
  } 
  self->GXFILE_tgxfileobj_DOT_ffile = ValueCast(
    GMSSTRM_tmibufferedstream,
    GMSSTRM_tmibufferedstream_DOT_createwithpath(ValueCast(
    GMSSTRM_tmibufferedstream,_P3alloc_object(&
    GMSSTRM_tmibufferedstream_CD)),filename,GMSSTRM_fmcreate,
    GXFILE_dllloadpath));
  *errnr = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
    GMSSTRM_txfilestream,self->GXFILE_tgxfileobj_DOT_ffile));
  if (*errnr != 0) {
    SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_ffile);
    if (*errnr == 5) 
      *errnr = GXFILE_err_zlib_not_found;
    self->GXFILE_tgxfileobj_DOT_lasterror = *errnr;
    return result;
  } 
  if (compr != 0) 
    if (self->GXFILE_tgxfileobj_DOT_ffile->
      GMSSTRM_tbufferedfilestream_DOT_fcancompress) { 
      compr = 1;
    } else 
      compr = 0;
  self->GXFILE_tgxfileobj_DOT_fcomprlev = compr;
  self->GXFILE_tgxfileobj_DOT_compressout = compr > 0;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_not_open;
  self->GXFILE_tgxfileobj_DOT_readptr = NULL;
  _P3strcpy(self->GXFILE_tgxfileobj_DOT_majcontext,63,_P3str1("\011OpenWrite"));
  self->GXFILE_tgxfileobj_DOT_tracelevel = GXFILE_trl_none;
  GXFILE_tgxfileobj_DOT_initerrors(self);
  self->GXFILE_tgxfileobj_DOT_namelist = ValueCast(
    STRHASH_txstrhashlist,STRHASH_txstrhashlist_DOT_create(ValueCast(
    STRHASH_txstrhashlist,_P3alloc_object(&STRHASH_txstrhashlist_CD))));
  self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_onebased = SYSTEM_true;
  self->GXFILE_tgxfileobj_DOT_ueltable = ValueCast(GXFILE_tueltable,
    GXFILE_tueltable_DOT_create(ValueCast(GXFILE_tueltable,
    _P3alloc_object(&GXFILE_tueltable_CD))));
  self->GXFILE_tgxfileobj_DOT_acronymlist = ValueCast(
    GXFILE_tacronymlist,GXFILE_tacronymlist_DOT_create(ValueCast(
    GXFILE_tacronymlist,_P3alloc_object(&GXFILE_tacronymlist_CD))));
  self->GXFILE_tgxfileobj_DOT_filterlist = ValueCast(
    GXFILE_tfilterlist,GXFILE_tfilterlist_DOT_create(ValueCast(
    GXFILE_tfilterlist,_P3alloc_object(&GXFILE_tfilterlist_CD))));
  GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),GXFILE_gdxheadernr);
  GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),GXFILE_gdxheaderid);
  self->GXFILE_tgxfileobj_DOT_versionread = GXFILE_version;
  GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),self->
    GXFILE_tgxfileobj_DOT_versionread);
  GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),compr);
  GDLAUDIT_gdlgetauditline(self->GXFILE_tgxfileobj_DOT_filesystemid,255);
  GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),self->
    GXFILE_tgxfileobj_DOT_filesystemid);
  _P3strcpy(self->GXFILE_tgxfileobj_DOT_fproducer,255,producer);
  _P3strclr(self->GXFILE_tgxfileobj_DOT_fproducer2);
  GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),self->GXFILE_tgxfileobj_DOT_fproducer);
  self->GXFILE_tgxfileobj_DOT_majorindexposition = VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
    GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile)));
  for (n = 1;n <= (SYSTEM_int32)10;++n) {
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),0);
  }
  self->GXFILE_tgxfileobj_DOT_settextlist = ValueCast(GMSOBJ_txstrpool,
    GMSOBJ_txhashedstringlist_DOT_create(ValueCast(
    GMSOBJ_txhashedstringlist,_P3alloc_object(&GMSOBJ_txstrpool_CD))));
  GMSOBJ_txhashedstringlist_DOT_add(ValueCast(
    GMSOBJ_txhashedstringlist,self->GXFILE_tgxfileobj_DOT_settextlist),_P3str1("\000"));
  GXFILE_tgxfileobj_DOT_gdxresetspecialvalues(self);
  self->GXFILE_tgxfileobj_DOT_nextwriteposition = VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
    GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile)));
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_init;
  self->GXFILE_tgxfileobj_DOT_fstatus = GXFILE_stat_write;
  result = GXFILE_mytrue;
  self->GXFILE_tgxfileobj_DOT_domainstrlist = ValueCast(
    STRHASH_txstrhashlist,STRHASH_txstrhashlist_DOT_create(ValueCast(
    STRHASH_txstrhashlist,_P3alloc_object(&STRHASH_txstrhashlist_CD))));
  self->GXFILE_tgxfileobj_DOT_domainstrlist->
    STRHASH_txstrhashlist_DOT_onebased = SYSTEM_true;
  return result;
}  /* gdxopenwriteex */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenwrite(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *producer,
  SYSTEM_integer *errnr)
{
  SYSTEM_integer result;

  result = GXFILE_tgxfileobj_DOT_gdxopenwriteex(self,filename,producer,
    GXFILE_getenvcompressflag(),errnr);
  return result;
}  /* gdxopenwrite */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenread(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  SYSTEM_integer *errnr)
{
  SYSTEM_integer result;

  result = GXFILE_tgxfileobj_DOT_gdxopenreadxx(self,filename,
    GMSSTRM_fmopenread,errnr);
  return result;
}  /* gdxopenread */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenreadxx(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *afn,
  SYSTEM_integer filemode,
  SYSTEM_integer *errnr)
{
  SYSTEM_integer result;
  SYSTEM_byte b;
  SYSTEM_shortstring s;
  SYSTEM_integer n;
  SYSTEM_integer nrelem;
  SYSTEM_int64 uelpos;
  SYSTEM_int64 symbpos;
  SYSTEM_int64 settextpos;
  SYSTEM_int64 acronympos;
  SYSTEM_int64 domstrpos;
  SYSTEM_integer compr;
  SYSTEM_integer synr;
  SYSTEM_integer d;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_fmode != GXFILE_f_not_open) {
    *errnr = GXFILE_err_filealreadyopen;
    return result;
  } 
  _P3strcpy(self->GXFILE_tgxfileobj_DOT_majcontext,63,_P3str1("\010OpenRead"));
  self->GXFILE_tgxfileobj_DOT_tracelevel = GXFILE_trl_none;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_not_open;
  self->GXFILE_tgxfileobj_DOT_readptr = NULL;
  GXFILE_tgxfileobj_DOT_initerrors(self);
  self->GXFILE_tgxfileobj_DOT_ffile = ValueCast(
    GMSSTRM_tmibufferedstream,
    GMSSTRM_tmibufferedstream_DOT_createwithpath(ValueCast(
    GMSSTRM_tmibufferedstream,_P3alloc_object(&
    GMSSTRM_tmibufferedstream_CD)),afn,filemode,GXFILE_dllloadpath));
  *errnr = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
    GMSSTRM_txfilestream,self->GXFILE_tgxfileobj_DOT_ffile));
  if (*errnr != 0) 
    goto _Lfilenogood_31;
  if (GMSSTRM_tmibufferedstream_DOT_goodbyteorder(self->
    GXFILE_tgxfileobj_DOT_ffile) != 0) {
    *errnr = GXFILE_err_baddataformat;
    goto _Lfilenogood_31;
  } 
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,
    GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile)) == GXFILE_gdxheadernr,
    GXFILE_err_open_fileheader)) 
    goto _Lfileerrornr_31;
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_gdxheaderid),GXFILE_err_open_filemarker)) 
      goto _Lfileerrornr_31;
  }
  self->GXFILE_tgxfileobj_DOT_versionread = VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
    GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile)));
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,self->
    GXFILE_tgxfileobj_DOT_versionread <= 7,
    GXFILE_err_open_fileversion)) 
    goto _Lfileerrornr_31;
  if (self->GXFILE_tgxfileobj_DOT_versionread <= 5) { 
    compr = 0;
  } else 
    compr = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
  self->GXFILE_tgxfileobj_DOT_douncompress = compr > 0;
  if (self->GXFILE_tgxfileobj_DOT_douncompress && !self->
    GXFILE_tgxfileobj_DOT_ffile->
    GMSSTRM_tbufferedfilestream_DOT_fcancompress) {
    *errnr = GXFILE_err_zlib_not_found;
    goto _Lfilenogood_31;
  } 
  self->GXFILE_tgxfileobj_DOT_fcomprlev = compr;
  GMSSTRM_txstream_DOT_readstring(self->
    GXFILE_tgxfileobj_DOT_filesystemid,255,ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile));
  GMSSTRM_txstream_DOT_readstring(self->
    GXFILE_tgxfileobj_DOT_fproducer,255,ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile));
  _P3strclr(self->GXFILE_tgxfileobj_DOT_fproducer2);
  self->GXFILE_tgxfileobj_DOT_majorindexposition = VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
    GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile)));
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
    GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile))) == GXFILE_mark_boi,
    GXFILE_err_open_boi)) 
    goto _Lfileerrornr_31;
  acronympos = 0;
  domstrpos = 0;
  if (self->GXFILE_tgxfileobj_DOT_versionread <= 5) {
    symbpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    uelpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    settextpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    self->GXFILE_tgxfileobj_DOT_nextwriteposition = VirtMethodCall(ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
      GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
  } else {
    symbpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    uelpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    settextpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    if (self->GXFILE_tgxfileobj_DOT_versionread >= 7) {
      acronympos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
        GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
      self->GXFILE_tgxfileobj_DOT_nextwriteposition = VirtMethodCall(ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
        GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
      domstrpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
        GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    } 
  } 
  GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
    GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
    self->GXFILE_tgxfileobj_DOT_douncompress);
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),symbpos));
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_mark_symb),GXFILE_err_open_symbolmarker1)) 
      goto _Lfileerrornr_31;
  }
  nrelem = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
  self->GXFILE_tgxfileobj_DOT_namelist = ValueCast(
    STRHASH_txstrhashlist,STRHASH_txstrhashlist_DOT_create(ValueCast(
    STRHASH_txstrhashlist,_P3alloc_object(&STRHASH_txstrhashlist_CD))));
  self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_onebased = SYSTEM_true;
  self->GXFILE_tgxfileobj_DOT_acronymlist = ValueCast(
    GXFILE_tacronymlist,GXFILE_tacronymlist_DOT_create(ValueCast(
    GXFILE_tacronymlist,_P3alloc_object(&GXFILE_tacronymlist_CD))));
  self->GXFILE_tgxfileobj_DOT_filterlist = ValueCast(
    GXFILE_tfilterlist,GXFILE_tfilterlist_DOT_create(ValueCast(
    GXFILE_tfilterlist,_P3alloc_object(&GXFILE_tfilterlist_CD))));
  { register SYSTEM_int32 _stop = nrelem;
    if ((n = 1) <=  _stop) do {
      GMSSTRM_txstream_DOT_readstring(s,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile));
      _P3new(self->GXFILE_tgxfileobj_DOT_cursyptr);
      { register GXFILE_tgdxsymbrecord *_W2=self->
        GXFILE_tgxfileobj_DOT_cursyptr;
        if (self->GXFILE_tgxfileobj_DOT_versionread <= 5) { 
          _W2->sposition = VirtMethodCall(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile), 
            GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
            GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        } else 
          _W2->sposition = VirtMethodCall(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile), 
            GMSSTRM_txstream_DOT_readint64_T, 9, (ValueCast(
            GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        _W2->sdim = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
          GXFILE_tgxfileobj_DOT_ffile), 
          GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        b = GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,
          self->GXFILE_tgxfileobj_DOT_ffile));
        _W2->sdatatype = ValueCast(GMSSPECS_tgdxdatatype,b);
        _W2->suserinfo = VirtMethodCall(ValueCast(GMSSTRM_txstream,
          self->GXFILE_tgxfileobj_DOT_ffile), 
          GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        _W2->sdatacount = VirtMethodCall(ValueCast(GMSSTRM_txstream,
          self->GXFILE_tgxfileobj_DOT_ffile), 
          GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        _W2->serrors = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
          GXFILE_tgxfileobj_DOT_ffile), 
          GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        b = GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,
          self->GXFILE_tgxfileobj_DOT_ffile));
        _W2->ssettext = b != 0;
        GMSSTRM_txstream_DOT_readstring(_W2->sexpltxt,255,ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile));
        if (self->GXFILE_tgxfileobj_DOT_versionread <= 5) { 
          _W2->siscompressed = SYSTEM_false;
        } else {
          b = GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile));
          _W2->siscompressed = b != 0;
        } 
        _W2->sdomsymbols = NULL;
        _W2->scommentslist = NULL;
        if (self->GXFILE_tgxfileobj_DOT_versionread >= 7) {
          if (GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile)) != 0) {
            _P3getmem(_W2->sdomsymbols,(_W2->sdim + 1) * sizeof(
              SYSTEM_longint));
            { register SYSTEM_int32 _stop = _W2->sdim;
              if ((d = 1) <=  _stop) do {
                (*_W2->sdomsymbols)[d] = VirtMethodCall(ValueCast(
                  GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
                  GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
                  GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
              } while (d++ !=  _stop);

            }
          } 
          nrelem = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
            GXFILE_tgxfileobj_DOT_ffile), 
            GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
            GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
          if (nrelem > 0) {
            _W2->scommentslist = ValueCast(GMSOBJ_txstrings,
              SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,
              _P3alloc_object(&GMSOBJ_txstrings_CD))));
            while (nrelem > 0) {
              {
                SYSTEM_shortstring _t1;

                GMSOBJ_txstrings_DOT_add(_W2->scommentslist,
                  GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
                  GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
              }
              nrelem = nrelem - 1;
            
}
          } 
        } 
        _W2->ssetbitmap = NULL;
        _W2->sdomstrings = NULL;

      }
      self->GXFILE_tgxfileobj_DOT_cursyptr->ssynr = 
        STRHASH_txstrhashlist_DOT_storeobject(self->
        GXFILE_tgxfileobj_DOT_namelist,s,ValueCast(SYSTEM_tobject,self->
        GXFILE_tgxfileobj_DOT_cursyptr));
    
    } while (n++ !=  _stop);

  }
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_mark_symb),GXFILE_err_open_symbolmarker2)) 
      goto _Lfileerrornr_31;
  }
  GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
    GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
    self->GXFILE_tgxfileobj_DOT_douncompress);
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),uelpos));
  self->GXFILE_tgxfileobj_DOT_ueltable = ValueCast(GXFILE_tueltable,
    GXFILE_tueltable_DOT_create(ValueCast(GXFILE_tueltable,
    _P3alloc_object(&GXFILE_tueltable_CD))));
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_mark_uel),GXFILE_err_open_uelmarker1)) 
      goto _Lfileerrornr_31;
  }
  { register GXFILE_tueltable_OD *_W2=self->
    GXFILE_tgxfileobj_DOT_ueltable;
    nrelem = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    {
      SYSTEM_shortstring _t1;

      if (_P3strcmpE(SYSTEM_copy(_t1,255,self->
        GXFILE_tgxfileobj_DOT_filesystemid,16,4),_P3str1("\0042001"))) 
        nrelem = nrelem - 1;
    }
    while (_W2->STRHASH_txstrhashlist_DOT_fcount < nrelem) {

      {
        SYSTEM_shortstring _t1;

        STRHASH_txstrhashlist_DOT_storeobject(ValueCast(
          STRHASH_txstrhashlist,_W2),GMSSTRM_txstream_DOT_readstring(
          _t1,255,ValueCast(GMSSTRM_txstream,self->
          GXFILE_tgxfileobj_DOT_ffile)),ValueCast(SYSTEM_tobject,
          GMSOBJ_copyint2ptr( -1)));
      }
}
    self->GXFILE_tgxfileobj_DOT_uelcntorig = _W2->
      STRHASH_txstrhashlist_DOT_fcount;

  }
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_mark_uel),GXFILE_err_open_uelmarker2)) 
      goto _Lfileerrornr_31;
  }
  GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
    GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
    self->GXFILE_tgxfileobj_DOT_douncompress);
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),settextpos));
  self->GXFILE_tgxfileobj_DOT_settextlist = ValueCast(GMSOBJ_txstrpool,
    GMSOBJ_txhashedstringlist_DOT_create(ValueCast(
    GMSOBJ_txhashedstringlist,_P3alloc_object(&GMSOBJ_txstrpool_CD))));
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_mark_sett),GXFILE_err_open_textmarker1)) 
      goto _Lfileerrornr_31;
  }
  nrelem = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
  GMSOBJ_txcustomstringlist_DOT_setcapacity(ValueCast(
    GMSOBJ_txcustomstringlist,self->GXFILE_tgxfileobj_DOT_settextlist),
    nrelem);
  { register SYSTEM_int32 _stop = nrelem;
    if ((n = 1) <=  _stop) do {
      {
        SYSTEM_shortstring _t1;

        GMSOBJ_txhashedstringlist_DOT_add(ValueCast(
          GMSOBJ_txhashedstringlist,self->
          GXFILE_tgxfileobj_DOT_settextlist),
          GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
      }
    } while (n++ !=  _stop);

  }
  {
    SYSTEM_shortstring _t1;

    if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
      GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
      GXFILE_mark_sett),GXFILE_err_open_textmarker2)) 
      goto _Lfileerrornr_31;
  }
  if (self->GXFILE_tgxfileobj_DOT_versionread >= 7) {
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_douncompress);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),acronympos));
    {
      SYSTEM_shortstring _t1;

      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
        GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
        GXFILE_mark_acro),GXFILE_err_open_acromarker1)) 
        goto _Lfileerrornr_31;
    }
    GXFILE_tacronymlist_DOT_loadfromstream(self->
      GXFILE_tgxfileobj_DOT_acronymlist,ValueCast(GMSSTRM_txstream,
      self->GXFILE_tgxfileobj_DOT_ffile));
    {
      SYSTEM_shortstring _t1;

      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
        GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
        GXFILE_mark_acro),GXFILE_err_open_acromarker2)) 
        goto _Lfileerrornr_31;
    }
  } 
  self->GXFILE_tgxfileobj_DOT_domainstrlist = ValueCast(
    STRHASH_txstrhashlist,STRHASH_txstrhashlist_DOT_create(ValueCast(
    STRHASH_txstrhashlist,_P3alloc_object(&STRHASH_txstrhashlist_CD))));
  self->GXFILE_tgxfileobj_DOT_domainstrlist->
    STRHASH_txstrhashlist_DOT_onebased = SYSTEM_true;
  if (self->GXFILE_tgxfileobj_DOT_versionread >= 7 && domstrpos != 0) {
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_douncompress);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),domstrpos));
    {
      SYSTEM_shortstring _t1;

      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
        GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
        GXFILE_mark_doms),GXFILE_err_open_domsmarker1)) 
        goto _Lfileerrornr_31;
    }
    STRHASH_txstrhashlist_DOT_loadfromstream(self->
      GXFILE_tgxfileobj_DOT_domainstrlist,ValueCast(GMSSTRM_txstream,
      self->GXFILE_tgxfileobj_DOT_ffile));
    {
      SYSTEM_shortstring _t1;

      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
        GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
        GXFILE_mark_doms),GXFILE_err_open_domsmarker2)) 
        goto _Lfileerrornr_31;
    }
    while (SYSTEM_true) {
      synr = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
        GXFILE_tgxfileobj_DOT_ffile), 
        GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
      if (synr <= 0) 
        SYSTEM_break(BRK_3);
      { register GXFILE_tgdxsymbrecord *_W2=ValueCast(
        GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,synr));
        _P3getmem(_W2->sdomstrings,(_W2->sdim + 1) * sizeof(
          SYSTEM_longint));
        { register SYSTEM_int32 _stop = _W2->sdim;
          if ((d = 1) <=  _stop) do {
            (*_W2->sdomstrings)[d] = VirtMethodCall(ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
              GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
          } while (d++ !=  _stop);

        }

      }
    
CNT_3:;
    }
BRK_3:;
    {
      SYSTEM_shortstring _t1;

      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
        GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
        GXFILE_mark_doms),GXFILE_err_open_domsmarker3)) 
        goto _Lfileerrornr_31;
    }
  } 
  self->GXFILE_tgxfileobj_DOT_lasterror = GXFILE_err_noerror;
  GXFILE_tgxfileobj_DOT_gdxresetspecialvalues(self);
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
  self->GXFILE_tgxfileobj_DOT_fstatus = GXFILE_stat_read;
  GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
    GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
    SYSTEM_false);
  result = GXFILE_mytrue;
  return result;
  _Lfilenogood_31:;
  self->GXFILE_tgxfileobj_DOT_lasterror = *errnr;
  _Lfileerrornr_31:;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_ffile));
  self->GXFILE_tgxfileobj_DOT_ffile = NULL;
  result = GXFILE_myfalse;
  return result;
}  /* gdxopenreadxx */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxclose(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  SYSTEM_integer n;
  SYSTEM_integer commcnt;
  SYSTEM_integer cnt;
  SYSTEM_int64 uelpos;
  SYSTEM_int64 symbpos;
  SYSTEM_int64 settextpos;
  SYSTEM_int64 acronympos;
  SYSTEM_int64 domstrpos;
  GXFILE_pgdxsymbrecord psy;
  SYSTEM_integer d;
  SYSTEM_shortstring fnconv;

  _P3strclr(fnconv);
  if (_P3SET_in_1(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fw_raw_data,
    _P3SET_in_1(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fw_map_data,
    _P3SET_equal(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fw_str_data)))) 
    GXFILE_tgxfileobj_DOT_gdxdatawritedone(self);
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fw_init) {
    _P3strcpy(fnconv,255,self->GXFILE_tgxfileobj_DOT_ffile->
      GMSSTRM_txfilestream_DOT_ffilename);
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_compressout);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),self->
      GXFILE_tgxfileobj_DOT_nextwriteposition));
    symbpos = self->GXFILE_tgxfileobj_DOT_nextwriteposition;
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_symb);
    GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),self->
      GXFILE_tgxfileobj_DOT_namelist->STRHASH_txstrhashlist_DOT_fcount);
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_namelist->
        STRHASH_txstrhashlist_DOT_fcount;
      if ((n = 1) <=  _stop) do {
        {
          SYSTEM_shortstring _t1;

          GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),
            STRHASH_txstrhashlist_DOT_getstring(_t1,255,self->
            GXFILE_tgxfileobj_DOT_namelist,n));
        }
        psy = ValueCast(GXFILE_pgdxsymbrecord,
          STRHASH_txstrhashlist_DOT_getobject(self->
          GXFILE_tgxfileobj_DOT_namelist,n));
        { register GXFILE_tgdxsymbrecord *_W2=psy;
          GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),_W2->sposition);
          GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),_W2->sdim);
          GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),SYSTEM_ord(_W2->
            sdatatype));
          GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),_W2->suserinfo);
          GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),_W2->sdatacount);
          GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),_W2->serrors);
          GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),SYSTEM_ord(_W2->
            ssettext));
          GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),_W2->sexpltxt);
          GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),SYSTEM_ord(_W2->
            siscompressed));
          GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),SYSTEM_ord(_W2->
            sdomsymbols != NULL));
          if (_W2->sdomsymbols != NULL) 
            { register SYSTEM_int32 _stop = _W2->sdim;
              if ((d = 1) <=  _stop) do {
                GMSSTRM_txstream_DOT_writeinteger(ValueCast(
                  GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),(*
                  _W2->sdomsymbols)[d]);
              } while (d++ !=  _stop);

            }
          if (_W2->scommentslist == NULL) { 
            commcnt = 0;
          } else 
            commcnt = _W2->scommentslist->GMSOBJ_txlist_DOT_fcount;
          GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),commcnt);
          { register SYSTEM_int32 _stop = commcnt - 1;
            if ((cnt = 0) <=  _stop) do {
              {
                SYSTEM_shortstring _t1;

                GMSSTRM_txstream_DOT_writestring(ValueCast(
                  GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),
                  GMSOBJ_txstrings_DOT_get(_t1,255,_W2->
                  scommentslist,cnt));
              }
            } while (cnt++ !=  _stop);

          }

        }
      
      } while (n++ !=  _stop);

    }
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_symb);
    settextpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_compressout);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_sett);
    GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),self->
      GXFILE_tgxfileobj_DOT_settextlist->
      GMSOBJ_txcustomstringlist_DOT_fcount);
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_settextlist->
        GMSOBJ_txcustomstringlist_DOT_fcount - 1;
      if ((n = 0) <=  _stop) do {
        {
          SYSTEM_shortstring _t1;

          GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),
            GMSOBJ_txcustomstringlist_DOT_getname(_t1,255,ValueCast(
            GMSOBJ_txcustomstringlist,self->
            GXFILE_tgxfileobj_DOT_settextlist),n));
        }
      } while (n++ !=  _stop);

    }
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_sett);
    uelpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_compressout);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_uel);
    STRHASH_txstrhashlist_DOT_savetostream(ValueCast(
      STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile));
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_uel);
    acronympos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_compressout);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_acro);
    GXFILE_tacronymlist_DOT_savetostream(self->
      GXFILE_tgxfileobj_DOT_acronymlist,ValueCast(GMSSTRM_txstream,
      self->GXFILE_tgxfileobj_DOT_ffile));
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_acro);
    domstrpos = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      self->GXFILE_tgxfileobj_DOT_compressout);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_doms);
    STRHASH_txstrhashlist_DOT_savetostream(self->
      GXFILE_tgxfileobj_DOT_domainstrlist,ValueCast(GMSSTRM_txstream,
      self->GXFILE_tgxfileobj_DOT_ffile));
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_doms);
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_namelist->
        STRHASH_txstrhashlist_DOT_fcount;
      if ((n = 1) <=  _stop) do {
        psy = ValueCast(GXFILE_pgdxsymbrecord,
          STRHASH_txstrhashlist_DOT_getobject(self->
          GXFILE_tgxfileobj_DOT_namelist,n));
        { register GXFILE_tgdxsymbrecord *_W2=psy;
          if (_W2->sdomstrings != NULL) {
            GMSSTRM_txstream_DOT_writeinteger(ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),n);
            { register SYSTEM_int32 _stop = _W2->sdim;
              if ((d = 1) <=  _stop) do {
                GMSSTRM_txstream_DOT_writeinteger(ValueCast(
                  GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),(*
                  _W2->sdomstrings)[d]);
              } while (d++ !=  _stop);

            }
          } 

        }
      
      } while (n++ !=  _stop);

    }
    GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), -1);
    GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_doms);
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),self->
      GXFILE_tgxfileobj_DOT_majorindexposition));
    GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
      GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
      SYSTEM_false);
    GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_boi);
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),symbpos);
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),uelpos);
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),settextpos);
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),acronympos);
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),self->
      GXFILE_tgxfileobj_DOT_nextwriteposition);
    GMSSTRM_txstream_DOT_writeint64(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),domstrpos);
  } 
  result = 0;
  if (self->GXFILE_tgxfileobj_DOT_ffile != NULL) {
    result = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
      GMSSTRM_txfilestream,self->GXFILE_tgxfileobj_DOT_ffile));
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
      GXFILE_tgxfileobj_DOT_ffile));
    self->GXFILE_tgxfileobj_DOT_ffile = NULL;
  } 
  if (self->GXFILE_tgxfileobj_DOT_namelist != NULL) {
    for (n = self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount;n >= (SYSTEM_int32)1;--n) {
      psy = ValueCast(GXFILE_pgdxsymbrecord,
        STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,n));
      { register GXFILE_tgdxsymbrecord *_W2=psy;
        if (_W2->sdomsymbols != NULL) {
          _P3freemem(_W2->sdomsymbols);
          _W2->sdomsymbols = NULL;
        } 
        if (_W2->scommentslist != NULL) 
          SYSUTILS_P3_freeandnil(&_W2->scommentslist);
        if (_W2->ssetbitmap != NULL) 
          SYSUTILS_P3_freeandnil(&_W2->ssetbitmap);
        if (_W2->sdomstrings != NULL) 
          _P3freemem(_W2->sdomstrings);

      }
      _P3freemem(psy);
    
    }
    SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_namelist);
  } 
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_errorlist);
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_settextlist);
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_ueltable);
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_sortlist);
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_filterlist);
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_acronymlist);
  SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_domainstrlist);
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_not_open;
  self->GXFILE_tgxfileobj_DOT_fstatus = GXFILE_stat_notopen;
  if (self->GXFILE_tgxfileobj_DOT_autoconvert != 0 && _P3strcmpN(
    fnconv,_P3str1("\000"))) {
    if (self->GXFILE_tgxfileobj_DOT_compressout) { 
      result = GXFILE_convertgdxfile(fnconv,_P3str1("\001C"));
    } else 
      result = GXFILE_convertgdxfile(fnconv,_P3str1("\001U"));
    if (result > 0) 
      result = result + 100;
  } 
  return result;
}  /* gdxclose */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxerrorcount(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;

  result = self->GXFILE_tgxfileobj_DOT_errcnttotal;
  return result;
}  /* gdxerrorcount */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetlasterror(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_ffile == NULL) {
    result = self->GXFILE_tgxfileobj_DOT_lasterror;
    self->GXFILE_tgxfileobj_DOT_lasterror = GXFILE_err_noerror;
  } else {
    result = GMSSTRM_txfilestream_DOT_getlastioresult(ValueCast(
      GMSSTRM_txfilestream,self->GXFILE_tgxfileobj_DOT_ffile));
    if (result == 0) {
      result = self->GXFILE_tgxfileobj_DOT_lasterror;
      self->GXFILE_tgxfileobj_DOT_lasterror = GXFILE_err_noerror;
    } 
  } 
  return result;
}  /* gdxgetlasterror */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsysteminfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *sycnt,
  SYSTEM_integer *uelcnt)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_ueltable == NULL) { 
    *uelcnt = 0;
  } else 
    *uelcnt = self->GXFILE_tgxfileobj_DOT_ueltable->
      STRHASH_txstrhashlist_DOT_fcount;
  if (self->GXFILE_tgxfileobj_DOT_namelist == NULL) { 
    *sycnt = 0;
  } else 
    *sycnt = self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount;
  result = GXFILE_mytrue;
  return result;
}  /* gdxsysteminfo */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymboldim(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr)
{
  SYSTEM_integer result;

  if (synr == 0) { 
    result = 1;
  } else 
    if (!(self->GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && 
      synr <= self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount)) { 
      result =  -1;
    } else 
      result = (ValueCast(GXFILE_pgdxsymbrecord,
        STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,synr)))->sdim;
  return result;
}  /* gdxsymboldim */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_ansichar *syid,
  SYSTEM_integer *dimen,
  SYSTEM_integer *typ)
{
  SYSTEM_integer result;

  if (synr == 0) {
    _P3strcpy(syid,255,_P3str1("\001*"));
    *dimen = 1;
    *typ = 0;
    result = GXFILE_mytrue;
  } else 
    if (!(self->GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && 
      synr <= self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount)) {
      result = GXFILE_myfalse;
      _P3strclr(syid);
      *dimen =  -1;
      *typ = 0;
    } else {
      result = GXFILE_mytrue;
      STRHASH_txstrhashlist_DOT_getstring(syid,255,self->
        GXFILE_tgxfileobj_DOT_namelist,synr);
      { register GXFILE_tgdxsymbrecord *_W2=ValueCast(
        GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,synr));
        *dimen = _W2->sdim;
        *typ = SYSTEM_ord(_W2->sdatatype);

      }
    } 
  return result;
}  /* gdxsymbolinfo */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolinfox(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *reccnt,
  SYSTEM_integer *userinfo,
  SYSTEM_ansichar *expltxt)
{
  SYSTEM_integer result;

  if (synr == 0) {
    *reccnt = self->GXFILE_tgxfileobj_DOT_uelcntorig;
    *userinfo = 0;
    _P3strcpy(expltxt,255,_P3str1("\010Universe"));
    result = GXFILE_mytrue;
  } else 
    if (!(self->GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && 
      synr <= self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount)) {
      result = GXFILE_myfalse;
      *reccnt = 0;
      *userinfo = 0;
      _P3strclr(expltxt);
    } else {
      result = GXFILE_mytrue;
      { register GXFILE_tgdxsymbrecord *_W2=ValueCast(
        GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,synr));
        if (_W2->sdim == 0) { 
          *reccnt = 1;
        } else 
          *reccnt = _W2->sdatacount;
        *userinfo = _W2->suserinfo;
        _P3strcpy(expltxt,255,_W2->sexpltxt);

      }
    } 
  return result;
}  /* gdxsymbolinfox */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfindsymbol(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  SYSTEM_integer *synr)
{
  SYSTEM_integer result;

  if (_P3stccmpE(syid,_P3char('*'))) {
    *synr = 0;
    result = GXFILE_mytrue;
  } else {
    result = GXFILE_myfalse;
    if (self->GXFILE_tgxfileobj_DOT_namelist != NULL) {
      *synr = STRHASH_txstrhashlist_DOT_indexof(self->
        GXFILE_tgxfileobj_DOT_namelist,syid);
      if (*synr >= 1) 
        result = GXFILE_mytrue;
    } 
  } 
  return result;
}  /* gdxfindsymbol */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetuel(
  GXFILE_tgxfileobj self,
  SYSTEM_integer uelnr,
  SYSTEM_ansichar *uel)
{
  SYSTEM_integer result;
  SYSTEM_integer en;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_ueltable == NULL) {
    _P3strclr(uel);
    return result;
  } 
  en = GXFILE_tintegermapping_DOT_getmapping(self->
    GXFILE_tgxfileobj_DOT_ueltable->GXFILE_tueltable_DOT_usruel2ent,
    uelnr);
  if (en >= 1) {
    result = GXFILE_mytrue;
    STRHASH_txstrhashlist_DOT_getstring(uel,255,ValueCast(
      STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),en);
  } else 
    {
      SYSTEM_shortstring _t1;

      _P3strcat(uel,255,GXFILE_baduel_prefix,STRUTILX_inttostr(_t1,255,
        uelnr));
    }
  return result;
}  /* gdxgetuel */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_preparesymbolwrite(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *caller,
  const SYSTEM_ansichar *aname,
  const SYSTEM_ansichar *atext,
  SYSTEM_integer adim,
  SYSTEM_integer atype,
  SYSTEM_integer auserinfo)
{
  SYSTEM_boolean result;
  static GXFILE_tgxmodeset allowedmodes = {4,0,0};
  SYSTEM_integer d;

  result = SYSTEM_false;
  self->GXFILE_tgxfileobj_DOT_cursyptr = NULL;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_errorlist));
  self->GXFILE_tgxfileobj_DOT_errorlist = NULL;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_sortlist));
  self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,caller,allowedmodes)) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_some) 
    {
      _P3STR_255 _t1;
      _P3STR_255 _t2;
      SYSTEM_shortstring _t3;
      _P3STR_255 _t4;

      GXFILE_tgxfileobj_DOT_writetrace(self,_P3strcat(_t4,255,
        _P3strcat(_t2,255,_P3strcat(_t1,255,_P3str1("\011Symbol = "),
        aname),_P3str1("\010, Dim = ")),STRUTILX_inttostr(_t3,255,
        adim)));
    }
  if (!GXFILE_tgxfileobj_DOT_isgoodnewsymbol(self,aname)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,adim >= 0 && adim <= 20,
    GXFILE_err_baddimension)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,atype >= 0 && 
    atype <= 3,GXFILE_err_baddatatype)) 
    return result;
  _P3new(self->GXFILE_tgxfileobj_DOT_cursyptr);
  { register GXFILE_tgdxsymbrecord *_W2=self->
    GXFILE_tgxfileobj_DOT_cursyptr;
    _W2->sposition = 0;
    _W2->sdim = adim;
    _W2->sdatatype = ValueCast(GMSSPECS_tgdxdatatype,atype);
    _W2->sdatacount = 0;
    _W2->serrors = 0;
    _W2->suserinfo = auserinfo;
    _W2->ssettext = SYSTEM_false;
    GXFILE_makegoodexpltext(_W2->sexpltxt,255,atext);
    _W2->siscompressed = self->GXFILE_tgxfileobj_DOT_compressout && 
      adim > 0;
    _W2->scommentslist = NULL;
    _W2->sdomsymbols = NULL;
    if (_P3SET_in_1(atype,0,_P3SET_equal(atype,4)) && adim == 1) { 
      _W2->ssetbitmap = ValueCast(GMSOBJ_tbooleanbitarray,
        GMSOBJ_tbooleanbitarray_DOT_create(ValueCast(
        GMSOBJ_tbooleanbitarray,_P3alloc_object(&
        GMSOBJ_tbooleanbitarray_CD))));
    } else 
      _W2->ssetbitmap = NULL;
    _W2->sdomstrings = NULL;

  }
  self->GXFILE_tgxfileobj_DOT_cursyptr->ssynr = 
    STRHASH_txstrhashlist_DOT_addobject(self->
    GXFILE_tgxfileobj_DOT_namelist,aname,ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_cursyptr));
  self->GXFILE_tgxfileobj_DOT_fcurrentdim = adim;
  if (SYSTEM_false) { 
    self->GXFILE_tgxfileobj_DOT_deltaforwrite = 244;
  } else 
    self->GXFILE_tgxfileobj_DOT_deltaforwrite = 255 - self->
      GXFILE_tgxfileobj_DOT_fcurrentdim - 1;
  self->GXFILE_tgxfileobj_DOT_datasize = GXDEFS_datatypsize[ValueCast(
    GMSSPECS_tgdxdatatype,atype)];
  if (self->GXFILE_tgxfileobj_DOT_datasize > 0) 
    self->GXFILE_tgxfileobj_DOT_lastdatafield = ValueCast(
      GMSSPECS_tvarvaltype,self->GXFILE_tgxfileobj_DOT_datasize - 1);
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = 
        GXFILE_index_initial;
      self->GXFILE_tgxfileobj_DOT_minelem[d - 1] = SYSTEM_maxint;
      self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] = 0;
      self->GXFILE_tgxfileobj_DOT_wrbitmaps[d - 1] = NULL;
    
    } while (d++ !=  _stop);

  }
  GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
    GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
    self->GXFILE_tgxfileobj_DOT_cursyptr->siscompressed);
  result = SYSTEM_true;
  return result;
}  /* preparesymbolwrite */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawriterawstart(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  const SYSTEM_ansichar *expltxt,
  SYSTEM_integer dimen,
  SYSTEM_integer typ,
  SYSTEM_integer userinfo)
{
  SYSTEM_integer result;
  SYSTEM_integer d;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_preparesymbolwrite(self,_P3str1("\021DataWriteRawStart"),
    syid,expltxt,dimen,typ,userinfo)) 
    return result;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      self->GXFILE_tgxfileobj_DOT_minelem[d - 1] = 0;
      self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] = SYSTEM_maxint;
    
    } while (d++ !=  _stop);

  }
  GXFILE_tgxfileobj_DOT_initdowrite(self, -1);
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_dom_raw;
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatawriterawstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritemapstart(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  const SYSTEM_ansichar *expltxt,
  SYSTEM_integer dimen,
  SYSTEM_integer typ,
  SYSTEM_integer userinfo)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_preparesymbolwrite(self,_P3str1("\021DataWriteMapStart"),
    syid,expltxt,dimen,typ,userinfo)) 
    return result;
  self->GXFILE_tgxfileobj_DOT_sortlist = ValueCast(
    DATASTORAGE_tlinkeddata,DATASTORAGE_tlinkeddata_DOT_create(ValueCast(
    DATASTORAGE_tlinkeddata,_P3alloc_object(&
    DATASTORAGE_tlinkeddata_CD)),self->
    GXFILE_tgxfileobj_DOT_fcurrentdim,self->
    GXFILE_tgxfileobj_DOT_datasize * sizeof(SYSTEM_double)));
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_dom_map;
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatawritemapstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritestrstart(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  const SYSTEM_ansichar *expltxt,
  SYSTEM_integer dimen,
  SYSTEM_integer typ,
  SYSTEM_integer userinfo)
{
  SYSTEM_integer result;
  SYSTEM_integer d;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_preparesymbolwrite(self,_P3str1("\021DataWriteStrStart"),
    syid,expltxt,dimen,typ,userinfo)) 
    return result;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      _P3strcpy(self->GXFILE_tgxfileobj_DOT_laststrelem[d - 1],255,_P3str1("\001\377"));
    } while (d++ !=  _stop);

  }
  self->GXFILE_tgxfileobj_DOT_sortlist = ValueCast(
    DATASTORAGE_tlinkeddata,DATASTORAGE_tlinkeddata_DOT_create(ValueCast(
    DATASTORAGE_tlinkeddata,_P3alloc_object(&
    DATASTORAGE_tlinkeddata_CD)),self->
    GXFILE_tgxfileobj_DOT_fcurrentdim,self->
    GXFILE_tgxfileobj_DOT_datasize * sizeof(SYSTEM_double)));
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_dom_str;
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatawritestrstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawriteraw(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *keyint,
  const SYSTEM_double *values)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_2 = {64,0,0};

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fw_dom_raw) 
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_raw_data;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_some || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_2)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\014DataWriteRaw"),
      allowedmodes_2)) 
      return result;
  if (GXFILE_tgxfileobj_DOT_dowrite(self,keyint,values)) 
    result = GXFILE_mytrue;
  return result;
}  /* gdxdatawriteraw */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritemap(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *keyint,
  const SYSTEM_double *values)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_3 = {128,0,0};
  SYSTEM_integer d;
  SYSTEM_integer kd;
  GMSSPECS_tindex keys;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fw_dom_map) 
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_map_data;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_3)) {
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\014DataWriteMap"),
      allowedmodes_3)) 
      return result;
    _Iplus_bgn();
    _P3write_s0(_P3str1("\012   Index ="));
    _Iplus_end();
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        _Iplus_bgn();
        _P3write_c0(_P3char(' '));
        _P3write_i0(keyint[d - 1]);
        _Iplus_end();
        if (d < self->GXFILE_tgxfileobj_DOT_fcurrentdim) {
          _Iplus_bgn();
          _P3write_c0(_P3char(','));
          _Iplus_end();
        } 
      
      } while (d++ !=  _stop);

    }
    _Iplus_bgn();
    _P3writeln();
    _Iplus_end();
  } 
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      kd = keyint[d - 1];
      kd = GXFILE_tintegermapping_DOT_getmapping(self->
        GXFILE_tgxfileobj_DOT_ueltable->
        GXFILE_tueltable_DOT_usruel2ent,kd);
      if (kd < 0) {
        GXFILE_tgxfileobj_DOT_reporterror(self,
          GXFILE_err_badelementindex);
        return result;
      } 
      keys[d - 1] = kd;
      if (kd < self->GXFILE_tgxfileobj_DOT_minelem[d - 1]) 
        self->GXFILE_tgxfileobj_DOT_minelem[d - 1] = kd;
      if (kd > self->GXFILE_tgxfileobj_DOT_maxelem[d - 1]) 
        self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] = kd;
    
    } while (d++ !=  _stop);

  }
  _P3_TRY {
    DATASTORAGE_tlinkeddata_DOT_additem(self->
      GXFILE_tgxfileobj_DOT_sortlist,ValueCast(GMSGEN_pintegerarrayone,&
      keys[0]),values);
    result = GXFILE_mytrue;
  } _P3_EXCEPT {
{
      GXFILE_tgxfileobj_DOT_seterror(self,GXFILE_err_out_of_memory);
      SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
        GXFILE_tgxfileobj_DOT_sortlist));
      self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
      self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_init;
    } 
  } _P3_END_TRY_EXCEPT;
  return result;
}  /* gdxdatawritemap */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritestr(
  GXFILE_tgxfileobj self,
  const SYSTEM_shortstring *keystr,
  const SYSTEM_double *values)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  SYSTEM_shortstring sv;
  SYSTEM_integer kd;
  static GXFILE_tgxmodeset allowedmodes_4 = {0,1,0};

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fw_dom_str) 
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_str_data;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_4)) {
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\014DataWriteStr"),
      allowedmodes_4)) 
      return result;
    _Iplus_bgn();
    _P3write_s0(_P3str1("\011  Index ="));
    _P3writeln();
    _Iplus_end();
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        _Iplus_bgn();
        _P3write_c0(_P3char(' '));
        _P3write_s0(keystr[d - 1]);
        _Iplus_end();
        if (d < self->GXFILE_tgxfileobj_DOT_fcurrentdim) {
          _Iplus_bgn();
          _P3write_c0(_P3char(','));
          _Iplus_end();
        } 
      
      } while (d++ !=  _stop);

    }
    _Iplus_bgn();
    _P3writeln();
    _Iplus_end();
  } 
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      SYSUTILS_P3_trimright(sv,255,keystr[d - 1]);
      if (_P3strcmpN(sv,self->GXFILE_tgxfileobj_DOT_laststrelem[d - 1])) {
        kd = STRHASH_txstrhashlist_DOT_indexof(ValueCast(
          STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
          sv);
        if (kd <= 0) {
          if (GXFILE_tgxfileobj_DOT_errorcondition(self,
            GXDEFS_gooduelstring(sv),GXFILE_err_baduelstr)) 
            return result;
          kd = STRHASH_txstrhashlist_DOT_addobject(ValueCast(
            STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
            sv,ValueCast(SYSTEM_tobject,GMSOBJ_copyint2ptr( -1)));
        } 
        self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = kd;
        _P3strcpy(self->GXFILE_tgxfileobj_DOT_laststrelem[d - 1],255,
          sv);
        if (kd < self->GXFILE_tgxfileobj_DOT_minelem[d - 1]) 
          self->GXFILE_tgxfileobj_DOT_minelem[d - 1] = kd;
        if (kd > self->GXFILE_tgxfileobj_DOT_maxelem[d - 1]) 
          self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] = kd;
      } 
    
    } while (d++ !=  _stop);

  }
  _P3_TRY {
    DATASTORAGE_tlinkeddata_DOT_additem(self->
      GXFILE_tgxfileobj_DOT_sortlist,ValueCast(GMSGEN_pintegerarrayone,&
      self->GXFILE_tgxfileobj_DOT_lastelem[0]),values);
    result = GXFILE_mytrue;
  } _P3_EXCEPT {
{
      GXFILE_tgxfileobj_DOT_seterror(self,GXFILE_err_out_of_memory);
      SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
        GXFILE_tgxfileobj_DOT_sortlist));
      self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
      self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_init;
    } 
  } _P3_END_TRY_EXCEPT;
  return result;
}  /* gdxdatawritestr */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritedone(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_5 = {248,1,0};
  GMSSPECS_tindex aelements;
  GXDEFS_tgdxvalues avals;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\015DataWriteDone"),
    allowedmodes_5)) 
    return result;
  if (!_P3SET_in_1(self->GXFILE_tgxfileobj_DOT_fmode,
    GXFILE_fw_raw_data,_P3SET_equal(self->GXFILE_tgxfileobj_DOT_fmode,
    GXFILE_fw_dom_raw))) {
    GXFILE_tgxfileobj_DOT_initdowrite(self,self->
      GXFILE_tgxfileobj_DOT_sortlist->
      DATASTORAGE_tlinkeddata_DOT_fcount);
    if (DATASTORAGE_tlinkeddata_DOT_startread(self->
      GXFILE_tgxfileobj_DOT_sortlist,&self->
      GXFILE_tgxfileobj_DOT_readptr,ValueCast(GMSGEN_pintegerarrayone,NULL))) 
      while (DATASTORAGE_tlinkeddata_DOT_getnextrecord(self->
        GXFILE_tgxfileobj_DOT_sortlist,&self->
        GXFILE_tgxfileobj_DOT_readptr,ValueCast(
        GMSGEN_pintegerarrayone,&aelements[0]),avals)) {

        GXFILE_tgxfileobj_DOT_dowrite(self,VariableCast(
          GMSSPECS_tindex,aelements,GMSSPECS_tindex),avals);
}
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
      GXFILE_tgxfileobj_DOT_sortlist));
    self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
  } 
  GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),255);
  self->GXFILE_tgxfileobj_DOT_nextwriteposition = VirtMethodCall(ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
    GMSSTRM_txstream_DOT_getposition_T, 1, (ValueCast(GMSSTRM_txstream,
    self->GXFILE_tgxfileobj_DOT_ffile)));
  self->GXFILE_tgxfileobj_DOT_cursyptr->sdatacount = self->
    GXFILE_tgxfileobj_DOT_datacount;
  self->GXFILE_tgxfileobj_DOT_cursyptr->serrors = self->
    GXFILE_tgxfileobj_DOT_errcnt;
  self->GXFILE_tgxfileobj_DOT_errcnt = 0;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_init;
  GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
    GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
    SYSTEM_false);
  self->GXFILE_tgxfileobj_DOT_cursyptr = NULL;
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatawritedone */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_preparesymbolread(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *caller,
  SYSTEM_integer synr,
  const SYSTEM_integer *adomainnrs,
  GXFILE_tgxfilemode newmode)
{
  SYSTEM_integer result;
  GXDEFS_tgdxvalues avals;
  static GXFILE_tgxmodeset allowedmodes_6 = {2,0,0};
  SYSTEM_integer n;
  SYSTEM_integer d;
  GMSSPECS_tgdxdimension d2;
  SYSTEM_integer v;
  SYSTEM_integer en;
  SYSTEM_boolean addnew;
  SYSTEM_boolean adderror;
  SYSTEM_integer fidim;
  SYSTEM_integer afdim;
  GXDEFS_tgdxuelindex aelements;
  SYSTEM_integer nrrecs;
  SYSTEM_boolean allocok;
  GXFILE_tintegermapping expndlist;

  if (_P3SET_in_1(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fr_str_data,
    _P3SET_in_1(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fr_map_data,
    _P3SET_in_1(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fr_mapr_data,
    _P3SET_equal(self->GXFILE_tgxfileobj_DOT_fmode,GXFILE_fr_raw_data))))) 
    GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  self->GXFILE_tgxfileobj_DOT_nrmappedadded = 0;
  expndlist = NULL;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_errorlist));
  self->GXFILE_tgxfileobj_DOT_errorlist = NULL;
  result =  -1;
  self->GXFILE_tgxfileobj_DOT_cursyptr = NULL;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_sortlist));
  self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,caller,allowedmodes_6)) {
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
    return result;
  } 
  self->GXFILE_tgxfileobj_DOT_readuniverse = synr == 0;
  if (!self->GXFILE_tgxfileobj_DOT_readuniverse) {
    if (GXFILE_tgxfileobj_DOT_errorcondition(self,synr >= 1 && 
      synr <= self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount,GXFILE_err_badsymbolindex)) 
      return result;
    self->GXFILE_tgxfileobj_DOT_cursyptr = ValueCast(
      GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
      GXFILE_tgxfileobj_DOT_namelist,synr));
    if (self->GXFILE_tgxfileobj_DOT_cursyptr->sdatatype == 
      GMSSPECS_dt_alias) {
      do {
        synr = self->GXFILE_tgxfileobj_DOT_cursyptr->suserinfo;
        if (synr == 0) {
          self->GXFILE_tgxfileobj_DOT_readuniverse = SYSTEM_true;
          SYSTEM_break(BRK_4);
        } 
        self->GXFILE_tgxfileobj_DOT_cursyptr = ValueCast(
          GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(
          self->GXFILE_tgxfileobj_DOT_namelist,synr));
      CNT_4:;
      } while (!(self->GXFILE_tgxfileobj_DOT_cursyptr->sdatatype != 
        GMSSPECS_dt_alias));
BRK_4:;
      if (!self->GXFILE_tgxfileobj_DOT_readuniverse) 
        SYSTEM_assert(self->GXFILE_tgxfileobj_DOT_cursyptr->sdatatype == 
          GMSSPECS_dt_set,_P3str1("\021Bad aliased set-1"));
    } 
  } 
  if (self->GXFILE_tgxfileobj_DOT_readuniverse) {
    self->GXFILE_tgxfileobj_DOT_fcurrentdim = 1;
    self->GXFILE_tgxfileobj_DOT_datasize = GXDEFS_datatypsize[
      GMSSPECS_dt_set];
    self->GXFILE_tgxfileobj_DOT_lastdatafield = ValueCast(
      GMSSPECS_tvarvaltype,self->GXFILE_tgxfileobj_DOT_datasize - 1);
    nrrecs = self->GXFILE_tgxfileobj_DOT_uelcntorig;
    self->GXFILE_tgxfileobj_DOT_universenr = 0;
    self->GXFILE_tgxfileobj_DOT_cursyptr = NULL;
  } else 
    { register GXFILE_tgdxsymbrecord *_W2=self->
      GXFILE_tgxfileobj_DOT_cursyptr;
      self->GXFILE_tgxfileobj_DOT_fcurrentdim = _W2->sdim;
      GMSSTRM_tbufferedfilestream_DOT_setcompression(ValueCast(
        GMSSTRM_tbufferedfilestream,self->GXFILE_tgxfileobj_DOT_ffile),
        _W2->siscompressed);
      VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
        GXFILE_tgxfileobj_DOT_ffile), 
        GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),_W2->
        sposition));
      self->GXFILE_tgxfileobj_DOT_datasize = GXDEFS_datatypsize[_W2->
        sdatatype];
      if (self->GXFILE_tgxfileobj_DOT_datasize > 0) 
        self->GXFILE_tgxfileobj_DOT_lastdatafield = ValueCast(
          GMSSPECS_tvarvaltype,self->GXFILE_tgxfileobj_DOT_datasize - 1);
      nrrecs = _W2->sdatacount;

    }
  if (self->GXFILE_tgxfileobj_DOT_versionread <= 6) { 
    self->GXFILE_tgxfileobj_DOT_deltaforread = GMSSPECS_maxdimv148;
  } else 
    self->GXFILE_tgxfileobj_DOT_deltaforread = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      { register GXFILE_tdomain *_W2= &self->
        GXFILE_tgxfileobj_DOT_domainlist[d - 1];
        _W2->dfilter = NULL;
        switch (adomainnrs[d - 1]) {
          case GXDEFS_domc_unmapped: 
            _W2->daction = GXFILE_dm_unmapped;
            break;
          case GXDEFS_domc_expand: 
            _W2->daction = GXFILE_dm_expand;
            if (expndlist == NULL) 
              expndlist = ValueCast(GXFILE_tintegermapping,
                GXFILE_tintegermapping_DOT_create(ValueCast(
                GXFILE_tintegermapping,_P3alloc_object(&
                GXFILE_tintegermapping_CD))));
            break;
          case GXDEFS_domc_strict: 
            _W2->daction = GXFILE_dm_strict;
            break;
          default:
            _W2->dfilter = GXFILE_tfilterlist_DOT_findfilter(self->
              GXFILE_tgxfileobj_DOT_filterlist,adomainnrs[d - 1]);
            if (_W2->dfilter != NULL) { 
              _W2->daction = GXFILE_dm_filter;
            } else {
              GXFILE_tgxfileobj_DOT_reporterror(self,
                GXFILE_err_unknownfilter);
              return result;
            } 
        }

      }
    } while (d++ !=  _stop);

  }
  if (!self->GXFILE_tgxfileobj_DOT_readuniverse) {
    {
      SYSTEM_shortstring _t1;

      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3strcmpE(
        GMSSTRM_txstream_DOT_readstring(_t1,255,ValueCast(
        GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)),
        GXFILE_mark_data),GXFILE_err_baddatamarker_data)) 
        return result;
    }
    if (GXFILE_tgxfileobj_DOT_errorcondition(self,ValueCast(
      SYSTEM_int32,GMSSTRM_txstream_DOT_readbyte(ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile))) == self->
      GXFILE_tgxfileobj_DOT_fcurrentdim,GXFILE_err_baddatamarker_dim)) 
      return result;
    VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
      GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
  } 
  if (self->GXFILE_tgxfileobj_DOT_fcurrentdim == 0 && nrrecs == 0) {
    self->GXFILE_tgxfileobj_DOT_cursyptr->sscalarfrst = SYSTEM_true;
    self->GXFILE_tgxfileobj_DOT_fmode = newmode;
    result = 1;
    return result;
  } 
  if (!self->GXFILE_tgxfileobj_DOT_readuniverse) {
    self->GXFILE_tgxfileobj_DOT_cursyptr->sscalarfrst = SYSTEM_false;
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = 
          GXFILE_index_initial;
        self->GXFILE_tgxfileobj_DOT_prevelem[d - 1] =  -1;
        self->GXFILE_tgxfileobj_DOT_minelem[d - 1] = VirtMethodCall(ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
          GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] = VirtMethodCall(ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
          GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
        self->GXFILE_tgxfileobj_DOT_elemtype[d - 1] = 
          GXFILE_getintegersize(self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] - 
          self->GXFILE_tgxfileobj_DOT_minelem[d - 1] + 1);
      
      } while (d++ !=  _stop);

    }
  } 
  allocok = SYSTEM_true;
  if (_P3SET_in_1(newmode,GXFILE_fr_raw_data,_P3SET_in_1(newmode,
    GXFILE_fr_str_data,_P3SET_equal(newmode,GXFILE_fr_slice)))) { 
    result = nrrecs;
  } else {
    SYSTEM_assert(newmode == GXFILE_fr_map_data,_P3str1("\032Expect to read mapped data"));
    if (GXFILE_tgxfileobj_DOT_resultwillbesorted(self,adomainnrs)) {
      result = nrrecs;
      newmode = GXFILE_fr_mapr_data;
    } else 
      _P3_TRY {
        self->GXFILE_tgxfileobj_DOT_sortlist = ValueCast(
          DATASTORAGE_tlinkeddata,DATASTORAGE_tlinkeddata_DOT_create(ValueCast(
          DATASTORAGE_tlinkeddata,_P3alloc_object(&
          DATASTORAGE_tlinkeddata_CD)),self->
          GXFILE_tgxfileobj_DOT_fcurrentdim,self->
          GXFILE_tgxfileobj_DOT_datasize * sizeof(SYSTEM_double)));
        fidim = self->GXFILE_tgxfileobj_DOT_fcurrentdim;
        addnew = SYSTEM_false;
        adderror = SYSTEM_false;
        while (GXFILE_tgxfileobj_DOT_doread(self,avals,&afdim)) {
          if (fidim < afdim) 
            afdim = fidim;
          fidim = self->GXFILE_tgxfileobj_DOT_fcurrentdim;
          { register SYSTEM_int32 _stop = self->
              GXFILE_tgxfileobj_DOT_fcurrentdim;
            if ((d = afdim) <=  _stop) do {
              { register GXFILE_tdomain *_W2= &self->
                GXFILE_tgxfileobj_DOT_domainlist[d - 1];
                switch (_W2->daction) {
                  case GXFILE_dm_unmapped: 
                    aelements[d - 1] = self->
                      GXFILE_tgxfileobj_DOT_lastelem[d - 1];
                    break;
                  case GXFILE_dm_filter: 
                    v = GXFILE_tueltable_DOT_getusermap(self->
                      GXFILE_tgxfileobj_DOT_ueltable,self->
                      GXFILE_tgxfileobj_DOT_lastelem[d - 1]);
                    if (GXFILE_tdfilter_DOT_infilter(_W2->dfilter,v)) { 
                      aelements[d - 1] = v;
                    } else {
                      adderror = SYSTEM_true;
                      fidim = d;
                      SYSTEM_break(BRK_5);
                    } 
                    break;
                  case GXFILE_dm_strict: 
                    v = GXFILE_tueltable_DOT_getusermap(self->
                      GXFILE_tgxfileobj_DOT_ueltable,self->
                      GXFILE_tgxfileobj_DOT_lastelem[d - 1]);
                    if (v >= 0) { 
                      aelements[d - 1] = v;
                    } else {
                      adderror = SYSTEM_true;
                      fidim = d;
                      SYSTEM_break(BRK_5);
                    } 
                    break;
                  case GXFILE_dm_expand: 
                    en = self->GXFILE_tgxfileobj_DOT_lastelem[d - 1];
                    v = GXFILE_tintegermapping_DOT_getmapping(
                      expndlist,en);
                    if (v >= 0) { 
                      aelements[d - 1] = v;
                    } else {
                      v = GXFILE_tueltable_DOT_getusermap(self->
                        GXFILE_tgxfileobj_DOT_ueltable,en);
                      if (v >= 0) {
                        aelements[d - 1] = v;
                        GXFILE_tintegermapping_DOT_setmapping(
                          expndlist,en,v);
                      } else {
                        aelements[d - 1] = -en;
                        addnew = SYSTEM_true;
                      } 
                    } 
                    break;
                  default: break;
                }

              }
CNT_5:;
            } while (d++ !=  _stop);
BRK_5:;

          }
          if (adderror) {
            GXFILE_tgxfileobj_DOT_addtoerrorlist(self,self->
              GXFILE_tgxfileobj_DOT_lastelem,avals);
            adderror = SYSTEM_false;
          } else {
            if (addnew) {
              { register SYSTEM_int32 _stop = self->
                  GXFILE_tgxfileobj_DOT_fcurrentdim;
                if ((d = 1) <=  _stop) do {
                  en = aelements[d - 1];
                  if (en < 0) {
                    v = GXFILE_tueltable_DOT_newusruel(self->
                      GXFILE_tgxfileobj_DOT_ueltable,-en);
                    aelements[d - 1] = v;
                    GXFILE_tintegermapping_DOT_setmapping(expndlist,-
                      en,v);
                    self->GXFILE_tgxfileobj_DOT_nrmappedadded = self->
                      GXFILE_tgxfileobj_DOT_nrmappedadded + 1;
                    { register SYSTEM_uint8 _stop = self->
                        GXFILE_tgxfileobj_DOT_fcurrentdim;
                      if ((d2 = d + 1) <=  _stop) do {
                        if (aelements[d2 - 1] == en) 
                          aelements[d2 - 1] = v;
                      } while (d2++ !=  _stop);

                    }
                  } 
                
                } while (d++ !=  _stop);

              }
              addnew = SYSTEM_false;
            } 
            DATASTORAGE_tlinkeddata_DOT_additem(self->
              GXFILE_tgxfileobj_DOT_sortlist,ValueCast(
              GMSGEN_pintegerarrayone,&aelements[0]),avals);
          } 
        
}
        SYSUTILS_P3_freeandnil(&expndlist);
        DATASTORAGE_tlinkeddata_DOT_startread(self->
          GXFILE_tgxfileobj_DOT_sortlist,&self->
          GXFILE_tgxfileobj_DOT_readptr,ValueCast(
          GMSGEN_pintegerarrayone,NULL));
        result = self->GXFILE_tgxfileobj_DOT_sortlist->
          DATASTORAGE_tlinkeddata_DOT_fcount;
      } _P3_EXCEPT {

          allocok = SYSTEM_false;
      } _P3_END_TRY_EXCEPT;
} 
  if (allocok) {
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] =  -1;
      } while (d++ !=  _stop);

    }
    self->GXFILE_tgxfileobj_DOT_fmode = newmode;
  } else {
    GXFILE_tgxfileobj_DOT_seterror(self,GXFILE_err_out_of_memory);
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
      GXFILE_tgxfileobj_DOT_sortlist));
    self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
    result =  -1;
  } 
  return result;
}  /* preparesymbolread */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadrawstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *nrrecs)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex xdomains;

  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_unmapped;
  }
  *nrrecs = GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\020DataReadRawStart"),
    synr,xdomains,GXFILE_fr_raw_data);
  if (*nrrecs >= 0) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxdatareadrawstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadmapstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *nrrecs)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex xdomains;

  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_strict;
  }
  *nrrecs = GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\020DataReadMapStart"),
    synr,xdomains,GXFILE_fr_map_data);
  if (*nrrecs >= 0) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxdatareadmapstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadstrstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *nrrecs)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex xdomains;

  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_unmapped;
  }
  *nrrecs = GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\020DataReadStrStart"),
    synr,xdomains,GXFILE_fr_str_data);
  if (*nrrecs >= 0) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxdatareadstrstart */

Function(SYSTEM_integer ) 
  GXFILE_tgxfileobj_DOT_gdxdatareadfilteredstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_integer *filteraction,
  SYSTEM_integer *nrrecs)
{
  SYSTEM_integer result;

  *nrrecs = GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\025DataReadStartFiltered"),
    synr,filteraction,GXFILE_fr_map_data);
  if (*nrrecs >= 0) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxdatareadfilteredstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadraw(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *keyint,
  SYSTEM_double *values,
  SYSTEM_integer *dimfrst)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_7 = {0,16,0};
  SYSTEM_integer d;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_7)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\013DataReadRaw"),
      allowedmodes_7)) 
      return result;
  if (!GXFILE_tgxfileobj_DOT_doread(self,values,dimfrst)) { 
    GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  } else {
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        keyint[d - 1] = self->GXFILE_tgxfileobj_DOT_lastelem[d - 1];
      } while (d++ !=  _stop);

    }
    result = GXFILE_mytrue;
  } 
  return result;
}  /* gdxdatareadraw */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadmap(
  GXFILE_tgxfileobj self,
  SYSTEM_integer recnr,
  SYSTEM_integer *keyint,
  SYSTEM_double *values,
  SYSTEM_integer *dimfrst)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_8 = {0,96,0};
  SYSTEM_integer d, d2;
  SYSTEM_integer en;
  SYSTEM_integer v;
  SYSTEM_integer fidim;
  SYSTEM_boolean addnew, adderror;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_8)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\013DataReadMap"),
      allowedmodes_8)) 
      return result;
  if (self->GXFILE_tgxfileobj_DOT_cursyptr != NULL && self->
    GXFILE_tgxfileobj_DOT_cursyptr->sscalarfrst) {
    self->GXFILE_tgxfileobj_DOT_cursyptr->sscalarfrst = SYSTEM_false;
    GXFILE_tgxfileobj_DOT_getdefaultrecord(self,values);
    *dimfrst = 0;
    result = GXFILE_mytrue;
    return result;
  } 
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fr_map_data) {
    *dimfrst = 0;
    if (self->GXFILE_tgxfileobj_DOT_readptr == NULL) 
      return result;
    DATASTORAGE_tlinkeddata_DOT_getnextrecord(self->
      GXFILE_tgxfileobj_DOT_sortlist,&self->
      GXFILE_tgxfileobj_DOT_readptr,ValueCast(GMSGEN_pintegerarrayone,&
      keyint[0]),values);
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        if (keyint[d - 1] != self->GXFILE_tgxfileobj_DOT_prevelem[
          d - 1]) {
          self->GXFILE_tgxfileobj_DOT_prevelem[d - 1] = keyint[d - 1];
          if (*dimfrst == 0) 
            *dimfrst = d;
        } 
      } while (d++ !=  _stop);

    }
    result = GXFILE_mytrue;
    return result;
  } 
  SYSTEM_assert(self->GXFILE_tgxfileobj_DOT_fmode == 
    GXFILE_fr_mapr_data,_P3str1("\025fr_MapR_data expected"));
  addnew = SYSTEM_false;
  adderror = SYSTEM_false;
  fidim = self->GXFILE_tgxfileobj_DOT_fcurrentdim;
  _Lagain_55:;
  if (!GXFILE_tgxfileobj_DOT_doread(self,values,dimfrst)) 
    return result;
  if (fidim < *dimfrst) 
    *dimfrst = fidim;
  fidim = self->GXFILE_tgxfileobj_DOT_fcurrentdim;
  if (*dimfrst > 0) 
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = *dimfrst) <=  _stop) do {
        { register GXFILE_tdomain *_W2= &self->
          GXFILE_tgxfileobj_DOT_domainlist[d - 1];
          switch (_W2->daction) {
            case GXFILE_dm_unmapped: 
              keyint[d - 1] = self->
                GXFILE_tgxfileobj_DOT_lastelem[d - 1];
              break;
            case GXFILE_dm_filter: 
              v = GXFILE_tueltable_DOT_getusermap(self->
                GXFILE_tgxfileobj_DOT_ueltable,self->
                GXFILE_tgxfileobj_DOT_lastelem[d - 1]);
              if (GXFILE_tdfilter_DOT_infilter(_W2->dfilter,v)) { 
                keyint[d - 1] = v;
              } else {
                adderror = SYSTEM_true;
                fidim = d;
                SYSTEM_break(BRK_6);
              } 
              break;
            case GXFILE_dm_strict: 
              v = GXFILE_tueltable_DOT_getusermap(self->
                GXFILE_tgxfileobj_DOT_ueltable,self->
                GXFILE_tgxfileobj_DOT_lastelem[d - 1]);
              if (v >= 0) { 
                keyint[d - 1] = v;
              } else {
                adderror = SYSTEM_true;
                fidim = d;
                SYSTEM_break(BRK_6);
              } 
              break;
            case GXFILE_dm_expand: 
              en = self->GXFILE_tgxfileobj_DOT_lastelem[d - 1];
              v = GXFILE_tueltable_DOT_getusermap(self->
                GXFILE_tgxfileobj_DOT_ueltable,en);
              if (v >= 0) { 
                keyint[d - 1] = v;
              } else {
                keyint[d - 1] = -en;
                addnew = SYSTEM_true;
              } 
              break;
            default: break;
          }

        }
CNT_6:;
      } while (d++ !=  _stop);
BRK_6:;

    }
  if (adderror) {
    GXFILE_tgxfileobj_DOT_addtoerrorlist(self,self->
      GXFILE_tgxfileobj_DOT_lastelem,values);
    adderror = SYSTEM_false;
    goto _Lagain_55;
  } 
  if (addnew) 
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        en = keyint[d - 1];
        if (en < 0) {
          v = GXFILE_tueltable_DOT_newusruel(self->
            GXFILE_tgxfileobj_DOT_ueltable,-en);
          keyint[d - 1] = v;
          self->GXFILE_tgxfileobj_DOT_nrmappedadded = self->
            GXFILE_tgxfileobj_DOT_nrmappedadded + 1;
          { register SYSTEM_int32 _stop = self->
              GXFILE_tgxfileobj_DOT_fcurrentdim;
            if ((d2 = d + 1) <=  _stop) do {
              if (keyint[d2 - 1] == en) 
                keyint[d2 - 1] = v;
            } while (d2++ !=  _stop);

          }
        } 
      
      } while (d++ !=  _stop);

    }
  *dimfrst = 0;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      if (self->GXFILE_tgxfileobj_DOT_prevelem[d - 1] != keyint[d - 1]) {
        self->GXFILE_tgxfileobj_DOT_prevelem[d - 1] = keyint[d - 1];
        if (*dimfrst == 0) 
          *dimfrst = d;
      } 
    } while (d++ !=  _stop);

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatareadmap */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadstr(
  GXFILE_tgxfileobj self,
  SYSTEM_shortstring *keystr,
  SYSTEM_double *values,
  SYSTEM_integer *dimfrst)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_9 = {0,128,0};
  SYSTEM_integer d;
  SYSTEM_integer led;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_9)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\013DataReadStr"),
      allowedmodes_9)) 
      return result;
  if (!GXFILE_tgxfileobj_DOT_doread(self,values,dimfrst)) { 
    GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  } else {
    result = GXFILE_mytrue;
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        led = self->GXFILE_tgxfileobj_DOT_lastelem[d - 1];
        if (led >= 1 && led <= self->
          GXFILE_tgxfileobj_DOT_ueltable->
          STRHASH_txstrhashlist_DOT_fcount) { 
          STRHASH_txstrhashlist_DOT_getstring(keystr[d - 1],255,ValueCast(
            STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
            led);
        } else 
          {
            SYSTEM_shortstring _t1;

            _P3strcat(keystr[d - 1],255,GXFILE_baduel_prefix,
              STRUTILX_inttostr(_t1,255,led));
          }
      
      } while (d++ !=  _stop);

    }
  } 
  return result;
}  /* gdxdatareadstr */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadslicestart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *elemcounts)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex xdomains;
  GXDEFS_tgdxvalues values;
  SYSTEM_integer _fdim;
  SYSTEM_integer cnt;
  SYSTEM_integer n;

  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_unmapped;
  }
  self->GXFILE_tgxfileobj_DOT_slicesynr = synr;
  GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\022DataReadSliceStart"),
    self->GXFILE_tgxfileobj_DOT_slicesynr,xdomains,GXFILE_fr_raw_data);
  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    elemcounts[d - 1] = 0;
  }
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      self->GXFILE_tgxfileobj_DOT_sliceindxs[d - 1] = ValueCast(
        GXFILE_tintegermapping,GXFILE_tintegermapping_DOT_create(ValueCast(
        GXFILE_tintegermapping,_P3alloc_object(&
        GXFILE_tintegermapping_CD))));
      self->GXFILE_tgxfileobj_DOT_slicerevmap[d - 1] = ValueCast(
        GXFILE_tintegermapping,GXFILE_tintegermapping_DOT_create(ValueCast(
        GXFILE_tintegermapping,_P3alloc_object(&
        GXFILE_tintegermapping_CD))));
    
    } while (d++ !=  _stop);

  }
  while (GXFILE_tgxfileobj_DOT_doread(self,values,&_fdim)) {

    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        GXFILE_tintegermapping_DOT_setmapping(self->
          GXFILE_tgxfileobj_DOT_sliceindxs[d - 1],self->
          GXFILE_tgxfileobj_DOT_lastelem[d - 1],1);
      } while (d++ !=  _stop);

    }
}
  GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      { register GXFILE_tintegermapping_OD *_W2=self->
        GXFILE_tgxfileobj_DOT_sliceindxs[d - 1];
        cnt = 0;
        { register SYSTEM_int32 _stop = _W2->
            GXFILE_tintegermapping_DOT_fhighestindex;
          if ((n = 0) <=  _stop) do {
            if (GXFILE_tintegermapping_DOT_getmapping(ValueCast(
              GXFILE_tintegermapping,_W2),n) >= 0) {
              GXFILE_tintegermapping_DOT_setmapping(ValueCast(
                GXFILE_tintegermapping,_W2),n,cnt);
              GXFILE_tintegermapping_DOT_setmapping(self->
                GXFILE_tgxfileobj_DOT_slicerevmap[d - 1],cnt,n);
              cnt = cnt + 1;
            } 
          } while (n++ !=  _stop);

        }

      }
      elemcounts[d - 1] = cnt;
    
    } while (d++ !=  _stop);

  }
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_slice;
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatareadslicestart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadslice(
  GXFILE_tgxfileobj self,
  const SYSTEM_shortstring *uelfilterstr,
  SYSTEM_integer *dimen,
  GXDEFS_tdatastoreproc dp)
{
  SYSTEM_integer result;
  GXDEFS_tgdxuelindex elemnrs;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex hisindx;
  SYSTEM_integer hisdim;
  GXDEFS_tgdxuelindex xdomains;
  GXDEFS_tgdxvalues values;
  SYSTEM_integer _fdim;
  SYSTEM_boolean goodindx;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\015DataReadSlice"),
    _P3set1("\0\0\002"))) 
    return result;
  goodindx = SYSTEM_true;
  *dimen = 0;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      _P3strcpy(self->GXFILE_tgxfileobj_DOT_sliceelems[d - 1],255,
        uelfilterstr[d - 1]);
      if (_P3strcmpE(uelfilterstr[d - 1],_P3str1("\000"))) {
        elemnrs[d - 1] =  -1;
        *dimen = *dimen + 1;
      } else {
        elemnrs[d - 1] = STRHASH_txstrhashlist_DOT_indexof(ValueCast(
          STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
          uelfilterstr[d - 1]);
        if (elemnrs[d - 1] < 0) 
          goodindx = SYSTEM_false;
      } 
    
    } while (d++ !=  _stop);

  }
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
  if (!goodindx) 
    return result;
  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_unmapped;
  }
  GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\015DataReadSlice"),
    self->GXFILE_tgxfileobj_DOT_slicesynr,xdomains,GXFILE_fr_slice);
  while (GXFILE_tgxfileobj_DOT_doread(self,values,&_fdim)) {
    goodindx = SYSTEM_true;
    hisdim = 0;
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        if (elemnrs[d - 1] ==  -1) {
          hisdim = hisdim + 1;
          hisindx[hisdim - 1] = 
            GXFILE_tintegermapping_DOT_getmapping(self->
            GXFILE_tgxfileobj_DOT_sliceindxs[d - 1],self->
            GXFILE_tgxfileobj_DOT_lastelem[d - 1]);
        } else 
          if (elemnrs[d - 1] != self->
            GXFILE_tgxfileobj_DOT_lastelem[d - 1]) 
            goodindx = SYSTEM_false;
      } while (d++ !=  _stop);

    }
    if (goodindx) 
      (*dp)(hisindx,values);
  
}
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatareadslice */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_shortstring *uelfilterstr,
  GXDEFS_tdatastorefiltproc dp)
{
  SYSTEM_integer result;
  GXDEFS_tgdxuelindex xdomains;
  SYSTEM_integer nrrecs;
  SYSTEM_integer filtdim;
  GXDEFS_tgdxuelindex elemdim;
  GXDEFS_tgdxuelindex elemnrs;
  GXDEFS_tgdxuelindex hisindx;
  SYSTEM_boolean goodindx;
  SYSTEM_integer d;
  GXDEFS_tgdxvalues values;
  SYSTEM_integer afdim;

  self->GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp = dp;
  result = GXFILE_myfalse;
  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_unmapped;
  }
  nrrecs = GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\026gdxDataReadRawFastFilt"),
    synr,xdomains,GXFILE_fr_raw_data);
  if (nrrecs >= 0) {
    goodindx = SYSTEM_true;
    filtdim = 0;
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        if (_P3strcmpN(uelfilterstr[d - 1],_P3str1("\000"))) {
          filtdim = filtdim + 1;
          elemdim[filtdim - 1] = d;
          elemnrs[filtdim - 1] = 
            STRHASH_txstrhashlist_DOT_indexof(ValueCast(
            STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
            uelfilterstr[d - 1]);
          if (elemnrs[filtdim - 1] < 0) 
            goodindx = SYSTEM_false;
        } 
      } while (d++ !=  _stop);

    }
    if (goodindx) {
      while (GXFILE_tgxfileobj_DOT_doread(self,values,&afdim)) {
        goodindx = SYSTEM_true;
        { register SYSTEM_int32 _stop = filtdim;
          if ((d = 1) <=  _stop) do {
            if (self->GXFILE_tgxfileobj_DOT_lastelem[elemdim[d - 1] - 1] != 
              elemnrs[d - 1]) {
              goodindx = SYSTEM_false;
              SYSTEM_break(BRK_7);
            } 
CNT_7:;
          } while (d++ !=  _stop);
BRK_7:;

        }
        if (goodindx) 
          if (GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_fc(self,
            self->GXFILE_tgxfileobj_DOT_lastelem,values,self) == 0) 
            SYSTEM_break(BRK_8);
      
CNT_8:;
      }
BRK_8:;
      result = GXFILE_mytrue;
    } 
    GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  } 
  return result;
}  /* gdxdatareadrawfastfilt */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatasliceuels(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *slicekeyint,
  SYSTEM_shortstring *keystr)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  SYSTEM_integer n;
  SYSTEM_integer hisdim;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\015DataSliceUELS"),
    _P3set1("\0\0\002"))) 
    return result;
  hisdim = 0;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      if (_P3strcmpN(self->GXFILE_tgxfileobj_DOT_sliceelems[d - 1],_P3str1("\000"))) { 
        _P3strcpy(keystr[d - 1],255,self->
          GXFILE_tgxfileobj_DOT_sliceelems[d - 1]);
      } else {
        hisdim = hisdim + 1;
        n = GXFILE_tintegermapping_DOT_getmapping(self->
          GXFILE_tgxfileobj_DOT_slicerevmap[d - 1],slicekeyint[
          hisdim - 1]);
        if (n < 0) { 
          _P3strcpy(keystr[d - 1],255,_P3str1("\001?"));
        } else 
          STRHASH_txstrhashlist_DOT_getstring(keystr[d - 1],255,ValueCast(
            STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
            n);
      } 
    } while (d++ !=  _stop);

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatasliceuels */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareaddone(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_10 = {2,240,2};
  SYSTEM_integer n;
  SYSTEM_integer en;
  SYSTEM_integer d;

  result = GXFILE_myfalse;
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tgxfileobj_DOT_sortlist));
  self->GXFILE_tgxfileobj_DOT_sortlist = NULL;
  self->GXFILE_tgxfileobj_DOT_cursyptr = NULL;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\014DataReadDone"),
    allowedmodes_10)) {
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
    return result;
  } 
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fr_slice) 
    for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
      SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_sliceindxs[d - 1]);
      SYSUTILS_P3_freeandnil(&self->GXFILE_tgxfileobj_DOT_slicerevmap[
        d - 1]);
    
    }
  if (self->GXFILE_tgxfileobj_DOT_nrmappedadded > 0) {
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_ueltable->
        GXFILE_tueltable_DOT_usruel2ent->
        GXFILE_tintegermapping_DOT_fhighestindex - self->
        GXFILE_tgxfileobj_DOT_nrmappedadded + 1;
      if ((n = self->GXFILE_tgxfileobj_DOT_ueltable->
        GXFILE_tueltable_DOT_usruel2ent->
        GXFILE_tintegermapping_DOT_fhighestindex) >=  _stop) do {
        en = GXFILE_tintegermapping_DOT_getmapping(self->
          GXFILE_tgxfileobj_DOT_ueltable->
          GXFILE_tueltable_DOT_usruel2ent,n);
        SYSTEM_assert(n >= 1,_P3str1("\016Wrong entry Nr"));
        SYSTEM_assert(GMSOBJ_copyptr2int(
          STRHASH_txstrhashlist_DOT_getobject(ValueCast(
          STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
          en)) ==  -1 || GMSOBJ_copyptr2int(
          STRHASH_txstrhashlist_DOT_getobject(ValueCast(
          STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
          en)) == n,_P3str1("\016Mapped already"));
        STRHASH_txstrhashlist_DOT_setobject(ValueCast(
          STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
          en,ValueCast(SYSTEM_tobject,GMSOBJ_copyint2ptr(n)));
      
      } while (n-- !=  _stop);

    }
    self->GXFILE_tgxfileobj_DOT_nrmappedadded = 0;
  } 
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
  result = GXFILE_mytrue;
  return result;
}  /* gdxdatareaddone */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdataerrorcount(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_errorlist == NULL) { 
    result = 0;
  } else 
    result = GMSDATA_ttblgamsdata_DOT_getcount(self->
      GXFILE_tgxfileobj_DOT_errorlist);
  return result;
}  /* gdxdataerrorcount */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdataerrorrecord(
  GXFILE_tgxfileobj self,
  SYSTEM_integer recnr,
  SYSTEM_integer *keyint,
  SYSTEM_double *values)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_11 = {198,97,0};

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_11)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\017DataErrorRecord"),
      allowedmodes_11)) 
      return result;
  if (self->GXFILE_tgxfileobj_DOT_errorlist != NULL) 
    if (recnr < 1 || recnr > GMSDATA_ttblgamsdata_DOT_getcount(
      self->GXFILE_tgxfileobj_DOT_errorlist)) { 
      GXFILE_tgxfileobj_DOT_reporterror(self,GXFILE_err_baderrorrecord);
    } else {
      GMSDATA_ttblgamsdata_DOT_getrecord(self->
        GXFILE_tgxfileobj_DOT_errorlist,recnr - 1,keyint,values);
      result = GXFILE_mytrue;
    } 
  return result;
}  /* gdxdataerrorrecord */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterregisterstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer filternr)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_12 = {2,0,0};

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\023FilterRegisterStart"),
    allowedmodes_12)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,filternr >= 1,
    GXFILE_err_bad_filter_nr)) 
    return result;
  self->GXFILE_tgxfileobj_DOT_curfilter = ValueCast(GXFILE_tdfilter,
    GXFILE_tdfilter_DOT_create(ValueCast(GXFILE_tdfilter,
    _P3alloc_object(&GXFILE_tdfilter_CD)),filternr,self->
    GXFILE_tgxfileobj_DOT_ueltable->GXFILE_tueltable_DOT_usruel2ent->
    GXFILE_tintegermapping_DOT_fhighestindex));
  GXFILE_tfilterlist_DOT_addfilter(self->
    GXFILE_tgxfileobj_DOT_filterlist,self->
    GXFILE_tgxfileobj_DOT_curfilter);
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_filter;
  result = GXFILE_mytrue;
  return result;
}  /* gdxfilterregisterstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterregister(
  GXFILE_tgxfileobj self,
  SYSTEM_integer uelmap)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_13 = {0,0,1};
  SYSTEM_integer en;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_13)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\016FilterRegister"),
      allowedmodes_13)) 
      return result;
  { register GXFILE_tdfilter_OD *_W2=self->
    GXFILE_tgxfileobj_DOT_curfilter;
    if (GXFILE_tgxfileobj_DOT_errorcondition(self,uelmap >= 1 && 
      uelmap <= _W2->GXFILE_tdfilter_DOT_filtmaxuel,
      GXFILE_err_bad_filter_indx)) 
      return result;
    en = GXFILE_tintegermapping_DOT_getmapping(self->
      GXFILE_tgxfileobj_DOT_ueltable->GXFILE_tueltable_DOT_usruel2ent,
      uelmap);
    if (en >= 1) { 
      GMSOBJ_tbooleanbitarray_DOT_setbit(_W2->
        GXFILE_tdfilter_DOT_filtmap,uelmap,SYSTEM_true);
    } else {
      GXFILE_tgxfileobj_DOT_reporterror(self,
        GXFILE_err_filter_unmapped);
      return result;
    } 

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxfilterregister */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterregisterdone(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_14 = {0,0,1};
  SYSTEM_integer lv;
  SYSTEM_integer n;
  SYSTEM_integer v;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\022FilterRegisterDone"),
    allowedmodes_14)) 
    return result;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
  result = GXFILE_mytrue;
  self->GXFILE_tgxfileobj_DOT_curfilter->
    GXFILE_tdfilter_DOT_filtsorted = SYSTEM_true;
  if (GXFILE_tueltable_DOT_getmaptouserstatus(self->
    GXFILE_tgxfileobj_DOT_ueltable) == GXFILE_map_unsorted) {
    lv =  -1;
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_ueltable->
        STRHASH_txstrhashlist_DOT_fcount;
      if ((n = 1) <=  _stop) do {
        v = GXFILE_tueltable_DOT_getusermap(self->
          GXFILE_tgxfileobj_DOT_ueltable,n);
        if (!GXFILE_tdfilter_DOT_infilter(self->
          GXFILE_tgxfileobj_DOT_curfilter,v)) 
          SYSTEM_continue(CNT_9);
        if (v <= lv) {
          self->GXFILE_tgxfileobj_DOT_curfilter->
            GXFILE_tdfilter_DOT_filtsorted = SYSTEM_false;
          SYSTEM_break(BRK_9);
        } 
        lv = v;
      
CNT_9:;
      } while (n++ !=  _stop);
BRK_9:;

    }
  } 
  self->GXFILE_tgxfileobj_DOT_curfilter = NULL;
  return result;
}  /* gdxfilterregisterdone */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterexists(
  GXFILE_tgxfileobj self,
  SYSTEM_integer filternr)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\014FilterExists"),
    GXFILE_anyreadmode)) 
    return result;
  if (GXFILE_tfilterlist_DOT_findfilter(self->
    GXFILE_tgxfileobj_DOT_filterlist,filternr) != NULL) 
    result = GXFILE_mytrue;
  return result;
}  /* gdxfilterexists */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxumuelinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *uelcnt,
  SYSTEM_integer *highmap)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_ueltable == NULL) {
    *uelcnt = 0;
    *highmap = 0;
    result = GXFILE_myfalse;
  } else {
    *uelcnt = self->GXFILE_tgxfileobj_DOT_ueltable->
      STRHASH_txstrhashlist_DOT_fcount;
    *highmap = self->GXFILE_tgxfileobj_DOT_ueltable->
      GXFILE_tueltable_DOT_usruel2ent->
      GXFILE_tintegermapping_DOT_fhighestindex;
    result = GXFILE_mytrue;
  } 
  return result;
}  /* gdxumuelinfo */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxumuelget(
  GXFILE_tgxfileobj self,
  SYSTEM_integer uelnr,
  SYSTEM_ansichar *uel,
  SYSTEM_integer *uelmap)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_ueltable != NULL && uelnr >= 1 && 
    uelnr <= self->GXFILE_tgxfileobj_DOT_ueltable->
    STRHASH_txstrhashlist_DOT_fcount) {
    result = GXFILE_mytrue;
    STRHASH_txstrhashlist_DOT_getstring(uel,255,ValueCast(
      STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
      uelnr);
    *uelmap = GMSOBJ_copyptr2int(STRHASH_txstrhashlist_DOT_getobject(ValueCast(
      STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
      uelnr));
  } else {
    result = GXFILE_myfalse;
    {
      SYSTEM_shortstring _t1;

      _P3strcat(uel,255,GXFILE_baduel_prefix,STRUTILX_inttostr(_t1,255,
        uelnr));
    }
    *uelmap =  -1;
  } 
  return result;
}  /* gdxumuelget */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxumfinduel(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *uel,
  SYSTEM_integer *uelnr,
  SYSTEM_integer *uelmap)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  *uelmap =  -1;
  if (self->GXFILE_tgxfileobj_DOT_ueltable == NULL) {
    *uelnr =  -1;
    return result;
  } 
  {
    SYSTEM_shortstring _t1;

    *uelnr = STRHASH_txstrhashlist_DOT_indexof(ValueCast(
      STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
      SYSUTILS_P3_trimright(_t1,255,uel));
  }
  if (*uelnr < 0) 
    return result;
  *uelmap = GMSOBJ_copyptr2int(STRHASH_txstrhashlist_DOT_getobject(ValueCast(
    STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),*uelnr));
  result = GXFILE_mytrue;
  return result;
}  /* gdxumfinduel */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxaddsettext(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *txt,
  SYSTEM_integer *txtnr)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  *txtnr = 0;
  if (self->GXFILE_tgxfileobj_DOT_settextlist == NULL) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\012AddSetText"),
      _P3empty_set)) 
      return result;
  {
    SYSTEM_shortstring _t1;

    *txtnr = GMSOBJ_txhashedstringlist_DOT_add(ValueCast(
      GMSOBJ_txhashedstringlist,self->
      GXFILE_tgxfileobj_DOT_settextlist),GXFILE_makegoodexpltext(_t1,255,
      txt));
  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxaddsettext */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetelemtext(
  GXFILE_tgxfileobj self,
  SYSTEM_integer txtnr,
  SYSTEM_ansichar *txt,
  SYSTEM_integer *node)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  *node = 0;
  if (self->GXFILE_tgxfileobj_DOT_settextlist == NULL) {
    _P3strclr(txt);
    return result;
  } 
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\013GetElemText"),
      _P3empty_set)) 
      return result;
  { register GMSOBJ_txstrpool_OD *_W2=self->
    GXFILE_tgxfileobj_DOT_settextlist;
    if (!(txtnr > 0 && txtnr < _W2->
      GMSOBJ_txcustomstringlist_DOT_fcount)) { 
      {
        SYSTEM_shortstring _t1;

        _P3strcat(txt,255,GXFILE_badstr_prefix,STRUTILX_inttostr(
          _t1,255,txtnr));
      }
    } else {
      result = GXFILE_mytrue;
      GMSOBJ_txcustomstringlist_DOT_getname(txt,255,ValueCast(
        GMSOBJ_txcustomstringlist,_W2),txtnr);
      *node = GMSOBJ_copyptr2int(
        GMSOBJ_txcustomstringlist_DOT_getobject(ValueCast(
        GMSOBJ_txcustomstringlist,_W2),txtnr));
    } 

  }
  return result;
}  /* gdxgetelemtext */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsettextnodenr(
  GXFILE_tgxfileobj self,
  SYSTEM_integer txtnr,
  SYSTEM_integer node)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_settextlist == NULL) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\015SetTextNodeNr"),
      _P3empty_set)) 
      return result;
  { register GMSOBJ_txstrpool_OD *_W2=self->
    GXFILE_tgxfileobj_DOT_settextlist;
    if (txtnr > 0 && txtnr < _W2->
      GMSOBJ_txcustomstringlist_DOT_fcount && 
      GMSOBJ_txcustomstringlist_DOT_getobject(ValueCast(
      GMSOBJ_txcustomstringlist,_W2),txtnr) == NULL) {
      result = GXFILE_mytrue;
      GMSOBJ_txcustomstringlist_DOT_putobject(ValueCast(
        GMSOBJ_txcustomstringlist,_W2),txtnr,ValueCast(SYSTEM_tobject,
        GMSOBJ_copyint2ptr(node)));
    } 

  }
  return result;
}  /* gdxsettextnodenr */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsethastext(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && 
    synr <= self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount && (ValueCast(
    GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
    GXFILE_tgxfileobj_DOT_namelist,synr)))->ssettext) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxsethastext */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxresetspecialvalues(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;

  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valund] = 
    GMSSPECS_valund;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valna] = 
    GMSSPECS_valna;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valpin] = 
    GMSSPECS_valpin;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valmin] = 
    GMSSPECS_valmin;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valeps] = 
    GMSSPECS_valeps;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_zero] = 0.0;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_one] = 1.0;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_mone] = -1.0;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_half] = 0.5;
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_two] = 2.0;
  _P3memcpy(self->GXFILE_tgxfileobj_DOT_readintvaluemap,sizeof(
    _arr_1GXFILE),self->GXFILE_tgxfileobj_DOT_intvaluemap);
  self->GXFILE_tgxfileobj_DOT_zvalacr = GMSSPECS_valacr;
  result = GXFILE_mytrue;
  return result;
}  /* gdxresetspecialvalues */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsetspecialvalues(
  GXFILE_tgxfileobj self,
  const SYSTEM_double *avals)
{
  SYSTEM_integer result;
  GMSSPECS_tgdxspecialvalue sv1;
  GMSSPECS_tgdxspecialvalue sv2;

  result = GXFILE_myfalse;
  { register GMSSPECS_tgdxspecialvalue _stop = GMSSPECS_sv_normal;
    if ((sv1 = GMSSPECS_sv_valund) <=  _stop) do {
      { register GMSSPECS_tgdxspecialvalue _stop = GMSSPECS_sv_valeps;
        if ((sv2 = SYSTEM_succ(sv1)) <=  _stop) do {
          if (avals[sv1] == avals[sv2]) {
            GXFILE_tgxfileobj_DOT_reporterror(self,
              GXFILE_err_duplicatespecval);
            return result;
          } 
        } while (sv2++ !=  _stop);

      }
    } while (sv1++ !=  _stop);

  }
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valund] = avals[
    GMSSPECS_sv_valund];
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valna] = avals[
    GMSSPECS_sv_valna];
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valpin] = avals[
    GMSSPECS_sv_valpin];
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valmin] = avals[
    GMSSPECS_sv_valmin];
  self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valeps] = avals[
    GMSSPECS_sv_valeps];
  _P3memcpy(self->GXFILE_tgxfileobj_DOT_readintvaluemap,sizeof(
    _arr_1GXFILE),self->GXFILE_tgxfileobj_DOT_intvaluemap);
  result = GXFILE_mytrue;
  return result;
}  /* gdxsetspecialvalues */

Function(SYSTEM_integer ) 
  GXFILE_tgxfileobj_DOT_gdxsetreadspecialvalues(
  GXFILE_tgxfileobj self,
  const SYSTEM_double *avals)
{
  SYSTEM_integer result;

  self->GXFILE_tgxfileobj_DOT_readintvaluemap[GXFILE_vm_valund] = 
    avals[GMSSPECS_sv_valund];
  self->GXFILE_tgxfileobj_DOT_readintvaluemap[GXFILE_vm_valna] = avals[
    GMSSPECS_sv_valna];
  self->GXFILE_tgxfileobj_DOT_readintvaluemap[GXFILE_vm_valpin] = 
    avals[GMSSPECS_sv_valpin];
  self->GXFILE_tgxfileobj_DOT_readintvaluemap[GXFILE_vm_valmin] = 
    avals[GMSSPECS_sv_valmin];
  self->GXFILE_tgxfileobj_DOT_readintvaluemap[GXFILE_vm_valeps] = 
    avals[GMSSPECS_sv_valeps];
  result = GXFILE_mytrue;
  return result;
}  /* gdxsetreadspecialvalues */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetspecialvalues(
  GXFILE_tgxfileobj self,
  SYSTEM_double *avals)
{
  SYSTEM_integer result;

  avals[GMSSPECS_sv_valund] = self->GXFILE_tgxfileobj_DOT_intvaluemap[
    GXFILE_vm_valund];
  avals[GMSSPECS_sv_valna] = self->GXFILE_tgxfileobj_DOT_intvaluemap[
    GXFILE_vm_valna];
  avals[GMSSPECS_sv_valpin] = self->GXFILE_tgxfileobj_DOT_intvaluemap[
    GXFILE_vm_valpin];
  avals[GMSSPECS_sv_valmin] = self->GXFILE_tgxfileobj_DOT_intvaluemap[
    GXFILE_vm_valmin];
  avals[GMSSPECS_sv_valeps] = self->GXFILE_tgxfileobj_DOT_intvaluemap[
    GXFILE_vm_valeps];
  avals[GMSSPECS_sv_acronym] = self->GXFILE_tgxfileobj_DOT_zvalacr;
  result = GXFILE_mytrue;
  return result;
}  /* gdxgetspecialvalues */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxmapvalue(
  GXFILE_tgxfileobj self,
  SYSTEM_double d,
  SYSTEM_integer *sv)
{
  SYSTEM_integer result;

  result = GXFILE_mytrue;
  if (d == self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valund]) { 
    *sv = 0;
  } else 
    if (d == self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_valna]) { 
      *sv = 1;
    } else 
      if (d == self->GXFILE_tgxfileobj_DOT_intvaluemap[
        GXFILE_vm_valpin]) { 
        *sv = 2;
      } else 
        if (d == self->GXFILE_tgxfileobj_DOT_intvaluemap[
          GXFILE_vm_valmin]) { 
          *sv = 3;
        } else 
          if (d == self->GXFILE_tgxfileobj_DOT_intvaluemap[
            GXFILE_vm_valeps]) { 
            *sv = 4;
          } else {
            *sv = 5;
            result = GXFILE_myfalse;
          } 
  return result;
}  /* gdxmapvalue */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_majorcheckmode(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *routine,
  const _P3set_elem *ms)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  _P3strcpy(self->GXFILE_tgxfileobj_DOT_majcontext,63,routine);
  self->GXFILE_tgxfileobj_DOT_lastreperror = GXFILE_err_noerror;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_some || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,ms)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,routine,ms)) 
      return result;
  result = SYSTEM_true;
  return result;
}  /* majorcheckmode */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_checkmode(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *routine,
  const _P3set_elem *ms)
{
  SYSTEM_boolean result;
  GXFILE_tgxfilemode m;
  SYSTEM_boolean f;

  result = _P3setcmpE(17,ms,_P3empty_set) || _P3SET_i(17,self->
    GXFILE_tgxfileobj_DOT_fmode,ms);
  if (result) { 
    GXFILE_tgxfileobj_DOT_writetrace(self,routine);
  } else {
    GXFILE_tgxfileobj_DOT_seterror(self,GXFILE_err_badmode);
    _Iplus_bgn();
    _P3write_s0(_P3str1("\014**** Error: "));
    _P3write_s0(routine);
    _P3write_s0(_P3str1("\026 called out of context"));
    _P3writeln();
    _Iplus_end();
    if (_P3strcmpN(self->GXFILE_tgxfileobj_DOT_majcontext,_P3str1("\000")) && 
      _P3strcmpN(self->GXFILE_tgxfileobj_DOT_majcontext,routine)) {
      _Iplus_bgn();
      _P3write_s0(_P3str1("\050     Previous major function called was "));
      _P3write_s0(self->GXFILE_tgxfileobj_DOT_majcontext);
      _P3writeln();
      _Iplus_end();
    } 
    _Iplus_bgn();
    _P3write_s0(_P3str1("\027     Current context = "));
    _P3write_s0(GXFILE_fmode_str[self->GXFILE_tgxfileobj_DOT_fmode]);
    _P3writeln();
    _Iplus_end();
    _Iplus_bgn();
    _P3write_s0(_P3str1("\020     Allowed = {"));
    _Iplus_end();
    f = SYSTEM_true;
    if ((m = GXFILE_f_not_open) <= (GXFILE_fr_slice)) do {
      if (_P3SET_i(17,m,ms)) {
        if (f) { 
          f = SYSTEM_false;
        } else {
          _Iplus_bgn();
          _P3write_c0(_P3char(','));
          _Iplus_end();
        } 
        _Iplus_bgn();
        _P3write_s0(GXFILE_fmode_str[m]);
        _Iplus_end();
      } 
    } while (m++ != (GXFILE_fr_slice));

    _Iplus_bgn();
    _P3write_c0(_P3char('}'));
    _P3writeln();
    _Iplus_end();
  } 
  return result;
}  /* checkmode */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsettracelevel(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  if (n <= 0) { 
    self->GXFILE_tgxfileobj_DOT_tracelevel = GXFILE_trl_none;
  } else {
    switch (n) {
      case 1: 
        self->GXFILE_tgxfileobj_DOT_tracelevel = GXFILE_trl_errors;
        break;
      case 2: 
        self->GXFILE_tgxfileobj_DOT_tracelevel = GXFILE_trl_some;
        break;
      default:
        self->GXFILE_tgxfileobj_DOT_tracelevel = GXFILE_trl_all;
    }
    _P3strcpy(self->GXFILE_tgxfileobj_DOT_tracestr,255,s);
  } 
  result = GXFILE_mytrue;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel > GXFILE_trl_errors) {
    _Iplus_bgn();
    _P3writeln();
    _Iplus_end();
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;

      GXFILE_tgxfileobj_DOT_writetrace(self,_P3strcat(_t2,255,_P3str1("\021Tracing at level "),
        STRUTILX_inttostr(_t1,255,SYSTEM_ord(self->
        GXFILE_tgxfileobj_DOT_tracelevel))));
    }
  } 
  return result;
}  /* gdxsettracelevel */

Procedure GXFILE_tgxfileobj_DOT_writetrace(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *s)
{
  _Iplus_bgn();
  _P3write_s0(_P3str1("\011gdxTrace "));
  _P3write_s0(self->GXFILE_tgxfileobj_DOT_tracestr);
  _P3write_s0(_P3str1("\002: "));
  _P3write_s0(s);
  _P3writeln();
  _Iplus_end();
}  /* writetrace */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfileversion(
  GXFILE_tgxfileobj self,
  SYSTEM_ansichar *filestr,
  SYSTEM_ansichar *producestr)
{
  SYSTEM_integer result;

  result = GXFILE_mytrue;
  _P3strcpy(filestr,255,self->GXFILE_tgxfileobj_DOT_filesystemid);
  _P3strcpy(producestr,255,self->GXFILE_tgxfileobj_DOT_fproducer);
  return result;
}  /* gdxfileversion */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterrawstart(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_15 = {4,0,0};

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\023UELRegisterRawStart"),
    allowedmodes_15)) 
    return result;
  self->GXFILE_tgxfileobj_DOT_fmode_aftreg = GXFILE_fw_init;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_raw_elem;
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregisterrawstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterraw(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *uel)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_16 = {0,2,0};
  SYSTEM_shortstring sv;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_16)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\016UELRegisterRaw"),
      allowedmodes_16)) 
      return result;
  SYSUTILS_P3_trimright(sv,255,uel);
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,GXDEFS_gooduelstring(
    sv),GXFILE_err_baduelstr)) 
    return result;
  STRHASH_txstrhashlist_DOT_addobject(ValueCast(STRHASH_txstrhashlist,
    self->GXFILE_tgxfileobj_DOT_ueltable),sv,ValueCast(SYSTEM_tobject,
    GMSOBJ_copyint2ptr( -1)));
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregisterraw */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregistermapstart(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_17 = {6,0,0};

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\023UELRegisterMapStart"),
    allowedmodes_17)) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fw_init) { 
    self->GXFILE_tgxfileobj_DOT_fmode_aftreg = GXFILE_fw_init;
  } else 
    self->GXFILE_tgxfileobj_DOT_fmode_aftreg = GXFILE_fr_init;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_map_elem;
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregistermapstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregistermap(
  GXFILE_tgxfileobj self,
  SYSTEM_integer umap,
  const SYSTEM_ansichar *uel)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_18 = {0,4,0};
  SYSTEM_shortstring sv;

  result = GXFILE_myfalse;
  SYSUTILS_P3_trimright(sv,255,uel);
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_18)) {
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\016UELRegisterMap"),
      allowedmodes_18)) 
      return result;
    _Iplus_bgn();
    _P3write_s0(_P3str1("\016   Enter UEL: "));
    _P3write_s0(sv);
    _P3write_s0(_P3str1("\015 with number "));
    _P3write_i0(umap);
    _P3writeln();
    _Iplus_end();
  } 
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,GXDEFS_gooduelstring(
    sv),GXFILE_err_baduelstr)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,
    GXFILE_tueltable_DOT_addusrindxnew(self->
    GXFILE_tgxfileobj_DOT_ueltable,sv,umap) >= 0,
    GXFILE_err_uelconflict)) 
    return result;
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregistermap */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterstrstart(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_19 = {6,0,0};

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\023UELRegisterStrStart"),
    allowedmodes_19)) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_fmode == GXFILE_fw_init) { 
    self->GXFILE_tgxfileobj_DOT_fmode_aftreg = GXFILE_fw_init;
  } else 
    self->GXFILE_tgxfileobj_DOT_fmode_aftreg = GXFILE_fr_init;
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_str_elem;
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregisterstrstart */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterstr(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *uel,
  SYSTEM_integer *uelnr)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_20 = {0,8,0};
  SYSTEM_shortstring sv;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_all || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_20)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\016UELRegisterStr"),
      allowedmodes_20)) 
      return result;
  SYSUTILS_P3_trimright(sv,255,uel);
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,GXDEFS_gooduelstring(
    sv),GXFILE_err_baduelstr)) 
    return result;
  *uelnr = GXFILE_tueltable_DOT_addusrnew(self->
    GXFILE_tgxfileobj_DOT_ueltable,sv);
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregisterstr */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterdone(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_21 = {0,14,0};

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\017UELRegisterDone"),
    allowedmodes_21)) 
    return result;
  self->GXFILE_tgxfileobj_DOT_fmode = self->
    GXFILE_tgxfileobj_DOT_fmode_aftreg;
  result = GXFILE_mytrue;
  return result;
}  /* gdxuelregisterdone */

Procedure GXFILE_tgxfileobj_DOT_initdowrite(
  GXFILE_tgxfileobj self,
  SYSTEM_integer nrrecs)
{
  SYSTEM_integer d;

  self->GXFILE_tgxfileobj_DOT_datacount = 0;
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),self->
    GXFILE_tgxfileobj_DOT_nextwriteposition));
  self->GXFILE_tgxfileobj_DOT_cursyptr->sposition = self->
    GXFILE_tgxfileobj_DOT_nextwriteposition;
  GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),GXFILE_mark_data);
  GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),self->
    GXFILE_tgxfileobj_DOT_fcurrentdim);
  GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile),nrrecs);
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = 
        GXFILE_index_initial;
      self->GXFILE_tgxfileobj_DOT_elemtype[d - 1] = 
        GXFILE_getintegersize(self->GXFILE_tgxfileobj_DOT_maxelem[d - 1] - 
        self->GXFILE_tgxfileobj_DOT_minelem[d - 1] + 1);
      GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
        self->GXFILE_tgxfileobj_DOT_ffile),self->
        GXFILE_tgxfileobj_DOT_minelem[d - 1]);
      GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,
        self->GXFILE_tgxfileobj_DOT_ffile),self->
        GXFILE_tgxfileobj_DOT_maxelem[d - 1]);
    
    } while (d++ !=  _stop);

  }
}  /* initdowrite */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_dowrite(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *aelements,
  const SYSTEM_double *avals)
{
  SYSTEM_boolean result;
  SYSTEM_integer d;
  SYSTEM_integer _fdim;
  SYSTEM_integer delta;
  GMSSPECS_tvarvaltype dv;
  SYSTEM_double x;
  GXFILE_tgdxintvaltyp xv;
  SYSTEM_byte bxv;

  result = SYSTEM_true;
  _fdim = self->GXFILE_tgxfileobj_DOT_fcurrentdim + 1;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      if (self->GXFILE_tgxfileobj_DOT_wrbitmaps[d - 1] != NULL && !
        GMSOBJ_tbooleanbitarray_DOT_getbit(self->
        GXFILE_tgxfileobj_DOT_wrbitmaps[d - 1],aelements[d - 1])) {
        GXFILE_tgxfileobj_DOT_reporterror(self,
          GXFILE_err_domainviolation);
        GXFILE_tgxfileobj_DOT_addtoerrorlist(self,aelements,avals);
        result = SYSTEM_false;
        return result;
      } 
      delta = aelements[d - 1] - self->
        GXFILE_tgxfileobj_DOT_lastelem[d - 1];
      if (delta != 0) {
        _fdim = d;
        SYSTEM_break(BRK_10);
      } 
    
CNT_10:;
    } while (d++ !=  _stop);
BRK_10:;

  }
  if (_fdim > self->GXFILE_tgxfileobj_DOT_fcurrentdim) {
    if (self->GXFILE_tgxfileobj_DOT_fcurrentdim > 0 || self->
      GXFILE_tgxfileobj_DOT_datacount >= 1) {
      GXFILE_tgxfileobj_DOT_reporterror(self,GXFILE_err_dataduplicate);
      GXFILE_tgxfileobj_DOT_addtoerrorlist(self,aelements,avals);
      result = SYSTEM_false;
      return result;
    } 
    GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
      GXFILE_tgxfileobj_DOT_ffile),1);
  } else {
    if (delta < 0) {
      GXFILE_tgxfileobj_DOT_reporterror(self,GXFILE_err_rawnotsorted);
      GXFILE_tgxfileobj_DOT_addtoerrorlist(self,aelements,avals);
      result = SYSTEM_false;
      return result;
    } 
    if (_fdim == self->GXFILE_tgxfileobj_DOT_fcurrentdim && delta <= 
      self->GXFILE_tgxfileobj_DOT_deltaforwrite) {
      GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
        GXFILE_tgxfileobj_DOT_ffile),self->
        GXFILE_tgxfileobj_DOT_fcurrentdim + delta);
      self->GXFILE_tgxfileobj_DOT_lastelem[self->
        GXFILE_tgxfileobj_DOT_fcurrentdim - 1] = aelements[self->
        GXFILE_tgxfileobj_DOT_fcurrentdim - 1];
    } else {
      GMSSTRM_txstream_DOT_writebyte(ValueCast(GMSSTRM_txstream,self->
        GXFILE_tgxfileobj_DOT_ffile),_fdim);
      { register SYSTEM_int32 _stop = self->
          GXFILE_tgxfileobj_DOT_fcurrentdim;
        if ((d = _fdim) <=  _stop) do {
          switch (self->GXFILE_tgxfileobj_DOT_elemtype[d - 1]) {
            case GXFILE_sz_integer: 
              GMSSTRM_txstream_DOT_writeinteger(ValueCast(
                GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),
                aelements[d - 1] - self->
                GXFILE_tgxfileobj_DOT_minelem[d - 1]);
              break;
            case GXFILE_sz_word: 
              GMSSTRM_txstream_DOT_writeword(ValueCast(
                GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),
                aelements[d - 1] - self->
                GXFILE_tgxfileobj_DOT_minelem[d - 1]);
              break;
            case GXFILE_sz_byte: 
              GMSSTRM_txstream_DOT_writebyte(ValueCast(
                GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),
                aelements[d - 1] - self->
                GXFILE_tgxfileobj_DOT_minelem[d - 1]);
              break;
            default: break;
          }
          self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = aelements[
            d - 1];
        
        } while (d++ !=  _stop);

      }
    } 
  } 
  if (self->GXFILE_tgxfileobj_DOT_datasize > 0) 
    { register GMSSPECS_tvarvaltype _stop = self->
        GXFILE_tgxfileobj_DOT_lastdatafield;
      if ((dv = GMSSPECS_vallevel) <=  _stop) do {
        x = avals[dv];
        if (!P3IEEEFP_p3isfinite(x)) { 
          switch (P3IEEEFP_fpclass(x)) {
            case P3IEEEFP_fp_ninf: 
              xv = GXFILE_vm_valmin;
              break;
            case P3IEEEFP_fp_pinf: 
              xv = GXFILE_vm_valpin;
              break;
            case P3IEEEFP_fp_snan: 
            case P3IEEEFP_fp_qnan: 
              xv = GXFILE_vm_valna;
              break;
            default:
              xv = GXFILE_vm_valund;
          }
        } else {
          self->GXFILE_tgxfileobj_DOT_intvaluemap[GXFILE_vm_normal] = 
            x;
          xv = GXFILE_vm_valund;
          while (x != self->GXFILE_tgxfileobj_DOT_intvaluemap[xv]) {

            _P3inc0(xv);
}
        } 
        bxv = SYSTEM_ord(xv);
        VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
          GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_write_T, 5, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),&bxv,1));
        if (xv == GXFILE_vm_normal) {
          GMSSTRM_txstream_DOT_writedouble(ValueCast(GMSSTRM_txstream,
            self->GXFILE_tgxfileobj_DOT_ffile),x);
          if (x >= self->GXFILE_tgxfileobj_DOT_zvalacr) 
            GXFILE_tacronymlist_DOT_checkentry(self->
              GXFILE_tgxfileobj_DOT_acronymlist,SYSTEM_round(x /  self->
              GXFILE_tgxfileobj_DOT_zvalacr));
        } 
      
      } while (dv++ !=  _stop);

    }
  _P3inc0(self->GXFILE_tgxfileobj_DOT_datacount);
  if (_P3SET_in_1(self->GXFILE_tgxfileobj_DOT_cursyptr->sdatatype,
    GMSSPECS_dt_set,_P3SET_equal(self->GXFILE_tgxfileobj_DOT_cursyptr->
    sdatatype,GMSSPECS_dt_alias))) {
    if (avals[GMSSPECS_vallevel] != 0.0) 
      self->GXFILE_tgxfileobj_DOT_cursyptr->ssettext = SYSTEM_true;
    if (self->GXFILE_tgxfileobj_DOT_fcurrentdim == 1) 
      GMSOBJ_tbooleanbitarray_DOT_setbit(self->
        GXFILE_tgxfileobj_DOT_cursyptr->ssetbitmap,self->
        GXFILE_tgxfileobj_DOT_lastelem[0],SYSTEM_true);
  } 
  return result;
}  /* dowrite */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_doread(
  GXFILE_tgxfileobj self,
  SYSTEM_double *avals,
  SYSTEM_integer *afdim)
{
  SYSTEM_boolean result;
  SYSTEM_byte b;
  SYSTEM_integer d;
  GMSSPECS_tvarvaltype dv;
  GXFILE_tgdxintvaltyp sv;
  SYSTEM_byte bsv;

  if (self->GXFILE_tgxfileobj_DOT_readuniverse) {
    self->GXFILE_tgxfileobj_DOT_universenr = self->
      GXFILE_tgxfileobj_DOT_universenr + 1;
    result = self->GXFILE_tgxfileobj_DOT_universenr <= self->
      GXFILE_tgxfileobj_DOT_uelcntorig;
    if (result) {
      self->GXFILE_tgxfileobj_DOT_lastelem[0] = self->
        GXFILE_tgxfileobj_DOT_universenr;
      avals[GMSSPECS_vallevel] = 0.0;
      *afdim = 1;
    } 
    return result;
  } 
  if (self->GXFILE_tgxfileobj_DOT_cursyptr->sscalarfrst) {
    self->GXFILE_tgxfileobj_DOT_cursyptr->sscalarfrst = SYSTEM_false;
    GXFILE_tgxfileobj_DOT_getdefaultrecord(self,avals);
    *afdim = 0;
    result = SYSTEM_true;
    return result;
  } 
  result = SYSTEM_true;
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),&b,1));
  if (ValueCast(SYSTEM_int32,b) > self->
    GXFILE_tgxfileobj_DOT_deltaforread) {
    if (b == 255) {
      result = SYSTEM_false;
      return result;
    } 
    *afdim = self->GXFILE_tgxfileobj_DOT_fcurrentdim;
    if (self->GXFILE_tgxfileobj_DOT_fcurrentdim > 0) 
      _P3inc1(self->GXFILE_tgxfileobj_DOT_lastelem[self->
        GXFILE_tgxfileobj_DOT_fcurrentdim - 1],b - self->
        GXFILE_tgxfileobj_DOT_deltaforread);
  } else {
    *afdim = b;
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = *afdim) <=  _stop) do {
        switch (self->GXFILE_tgxfileobj_DOT_elemtype[d - 1]) {
          case GXFILE_sz_integer: 
            self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = VirtMethodCall(ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
              GMSSTRM_txstream_DOT_readinteger_T, 7, (ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile))) + 
              self->GXFILE_tgxfileobj_DOT_minelem[d - 1];
            break;
          case GXFILE_sz_word: 
            self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = VirtMethodCall(ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile), 
              GMSSTRM_txstream_DOT_readword_T, 8, (ValueCast(
              GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile))) + 
              self->GXFILE_tgxfileobj_DOT_minelem[d - 1];
            break;
          case GXFILE_sz_byte: 
            self->GXFILE_tgxfileobj_DOT_lastelem[d - 1] = 
              GMSSTRM_txstream_DOT_readbyte(ValueCast(GMSSTRM_txstream,
              self->GXFILE_tgxfileobj_DOT_ffile)) + self->
              GXFILE_tgxfileobj_DOT_minelem[d - 1];
            break;
          default: break;
        }
      } while (d++ !=  _stop);

    }
  } 
  if (self->GXFILE_tgxfileobj_DOT_datasize > 0) 
    { register GMSSPECS_tvarvaltype _stop = self->
        GXFILE_tgxfileobj_DOT_lastdatafield;
      if ((dv = GMSSPECS_vallevel) <=  _stop) do {
        VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
          GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_read_T, 4, (ValueCast(
          GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),&bsv,1));
        sv = ValueCast(GXFILE_tgdxintvaltyp,bsv);
        if (sv != GXFILE_vm_normal) { 
          avals[dv] = self->GXFILE_tgxfileobj_DOT_readintvaluemap[sv];
        } else {
          avals[dv] = VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
            GXFILE_tgxfileobj_DOT_ffile), 
            GMSSTRM_txstream_DOT_readdouble_T, 6, (ValueCast(
            GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile)));
          if (avals[dv] >= self->GXFILE_tgxfileobj_DOT_zvalacr) 
            avals[dv] = GXFILE_tgxfileobj_DOT_acronymremap(self,avals[
              dv]);
        } 
      
      } while (dv++ !=  _stop);

    }
  return result;
}  /* doread */

Procedure GXFILE_tgxfileobj_DOT_addtoerrorlist(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *aelements,
  const SYSTEM_double *avals)
{
  if (self->GXFILE_tgxfileobj_DOT_errorlist == NULL) { 
    self->GXFILE_tgxfileobj_DOT_errorlist = ValueCast(
      GMSDATA_ttblgamsdata,GMSDATA_ttblgamsdata_DOT_create(ValueCast(
      GMSDATA_ttblgamsdata,_P3alloc_object(&GMSDATA_ttblgamsdata_CD)),
      self->GXFILE_tgxfileobj_DOT_fcurrentdim,self->
      GXFILE_tgxfileobj_DOT_datasize * sizeof(SYSTEM_double)));
  } else 
    if (GMSDATA_ttblgamsdata_DOT_getcount(self->
      GXFILE_tgxfileobj_DOT_errorlist) >= 11) 
      return;
  GMSDATA_ttblgamsdata_DOT_addrecord(self->
    GXFILE_tgxfileobj_DOT_errorlist,aelements,avals);
}  /* addtoerrorlist */

Procedure GXFILE_tgxfileobj_DOT_getdefaultrecord(
  GXFILE_tgxfileobj self,
  SYSTEM_double *avals)
{
  SYSTEM_integer ui;

  switch (self->GXFILE_tgxfileobj_DOT_cursyptr->sdatatype) {
    case GMSSPECS_dt_set: 
    case GMSSPECS_dt_alias: 
      avals[GMSSPECS_vallevel] = 0.0;
      break;
    case GMSSPECS_dt_par: 
      avals[GMSSPECS_vallevel] = 0.0;
      break;
    case GMSSPECS_dt_var: 
      ui = self->GXFILE_tgxfileobj_DOT_cursyptr->suserinfo;
      if (ui >= 0 && ui <= 9) { 
        _P3memcpy(avals,sizeof(GMSSPECS_tvarreca),GMSGLOB_defrecvar[ValueCast(
          GMSSPECS_tvarstyp,ui)]);
      } else 
        _P3memcpy(avals,sizeof(GMSSPECS_tvarreca),GMSGLOB_defrecvar[
          GMSSPECS_stypunknwn]);
      break;
    case GMSSPECS_dt_equ: 
      ui = self->GXFILE_tgxfileobj_DOT_cursyptr->suserinfo;
      if (ui >= 0 && ui <= 85) { 
        _P3memcpy(avals,sizeof(GMSSPECS_tvarreca),GMSGLOB_defrecequ[ValueCast(
          GMSGLOB_tssymbol,ui) - 53]);
      } else 
        _P3memcpy(avals,sizeof(GMSSPECS_tvarreca),GMSGLOB_defrecequ[0]);
      break;
    default:
      SYSTEM_assert(SYSTEM_false,_P3str1("\022GetDefaultRecord-2"));
  }
}  /* getdefaultrecord */

static Function(SYSTEM_double ) getasacronym(
  SYSTEM_double v,
  GXFILE_tgxfileobj *_2self)
{
  SYSTEM_double result;
  SYSTEM_integer n;
  SYSTEM_integer newindx;
  SYSTEM_integer orgindx;

  orgindx = SYSTEM_round(v /  (*_2self)->GXFILE_tgxfileobj_DOT_zvalacr);
  n = GXFILE_tacronymlist_DOT_findentry((*_2self)->
    GXFILE_tgxfileobj_DOT_acronymlist,orgindx);
  if (n < 0) { 
    if ((*_2self)->GXFILE_tgxfileobj_DOT_nextautoacronym <= 0) { 
      newindx = orgindx;
    } else {
      newindx = (*_2self)->GXFILE_tgxfileobj_DOT_nextautoacronym;
      (*_2self)->GXFILE_tgxfileobj_DOT_nextautoacronym = (*_2self)->
        GXFILE_tgxfileobj_DOT_nextautoacronym + 1;
      n = GXFILE_tacronymlist_DOT_addentry((*_2self)->
        GXFILE_tgxfileobj_DOT_acronymlist,_P3str1("\000"),_P3str1("\000"),
        orgindx);
      { register GXFILE_tacronym_OD *_W3=ValueCast(GXFILE_tacronym,
        GMSOBJ_txlist_DOT_get((*_2self)->
        GXFILE_tgxfileobj_DOT_acronymlist->
        GXFILE_tacronymlist_DOT_flist,n));
        _W3->GXFILE_tacronym_DOT_acrreadmap = newindx;
        _W3->GXFILE_tacronym_DOT_acrautogen = SYSTEM_true;

      }
    } 
  } else {
    newindx = (ValueCast(GXFILE_tacronym,GMSOBJ_txlist_DOT_get((*
      _2self)->GXFILE_tgxfileobj_DOT_acronymlist->
      GXFILE_tacronymlist_DOT_flist,n)))->
      GXFILE_tacronym_DOT_acrreadmap;
    if (newindx <= 0) 
      if ((*_2self)->GXFILE_tgxfileobj_DOT_nextautoacronym <= 0) { 
        newindx = orgindx;
      } else {
        newindx = (*_2self)->GXFILE_tgxfileobj_DOT_nextautoacronym;
        (*_2self)->GXFILE_tgxfileobj_DOT_nextautoacronym = (*_2self)->
          GXFILE_tgxfileobj_DOT_nextautoacronym + 1;
        { register GXFILE_tacronym_OD *_W3=ValueCast(GXFILE_tacronym,
          GMSOBJ_txlist_DOT_get((*_2self)->
          GXFILE_tgxfileobj_DOT_acronymlist->
          GXFILE_tacronymlist_DOT_flist,n));
          _W3->GXFILE_tacronym_DOT_acrreadmap = newindx;
          _W3->GXFILE_tacronym_DOT_acrautogen = SYSTEM_true;

        }
      } 
  } 
  result = (*_2self)->GXFILE_tgxfileobj_DOT_zvalacr * newindx;
  return result;
}  /* getasacronym */

Function(SYSTEM_double ) GXFILE_tgxfileobj_DOT_acronymremap(
  GXFILE_tgxfileobj self,
  SYSTEM_double v)
{
  SYSTEM_double result;

  if (v < self->GXFILE_tgxfileobj_DOT_zvalacr) { 
    result = v;
  } else 
    switch (P3IEEEFP_fpclass(v)) {
      case P3IEEEFP_fp_snan: 
      case P3IEEEFP_fp_qnan: 
      case P3IEEEFP_fp_ndenorm: 
      case P3IEEEFP_fp_pdenorm: 
        result = self->GXFILE_tgxfileobj_DOT_intvaluemap[
          GXFILE_vm_valna];
        break;
      case P3IEEEFP_fp_ninf: 
        result = self->GXFILE_tgxfileobj_DOT_intvaluemap[
          GXFILE_vm_valmin];
        break;
      case P3IEEEFP_fp_pinf: 
        result = self->GXFILE_tgxfileobj_DOT_intvaluemap[
          GXFILE_vm_valpin];
        break;
      case P3IEEEFP_fp_nzero: 
      case P3IEEEFP_fp_pzero: 
        result = 0.0;
        break;
      case P3IEEEFP_fp_nnorm: 
        result = v;
        break;
      case P3IEEEFP_fp_pnorm: 
        result = getasacronym(v,&self);
        break;
      default:
        result = self->GXFILE_tgxfileobj_DOT_intvaluemap[
          GXFILE_vm_valna];
    }
  return result;
}  /* acronymremap */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_isgoodnewsymbol(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_boolean result;

  result = SYSTEM_false;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,
    STRHASH_txstrhashlist_DOT_indexof(self->
    GXFILE_tgxfileobj_DOT_namelist,s) < 1,
    GXFILE_err_duplicatesymbol)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,
    GXFILE_tacronymlist_DOT_findname(self->
    GXFILE_tgxfileobj_DOT_acronymlist,s) < 0,
    GXFILE_err_duplicatesymbol)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,GXFILE_isgoodident(s),
    GXFILE_err_badidentformat)) 
    return result;
  result = SYSTEM_true;
  return result;
}  /* isgoodnewsymbol */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbmaxlength(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  SYSTEM_integer n;
  SYSTEM_integer l;

  result = 0;
  { register SYSTEM_int32 _stop = self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount;
    if ((n = 1) <=  _stop) do {
      l = STRHASH_txstrhashlist_DOT_getstringlength(self->
        GXFILE_tgxfileobj_DOT_namelist,n);
      if (l > result) 
        result = l;
    
    } while (n++ !=  _stop);

  }
  return result;
}  /* gdxsymbmaxlength */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelmaxlength(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;
  SYSTEM_integer n;
  SYSTEM_integer l;

  result = 0;
  { register SYSTEM_int32 _stop = self->GXFILE_tgxfileobj_DOT_ueltable->
      STRHASH_txstrhashlist_DOT_fcount;
    if ((n = 1) <=  _stop) do {
      l = STRHASH_txstrhashlist_DOT_getstringlength(ValueCast(
        STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),n);
      if (l > result) 
        result = l;
    
    } while (n++ !=  _stop);

  }
  return result;
}  /* gdxuelmaxlength */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbindxmaxlength(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *lengthinfo)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_22 = {2,0,0};
  SYSTEM_integer d;
  SYSTEM_integer nrrecs;
  GXDEFS_tgdxvalues avals;
  SYSTEM_integer afdim;
  SYSTEM_integer uel;
  SYSTEM_integer ueltablecount;
  SYSTEM_integer l;

  result = 0;
  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    lengthinfo[d - 1] = 0;
  }
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_some || !
    _P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,allowedmodes_22)) 
    if (!GXFILE_tgxfileobj_DOT_checkmode(self,_P3str1("\021SymbIndxMaxLength"),
      allowedmodes_22)) 
      return result;
  if (!(synr >= 0 && synr <= self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount)) 
    return result;
  if (GXFILE_tgxfileobj_DOT_gdxdatareadrawstart(self,synr,&nrrecs) == 0) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_fcurrentdim > 0) {
    ueltablecount = self->GXFILE_tgxfileobj_DOT_ueltable->
      STRHASH_txstrhashlist_DOT_fcount;
    while (GXFILE_tgxfileobj_DOT_doread(self,avals,&afdim)) {

      { register SYSTEM_int32 _stop = self->
          GXFILE_tgxfileobj_DOT_fcurrentdim;
        if ((d = afdim) <=  _stop) do {
          uel = self->GXFILE_tgxfileobj_DOT_lastelem[d - 1];
          if (uel >= 1 && uel <= ueltablecount) {
            l = STRHASH_txstrhashlist_DOT_getstringlength(ValueCast(
              STRHASH_txstrhashlist,self->
              GXFILE_tgxfileobj_DOT_ueltable),uel);
            if (l > lengthinfo[d - 1]) 
              lengthinfo[d - 1] = l;
          } 
        
        } while (d++ !=  _stop);

      }
}
    { register SYSTEM_int32 _stop = self->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((d = 1) <=  _stop) do {
        if (lengthinfo[d - 1] > result) 
          result = lengthinfo[d - 1];
      } while (d++ !=  _stop);

    }
  } 
  GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  return result;
}  /* gdxsymbindxmaxlength */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymindex(
  GXFILE_tgxfileobj self,
  SYSTEM_double v)
{
  SYSTEM_integer result;

  if (v < self->GXFILE_tgxfileobj_DOT_zvalacr) { 
    result = 0;
  } else 
    result = SYSTEM_round(v /  self->GXFILE_tgxfileobj_DOT_zvalacr);
  return result;
}  /* gdxacronymindex */

Function(SYSTEM_double ) GXFILE_tgxfileobj_DOT_gdxacronymvalue(
  GXFILE_tgxfileobj self,
  SYSTEM_integer aindx)
{
  SYSTEM_double result;

  if (aindx <= 0) { 
    result = 0.0;
  } else 
    result = self->GXFILE_tgxfileobj_DOT_zvalacr * aindx;
  return result;
}  /* gdxacronymvalue */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymname(
  GXFILE_tgxfileobj self,
  SYSTEM_double v,
  SYSTEM_ansichar *aname)
{
  SYSTEM_integer result;
  SYSTEM_integer n;
  SYSTEM_integer indx;

  result = GXFILE_myfalse;
  indx = GXFILE_tgxfileobj_DOT_gdxacronymindex(self,v);
  if (indx <= 0) { 
    _P3strclr(aname);
  } else {
    n = GXFILE_tacronymlist_DOT_findentry(self->
      GXFILE_tgxfileobj_DOT_acronymlist,indx);
    if (n < 0) { 
      {
        SYSTEM_shortstring _t1;

        _P3strcat(aname,255,_P3str1("\016UnknownAcronym"),
          STRUTILX_inttostr(_t1,255,indx));
      }
    } else 
      STRUTILX_getstring(aname,255,(ValueCast(GXFILE_tacronym,
        GMSOBJ_txlist_DOT_get(self->GXFILE_tgxfileobj_DOT_acronymlist->
        GXFILE_tacronymlist_DOT_flist,n)))->
        GXFILE_tacronym_DOT_acrname);
    result = GXFILE_mytrue;
  } 
  return result;
}  /* gdxacronymname */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymgetinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  SYSTEM_ansichar *aname,
  SYSTEM_ansichar *txt,
  SYSTEM_integer *aindx)
{
  SYSTEM_integer result;

  if (!(n >= 1 && n <= self->GXFILE_tgxfileobj_DOT_acronymlist->
    GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount)) {
    result = GXFILE_myfalse;
    _P3strclr(aname);
    _P3strclr(txt);
    *aindx = 0;
  } else {
    result = GXFILE_mytrue;
    { register GXFILE_tacronym_OD *_W2=ValueCast(GXFILE_tacronym,
      GMSOBJ_txlist_DOT_get(self->GXFILE_tgxfileobj_DOT_acronymlist->
      GXFILE_tacronymlist_DOT_flist,n - 1));
      STRUTILX_getstring(aname,255,_W2->GXFILE_tacronym_DOT_acrname);
      STRUTILX_getstring(txt,255,_W2->GXFILE_tacronym_DOT_acrtext);
      *aindx = _W2->GXFILE_tacronym_DOT_acrmap;

    }
  } 
  return result;
}  /* gdxacronymgetinfo */

static Function(SYSTEM_boolean ) mapisunique(
  SYSTEM_integer indx,
  GXFILE_tgxfileobj *_2self)
{
  SYSTEM_boolean result;
  SYSTEM_integer n;

  result = SYSTEM_false;
  { register SYSTEM_int32 _stop = (*_2self)->
      GXFILE_tgxfileobj_DOT_acronymlist->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if ((ValueCast(GXFILE_tacronym,GMSOBJ_txlist_DOT_get((*_2self)->
        GXFILE_tgxfileobj_DOT_acronymlist->
        GXFILE_tacronymlist_DOT_flist,n)))->
        GXFILE_tacronym_DOT_acrreadmap == indx) 
        return result;
    } while (n++ !=  _stop);

  }
  result = SYSTEM_true;
  return result;
}  /* mapisunique */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymsetinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  const SYSTEM_ansichar *aname,
  const SYSTEM_ansichar *txt,
  SYSTEM_integer aindx)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_some) 
    {
      _P3STR_255 _t1;
      _P3STR_255 _t2;
      SYSTEM_shortstring _t3;
      _P3STR_255 _t4;

      GXFILE_tgxfileobj_DOT_writetrace(self,_P3strcat(_t4,255,
        _P3strcat(_t2,255,_P3strcat(_t1,255,_P3str1("\020AcronymSetInfo: "),
        aname),_P3str1("\011 index = ")),STRUTILX_inttostr(_t3,255,
        aindx)));
    }
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,n >= 1 || n <= 
    self->GXFILE_tgxfileobj_DOT_acronymlist->
    GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount,
    GXFILE_err_badacronumber)) 
    return result;
  { register GXFILE_tacronym_OD *_W2=ValueCast(GXFILE_tacronym,
    GMSOBJ_txlist_DOT_get(self->GXFILE_tgxfileobj_DOT_acronymlist->
    GXFILE_tacronymlist_DOT_flist,n - 1));
    if (_P3SET_i(17,self->GXFILE_tgxfileobj_DOT_fmode,
      GXFILE_anywritemode) || _W2->GXFILE_tacronym_DOT_acrautogen) {
      if (GXFILE_tgxfileobj_DOT_errorcondition(self,
        GXFILE_tgxfileobj_DOT_isgoodnewsymbol(self,aname),
        GXFILE_err_badacroname)) 
        return result;
      if (_W2->GXFILE_tacronym_DOT_acrautogen) {
        SYSTEM_assert(_W2->GXFILE_tacronym_DOT_acrreadmap == aindx,_P3str1("\021gdxAcronymSetInfo"));
        _W2->GXFILE_tacronym_DOT_acrautogen = SYSTEM_false;
      } else 
        if (GXFILE_tgxfileobj_DOT_errorcondition(self,aindx == _W2->
          GXFILE_tacronym_DOT_acrmap,GXFILE_err_badacroindex)) 
          return result;
      STRUTILX_strassignm(&_W2->GXFILE_tacronym_DOT_acrname,aname,&_W2->
        GXFILE_tacronym_DOT_fstrmemory);
      {
        SYSTEM_shortstring _t1;

        STRUTILX_strassignm(&_W2->GXFILE_tacronym_DOT_acrtext,
          GXFILE_makegoodexpltext(_t1,255,txt),&_W2->
          GXFILE_tacronym_DOT_fstrmemory);
      }
    } else 
      if (_W2->GXFILE_tacronym_DOT_acrreadmap != aindx) {
        if (GXFILE_tgxfileobj_DOT_errorcondition(self,
          STRUTILX_pstruequal(ValueCast(SYSTEM_P3_pshortstring,aname),
          _W2->GXFILE_tacronym_DOT_acrname),GXFILE_err_badacroname)) 
          return result;
        if (GXFILE_tgxfileobj_DOT_errorcondition(self,mapisunique(
          aindx,&self),GXFILE_err_acrodupemap)) 
          return result;
        _W2->GXFILE_tacronym_DOT_acrreadmap = aindx;
      } 

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxacronymsetinfo */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymcount(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;

  if (self->GXFILE_tgxfileobj_DOT_acronymlist == NULL) { 
    result = 0;
  } else 
    result = self->GXFILE_tgxfileobj_DOT_acronymlist->
      GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount;
  return result;
}  /* gdxacronymcount */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymnextnr(
  GXFILE_tgxfileobj self,
  SYSTEM_integer nv)
{
  SYSTEM_integer result;

  result = self->GXFILE_tgxfileobj_DOT_nextautoacronym;
  if (nv >= 0) 
    self->GXFILE_tgxfileobj_DOT_nextautoacronym = nv;
  return result;
}  /* gdxacronymnextnr */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymgetmapping(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  SYSTEM_integer *orgindx,
  SYSTEM_integer *newindx,
  SYSTEM_integer *autoindex)
{
  SYSTEM_integer result;

  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_tracelevel >= GXFILE_trl_some) 
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;

      GXFILE_tgxfileobj_DOT_writetrace(self,_P3strcat(_t2,255,_P3str1("\027AcronymGetMapping: N = "),
        STRUTILX_inttostr(_t1,255,n)));
    }
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,n >= 1 || n <= 
    self->GXFILE_tgxfileobj_DOT_acronymlist->
    GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount,
    GXFILE_err_badacronumber)) 
    return result;
  { register GXFILE_tacronym_OD *_W2=ValueCast(GXFILE_tacronym,
    GMSOBJ_txlist_DOT_get(self->GXFILE_tgxfileobj_DOT_acronymlist->
    GXFILE_tacronymlist_DOT_flist,n - 1));
    *orgindx = _W2->GXFILE_tacronym_DOT_acrmap;
    *newindx = _W2->GXFILE_tacronym_DOT_acrreadmap;
    if (_W2->GXFILE_tacronym_DOT_acrautogen) { 
      *autoindex = GXFILE_mytrue;
    } else 
      *autoindex = GXFILE_myfalse;

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxacronymgetmapping */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymadd(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *aname,
  const SYSTEM_ansichar *txt,
  SYSTEM_integer aindx)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  result =  -1;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_acronymlist->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      { register GXFILE_tacronym_OD *_W2=ValueCast(GXFILE_tacronym,
        GMSOBJ_txlist_DOT_get(self->GXFILE_tgxfileobj_DOT_acronymlist->
        GXFILE_tacronymlist_DOT_flist,n));
        if (STRUTILX_pstruequal(_W2->GXFILE_tacronym_DOT_acrname,ValueCast(
          SYSTEM_P3_pshortstring,aname))) {
          if (GXFILE_tgxfileobj_DOT_errorcondition(self,_W2->
            GXFILE_tacronym_DOT_acrmap == aindx,
            GXFILE_err_acrobadaddition)) 
            return result;
          result = n;
          return result;
        } 
        if (GXFILE_tgxfileobj_DOT_errorcondition(self,_W2->
          GXFILE_tacronym_DOT_acrmap != aindx,
          GXFILE_err_acrobadaddition)) 
          return result;

      }
    } while (n++ !=  _stop);

  }
  result = GXFILE_tacronymlist_DOT_addentry(self->
    GXFILE_tgxfileobj_DOT_acronymlist,aname,txt,aindx);
  (ValueCast(GXFILE_tacronym,GMSOBJ_txlist_DOT_get(self->
    GXFILE_tgxfileobj_DOT_acronymlist->GXFILE_tacronymlist_DOT_flist,
    result)))->GXFILE_tacronym_DOT_acrreadmap = aindx;
  result = result + 1;
  return result;
}  /* gdxacronymadd */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymboladdcomment(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_ansichar *txt)
{
  SYSTEM_integer result;
  GXFILE_pgdxsymbrecord syptr;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\020SymbolAddComment"),
    GXFILE_anywritemode)) 
    return result;
  if (synr <= 0) { 
    syptr = self->GXFILE_tgxfileobj_DOT_cursyptr;
  } else 
    if (self->GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && 
      synr <= self->GXFILE_tgxfileobj_DOT_namelist->
      STRHASH_txstrhashlist_DOT_fcount) { 
      syptr = ValueCast(GXFILE_pgdxsymbrecord,
        STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,synr));
    } else 
      syptr = NULL;
  if (syptr == NULL) {
    GXFILE_tgxfileobj_DOT_reporterror(self,
      GXFILE_err_nosymbolforcomment);
    return result;
  } 
  { register GXFILE_tgdxsymbrecord *_W2=syptr;
    if (_W2->scommentslist == NULL) 
      _W2->scommentslist = ValueCast(GMSOBJ_txstrings,
        SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,
        _P3alloc_object(&GMSOBJ_txstrings_CD))));
    GMSOBJ_txstrings_DOT_add(_W2->scommentslist,txt);

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxsymboladdcomment */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolgetcomment(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer n,
  SYSTEM_ansichar *txt)
{
  SYSTEM_integer result;

  _P3strclr(txt);
  result = GXFILE_myfalse;
  if (self->GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && 
    synr <= self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount) 
    { register GXFILE_tgdxsymbrecord *_W2=ValueCast(
      GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
      GXFILE_tgxfileobj_DOT_namelist,synr));
      if (_W2->scommentslist != NULL) 
        if (n >= 1 && n <= _W2->scommentslist->
          GMSOBJ_txlist_DOT_fcount) {
          result = GXFILE_mytrue;
          GMSOBJ_txstrings_DOT_get(txt,255,_W2->scommentslist,n - 1);
        } 

    }
  return result;
}  /* gdxsymbolgetcomment */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolsetdomain(
  GXFILE_tgxfileobj self,
  const SYSTEM_shortstring *domainids)
{
  SYSTEM_integer result;
  static GXFILE_tgxmodeset allowedmodes_23 = {56,0,0};
  SYSTEM_integer d;
  SYSTEM_integer domsy;
  SYSTEM_integer synr;
  SYSTEM_boolean domap;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\017SymbolSetDomain"),
    allowedmodes_23)) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_cursyptr == NULL) 
    return result;
  result = GXFILE_mytrue;
  { register GXFILE_tgdxsymbrecord *_W2=self->
    GXFILE_tgxfileobj_DOT_cursyptr;
    SYSTEM_assert(_W2->sdomsymbols == NULL,_P3str1("\017SymbolSetDomain"));
    _P3getmem(_W2->sdomsymbols,(_W2->sdim + 1) * sizeof(
      SYSTEM_longint));
    { register SYSTEM_int32 _stop = _W2->sdim;
      if ((d = 1) <=  _stop) do {
        domap = SYSTEM_true;
        if (_P3stccmpE(domainids[d - 1],_P3char('*'))) { 
          domsy = 0;
        } else {
          domsy = STRHASH_txstrhashlist_DOT_indexof(self->
            GXFILE_tgxfileobj_DOT_namelist,domainids[d - 1]);
          if (domsy <= 0) {
            GXFILE_tgxfileobj_DOT_reporterror(self,
              GXFILE_err_unknowndomain);
            result = GXFILE_myfalse;
            domsy = 0;
          } 
        } 
        if (domsy > 0) {
          synr = domsy;
          do {
            { register GXFILE_tgdxsymbrecord *_W3=ValueCast(
              GXFILE_pgdxsymbrecord,
              STRHASH_txstrhashlist_DOT_getobject(self->
              GXFILE_tgxfileobj_DOT_namelist,synr));
              if (_W3->sdatatype == GMSSPECS_dt_set) 
                SYSTEM_break(BRK_11);
              if (_W3->sdatatype == GMSSPECS_dt_alias) {
                synr = _W3->suserinfo;
                if (synr > 0) 
                  SYSTEM_continue(CNT_11);
                domap = SYSTEM_false;
                SYSTEM_break(BRK_11);
              } 
              GXFILE_tgxfileobj_DOT_reporterror(self,
                GXFILE_err_aliassetexpected);
              result = GXFILE_myfalse;
              domsy = 0;
              SYSTEM_break(BRK_11);

            }
          CNT_11:;
          } while (SYSTEM_true);
BRK_11:;
        } 
        (*_W2->sdomsymbols)[d] = domsy;
        if (domap && domsy > 0) 
          if (!(_W2->sdim == 1 && self->
            GXFILE_tgxfileobj_DOT_cursyptr == ValueCast(
            GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(
            self->GXFILE_tgxfileobj_DOT_namelist,domsy)))) 
            self->GXFILE_tgxfileobj_DOT_wrbitmaps[d - 1] = (ValueCast(
              GXFILE_pgdxsymbrecord,
              STRHASH_txstrhashlist_DOT_getobject(self->
              GXFILE_tgxfileobj_DOT_namelist,synr)))->ssetbitmap;
      
      } while (d++ !=  _stop);

    }

  }
  switch (self->GXFILE_tgxfileobj_DOT_fmode) {
    case GXFILE_fw_dom_raw: 
      self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_raw_data;
      break;
    case GXFILE_fw_dom_map: 
      self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_map_data;
      break;
    case GXFILE_fw_dom_str: 
      self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_str_data;
      break;
    default: break;
  }
  return result;
}  /* gdxsymbolsetdomain */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolsetdomainx(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_shortstring *domainids)
{
  SYSTEM_integer result;
  GXFILE_pgdxsymbrecord syptr;
  SYSTEM_integer d;
  SYSTEM_shortstring s;

  result = GXFILE_myfalse;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,synr >= 1 && synr <= 
    self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount,GXFILE_err_badsymbolindex)) 
    return result;
  syptr = ValueCast(GXFILE_pgdxsymbrecord,
    STRHASH_txstrhashlist_DOT_getobject(self->
    GXFILE_tgxfileobj_DOT_namelist,synr));
  { register GXFILE_tgdxsymbrecord *_W2=syptr;
    if (_W2->sdim > 0) {
      if (_W2->sdomstrings == NULL) 
        _P3getmem(_W2->sdomstrings,(_W2->sdim + 1) * sizeof(
          SYSTEM_longint));
      { register SYSTEM_int32 _stop = _W2->sdim;
        if ((d = 1) <=  _stop) do {
          _P3strcpy(s,255,domainids[d - 1]);
          if (_P3strcmpE(s,_P3str1("\000")) || _P3stccmpE(s,_P3char('*'))) { 
            (*_W2->sdomstrings)[d] = 0;
          } else 
            if (!GXFILE_isgoodident(s)) { 
              (*_W2->sdomstrings)[d] = 0;
            } else {
              (*_W2->sdomstrings)[d] = 
                STRHASH_txstrhashlist_DOT_indexof(self->
                GXFILE_tgxfileobj_DOT_domainstrlist,s);
              if ((*_W2->sdomstrings)[d] <= 0) 
                (*_W2->sdomstrings)[d] = STRHASH_txstrhashlist_DOT_add(
                  self->GXFILE_tgxfileobj_DOT_domainstrlist,s);
            } 
        
        } while (d++ !=  _stop);

      }
    } 

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxsymbolsetdomainx */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolgetdomain(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *domainsynrs)
{
  SYSTEM_integer result;
  GXFILE_pgdxsymbrecord syptr;
  SYSTEM_integer d;

  result = GXFILE_myfalse;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,synr >= 1 && synr <= 
    self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount,GXFILE_err_badsymbolindex)) 
    return result;
  syptr = ValueCast(GXFILE_pgdxsymbrecord,
    STRHASH_txstrhashlist_DOT_getobject(self->
    GXFILE_tgxfileobj_DOT_namelist,synr));
  { register GXFILE_tgdxsymbrecord *_W2=syptr;
    { register SYSTEM_int32 _stop = _W2->sdim;
      if ((d = 1) <=  _stop) do {
        if (_W2->sdomsymbols == NULL || (*_W2->sdomsymbols)[d] == 0) { 
          domainsynrs[d - 1] = 0;
        } else 
          domainsynrs[d - 1] = (*_W2->sdomsymbols)[d];
      } while (d++ !=  _stop);

    }

  }
  result = GXFILE_mytrue;
  return result;
}  /* gdxsymbolgetdomain */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolgetdomainx(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_shortstring *domainids)
{
  SYSTEM_integer result;
  GXFILE_pgdxsymbrecord syptr;
  SYSTEM_integer d;

  result = 0;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,self->
    GXFILE_tgxfileobj_DOT_namelist != NULL && synr >= 1 && synr <= 
    self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount,GXFILE_err_badsymbolindex)) 
    return result;
  syptr = ValueCast(GXFILE_pgdxsymbrecord,
    STRHASH_txstrhashlist_DOT_getobject(self->
    GXFILE_tgxfileobj_DOT_namelist,synr));
  { register GXFILE_tgdxsymbrecord *_W2=syptr;
    { register SYSTEM_int32 _stop = _W2->sdim;
      if ((d = 1) <=  _stop) do {
        _P3strcpy(domainids[d - 1],255,_P3str1("\001*"));
      } while (d++ !=  _stop);

    }
    if (_W2->sdomstrings != NULL) {
      { register SYSTEM_int32 _stop = _W2->sdim;
        if ((d = 1) <=  _stop) do {
          if ((*_W2->sdomstrings)[d] != 0) 
            STRHASH_txstrhashlist_DOT_getstring(domainids[d - 1],255,
              self->GXFILE_tgxfileobj_DOT_domainstrlist,(*_W2->
              sdomstrings)[d]);
        } while (d++ !=  _stop);

      }
      result = 2;
    } else 
      if (_W2->sdomsymbols == NULL) { 
        result = 1;
      } else {
        { register SYSTEM_int32 _stop = _W2->sdim;
          if ((d = 1) <=  _stop) do {
            if ((*_W2->sdomsymbols)[d] != 0) 
              STRHASH_txstrhashlist_DOT_getstring(domainids[d - 1],255,
                self->GXFILE_tgxfileobj_DOT_namelist,(*_W2->
                sdomsymbols)[d]);
          } while (d++ !=  _stop);

        }
        result = 3;
      } 

  }
  return result;
}  /* gdxsymbolgetdomainx */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxaddalias(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *id1,
  const SYSTEM_ansichar *id2)
{
  SYSTEM_integer result;
  GXFILE_pgdxsymbrecord syptr;
  SYSTEM_integer synr1, synr2, synr;
  SYSTEM_shortstring aname;

  result = GXFILE_myfalse;
  if (!GXFILE_tgxfileobj_DOT_majorcheckmode(self,_P3str1("\010AddAlias"),
    GXFILE_anywritemode)) 
    return result;
  if (_P3stccmpE(id1,_P3char('*'))) { 
    synr1 = SYSTEM_maxint;
  } else 
    synr1 = STRHASH_txstrhashlist_DOT_indexof(self->
      GXFILE_tgxfileobj_DOT_namelist,id1);
  if (_P3stccmpE(id2,_P3char('*'))) { 
    synr2 = SYSTEM_maxint;
  } else 
    synr2 = STRHASH_txstrhashlist_DOT_indexof(self->
      GXFILE_tgxfileobj_DOT_namelist,id2);
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,synr1 >= 0 != 
    synr2 >= 0,GXFILE_err_aliassetexpected)) 
    return result;
  if (synr1 > 0) {
    synr = synr1;
    _P3strcpy(aname,255,id2);
  } else {
    synr = synr2;
    _P3strcpy(aname,255,id1);
  } 
  if (synr == SYSTEM_maxint) { 
    synr = 0;
  } else 
    { register GXFILE_tgdxsymbrecord *_W2=ValueCast(
      GXFILE_pgdxsymbrecord,STRHASH_txstrhashlist_DOT_getobject(self->
      GXFILE_tgxfileobj_DOT_namelist,synr));
      if (GXFILE_tgxfileobj_DOT_errorcondition(self,_P3SET_in_1(_W2->
        sdatatype,GMSSPECS_dt_set,_P3SET_equal(_W2->sdatatype,
        GMSSPECS_dt_alias)),GXFILE_err_aliassetexpected)) 
        return result;

    }
  if (!GXFILE_tgxfileobj_DOT_isgoodnewsymbol(self,aname)) 
    return result;
  _P3new(syptr);
  SYSTEM_P3_fillchar(syptr,sizeof(GXFILE_tgdxsymbrecord),0);
  { register GXFILE_tgdxsymbrecord *_W2=syptr;
    _W2->sdatatype = GMSSPECS_dt_alias;
    _W2->suserinfo = synr;
    if (synr == 0) {
      _W2->sdim = 1;
      _P3strcpy(_W2->sexpltxt,255,_P3str1("\016Aliased with *"));
    } else {
      _W2->sdim = (ValueCast(GXFILE_pgdxsymbrecord,
        STRHASH_txstrhashlist_DOT_getobject(self->
        GXFILE_tgxfileobj_DOT_namelist,synr)))->sdim;
      {
        SYSTEM_shortstring _t1;

        _P3strcat(_W2->sexpltxt,255,_P3str1("\015Aliased with "),
          STRHASH_txstrhashlist_DOT_getstring(_t1,255,self->
          GXFILE_tgxfileobj_DOT_namelist,synr));
      }
    } 

  }
  STRHASH_txstrhashlist_DOT_addobject(self->
    GXFILE_tgxfileobj_DOT_namelist,aname,ValueCast(SYSTEM_tobject,
    syptr));
  result = GXFILE_mytrue;
  return result;
}  /* gdxaddalias */

Constructor(GXFILE_tgxfileobj ) GXFILE_tgxfileobj_DOT_create(
  GXFILE_tgxfileobj self,
  SYSTEM_ansichar *errmsg)
{
  _P3strclr(errmsg);
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_f_not_open;
  self->GXFILE_tgxfileobj_DOT_fstatus = GXFILE_stat_notopen;
  self->GXFILE_tgxfileobj_DOT_autoconvert = 1;
  GXFILE_tgxfileobj_DOT_gdxresetspecialvalues(self);
  self->GXFILE_tgxfileobj_DOT_nextautoacronym = 0;
  self->GXFILE_tgxfileobj_DOT_appendactive = SYSTEM_false;
  return self;
}  /* create */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetdllversion(
  GXFILE_tgxfileobj self,
  SYSTEM_ansichar *v)
{
  SYSTEM_integer result;

  result = GXFILE_mytrue;
  GDLAUDIT_gdlgetauditline(v,255);
  return result;
}  /* gdxgetdllversion */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfileinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *filever,
  SYSTEM_integer *comprlev)
{
  SYSTEM_integer result;

  result = GXFILE_mytrue;
  switch (self->GXFILE_tgxfileobj_DOT_fstatus) {
    case GXFILE_stat_notopen: 
      *filever = 0;
      *comprlev = 0;
      break;
    case GXFILE_stat_read: 
      *filever = self->GXFILE_tgxfileobj_DOT_versionread;
      *comprlev = self->GXFILE_tgxfileobj_DOT_fcomprlev;
      break;
    case GXFILE_stat_write: 
      *filever = GXFILE_version;
      *comprlev = self->GXFILE_tgxfileobj_DOT_fcomprlev;
      break;
    default: break;
  }
  return result;
}  /* gdxfileinfo */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxerrorstr(
  GXFILE_tgxfileobj self,
  SYSTEM_integer errnr,
  SYSTEM_ansichar *errmsg)
{
  SYSTEM_integer result;

  result = GXFILE_mytrue;
  switch (errnr) {
    case GXFILE_err_nofile: 
      _P3strcpy(errmsg,255,_P3str1("\022File name is empty"));
      break;
    case GXFILE_err_fileerror: 
      _P3strcpy(errmsg,255,_P3str1("\016File I/O error"));
      break;
    case GXFILE_err_noerror: 
      _P3strcpy(errmsg,255,_P3str1("\010No error"));
      break;
    case GXFILE_err_badmode: 
      _P3strcpy(errmsg,255,_P3str1("\010Bad mode"));
      break;
    case GXFILE_err_baddimension: 
      _P3strcpy(errmsg,255,_P3str1("\015Bad dimension"));
      break;
    case GXFILE_err_bad_alias_dim: 
      _P3strcpy(errmsg,255,_P3str1("\035Bad dimension for aliased set"));
      break;
    case GXFILE_err_badelementindex: 
      _P3strcpy(errmsg,255,_P3str1("\012Bad UEL Nr"));
      break;
    case GXFILE_err_badsymbolindex: 
      _P3strcpy(errmsg,255,_P3str1("\021Bad symbol number"));
      break;
    case GXFILE_err_elementsequence: 
      _P3strcpy(errmsg,255,_P3str1("\027Element out of sequence"));
      break;
    case GXFILE_err_duplicatesymbol: 
      _P3strcpy(errmsg,255,_P3str1("\020Duplicate symbol"));
      break;
    case GXFILE_err_datanotsorted: 
      _P3strcpy(errmsg,255,_P3str1("\022Data is not sorted"));
      break;
    case GXFILE_err_dataduplicate: 
      _P3strcpy(errmsg,255,_P3str1("\016Duplicate keys"));
      break;
    case GXFILE_err_unknownfilter: 
      _P3strcpy(errmsg,255,_P3str1("\016Unknown filter"));
      break;
    case GXFILE_err_badstringformat: 
      _P3strcpy(errmsg,255,_P3str1("\012Bad quotes"));
      break;
    case GXFILE_err_badidentformat: 
      _P3strcpy(errmsg,255,_P3str1("\022Illegal identifier"));
      break;
    case GXFILE_err_uelconflict: 
      _P3strcpy(errmsg,255,_P3str1("\037UEL string with different index"));
      break;
    case GXFILE_err_duplicatespecval: 
      _P3strcpy(errmsg,255,_P3str1("\027Duplicate special value"));
      break;
    case GXFILE_err_baderrorrecord: 
      _P3strcpy(errmsg,255,_P3str1("\027Bad Error record number"));
      break;
    case GXFILE_err_duplicateuel: 
      _P3strcpy(errmsg,255,_P3str1("\015Duplicate UEL"));
      break;
    case GXFILE_err_baduelstr: 
      _P3strcpy(errmsg,255,_P3str1("\016Bad UEL string"));
      break;
    case GXFILE_err_undefuel: 
      _P3strcpy(errmsg,255,_P3str1("\013Unknown UEL"));
      break;
    case GXFILE_err_uelsecondwrite: 
      _P3strcpy(errmsg,255,_P3str1("\036gdx file has UEL table already"));
      break;
    case GXFILE_err_uelnotempty: 
      _P3strcpy(errmsg,255,_P3str1("\026UEL table is not empty"));
      break;
    case GXFILE_err_bad_filter_nr: 
      _P3strcpy(errmsg,255,_P3str1("\021Bad filter number"));
      break;
    case GXFILE_err_bad_filter_indx: 
      _P3strcpy(errmsg,255,_P3str1("\023Bad index in filter"));
      break;
    case GXFILE_err_filter_unmapped: 
      _P3strcpy(errmsg,255,_P3str1("\030Unmapped index in filter"));
      break;
    case GXFILE_err_obsolete_function: 
      _P3strcpy(errmsg,255,_P3str1("\030Use of obsolete function"));
      break;
    case GXFILE_err_rawnotsorted: 
      _P3strcpy(errmsg,255,_P3str1("\040Data not sorted when writing raw"));
      break;
    case GXFILE_err_badacroindex: 
      _P3strcpy(errmsg,255,_P3str1("\025Bad index for acronym"));
      break;
    case GXFILE_err_badacronumber: 
      _P3strcpy(errmsg,255,_P3str1("\031Bad acronym record number"));
      break;
    case GXFILE_err_badacroname: 
      _P3strcpy(errmsg,255,_P3str1("\033Bad acronym name for update"));
      break;
    case GXFILE_err_acrodupemap: 
      _P3strcpy(errmsg,255,_P3str1("\034Bad acronym index for update"));
      break;
    case GXFILE_err_acrobadaddition: 
      _P3strcpy(errmsg,255,_P3str1("\035Bad addition to acronym table"));
      break;
    case GXFILE_err_unknowndomain: 
      _P3strcpy(errmsg,255,_P3str1("\016Unknown domain"));
      break;
    case GXFILE_err_baddomain: 
      _P3strcpy(errmsg,255,_P3str1("\031Domain not set with dim=1"));
      break;
    case GXFILE_err_nodomaindata: 
      _P3strcpy(errmsg,255,_P3str1("\017Set has no data"));
      break;
    case GXFILE_err_aliassetexpected: 
      _P3strcpy(errmsg,255,_P3str1("\027Set expected for domain"));
      break;
    case GXFILE_err_baddatatype: 
      _P3strcpy(errmsg,255,_P3str1("\015Bad data type"));
      break;
    case GXFILE_err_nosymbolforcomment: 
      _P3strcpy(errmsg,255,_P3str1("\033No symbol to add comment to"));
      break;
    case GXFILE_err_domainviolation: 
      _P3strcpy(errmsg,255,_P3str1("\020Domain violation"));
      break;
    case GXFILE_err_filealreadyopen: 
      _P3strcpy(errmsg,255,_P3str1("\024File is already open"));
      break;
    case GXFILE_err_filetooldforappend: 
      _P3strcpy(errmsg,255,_P3str1("\036File version to old for append"));
      break;
    case GXFILE_err_open_domsmarker1: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (DOMS_1) not found in GDX file"));
      break;
    case GXFILE_err_open_domsmarker2: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (DOMS_2) not found in GDX file"));
      break;
    case GXFILE_err_open_domsmarker3: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (DOMS_3) not found in GDX file"));
      break;
    case GXFILE_err_baddatamarker_data: 
      _P3strcpy(errmsg,255,_P3str1("\061Expected data marker (DATA) not found in GDX file"));
      break;
    case GXFILE_err_baddatamarker_dim: 
      _P3strcpy(errmsg,255,_P3str1("\060Expected data marker (DIM) not found in GDX file"));
      break;
    case GXFILE_err_open_boi: 
      _P3strcpy(errmsg,255,_P3str1("\060Expected data marker (BOI) not found in GDX file"));
      break;
    case GXFILE_err_open_fileheader: 
      _P3strcpy(errmsg,255,_P3str1("\067Expected data marker (FILEHEADER) not found in GDX file"));
      break;
    case GXFILE_err_open_filemarker: 
      _P3strcpy(errmsg,255,_P3str1("\067Expected data marker (FILEMARKER) not found in GDX file"));
      break;
    case GXFILE_err_open_symbolmarker1: 
      _P3strcpy(errmsg,255,_P3str1("\065Expected data marker (SYMBOL_1) not found in GDX file"));
      break;
    case GXFILE_err_open_symbolmarker2: 
      _P3strcpy(errmsg,255,_P3str1("\065Expected data marker (SYMBOL_2) not found in GDX file"));
      break;
    case GXFILE_err_open_uelmarker1: 
      _P3strcpy(errmsg,255,_P3str1("\062Expected data marker (UEL_1) not found in GDX file"));
      break;
    case GXFILE_err_open_uelmarker2: 
      _P3strcpy(errmsg,255,_P3str1("\062Expected data marker (UEL_2) not found in GDX file"));
      break;
    case GXFILE_err_open_textmarker1: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (TEXT_1) not found in GDX file"));
      break;
    case GXFILE_err_open_textmarker2: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (TEXT_2) not found in GDX file"));
      break;
    case GXFILE_err_open_acromarker1: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (ACRO_1) not found in GDX file"));
      break;
    case GXFILE_err_open_acromarker2: 
      _P3strcpy(errmsg,255,_P3str1("\063Expected data marker (ACRO_2) not found in GDX file"));
      break;
    case GXFILE_err_open_fileversion: 
      _P3strcpy(errmsg,255,_P3str1("\036GDX file version not supported"));
      break;
    case GXFILE_err_baddataformat: 
      _P3strcpy(errmsg,255,_P3str1("\041File not recognized as a GDX file"));
      break;
    case GXFILE_err_out_of_memory: 
      _P3strcpy(errmsg,255,_P3str1("\015Out of memory"));
      break;
    case GXFILE_err_zlib_not_found: 
      _P3strcpy(errmsg,255,_P3str1("\035Compression library not found"));
      break;
    case GXFILE_err_gdxcopy: 
      _P3strcpy(errmsg,255,_P3str1("\026GDXCOPY: Unknown error"));
      break;
    case GXFILE_err_parameter: 
      _P3strcpy(errmsg,255,_P3str1("\030GDXCOPY: Parameter error"));
      break;
    case GXFILE_err_dll_not_found: 
      _P3strcpy(errmsg,255,_P3str1("\026GDXCOPY: DLL not found"));
      break;
    case GXFILE_err_create_dir: 
      _P3strcpy(errmsg,255,_P3str1("\040GDXCOPY: Cannot create directory"));
      break;
    case GXFILE_err_file_open: 
      _P3strcpy(errmsg,255,_P3str1("\031GDXCOPY: File open failed"));
      break;
    case GXFILE_err_file_write: 
      _P3strcpy(errmsg,255,_P3str1("\043GDXCOPY: Cannot open file for write"));
      break;
    case GXFILE_err_uel_length: 
      _P3strcpy(errmsg,255,_P3str1("\043GDXCOPY: UEL length exceeds maximum"));
      break;
    case GXFILE_err_uel_register: 
      _P3strcpy(errmsg,255,_P3str1("\035GDXCOPY: Cannot register UELs"));
      break;
    case GXFILE_err_expl_text: 
      _P3strcpy(errmsg,255,_P3str1("\045GDXCOPY: Cannot save explanatory text"));
      break;
    case GXFILE_err_dimension: 
      _P3strcpy(errmsg,255,_P3str1("\042GDXCOPY: Dimension exceeds maximum"));
      break;
    case GXFILE_err_write_symbol: 
      _P3strcpy(errmsg,255,_P3str1("\035GDXCOPY: Error writing symbol"));
      break;
    case GXFILE_err_close_file: 
      _P3strcpy(errmsg,255,_P3str1("\033GDXCOPY: Error closing file"));
      break;
    case GXFILE_err_cannot_delete: 
      _P3strcpy(errmsg,255,_P3str1("\033GDXCOPY: Cannot delete file"));
      break;
    case GXFILE_err_cannot_rename: 
      _P3strcpy(errmsg,255,_P3str1("\033GDXCOPY: Cannot rename file"));
      break;
    default:
      SYSUTILS_P3_syserrormessage(errmsg,255,errnr);
  }
  return result;
}  /* gdxerrorstr */

Function(SYSTEM_int64 ) GXFILE_tueltable_DOT_memoryused(
  GXFILE_tueltable self)
{
  SYSTEM_int64 result;

  result = STRHASH_txstrhashlist_DOT_memoryused(ValueCast(
    STRHASH_txstrhashlist,self)) + 
    GXFILE_tintegermapping_DOT_memoryused(self->
    GXFILE_tueltable_DOT_usruel2ent);
  return result;
}  /* memoryused */

Constructor(GXFILE_tacronym ) GXFILE_tacronym_DOT_create(
  GXFILE_tacronym self,
  const SYSTEM_ansichar *name,
  const SYSTEM_ansichar *text,
  SYSTEM_integer map)
{
  ValueCast(GXFILE_tacronym,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GXFILE_tacronym_DOT_acrname = STRUTILX_newstringm(name,&self->
    GXFILE_tacronym_DOT_fstrmemory);
  {
    SYSTEM_shortstring _t1;

    self->GXFILE_tacronym_DOT_acrtext = STRUTILX_newstringm(
      GXFILE_makegoodexpltext(_t1,255,text),&self->
      GXFILE_tacronym_DOT_fstrmemory);
  }
  self->GXFILE_tacronym_DOT_acrmap = map;
  self->GXFILE_tacronym_DOT_acrreadmap =  -1;
  self->GXFILE_tacronym_DOT_acrautogen = SYSTEM_false;
  return self;
}  /* create */

Constructor(GXFILE_tacronym ) GXFILE_tacronym_DOT_createfromstream(
  GXFILE_tacronym self,
  GMSSTRM_txstream s)
{
  SYSTEM_shortstring s1, s2;

  GMSSTRM_txstream_DOT_readstring(s1,255,s);
  GMSSTRM_txstream_DOT_readstring(s2,255,s);
  ValueCast(GXFILE_tacronym,GXFILE_tacronym_DOT_create(self,s1,s2,VirtMethodCall(
    s, GMSSTRM_txstream_DOT_readinteger_T, 7, (s))));
  return self;
}  /* createfromstream */

Destructor(GXFILE_tacronym ) GXFILE_tacronym_DOT_destroy(
  GXFILE_tacronym self)
{
  STRUTILX_disposestring(self->GXFILE_tacronym_DOT_acrname);
  STRUTILX_disposestring(self->GXFILE_tacronym_DOT_acrtext);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_int64 ) GXFILE_tacronym_DOT_memoryused(
  GXFILE_tacronym self)
{
  SYSTEM_int64 result;

  result = self->GXFILE_tacronym_DOT_fstrmemory;
  return result;
}  /* memoryused */

Procedure GXFILE_tacronym_DOT_savetostream(
  GXFILE_tacronym self,
  GMSSTRM_txstream s)
{
  { register GMSSTRM_txstream_OD *_W2=s;
    if (self->GXFILE_tacronym_DOT_acrname == NULL) { 
      {
        SYSTEM_shortstring _t1;
        _P3STR_255 _t2;

        GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,
          _W2),_P3strcat(_t2,255,_P3str1("\013UnknownACRO"),
          STRUTILX_inttostr(_t1,255,self->
          GXFILE_tacronym_DOT_acrmap)));
      }
    } else 
      {
        SYSTEM_shortstring _t1;

        GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,
          _W2),STRUTILX_getstring(_t1,255,self->
          GXFILE_tacronym_DOT_acrname));
      }
    {
      SYSTEM_shortstring _t1;

      GMSSTRM_txstream_DOT_writestring(ValueCast(GMSSTRM_txstream,_W2),
        STRUTILX_getstring(_t1,255,self->
        GXFILE_tacronym_DOT_acrtext));
    }
    GMSSTRM_txstream_DOT_writeinteger(ValueCast(GMSSTRM_txstream,_W2),
      self->GXFILE_tacronym_DOT_acrmap);

  }
}  /* savetostream */

Function(SYSTEM_integer ) GXFILE_tacronymlist_DOT_addentry(
  GXFILE_tacronymlist self,
  const SYSTEM_ansichar *name,
  const SYSTEM_ansichar *text,
  SYSTEM_integer map)
{
  SYSTEM_integer result;
  GXFILE_tacronym a;

  a = ValueCast(GXFILE_tacronym,GXFILE_tacronym_DOT_create(ValueCast(
    GXFILE_tacronym,_P3alloc_object(&GXFILE_tacronym_CD)),name,text,
    map));
  result = GMSOBJ_txlist_DOT_add(self->GXFILE_tacronymlist_DOT_flist,a);
  return result;
}  /* addentry */

Procedure GXFILE_tacronymlist_DOT_checkentry(
  GXFILE_tacronymlist self,
  SYSTEM_integer map)
{
  if (GXFILE_tacronymlist_DOT_findentry(self,map) < 0) 
    GXFILE_tacronymlist_DOT_addentry(self,_P3str1("\000"),_P3str1("\000"),
      map);
}  /* checkentry */

Constructor(GXFILE_tacronymlist ) GXFILE_tacronymlist_DOT_create(
  GXFILE_tacronymlist self)
{
  self->GXFILE_tacronymlist_DOT_flist = ValueCast(GMSOBJ_txlist,
    SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,_P3alloc_object(&
    GMSOBJ_txlist_CD))));
  return self;
}  /* create */

Destructor(GXFILE_tacronymlist ) GXFILE_tacronymlist_DOT_destroy(
  GXFILE_tacronymlist self)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,
        GMSOBJ_txlist_DOT_get(self->GXFILE_tacronymlist_DOT_flist,n)));
    } while (n++ !=  _stop);

  }
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tacronymlist_DOT_flist));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_integer ) GXFILE_tacronymlist_DOT_findentry(
  GXFILE_tacronymlist self,
  SYSTEM_integer map)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      { register GXFILE_tacronym_OD *_W2=ValueCast(GXFILE_tacronym,
        GMSOBJ_txlist_DOT_get(self->GXFILE_tacronymlist_DOT_flist,n));
        if (_W2->GXFILE_tacronym_DOT_acrmap == map) {
          result = n;
          return result;
        } 

      }
    } while (n++ !=  _stop);

  }
  result =  -1;
  return result;
}  /* findentry */

Function(SYSTEM_integer ) GXFILE_tacronymlist_DOT_findname(
  GXFILE_tacronymlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if (STRUTILX_pstruequal((ValueCast(GXFILE_tacronym,
        GMSOBJ_txlist_DOT_get(self->GXFILE_tacronymlist_DOT_flist,n)))->
        GXFILE_tacronym_DOT_acrname,ValueCast(SYSTEM_P3_pshortstring,s))) {
        result = n;
        return result;
      } 
    } while (n++ !=  _stop);

  }
  result =  -1;
  return result;
}  /* findname */

Procedure GXFILE_tacronymlist_DOT_loadfromstream(
  GXFILE_tacronymlist self,
  GMSSTRM_txstream s)
{
  SYSTEM_integer cnt;
  GXFILE_tacronym a;

  cnt = VirtMethodCall(s, GMSSTRM_txstream_DOT_readinteger_T, 7, (s));
  GMSOBJ_txlist_DOT_clear(self->GXFILE_tacronymlist_DOT_flist);
  GMSOBJ_txlist_DOT_setcapacity(self->GXFILE_tacronymlist_DOT_flist,
    cnt);
  while (self->GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount < 
    cnt) {
    a = ValueCast(GXFILE_tacronym,GXFILE_tacronym_DOT_createfromstream(ValueCast(
      GXFILE_tacronym,_P3alloc_object(&GXFILE_tacronym_CD)),s));
    GMSOBJ_txlist_DOT_add(self->GXFILE_tacronymlist_DOT_flist,a);
  
}
}  /* loadfromstream */

Function(SYSTEM_int64 ) GXFILE_tacronymlist_DOT_memoryused(
  GXFILE_tacronymlist self)
{
  SYSTEM_int64 result;
  SYSTEM_integer n;

  result = GMSOBJ_txlist_DOT_memoryused(self->
    GXFILE_tacronymlist_DOT_flist) + self->
    GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount * sizeof(
    GXFILE_tacronym);
  { register SYSTEM_int32 _stop = self->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      result = result + GXFILE_tacronym_DOT_memoryused(ValueCast(
        GXFILE_tacronym,GMSOBJ_txlist_DOT_get(self->
        GXFILE_tacronymlist_DOT_flist,n)));
    } while (n++ !=  _stop);

  }
  return result;
}  /* memoryused */

Procedure GXFILE_tacronymlist_DOT_savetostream(
  GXFILE_tacronymlist self,
  GMSSTRM_txstream s)
{
  SYSTEM_integer n;

  GMSSTRM_txstream_DOT_writeinteger(s,self->
    GXFILE_tacronymlist_DOT_flist->GMSOBJ_txlist_DOT_fcount);
  { register SYSTEM_int32 _stop = self->GXFILE_tacronymlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      GXFILE_tacronym_DOT_savetostream(ValueCast(GXFILE_tacronym,
        GMSOBJ_txlist_DOT_get(self->GXFILE_tacronymlist_DOT_flist,n)),
        s);
    } while (n++ !=  _stop);

  }
}  /* savetostream */

Destructor(GXFILE_tgxfileobj ) GXFILE_tgxfileobj_DOT_destroy(
  GXFILE_tgxfileobj self)
{
  if (self->GXFILE_tgxfileobj_DOT_fmode != GXFILE_f_not_open) {
    self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fr_init;
    GXFILE_tgxfileobj_DOT_gdxclose(self);
  } 
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxautoconvert(
  GXFILE_tgxfileobj self,
  SYSTEM_integer nv)
{
  SYSTEM_integer result;

  result = self->GXFILE_tgxfileobj_DOT_autoconvert;
  self->GXFILE_tgxfileobj_DOT_autoconvert = nv;
  return result;
}  /* gdxautoconvert */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxcurrentdim(
  GXFILE_tgxfileobj self)
{
  SYSTEM_integer result;

  result = self->GXFILE_tgxfileobj_DOT_fcurrentdim;
  return result;
}  /* gdxcurrentdim */

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_resultwillbesorted(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *adomainnrs)
{
  SYSTEM_boolean result;
  SYSTEM_integer d;
  SYSTEM_integer n;
  GXFILE_tdfilter pfilter;

  result = SYSTEM_false;
  { register SYSTEM_int32 _stop = self->
      GXFILE_tgxfileobj_DOT_fcurrentdim;
    if ((d = 1) <=  _stop) do {
      switch (adomainnrs[d - 1]) {
        case GXDEFS_domc_unmapped: 
          SYSTEM_continue(CNT_12);
          break;
        case GXDEFS_domc_expand: 
          if (GXFILE_tueltable_DOT_getmaptouserstatus(self->
            GXFILE_tgxfileobj_DOT_ueltable) == GXFILE_map_unsorted) 
            return result;
          if (d == 1) { 
            if (GXFILE_tueltable_DOT_getmaptouserstatus(self->
              GXFILE_tgxfileobj_DOT_ueltable) >= GXFILE_map_sortgrow) { 
              SYSTEM_continue(CNT_12);
            } else 
              return result;
          } else 
            if (GXFILE_tueltable_DOT_getmaptouserstatus(self->
              GXFILE_tgxfileobj_DOT_ueltable) == GXFILE_map_sortfull) { 
              SYSTEM_continue(CNT_12);
            } else 
              return result;
          break;
        case GXDEFS_domc_strict: 
          if (GXFILE_tueltable_DOT_getmaptouserstatus(self->
            GXFILE_tgxfileobj_DOT_ueltable) == GXFILE_map_unsorted) 
            return result;
          break;
        default:
          if (GXFILE_tueltable_DOT_getmaptouserstatus(self->
            GXFILE_tgxfileobj_DOT_ueltable) >= GXFILE_map_sorted) 
            SYSTEM_continue(CNT_12);
          pfilter = GXFILE_tfilterlist_DOT_findfilter(self->
            GXFILE_tgxfileobj_DOT_filterlist,adomainnrs[d - 1]);
          if (!pfilter->GXFILE_tdfilter_DOT_filtsorted) 
            return result;
      }
CNT_12:;
    } while (d++ !=  _stop);
BRK_12:;

  }
  result = SYSTEM_true;
  return result;
}  /* resultwillbesorted */

Constructor(GXFILE_tdfilter ) GXFILE_tdfilter_DOT_create(
  GXFILE_tdfilter self,
  SYSTEM_integer nr,
  SYSTEM_integer userhigh)
{
  ValueCast(GXFILE_tdfilter,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GXFILE_tdfilter_DOT_filtnumber = nr;
  self->GXFILE_tdfilter_DOT_filtmap = ValueCast(
    GMSOBJ_tbooleanbitarray,GMSOBJ_tbooleanbitarray_DOT_create(ValueCast(
    GMSOBJ_tbooleanbitarray,_P3alloc_object(&
    GMSOBJ_tbooleanbitarray_CD))));
  self->GXFILE_tdfilter_DOT_filtmaxuel = userhigh;
  GMSOBJ_tbooleanbitarray_DOT_sethighindex(self->
    GXFILE_tdfilter_DOT_filtmap,self->GXFILE_tdfilter_DOT_filtmaxuel);
  self->GXFILE_tdfilter_DOT_filtsorted = SYSTEM_false;
  return self;
}  /* create */

Destructor(GXFILE_tdfilter ) GXFILE_tdfilter_DOT_destroy(
  GXFILE_tdfilter self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tdfilter_DOT_filtmap));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_boolean ) GXFILE_tdfilter_DOT_infilter(
  GXFILE_tdfilter self,
  SYSTEM_integer v)
{
  SYSTEM_boolean result;

  result = v >= 0 && v <= self->GXFILE_tdfilter_DOT_filtmaxuel && 
    GMSOBJ_tbooleanbitarray_DOT_getbit(self->
    GXFILE_tdfilter_DOT_filtmap,v);
  return result;
}  /* infilter */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadrawfast(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  GXDEFS_tdatastoreproc dp,
  SYSTEM_integer *nrrecs)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex xdomains;
  GMSSPECS_tvarreca avals;
  SYSTEM_integer afdim;

  for (d = 1;d <= (SYSTEM_int32)GMSSPECS_maxdim;++d) {
    xdomains[d - 1] = GXDEFS_domc_unmapped;
  }
  *nrrecs = GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\022gdxDataReadRawFast"),
    synr,xdomains,GXFILE_fr_raw_data);
  while (GXFILE_tgxfileobj_DOT_doread(self,avals,&afdim)) {

    (*dp)(self->GXFILE_tgxfileobj_DOT_lastelem,avals);
}
  GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  if (*nrrecs >= 0) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxdatareadrawfast */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenappend(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *producer,
  SYSTEM_integer *errnr)
{
  SYSTEM_integer result;

  _P3strcpy(self->GXFILE_tgxfileobj_DOT_fproducer2,255,producer);
  self->GXFILE_tgxfileobj_DOT_appendactive = SYSTEM_true;
  result = GXFILE_tgxfileobj_DOT_gdxopenreadxx(self,filename,
    GMSSTRM_fmopenreadwrite,errnr);
  if (result == 0 || *errnr != 0) 
    return result;
  if (self->GXFILE_tgxfileobj_DOT_versionread < 7) {
    GXFILE_tgxfileobj_DOT_reporterror(self,
      GXFILE_err_filetooldforappend);
    GXFILE_tgxfileobj_DOT_gdxclose(self);
    return result;
  } 
  self->GXFILE_tgxfileobj_DOT_fmode = GXFILE_fw_init;
  self->GXFILE_tgxfileobj_DOT_fstatus = GXFILE_stat_write;
  VirtMethodCall(ValueCast(GMSSTRM_txstream,self->
    GXFILE_tgxfileobj_DOT_ffile), GMSSTRM_txstream_DOT_setposition_T, 2, (ValueCast(
    GMSSTRM_txstream,self->GXFILE_tgxfileobj_DOT_ffile),self->
    GXFILE_tgxfileobj_DOT_nextwriteposition));
  self->GXFILE_tgxfileobj_DOT_compressout = self->
    GXFILE_tgxfileobj_DOT_douncompress;
  return result;
}  /* gdxopenappend */

Function(SYSTEM_int64 ) GXFILE_tgxfileobj_DOT_gdxgetmemoryused(
  GXFILE_tgxfileobj self)
{
  SYSTEM_int64 result;

  result = 0;
  if (self->GXFILE_tgxfileobj_DOT_ueltable != NULL) 
    result = result + GXFILE_tueltable_DOT_memoryused(self->
      GXFILE_tgxfileobj_DOT_ueltable);
  if (self->GXFILE_tgxfileobj_DOT_settextlist != NULL) 
    result = result + GMSOBJ_txcustomstringlist_DOT_memoryused(ValueCast(
      GMSOBJ_txcustomstringlist,self->
      GXFILE_tgxfileobj_DOT_settextlist));
  if (self->GXFILE_tgxfileobj_DOT_namelist != NULL) 
    result = result + STRHASH_txstrhashlist_DOT_memoryused(self->
      GXFILE_tgxfileobj_DOT_namelist);
  if (self->GXFILE_tgxfileobj_DOT_domainstrlist != NULL) 
    result = result + STRHASH_txstrhashlist_DOT_memoryused(self->
      GXFILE_tgxfileobj_DOT_domainstrlist);
  if (self->GXFILE_tgxfileobj_DOT_sortlist != NULL) 
    result = result + DATASTORAGE_tlinkeddata_DOT_memoryused(self->
      GXFILE_tgxfileobj_DOT_sortlist);
  if (self->GXFILE_tgxfileobj_DOT_errorlist != NULL) 
    result = result + GMSDATA_ttblgamsdata_DOT_memoryused(self->
      GXFILE_tgxfileobj_DOT_errorlist);
  if (self->GXFILE_tgxfileobj_DOT_filterlist != NULL) 
    result = result + GXFILE_tfilterlist_DOT_memoryused(self->
      GXFILE_tgxfileobj_DOT_filterlist);
  return result;
}  /* gdxgetmemoryused */

Function(SYSTEM_int64 ) GXFILE_tdfilter_DOT_memoryused(
  GXFILE_tdfilter self)
{
  SYSTEM_int64 result;

  result = GMSOBJ_tbooleanbitarray_DOT_memoryused(self->
    GXFILE_tdfilter_DOT_filtmap);
  return result;
}  /* memoryused */

Procedure GXFILE_tfilterlist_DOT_addfilter(
  GXFILE_tfilterlist self,
  GXFILE_tdfilter f)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GXFILE_tfilterlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if ((ValueCast(GXFILE_tdfilter,GMSOBJ_txlist_DOT_get(self->
        GXFILE_tfilterlist_DOT_flist,n)))->
        GXFILE_tdfilter_DOT_filtnumber == f->
        GXFILE_tdfilter_DOT_filtnumber) {
        GXFILE_tfilterlist_DOT_deletefilter(self,n);
        SYSTEM_break(BRK_13);
      } 
CNT_13:;
    } while (n++ !=  _stop);
BRK_13:;

  }
  GMSOBJ_txlist_DOT_add(self->GXFILE_tfilterlist_DOT_flist,f);
}  /* addfilter */

Constructor(GXFILE_tfilterlist ) GXFILE_tfilterlist_DOT_create(
  GXFILE_tfilterlist self)
{
  ValueCast(GXFILE_tfilterlist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GXFILE_tfilterlist_DOT_flist = ValueCast(GMSOBJ_txlist,
    SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,_P3alloc_object(&
    GMSOBJ_txlist_CD))));
  return self;
}  /* create */

Procedure GXFILE_tfilterlist_DOT_deletefilter(
  GXFILE_tfilterlist self,
  SYSTEM_integer ix)
{
  GXFILE_tdfilter p;

  p = ValueCast(GXFILE_tdfilter,GMSOBJ_txlist_DOT_get(self->
    GXFILE_tfilterlist_DOT_flist,ix));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,p));
  GMSOBJ_txlist_DOT_delete(self->GXFILE_tfilterlist_DOT_flist,ix);
}  /* deletefilter */

Destructor(GXFILE_tfilterlist ) GXFILE_tfilterlist_DOT_destroy(
  GXFILE_tfilterlist self)
{
  while (self->GXFILE_tfilterlist_DOT_flist->GMSOBJ_txlist_DOT_fcount > 0) {

    GXFILE_tfilterlist_DOT_deletefilter(self,self->
      GXFILE_tfilterlist_DOT_flist->GMSOBJ_txlist_DOT_fcount - 1);
}
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GXFILE_tfilterlist_DOT_flist));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(GXFILE_tdfilter ) GXFILE_tfilterlist_DOT_findfilter(
  GXFILE_tfilterlist self,
  SYSTEM_integer nr)
{
  GXFILE_tdfilter result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GXFILE_tfilterlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      result = ValueCast(GXFILE_tdfilter,GMSOBJ_txlist_DOT_get(self->
        GXFILE_tfilterlist_DOT_flist,n));
      if (result->GXFILE_tdfilter_DOT_filtnumber == nr) 
        return result;
    
    } while (n++ !=  _stop);

  }
  result = NULL;
  return result;
}  /* findfilter */

Function(SYSTEM_int64 ) GXFILE_tfilterlist_DOT_memoryused(
  GXFILE_tfilterlist self)
{
  SYSTEM_int64 result;
  SYSTEM_integer n;

  result = GMSOBJ_txlist_DOT_memoryused(self->
    GXFILE_tfilterlist_DOT_flist) + self->GXFILE_tfilterlist_DOT_flist->
    GMSOBJ_txlist_DOT_fcount * sizeof(GXFILE_tdfilter);
  { register SYSTEM_int32 _stop = self->GXFILE_tfilterlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      result = result + GXFILE_tdfilter_DOT_memoryused(ValueCast(
        GXFILE_tdfilter,GMSOBJ_txlist_DOT_get(self->
        GXFILE_tfilterlist_DOT_flist,n)));
    } while (n++ !=  _stop);

  }
  return result;
}  /* memoryused */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetdomainelements(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer dimpos,
  SYSTEM_integer filternr,
  GXDEFS_tdomainindexproc dp,
  SYSTEM_integer *nrelem,
  SYSTEM_pointer uptr)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GXDEFS_tgdxuelindex xdomains;
  GMSSPECS_tvarreca avals;
  SYSTEM_integer afdim;
  SYSTEM_integer dim;
  GXFILE_tintegermapping domainindxs;
  SYSTEM_integer n;
  GXFILE_tdfilter dfilter;
  SYSTEM_integer rawnr;
  SYSTEM_integer mapnr;
  GMSDATA_ttblgamsdata sortl;
  GMSSPECS_tindex index;

  self->GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp = dp;
  result = GXFILE_myfalse;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,synr >= 1 && synr <= 
    self->GXFILE_tgxfileobj_DOT_namelist->
    STRHASH_txstrhashlist_DOT_fcount,GXFILE_err_badsymbolindex)) 
    return result;
  dim = (ValueCast(GXFILE_pgdxsymbrecord,
    STRHASH_txstrhashlist_DOT_getobject(self->
    GXFILE_tgxfileobj_DOT_namelist,synr)))->sdim;
  if (dim == 0) 
    return result;
  if (GXFILE_tgxfileobj_DOT_errorcondition(self,dimpos >= 1 && 
    dimpos <= dim,GXFILE_err_baddimension)) 
    return result;
  if (filternr == GXDEFS_domc_expand) { 
    dfilter = NULL;
  } else {
    dfilter = GXFILE_tfilterlist_DOT_findfilter(self->
      GXFILE_tgxfileobj_DOT_filterlist,filternr);
    if (dfilter == NULL) {
      GXFILE_tgxfileobj_DOT_reporterror(self,GXFILE_err_unknownfilter);
      return result;
    } 
  } 
  domainindxs = ValueCast(GXFILE_tintegermapping,
    GXFILE_tintegermapping_DOT_create(ValueCast(GXFILE_tintegermapping,
    _P3alloc_object(&GXFILE_tintegermapping_CD))));
  { register SYSTEM_int32 _stop = dim;
    if ((d = 1) <=  _stop) do {
      xdomains[d - 1] = GXDEFS_domc_unmapped;
    } while (d++ !=  _stop);

  }
  GXFILE_tgxfileobj_DOT_preparesymbolread(self,_P3str1("\014gdxGetDomain"),
    synr,xdomains,GXFILE_fr_raw_data);
  while (GXFILE_tgxfileobj_DOT_doread(self,avals,&afdim)) {
    rawnr = self->GXFILE_tgxfileobj_DOT_lastelem[dimpos - 1];
    if (dfilter != NULL) {
      mapnr = GXFILE_tueltable_DOT_getusermap(self->
        GXFILE_tgxfileobj_DOT_ueltable,rawnr);
      if (!GXFILE_tdfilter_DOT_infilter(dfilter,mapnr)) {
        GXFILE_tgxfileobj_DOT_addtoerrorlist(self,self->
          GXFILE_tgxfileobj_DOT_lastelem,avals);
        SYSTEM_continue(CNT_14);
      } 
    } 
    GXFILE_tintegermapping_DOT_setmapping(domainindxs,rawnr,1);
  
CNT_14:;
  }
BRK_14:;
  GXFILE_tgxfileobj_DOT_gdxdatareaddone(self);
  *nrelem = 0;
  if (!_P3assigned(dp)) { 
    { register SYSTEM_int32 _stop = domainindxs->
        GXFILE_tintegermapping_DOT_fhighestindex;
      if ((n = 1) <=  _stop) do {
        if (GXFILE_tintegermapping_DOT_getmapping(domainindxs,n) == 1) 
          *nrelem = *nrelem + 1;
      } while (n++ !=  _stop);

    }
  } else {
    sortl = ValueCast(GMSDATA_ttblgamsdata,
      GMSDATA_ttblgamsdata_DOT_create(ValueCast(GMSDATA_ttblgamsdata,
      _P3alloc_object(&GMSDATA_ttblgamsdata_CD)),1,sizeof(
      SYSTEM_longint)));
    { register SYSTEM_int32 _stop = domainindxs->
        GXFILE_tintegermapping_DOT_fhighestindex;
      if ((n = 1) <=  _stop) do {
        if (GXFILE_tintegermapping_DOT_getmapping(domainindxs,n) == 1) {
          *nrelem = *nrelem + 1;
          index[0] = GXFILE_tueltable_DOT_newusruel(self->
            GXFILE_tgxfileobj_DOT_ueltable,n);
          GMSDATA_ttblgamsdata_DOT_addrecord(sortl,index,&n);
        } 
      } while (n++ !=  _stop);

    }
    GMSDATA_ttblgamsdata_DOT_sort(sortl);
    { register SYSTEM_int32 _stop = GMSDATA_ttblgamsdata_DOT_getcount(
        sortl) - 1;
      if ((n = 0) <=  _stop) do {
        GMSDATA_ttblgamsdata_DOT_getrecord(sortl,n,index,&rawnr);
        GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_fc(self,rawnr,
          index[0],uptr);
      
      } while (n++ !=  _stop);

    }
    SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,sortl));
  } 
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,domainindxs));
  if (*nrelem >= 0) { 
    result = GXFILE_mytrue;
  } else 
    result = GXFILE_myfalse;
  return result;
}  /* gdxgetdomainelements */

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxrenameuel(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *oldname,
  const SYSTEM_ansichar *newname)
{
  SYSTEM_integer result;
  SYSTEM_integer n;
  SYSTEM_shortstring s;

  if (self->GXFILE_tgxfileobj_DOT_ueltable == NULL) {
    result =  -1;
    return result;
  } 
  SYSUTILS_P3_trimright(s,255,newname);
  if (!GXDEFS_gooduelstring(s)) {
    result = GXFILE_err_baduelstr;
    return result;
  } 
  {
    SYSTEM_shortstring _t1;

    n = STRHASH_txstrhashlist_DOT_indexof(ValueCast(
      STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),
      SYSUTILS_P3_trimright(_t1,255,oldname));
  }
  if (n < 0) {
    result = 2;
    return result;
  } 
  if (STRHASH_txstrhashlist_DOT_indexof(ValueCast(
    STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),s) >= 0) {
    result = 3;
    return result;
  } 
  STRHASH_txstrhashlist_DOT_renameentry(ValueCast(
    STRHASH_txstrhashlist,self->GXFILE_tgxfileobj_DOT_ueltable),n,s);
  result = 0;
  return result;
}  /* gdxrenameuel */

/* unit gxfile */
void _Init_Module_gxfile(void)
{
  GDLAUDIT_gdlsetauditlinelib(_P3str1("\124_GAMS_GDX Library      24.5.0 r50613 ALFA Released 27Jan15 LEG x86 64bit/Linux_SMAG_"));
  _P3strclr(GXFILE_dllloadpath);
} /* _Init_Module_gxfile */

void _Final_Module_gxfile(void)
{
} /* _Final_Module_gxfile */

