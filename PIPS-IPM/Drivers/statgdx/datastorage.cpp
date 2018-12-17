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
#include "gmsheapnew.h"
#include "datastorage.h"

/* NEEDED? (g++ doesn't like) -> 
SYSTEM_classdescriptor ** from sy_fwd_class(1) ** 
  DATASTORAGE_tgamsdatasearcher_CD;
  */

void * const DATASTORAGE_tgamsdatastore_VT[] = {(void*)&
  DATASTORAGE_tgamsdatastore_DOT_destroy, (void*)&_P3_abstract_call1, (void*)&
  _P3_abstract_call1, (void*)&_P3_abstract_call1, (void*)&
  _P3_abstract_call1, (void*)&DATASTORAGE_tgamsdatastore_DOT_startread, (void*)&
  _P3_abstract_call1, (void*)&DATASTORAGE_tgamsdatastore_DOT_endread, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startassign, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_endassign, (void*)&_P3_abstract_call1, (void*)&
  _P3_abstract_call1};

/* Class descriptor for 'tgamsdatastore' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatastore_CD = {
  _P3str1("\016tgamsdatastore"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatastore_OD), DATASTORAGE_tgamsdatastore_VT, NULL};


void * const DATASTORAGE_tgamsdatasearcher_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatasearcher_DOT_startsearch, (void*)&
  _P3_abstract_call1};

/* Class descriptor for 'tgamsdatasearcher' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasearcher_CD = {
  _P3str1("\021tgamsdatasearcher"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatasearcher_OD), 
    DATASTORAGE_tgamsdatasearcher_VT, NULL};


void * const DATASTORAGE_treclist_VT[] = {(void*)&
  DATASTORAGE_treclist_DOT_destroy, (void*)&
  DATASTORAGE_treclist_DOT_grow};

/* Class descriptor for 'treclist' */
const SYSTEM_classdescriptor_t DATASTORAGE_treclist_CD = {
  _P3str1("\010treclist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(DATASTORAGE_treclist_OD), DATASTORAGE_treclist_VT, NULL};

/* NEEDED? (g++ doesn't like) -> 
SYSTEM_classdescriptor ** from sy_fwd_class(1) ** 
  DATASTORAGE_tgamsdatatablesearcher_CD;
  */

void * const DATASTORAGE_tgamsdatatable_VT[] = {(void*)&
  DATASTORAGE_tgamsdatatable_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_getcount, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_clear, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_insertrecord, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_loadrecord, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startread, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_getnextkey, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_endread, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startassign, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_endassign, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_memoryused, (void*)&
  DATASTORAGE_tgamsdatatable_DOT_createsearcher};

/* Class descriptor for 'tgamsdatatable' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatatable_CD = {
  _P3str1("\016tgamsdatatable"), 
  &DATASTORAGE_tgamsdatastore_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatatable_OD), DATASTORAGE_tgamsdatatable_VT, NULL};


void * const DATASTORAGE_tgamsdatatablesearcher_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatatablesearcher_DOT_startsearch, (void*)&
  DATASTORAGE_tgamsdatatablesearcher_DOT_search};

/* Class descriptor for 'tgamsdatatablesearcher' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatatablesearcher_CD = {
  _P3str1("\026tgamsdatatablesearcher"), 
  &DATASTORAGE_tgamsdatasearcher_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatatablesearcher_OD), 
    DATASTORAGE_tgamsdatatablesearcher_VT, NULL};


void * const DATASTORAGE_tintegermapping_VT[] = {(void*)&
  DATASTORAGE_tintegermapping_DOT_destroy};

/* Class descriptor for 'tintegermapping' */
const SYSTEM_classdescriptor_t DATASTORAGE_tintegermapping_CD = {
  _P3str1("\017tintegermapping"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(DATASTORAGE_tintegermapping_OD), 
    DATASTORAGE_tintegermapping_VT, NULL};


void * const DATASTORAGE_tgamsdatafull_VT[] = {(void*)&
  DATASTORAGE_tgamsdatafull_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatafull_DOT_getcount, (void*)&
  DATASTORAGE_tgamsdatafull_DOT_clear, (void*)&_P3_abstract_call1, (void*)&
  DATASTORAGE_tgamsdatafull_DOT_loadrecord, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startread, (void*)&_P3_abstract_call1, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_endread, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startassign, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_endassign, (void*)&
  DATASTORAGE_tgamsdatafull_DOT_memoryused, (void*)&
  DATASTORAGE_tgamsdatafull_DOT_createsearcher};

/* Class descriptor for 'tgamsdatafull' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatafull_CD = {
  _P3str1("\015tgamsdatafull"), 
  &DATASTORAGE_tgamsdatastore_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatafull_OD), DATASTORAGE_tgamsdatafull_VT, NULL};


void * const DATASTORAGE_tgamsdatafullsearcher_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatasearcher_DOT_startsearch, (void*)&
  DATASTORAGE_tgamsdatafullsearcher_DOT_search};

/* Class descriptor for 'tgamsdatafullsearcher' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatafullsearcher_CD = {
  _P3str1("\025tgamsdatafullsearcher"), 
  &DATASTORAGE_tgamsdatasearcher_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatafullsearcher_OD), 
    DATASTORAGE_tgamsdatafullsearcher_VT, NULL};


void * const DATASTORAGE_tgamsdatasparse_VT[] = {(void*)&
  DATASTORAGE_tgamsdatasparse_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_getcount, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_clear, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_insertrecord, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_loadrecord, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_startread, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_getnextkey, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_endread, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startassign, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_endassign, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_memoryused, (void*)&
  DATASTORAGE_tgamsdatasparse_DOT_createsearcher};

/* Class descriptor for 'tgamsdatasparse' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasparse_CD = {
  _P3str1("\017tgamsdatasparse"), 
  &DATASTORAGE_tgamsdatastore_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatasparse_OD), 
    DATASTORAGE_tgamsdatasparse_VT, NULL};


void * const DATASTORAGE_tgamsdatasparsesearcher_VT[] = {(void*)&
  DATASTORAGE_tgamsdatasparsesearcher_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatasparsesearcher_DOT_startsearch, (void*)&
  DATASTORAGE_tgamsdatasparsesearcher_DOT_search};

/* Class descriptor for 'tgamsdatasparsesearcher' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasparsesearcher_CD = {
  _P3str1("\027tgamsdatasparsesearcher"), 
  &DATASTORAGE_tgamsdatasearcher_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatasparsesearcher_OD), 
    DATASTORAGE_tgamsdatasparsesearcher_VT, NULL};


void * const DATASTORAGE_tlinkeddata_VT[] = {(void*)&
  DATASTORAGE_tlinkeddata_DOT_destroy};

/* Class descriptor for 'tlinkeddata' */
const SYSTEM_classdescriptor_t DATASTORAGE_tlinkeddata_CD = {
  _P3str1("\013tlinkeddata"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(DATASTORAGE_tlinkeddata_OD), DATASTORAGE_tlinkeddata_VT, NULL};


void * const DATASTORAGE_tgamshashlist_VT[] = {(void*)&
  DATASTORAGE_tgamshashlist_DOT_destroy};

/* Class descriptor for 'tgamshashlist' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamshashlist_CD = {
  _P3str1("\015tgamshashlist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamshashlist_OD), DATASTORAGE_tgamshashlist_VT, NULL};

/* NEEDED? (g++ doesn't like) -> 
SYSTEM_classdescriptor ** from sy_fwd_class(1) ** 
  DATASTORAGE_tgamsdatahashedsearcher_CD;
  */

void * const DATASTORAGE_tgamsdatahashed_VT[] = {(void*)&
  DATASTORAGE_tgamsdatahashed_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_getcount, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_clear, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_insertrecord, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_loadrecord, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_startread, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_getnextkey, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_endread, (void*)&
  DATASTORAGE_tgamsdatastore_DOT_startassign, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_endassign, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_memoryused, (void*)&
  DATASTORAGE_tgamsdatahashed_DOT_createsearcher};

/* Class descriptor for 'tgamsdatahashed' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatahashed_CD = {
  _P3str1("\017tgamsdatahashed"), 
  &DATASTORAGE_tgamsdatastore_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatahashed_OD), 
    DATASTORAGE_tgamsdatahashed_VT, NULL};


void * const DATASTORAGE_tgamsdatahashedsearcher_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy, (void*)&
  DATASTORAGE_tgamsdatasearcher_DOT_startsearch, (void*)&
  DATASTORAGE_tgamsdatahashedsearcher_DOT_search};

/* Class descriptor for 'tgamsdatahashedsearcher' */
const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatahashedsearcher_CD = {
  _P3str1("\027tgamsdatahashedsearcher"), 
  &DATASTORAGE_tgamsdatasearcher_CD, NULL, 0, 
  sizeof(DATASTORAGE_tgamsdatahashedsearcher_OD), 
    DATASTORAGE_tgamsdatahashedsearcher_VT, NULL};


static Function(SYSTEM_boolean ) DATASTORAGE_dataequal(
  const SYSTEM_untyped *adata,
  const SYSTEM_untyped *adef,
  SYSTEM_integer asize)
{
  SYSTEM_boolean result;
  SYSTEM_integer k;
  GMSGEN_pbytedataarray s1, s2;

  result = SYSTEM_false;
  s1 = ValueCast(GMSGEN_pbytedataarray,adata);
  s2 = ValueCast(GMSGEN_pbytedataarray,adef);
  { register SYSTEM_int32 _stop = asize - 1;
    if ((k = 0) <=  _stop) do {
      if ((*s1)[k] != (*s2)[k]) 
        return result;
    } while (k++ !=  _stop);

  }
  result = SYSTEM_true;
  return result;
}  /* dataequal */

Function(GMSGEN_pintegerarrayone ) 
  DATASTORAGE_tgamsdatastore_DOT_allocindex(
  DATASTORAGE_tgamsdatastore self)
{
  GMSGEN_pintegerarrayone result;

  result = ValueCast(GMSGEN_pintegerarrayone,
    GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize));
  return result;
}  /* allocindex */

Procedure DATASTORAGE_tgamsdatastore_DOT_freeindex(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone p)
{
  GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,p,self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize);
}  /* freeindex */

Function(SYSTEM_P3_ppointerarray ) 
  DATASTORAGE_tgamsdatastore_DOT_allocptrs(
  DATASTORAGE_tgamsdatastore self)
{
  SYSTEM_P3_ppointerarray result;

  result = ValueCast(SYSTEM_P3_ppointerarray,
    GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,(self->
    DATASTORAGE_tgamsdatastore_DOT_fdimension + 1) * sizeof(
    SYSTEM_pointer)));
  return result;
}  /* allocptrs */

Procedure DATASTORAGE_tgamsdatastore_DOT_freeptrs(
  DATASTORAGE_tgamsdatastore self,
  SYSTEM_P3_ppointerarray p)
{
  GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,p,(self->
    DATASTORAGE_tgamsdatastore_DOT_fdimension + 1) * sizeof(
    SYSTEM_pointer));
}  /* freeptrs */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatastore_DOT_comparekeys(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone key1,
  GMSGEN_pintegerarrayone key2)
{
  SYSTEM_integer result;
  SYSTEM_integer d;

  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      result = (*key1)[d - 1] - (*key2)[d - 1];
      if (result != 0) 
        return result;
    
    } while (d++ !=  _stop);

  }
  result = 0;
  return result;
}  /* comparekeys */

Constructor(DATASTORAGE_tgamsdatastore ) 
  DATASTORAGE_tgamsdatastore_DOT_create(
  DATASTORAGE_tgamsdatastore self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec)
{
  ValueCast(DATASTORAGE_tgamsdatastore,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->DATASTORAGE_tgamsdatastore_DOT_fdimension = adimension;
  self->DATASTORAGE_tgamsdatastore_DOT_fkeysize = adimension * sizeof(
    SYSTEM_longint);
  self->DATASTORAGE_tgamsdatastore_DOT_fdatasize = adatasize;
  self->DATASTORAGE_tgamsdatastore_DOT_ftotalsize = self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize + self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  if (adatasize == 0) { 
    self->DATASTORAGE_tgamsdatastore_DOT_pdefrec = NULL;
  } else {
    self->DATASTORAGE_tgamsdatastore_DOT_pdefrec = ValueCast(
      GMSGEN_pbytedataarray,GMSHEAPNEW_theapmgr_DOT_xgetmem(
      GMSHEAPNEW_gheap,self->DATASTORAGE_tgamsdatastore_DOT_fdatasize));
    GMSOBJ_cmove(adefrec,&(*self->
      DATASTORAGE_tgamsdatastore_DOT_pdefrec)[0],self->
      DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  } 
  return self;
}  /* create */

Destructor(DATASTORAGE_tgamsdatastore ) 
  DATASTORAGE_tgamsdatastore_DOT_destroy(
  DATASTORAGE_tgamsdatastore self)
{
  if (self->DATASTORAGE_tgamsdatastore_DOT_fdatasize > 0) 
    GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,self->
      DATASTORAGE_tgamsdatastore_DOT_pdefrec,self->
      DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatastore_DOT_getnextrecord(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_untyped *adata)
{
  SYSTEM_boolean result;
  SYSTEM_P3_pbyte p;

  p = VirtMethodCall(self, DATASTORAGE_tgamsdatastore_DOT_getnextkey_T, 6, (
    self,akey));
  result = p != NULL;
  if (result) 
    GMSOBJ_cmove(p,adata,self->
      DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  return result;
}  /* getnextrecord */

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatastore_DOT_startread(
  DATASTORAGE_tgamsdatastore self)
{
  SYSTEM_boolean result;

  self->DATASTORAGE_tgamsdatastore_DOT_frecnr = 0;
  result = VirtMethodCall(self, 
    DATASTORAGE_tgamsdatastore_DOT_getcount_T, 1, (self)) > 0;
  return result;
}  /* startread */

Procedure DATASTORAGE_tgamsdatastore_DOT_endread(
  DATASTORAGE_tgamsdatastore self)
{
}  /* endread */

Procedure DATASTORAGE_tgamsdatastore_DOT_startassign(
  DATASTORAGE_tgamsdatastore self)
{
  self->DATASTORAGE_tgamsdatastore_DOT_fasrch = VirtMethodCall(self, 
    DATASTORAGE_tgamsdatastore_DOT_createsearcher_T, 11, (self));
  VirtMethodCall(self->DATASTORAGE_tgamsdatastore_DOT_fasrch, 
    DATASTORAGE_tgamsdatasearcher_DOT_startsearch_T, 1, (self->
    DATASTORAGE_tgamsdatastore_DOT_fasrch));
  self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec = SYSTEM_false;
}  /* startassign */

Procedure DATASTORAGE_tgamsdatastore_DOT_endassign(
  DATASTORAGE_tgamsdatastore self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    DATASTORAGE_tgamsdatastore_DOT_fasrch));
  self->DATASTORAGE_tgamsdatastore_DOT_fasrch = NULL;
}  /* endassign */

Procedure DATASTORAGE_tgamsdatastore_DOT_assignrecord(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_pointer pdata;

  if (VirtMethodCall(self->DATASTORAGE_tgamsdatastore_DOT_fasrch, 
    DATASTORAGE_tgamsdatasearcher_DOT_search_T, 2, (self->
    DATASTORAGE_tgamsdatastore_DOT_fasrch,akey,&pdata))) {
    if (self->DATASTORAGE_tgamsdatastore_DOT_fdatasize > 0 && !
      self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec) 
      self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec = 
        DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(self,adata);
    GMSOBJ_cmove(adata,ValueCast(SYSTEM_P3_pbyte,pdata),self->
      DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  } else 
    if (self->DATASTORAGE_tgamsdatastore_DOT_fdatasize == 0 || !
      DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(self,adata)) 
      VirtMethodCall(self, 
        DATASTORAGE_tgamsdatastore_DOT_insertrecord_T, 3, (self,self->
        DATASTORAGE_tgamsdatastore_DOT_fasrch,akey,adata));
}  /* assignrecord */

static Procedure printkey(
  GMSGEN_pintegerarrayone xkey,
  GMSGEN_pintegerarrayone *_2keys,
  DATASTORAGE_tgamsdatastore *_2self)
{
  SYSTEM_integer d;

  { register SYSTEM_int32 _stop = (*_2self)->
      DATASTORAGE_tgamsdatastore_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      _Iplus_bgn();
      _P3write_i1((**_2keys)[d - 1],6);
      _Iplus_end();
    
    } while (d++ !=  _stop);

  }
  _Iplus_bgn();
  _P3writeln();
  _Iplus_end();
}  /* printkey */

Procedure DATASTORAGE_tgamsdatastore_DOT_verify(
  DATASTORAGE_tgamsdatastore self,
  SYSTEM_boolean print)
{
  GMSGEN_pintegerarrayone keys;
  GMSGEN_pintegerarrayone keysx;
  SYSTEM_integer d;

  keys = DATASTORAGE_tgamsdatastore_DOT_allocindex(self);
  keysx = DATASTORAGE_tgamsdatastore_DOT_allocindex(self);
  if (VirtMethodCall(self, DATASTORAGE_tgamsdatastore_DOT_startread_T, 5, (
    self))) {
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        (*keysx)[d - 1] =  -1;
      } while (d++ !=  _stop);

    }
    while (VirtMethodCall(self, 
      DATASTORAGE_tgamsdatastore_DOT_getnextkey_T, 6, (self,keys)) != NULL) {
      if (print) 
        printkey(keys,&keys,&self);
      if (DATASTORAGE_tgamsdatastore_DOT_comparekeys(self,keysx,keys) >= 0) {
        if (!print) 
          printkey(keys,&keys,&self);
        printkey(keysx,&keys,&self);
        _Iplus_bgn();
        _P3write_s0(_P3str1("\021Keys out of order"));
        _P3writeln();
        _Iplus_end();
      } 
      { register SYSTEM_int32 _stop = self->
          DATASTORAGE_tgamsdatastore_DOT_fdimension;
        if ((d = 1) <=  _stop) do {
          (*keysx)[d - 1] = (*keys)[d - 1];
        } while (d++ !=  _stop);

      }
    
}
  } 
  DATASTORAGE_tgamsdatastore_DOT_freeindex(self,keys);
  DATASTORAGE_tgamsdatastore_DOT_freeindex(self,keysx);
}  /* verify */

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(
  DATASTORAGE_tgamsdatastore self,
  const SYSTEM_untyped *adata)
{
  SYSTEM_boolean result;

  result = DATASTORAGE_dataequal(adata,&(*self->
    DATASTORAGE_tgamsdatastore_DOT_pdefrec)[0],self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  return result;
}  /* isdefaultdata */

Constructor(DATASTORAGE_tgamsdatatable ) 
  DATASTORAGE_tgamsdatatable_DOT_create(
  DATASTORAGE_tgamsdatatable self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec)
{
  ValueCast(DATASTORAGE_tgamsdatatable,
    DATASTORAGE_tgamsdatastore_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatastore,self),adimension,adatasize,adefrec));
  self->DATASTORAGE_tgamsdatatable_DOT_xlist = ValueCast(
    DATASTORAGE_treclist,DATASTORAGE_treclist_DOT_create(ValueCast(
    DATASTORAGE_treclist,_P3alloc_object(&DATASTORAGE_treclist_CD)),
    self->DATASTORAGE_tgamsdatastore_DOT_ftotalsize));
  return self;
}  /* create */

Destructor(DATASTORAGE_tgamsdatatable ) 
  DATASTORAGE_tgamsdatatable_DOT_destroy(
  DATASTORAGE_tgamsdatatable self)
{
  VirtMethodCall(ValueCast(DATASTORAGE_tgamsdatastore,self), 
    DATASTORAGE_tgamsdatastore_DOT_clear_T, 2, (ValueCast(
    DATASTORAGE_tgamsdatastore,self)));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    DATASTORAGE_tgamsdatatable_DOT_xlist));
  DATASTORAGE_tgamsdatastore_DOT_destroy(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  return self;
}  /* destroy */

Procedure DATASTORAGE_tgamsdatatable_DOT_clear(
  DATASTORAGE_tgamsdatatable self)
{
  DATASTORAGE_treclist_DOT_clear(self->
    DATASTORAGE_tgamsdatatable_DOT_xlist);
}  /* clear */

Procedure DATASTORAGE_tgamsdatatable_DOT_loadrecord(
  DATASTORAGE_tgamsdatatable self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  GMSGEN_pbytedataarray p;

  p = ValueCast(GMSGEN_pbytedataarray,DATASTORAGE_treclist_DOT_additem(
    self->DATASTORAGE_tgamsdatatable_DOT_xlist));
  GMSOBJ_cmove(&(*akey)[0],&(*p)[0],self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize);
  GMSOBJ_cmove(adata,&(*p)[self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize],self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize);
}  /* loadrecord */

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatatable_DOT_getnextkey(
  DATASTORAGE_tgamsdatatable self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_P3_pbyte result;
  GMSGEN_pbytedataarray p;

  if (self->DATASTORAGE_tgamsdatastore_DOT_frecnr >= self->
    DATASTORAGE_tgamsdatatable_DOT_xlist->
    DATASTORAGE_treclist_DOT_fcount) { 
    result = NULL;
  } else {
    p = ValueCast(GMSGEN_pbytedataarray,
      DATASTORAGE_treclist_DOT_getitem(self->
      DATASTORAGE_tgamsdatatable_DOT_xlist,self->
      DATASTORAGE_tgamsdatastore_DOT_frecnr));
    self->DATASTORAGE_tgamsdatastore_DOT_frecnr = self->
      DATASTORAGE_tgamsdatastore_DOT_frecnr + 1;
    GMSOBJ_cmove(&(*p)[0],&(*akey)[0],self->
      DATASTORAGE_tgamsdatastore_DOT_fkeysize);
    result = ValueCast(SYSTEM_P3_pbyte,&(*p)[self->
      DATASTORAGE_tgamsdatastore_DOT_fkeysize]);
  } 
  return result;
}  /* getnextkey */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatatable_DOT_memoryused(
  DATASTORAGE_tgamsdatatable self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamsdatatable_DOT_xlist->
    DATASTORAGE_treclist_DOT_fcapacity * sizeof(SYSTEM_longint) + self->
    DATASTORAGE_tgamsdatatable_DOT_xlist->
    DATASTORAGE_treclist_DOT_fcount * self->
    DATASTORAGE_tgamsdatastore_DOT_ftotalsize;
  return result;
}  /* memoryused */

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatatable_DOT_createsearcher(
  DATASTORAGE_tgamsdatatable self)
{
  DATASTORAGE_tgamsdatasearcher result;

  result = ValueCast(DATASTORAGE_tgamsdatasearcher,
    DATASTORAGE_tgamsdatasearcher_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatasearcher,_P3alloc_object(&
    DATASTORAGE_tgamsdatatablesearcher_CD)),ValueCast(
    DATASTORAGE_tgamsdatastore,self)));
  return result;
}  /* createsearcher */

Procedure DATASTORAGE_tgamsdatatable_DOT_insertrecord(
  DATASTORAGE_tgamsdatatable self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_integer l;
  GMSGEN_pbytedataarray p;

  l = (ValueCast(DATASTORAGE_tgamsdatatablesearcher,ads))->
    DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex;
  p = ValueCast(GMSGEN_pbytedataarray,DATASTORAGE_treclist_DOT_insert(
    self->DATASTORAGE_tgamsdatatable_DOT_xlist,l));
  GMSOBJ_cmove(&(*akey)[0],&(*p)[0],self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize);
  GMSOBJ_cmove(adata,&(*p)[self->
    DATASTORAGE_tgamsdatastore_DOT_fkeysize],self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize);
}  /* insertrecord */

Procedure DATASTORAGE_tgamsdatatable_DOT_endassign(
  DATASTORAGE_tgamsdatatable self)
{
  SYSTEM_integer n;
  GMSGEN_pbytedataarray pr;

  DATASTORAGE_tgamsdatastore_DOT_endassign(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  if (self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec) {
    self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec = SYSTEM_false;
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatatable_DOT_xlist->
        DATASTORAGE_treclist_DOT_fcount - 1;
      if ((n = 0) <=  _stop) do {
        pr = ValueCast(GMSGEN_pbytedataarray,
          DATASTORAGE_treclist_DOT_getitem(self->
          DATASTORAGE_tgamsdatatable_DOT_xlist,n));
        if (DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(ValueCast(
          DATASTORAGE_tgamsdatastore,self),&(*pr)[self->
          DATASTORAGE_tgamsdatastore_DOT_fkeysize])) {
          DATASTORAGE_treclist_DOT_remove(self->
            DATASTORAGE_tgamsdatatable_DOT_xlist,n);
          self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec = 
            SYSTEM_true;
        } 
      
      } while (n++ !=  _stop);

    }
    if (self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec) 
      DATASTORAGE_treclist_DOT_cleanup(self->
        DATASTORAGE_tgamsdatatable_DOT_xlist);
  } 
}  /* endassign */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatatable_DOT_getcount(
  DATASTORAGE_tgamsdatatable self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamsdatatable_DOT_xlist->
    DATASTORAGE_treclist_DOT_fcount;
  return result;
}  /* getcount */

static Procedure quicksort(
  SYSTEM_integer l,
  SYSTEM_integer r,
  DATASTORAGE_tgamsdatatable *_2self)
{
  SYSTEM_integer i, j, p;
  GMSGEN_pintegerarrayone prec;

  do {
    i = l;
    j = r;
    p = ValueCast(SYSTEM_uint32,l + r) >> 1;
    do {
      prec = ValueCast(GMSGEN_pintegerarrayone,
        DATASTORAGE_treclist_DOT_getitem((*_2self)->
        DATASTORAGE_tgamsdatatable_DOT_xlist,p));
      while (DATASTORAGE_tgamsdatastore_DOT_comparekeys(ValueCast(
        DATASTORAGE_tgamsdatastore,*_2self),ValueCast(
        GMSGEN_pintegerarrayone,DATASTORAGE_treclist_DOT_getitem((*
        _2self)->DATASTORAGE_tgamsdatatable_DOT_xlist,i)),prec) < 0) {

        _P3inc0(i);
}
      while (DATASTORAGE_tgamsdatastore_DOT_comparekeys(ValueCast(
        DATASTORAGE_tgamsdatastore,*_2self),ValueCast(
        GMSGEN_pintegerarrayone,DATASTORAGE_treclist_DOT_getitem((*
        _2self)->DATASTORAGE_tgamsdatatable_DOT_xlist,j)),prec) > 0) {

        _P3dec0(j);
}
      if (i <= j) {
        DATASTORAGE_treclist_DOT_exchange((*_2self)->
          DATASTORAGE_tgamsdatatable_DOT_xlist,i,j);
        if (p == i) { 
          p = j;
        } else 
          if (p == j) 
            p = i;
        _P3inc0(i);
        _P3dec0(j);
      } 
    } while (!(i > j));
    if (l < j) 
      quicksort(l,j,_2self);
    l = i;
  } while (!(i >= r));
}  /* quicksort */

Procedure DATASTORAGE_tgamsdatatable_DOT_sort(
  DATASTORAGE_tgamsdatatable self)
{
  if (VirtMethodCall(ValueCast(DATASTORAGE_tgamsdatastore,self), 
    DATASTORAGE_tgamsdatastore_DOT_getcount_T, 1, (ValueCast(
    DATASTORAGE_tgamsdatastore,self))) >= 2) 
    quicksort(0,VirtMethodCall(ValueCast(
      DATASTORAGE_tgamsdatastore,self), 
      DATASTORAGE_tgamsdatastore_DOT_getcount_T, 1, (ValueCast(
      DATASTORAGE_tgamsdatastore,self))) - 1,&self);
}  /* sort */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatatable_DOT_getcapacity(
  DATASTORAGE_tgamsdatatable self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamsdatatable_DOT_xlist->
    DATASTORAGE_treclist_DOT_fcapacity;
  return result;
}  /* getcapacity */

Procedure DATASTORAGE_tgamsdatatable_DOT_setcapacity(
  DATASTORAGE_tgamsdatatable self,
  SYSTEM_integer n)
{
  DATASTORAGE_treclist_DOT_setcapacity(self->
    DATASTORAGE_tgamsdatatable_DOT_xlist,n);
}  /* setcapacity */

Constructor(DATASTORAGE_tgamsdatafull ) 
  DATASTORAGE_tgamsdatafull_DOT_create(
  DATASTORAGE_tgamsdatafull self,
  DATASTORAGE_tgamsdatastore ds)
{
  SYSTEM_integer d;
  GMSGEN_pintegerarrayone key;
  SYSTEM_integer n;
  SYSTEM_integer nr;

  ValueCast(DATASTORAGE_tgamsdatafull,
    DATASTORAGE_tgamsdatastore_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatastore,self),ds->
    DATASTORAGE_tgamsdatastore_DOT_fdimension,ds->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize,&(*ds->
    DATASTORAGE_tgamsdatastore_DOT_pdefrec)[0]));
  self->DATASTORAGE_tgamsdatafull_DOT_xmin = 
    DATASTORAGE_tgamsdatastore_DOT_allocindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  self->DATASTORAGE_tgamsdatafull_DOT_xmult = 
    DATASTORAGE_tgamsdatastore_DOT_allocindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  self->DATASTORAGE_tgamsdatafull_DOT_xcnt = 
    DATASTORAGE_tgamsdatastore_DOT_allocindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  self->DATASTORAGE_tgamsdatafull_DOT_finxmap = ValueCast(
    DATASTORAGE_pintegermappingarray,GMSHEAPNEW_theapmgr_DOT_xgetmem(
    GMSHEAPNEW_gheap,self->DATASTORAGE_tgamsdatastore_DOT_fdimension * sizeof(
    DATASTORAGE_tintegermappingarray)));
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      (*self->DATASTORAGE_tgamsdatafull_DOT_finxmap)[d - 1] = ValueCast(
        DATASTORAGE_tintegermapping,
        DATASTORAGE_tintegermapping_DOT_create(ValueCast(
        DATASTORAGE_tintegermapping,_P3alloc_object(&
        DATASTORAGE_tintegermapping_CD))));
    } while (d++ !=  _stop);

  }
  key = DATASTORAGE_tgamsdatastore_DOT_allocindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  if (VirtMethodCall(ds, DATASTORAGE_tgamsdatastore_DOT_startread_T, 5, (
    ds))) {
    while (VirtMethodCall(ds, 
      DATASTORAGE_tgamsdatastore_DOT_getnextkey_T, 6, (ds,key)) != NULL) {

      { register SYSTEM_int32 _stop = self->
          DATASTORAGE_tgamsdatastore_DOT_fdimension;
        if ((d = 1) <=  _stop) do {
          DATASTORAGE_tintegermapping_DOT_setmapping((*self->
            DATASTORAGE_tgamsdatafull_DOT_finxmap)[d - 1],(*key)[
            d - 1],0);
        } while (d++ !=  _stop);

      }
}
    VirtMethodCall(ds, DATASTORAGE_tgamsdatastore_DOT_endread_T, 7, (ds));
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        nr =  -1;
        { register DATASTORAGE_tintegermapping_OD *_W2=(*self->
          DATASTORAGE_tgamsdatafull_DOT_finxmap)[d - 1];
          { register SYSTEM_int32 _stop = _W2->
              DATASTORAGE_tintegermapping_DOT_fhighestindex;
            if ((n = 0) <=  _stop) do {
              if ((*_W2->DATASTORAGE_tintegermapping_DOT_pmap)[n] == 0) {
                if (nr ==  -1) 
                  (*self->DATASTORAGE_tgamsdatafull_DOT_xmin)[d - 1] = 
                    n;
                nr = nr + 1;
                (*_W2->DATASTORAGE_tintegermapping_DOT_pmap)[n] = nr;
              } 
            } while (n++ !=  _stop);

          }

        }
        if (nr == 0) { 
          (*self->DATASTORAGE_tgamsdatafull_DOT_xcnt)[d - 1] = 0;
        } else 
          (*self->DATASTORAGE_tgamsdatafull_DOT_xcnt)[d - 1] = nr + 1;
      
      } while (d++ !=  _stop);

    }
  } 
  self->DATASTORAGE_tgamsdatafull_DOT_fallocsize = self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  (*self->DATASTORAGE_tgamsdatafull_DOT_xmult)[self->
    DATASTORAGE_tgamsdatastore_DOT_fdimension - 1] = self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  for (d = self->DATASTORAGE_tgamsdatastore_DOT_fdimension - 1;d >= (
    SYSTEM_int32)1;--d) {
    if ((*self->DATASTORAGE_tgamsdatafull_DOT_xcnt)[d + 1 - 1] != 0) 
      self->DATASTORAGE_tgamsdatafull_DOT_fallocsize = self->
        DATASTORAGE_tgamsdatafull_DOT_fallocsize * (*self->
        DATASTORAGE_tgamsdatafull_DOT_xcnt)[d + 1 - 1];
    (*self->DATASTORAGE_tgamsdatafull_DOT_xmult)[d - 1] = self->
      DATASTORAGE_tgamsdatafull_DOT_fallocsize;
  
  }
  if ((*self->DATASTORAGE_tgamsdatafull_DOT_xcnt)[0] != 0) 
    self->DATASTORAGE_tgamsdatafull_DOT_fallocsize = self->
      DATASTORAGE_tgamsdatafull_DOT_fallocsize * (*self->
      DATASTORAGE_tgamsdatafull_DOT_xcnt)[0];
  DATASTORAGE_tgamsdatastore_DOT_freeindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self),key);
  self->DATASTORAGE_tgamsdatafull_DOT_pdata = NULL;
  return self;
}  /* create */

Destructor(DATASTORAGE_tgamsdatafull ) 
  DATASTORAGE_tgamsdatafull_DOT_destroy(
  DATASTORAGE_tgamsdatafull self)
{
  SYSTEM_integer d;

  GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,self->
    DATASTORAGE_tgamsdatafull_DOT_pdata,self->
    DATASTORAGE_tgamsdatafull_DOT_fallocsize);
  DATASTORAGE_tgamsdatastore_DOT_freeindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self),self->
    DATASTORAGE_tgamsdatafull_DOT_xcnt);
  DATASTORAGE_tgamsdatastore_DOT_freeindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self),self->
    DATASTORAGE_tgamsdatafull_DOT_xmult);
  DATASTORAGE_tgamsdatastore_DOT_freeindex(ValueCast(
    DATASTORAGE_tgamsdatastore,self),self->
    DATASTORAGE_tgamsdatafull_DOT_xmin);
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,(*self->
        DATASTORAGE_tgamsdatafull_DOT_finxmap)[d - 1]));
    } while (d++ !=  _stop);

  }
  DATASTORAGE_tgamsdatastore_DOT_destroy(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  return self;
}  /* destroy */

Procedure DATASTORAGE_tgamsdatafull_DOT_loadrecord(
  DATASTORAGE_tgamsdatafull self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_P3_pbyte p;

  p = DATASTORAGE_tgamsdatafull_DOT_getoffset(self,akey);
  if (p != NULL) {
    GMSOBJ_cmove(adata,DATASTORAGE_tgamsdatafull_DOT_getoffset(self,
      akey),self->DATASTORAGE_tgamsdatastore_DOT_fdatasize);
    _P3inc0(self->DATASTORAGE_tgamsdatafull_DOT_fcount);
  } 
}  /* loadrecord */

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatafull_DOT_getoffset(
  DATASTORAGE_tgamsdatafull self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_P3_pbyte result;
  SYSTEM_integer d;
  SYSTEM_integer k;

  result = self->DATASTORAGE_tgamsdatafull_DOT_pdata;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      k = DATASTORAGE_tintegermapping_DOT_getmapping((*self->
        DATASTORAGE_tgamsdatafull_DOT_finxmap)[d - 1],(*akey)[d - 1]);
      if (k < 0) {
        result = NULL;
        return result;
      } 
      _P3inc1(result,k * (*self->DATASTORAGE_tgamsdatafull_DOT_xmult)[
        d - 1]);
    
    } while (d++ !=  _stop);

  }
  return result;
}  /* getoffset */

Procedure DATASTORAGE_tgamsdatafull_DOT_allocatememory(
  DATASTORAGE_tgamsdatafull self,
  DATASTORAGE_tgamsdatastore ds)
{
  GMSGEN_pintegerarrayone key;
  GMSGEN_pbytedataarray data;

  if (self->DATASTORAGE_tgamsdatafull_DOT_pdata == NULL) {
    self->DATASTORAGE_tgamsdatafull_DOT_pdata = ValueCast(
      SYSTEM_P3_pbyte,GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,
      self->DATASTORAGE_tgamsdatafull_DOT_fallocsize));
    VirtMethodCall(ValueCast(DATASTORAGE_tgamsdatastore,self), 
      DATASTORAGE_tgamsdatastore_DOT_clear_T, 2, (ValueCast(
      DATASTORAGE_tgamsdatastore,self)));
    if (VirtMethodCall(ds, DATASTORAGE_tgamsdatastore_DOT_startread_T, 5, (
      ds))) {
      key = DATASTORAGE_tgamsdatastore_DOT_allocindex(ValueCast(
        DATASTORAGE_tgamsdatastore,self));
      data = ValueCast(GMSGEN_pbytedataarray,
        GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,self->
        DATASTORAGE_tgamsdatastore_DOT_fdatasize));
      while (DATASTORAGE_tgamsdatastore_DOT_getnextrecord(ds,key,&(*
        data)[0])) {

        VirtMethodCall(ValueCast(DATASTORAGE_tgamsdatastore,self), 
          DATASTORAGE_tgamsdatastore_DOT_loadrecord_T, 4, (ValueCast(
          DATASTORAGE_tgamsdatastore,self),key,&(*data)[0]));
}
      VirtMethodCall(ds, DATASTORAGE_tgamsdatastore_DOT_endread_T, 7, (
        ds));
      DATASTORAGE_tgamsdatastore_DOT_freeindex(ValueCast(
        DATASTORAGE_tgamsdatastore,self),key);
      GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,data,self->
        DATASTORAGE_tgamsdatastore_DOT_fdatasize);
    } 
  } 
}  /* allocatememory */

Procedure DATASTORAGE_tgamsdatafull_DOT_clear(
  DATASTORAGE_tgamsdatafull self)
{
  SYSTEM_integer cnt;
  SYSTEM_P3_pbyte p;

  p = self->DATASTORAGE_tgamsdatafull_DOT_pdata;
  cnt = 0;
  while (cnt < self->DATASTORAGE_tgamsdatafull_DOT_fallocsize) {
    GMSOBJ_cmove(&(*self->DATASTORAGE_tgamsdatastore_DOT_pdefrec)[0],
      p,self->DATASTORAGE_tgamsdatastore_DOT_fdatasize);
    _P3inc1(p,self->DATASTORAGE_tgamsdatastore_DOT_fdatasize);
    _P3inc1(cnt,self->DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  
}
}  /* clear */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatafull_DOT_memoryused(
  DATASTORAGE_tgamsdatafull self)
{
  SYSTEM_integer result;
  SYSTEM_integer d;

  result = self->DATASTORAGE_tgamsdatafull_DOT_fallocsize;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      _P3inc1(result,DATASTORAGE_tintegermapping_DOT_memoryused((*self->
        DATASTORAGE_tgamsdatafull_DOT_finxmap)[d - 1]));
    } while (d++ !=  _stop);

  }
  return result;
}  /* memoryused */

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatafull_DOT_createsearcher(
  DATASTORAGE_tgamsdatafull self)
{
  DATASTORAGE_tgamsdatasearcher result;

  result = ValueCast(DATASTORAGE_tgamsdatasearcher,
    DATASTORAGE_tgamsdatasearcher_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatasearcher,_P3alloc_object(&
    DATASTORAGE_tgamsdatafullsearcher_CD)),ValueCast(
    DATASTORAGE_tgamsdatastore,self)));
  return result;
}  /* createsearcher */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatafull_DOT_getcount(
  DATASTORAGE_tgamsdatafull self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamsdatafull_DOT_fcount;
  return result;
}  /* getcount */

Constructor(DATASTORAGE_tintegermapping ) 
  DATASTORAGE_tintegermapping_DOT_create(
  DATASTORAGE_tintegermapping self)
{
  ValueCast(DATASTORAGE_tintegermapping,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->DATASTORAGE_tintegermapping_DOT_fcapacity = 0;
  self->DATASTORAGE_tintegermapping_DOT_fhighestindex = 0;
  self->DATASTORAGE_tintegermapping_DOT_pmap = NULL;
  return self;
}  /* create */

Destructor(DATASTORAGE_tintegermapping ) 
  DATASTORAGE_tintegermapping_DOT_destroy(
  DATASTORAGE_tintegermapping self)
{
  if (self->DATASTORAGE_tintegermapping_DOT_pmap != NULL) 
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      DATASTORAGE_tintegermapping_DOT_pmap),0);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_integer ) DATASTORAGE_tintegermapping_DOT_getmapping(
  DATASTORAGE_tintegermapping self,
  SYSTEM_integer f)
{
  SYSTEM_integer result;

  if (f >= 0 && f < self->
    DATASTORAGE_tintegermapping_DOT_fcapacity) { 
    result = (*self->DATASTORAGE_tintegermapping_DOT_pmap)[f];
  } else 
    result =  -1;
  return result;
}  /* getmapping */

Function(SYSTEM_integer ) DATASTORAGE_tintegermapping_DOT_memoryused(
  DATASTORAGE_tintegermapping self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tintegermapping_DOT_fcapacity * sizeof(
    SYSTEM_longint);
  return result;
}  /* memoryused */

Procedure DATASTORAGE_tintegermapping_DOT_setmapping(
  DATASTORAGE_tintegermapping self,
  SYSTEM_integer f,
  SYSTEM_integer t)
{
  SYSTEM_integer n;
  SYSTEM_integer delta;

  if (f >= self->DATASTORAGE_tintegermapping_DOT_fcapacity) {
    delta = 0;
    do {
      if (self->DATASTORAGE_tintegermapping_DOT_fcapacity == 0) { 
        _P3inc1(delta,1024);
      } else 
        if (self->DATASTORAGE_tintegermapping_DOT_fcapacity <= 32768) { 
          _P3inc1(delta,self->
            DATASTORAGE_tintegermapping_DOT_fcapacity);
        } else 
          _P3inc1(delta,self->
            DATASTORAGE_tintegermapping_DOT_fcapacity /  4);
    } while (!(f < self->DATASTORAGE_tintegermapping_DOT_fcapacity + 
      delta));
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      DATASTORAGE_tintegermapping_DOT_pmap),(self->
      DATASTORAGE_tintegermapping_DOT_fcapacity + delta) * sizeof(
      SYSTEM_longint));
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tintegermapping_DOT_fcapacity + delta - 1;
      if ((n = self->DATASTORAGE_tintegermapping_DOT_fcapacity) <=  _stop) do {
        (*self->DATASTORAGE_tintegermapping_DOT_pmap)[n] =  -1;
      } while (n++ !=  _stop);

    }
    self->DATASTORAGE_tintegermapping_DOT_fcapacity = self->
      DATASTORAGE_tintegermapping_DOT_fcapacity + delta;
  } 
  (*self->DATASTORAGE_tintegermapping_DOT_pmap)[f] = t;
  if (f > self->DATASTORAGE_tintegermapping_DOT_fhighestindex) 
    self->DATASTORAGE_tintegermapping_DOT_fhighestindex = f;
}  /* setmapping */
typedef struct DATASTORAGE_tsparsecell_S *DATASTORAGE_psparsecell;
typedef struct DATASTORAGE_tsparsecell_S {
  DATASTORAGE_psparsecell pcelldown;
  SYSTEM_integer cellkey;
  DATASTORAGE_psparsecell pcellright;
} DATASTORAGE_tsparsecell;

typedef struct DATASTORAGE_tsparsedatacell_S *
  DATASTORAGE_psparsedatacell;
typedef struct DATASTORAGE_tsparsedatacell_S {
  DATASTORAGE_psparsecell pcelldown;
  SYSTEM_integer cellkey;
  SYSTEM_byte celldata;
} DATASTORAGE_tsparsedatacell;

cnstdef {DATASTORAGE_szcell = sizeof(DATASTORAGE_tsparsecell)};
cnstdef {DATASTORAGE_szdatacell = sizeof(DATASTORAGE_tsparsecell) - sizeof(
  SYSTEM_pointer)};

Constructor(DATASTORAGE_tgamsdatasparse ) 
  DATASTORAGE_tgamsdatasparse_DOT_create(
  DATASTORAGE_tgamsdatasparse self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec)
{
  ValueCast(DATASTORAGE_tgamsdatasparse,
    DATASTORAGE_tgamsdatastore_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatastore,self),adimension,adatasize,adefrec));
  self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs = 
    DATASTORAGE_tgamsdatastore_DOT_allocptrs(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0] = 
    DATASTORAGE_tgamsdatasparse_DOT_getcell(self,0);
  (ValueCast(DATASTORAGE_psparsecell,(*self->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0]))->pcellright = NULL;
  self->DATASTORAGE_tgamsdatasparse_DOT_fcellcount = 0;
  return self;
}  /* create */

Destructor(DATASTORAGE_tgamsdatasparse ) 
  DATASTORAGE_tgamsdatasparse_DOT_destroy(
  DATASTORAGE_tgamsdatasparse self)
{
  VirtMethodCall(ValueCast(DATASTORAGE_tgamsdatastore,self), 
    DATASTORAGE_tgamsdatastore_DOT_clear_T, 2, (ValueCast(
    DATASTORAGE_tgamsdatastore,self)));
  DATASTORAGE_tgamsdatasparse_DOT_freecell(self,(*self->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0],0);
  DATASTORAGE_tgamsdatastore_DOT_freeptrs(ValueCast(
    DATASTORAGE_tgamsdatastore,self),self->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs);
  DATASTORAGE_tgamsdatastore_DOT_destroy(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  return self;
}  /* destroy */

static Procedure freenode(
  DATASTORAGE_psparsecell p,
  SYSTEM_integer d,
  DATASTORAGE_tgamsdatasparse *_2self)
{
  DATASTORAGE_psparsecell pn;

  while (p != NULL) {
    if (d < (*_2self)->DATASTORAGE_tgamsdatastore_DOT_fdimension) 
      freenode(p->pcellright,d + 1,_2self);
    pn = p->pcelldown;
    DATASTORAGE_tgamsdatasparse_DOT_freecell(*_2self,p,d);
    p = pn;
  
}
}  /* freenode */

Procedure DATASTORAGE_tgamsdatasparse_DOT_clear(
  DATASTORAGE_tgamsdatasparse self)
{
  freenode((ValueCast(DATASTORAGE_psparsecell,(*self->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0]))->pcellright,1,&
    self);
  (ValueCast(DATASTORAGE_psparsecell,(*self->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0]))->pcellright = NULL;
}  /* clear */

Procedure DATASTORAGE_tgamsdatasparse_DOT_loadrecord(
  DATASTORAGE_tgamsdatasparse self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_integer d, d2;
  DATASTORAGE_psparsecell cp;

  if (self->DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount == 0) { 
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d] = 
          DATASTORAGE_tgamsdatasparse_DOT_getcell(self,d);
        (ValueCast(DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->cellkey = (*
          akey)[d - 1];
        if (d < self->DATASTORAGE_tgamsdatastore_DOT_fdimension) 
          (ValueCast(DATASTORAGE_psparsecell,(*self->
            DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->pcelldown = NULL;
        (ValueCast(DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d - 1]))->
          pcellright = ValueCast(DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]);
      
      } while (d++ !=  _stop);

    }
  } else 
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        if ((ValueCast(DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->cellkey != (*
          akey)[d - 1]) {
          cp = ValueCast(DATASTORAGE_psparsecell,
            DATASTORAGE_tgamsdatasparse_DOT_getcell(self,d));
          cp->cellkey = (*akey)[d - 1];
          cp->pcelldown = NULL;
          (ValueCast(DATASTORAGE_psparsecell,(*self->
            DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->pcelldown = 
            cp;
          (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d] = cp;
          { register SYSTEM_int32 _stop = self->
              DATASTORAGE_tgamsdatastore_DOT_fdimension;
            if ((d2 = d + 1) <=  _stop) do {
              cp = ValueCast(DATASTORAGE_psparsecell,
                DATASTORAGE_tgamsdatasparse_DOT_getcell(self,d2));
              cp->cellkey = (*akey)[d2 - 1];
              cp->pcelldown = NULL;
              (ValueCast(DATASTORAGE_psparsecell,(*self->
                DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d2 - 1]))->
                pcellright = cp;
              (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d2] = 
                cp;
            
            } while (d2++ !=  _stop);

          }
          SYSTEM_break(BRK_1);
        } 
CNT_1:;
      } while (d++ !=  _stop);
BRK_1:;

    }
  GMSOBJ_cmove(adata,&(ValueCast(DATASTORAGE_psparsedatacell,(*self->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[self->
    DATASTORAGE_tgamsdatastore_DOT_fdimension]))->celldata,self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize);
}  /* loadrecord */

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatasparse_DOT_startread(
  DATASTORAGE_tgamsdatasparse self)
{
  SYSTEM_boolean result;
  SYSTEM_integer d;

  result = DATASTORAGE_tgamsdatastore_DOT_startread(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  if (result) 
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d] = (ValueCast(
          DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d - 1]))->
          pcellright;
      } while (d++ !=  _stop);

    }
  return result;
}  /* startread */

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatasparse_DOT_getnextkey(
  DATASTORAGE_tgamsdatasparse self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_P3_pbyte result;
  SYSTEM_integer d;
  SYSTEM_integer d2;

  if ((*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[1] == NULL) { 
    result = NULL;
  } else {
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        (*akey)[d - 1] = (ValueCast(DATASTORAGE_psparsecell,(*
          self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->
          cellkey;
      } while (d++ !=  _stop);

    }
    result = ValueCast(SYSTEM_P3_pbyte,&(ValueCast(
      DATASTORAGE_psparsedatacell,(*self->
      DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension]))->celldata);
    (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension] = (ValueCast(
      DATASTORAGE_psparsecell,(*self->
      DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension]))->pcelldown;
    if ((*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[self->
      DATASTORAGE_tgamsdatastore_DOT_fdimension] != NULL) 
      return result;
    d2 = 0;
    for (d = self->DATASTORAGE_tgamsdatastore_DOT_fdimension - 1;
      d >= (SYSTEM_int32)1;--d) {
      (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d] = (ValueCast(
        DATASTORAGE_psparsecell,(*self->
        DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->pcelldown;
      if ((*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d] != NULL) {
        d2 = d;
        SYSTEM_break(BRK_2);
      } 
    
    CNT_2:;
}
BRK_2:;
    if (d2 >= 1) 
      { register SYSTEM_int32 _stop = self->
          DATASTORAGE_tgamsdatastore_DOT_fdimension - 1;
        if ((d = d2) <=  _stop) do {
          (*self->DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d + 1] = (ValueCast(
            DATASTORAGE_psparsecell,(*self->
            DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[d]))->pcellright;
        } while (d++ !=  _stop);

      }
  } 
  return result;
}  /* getnextkey */

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatasparse_DOT_createsearcher(
  DATASTORAGE_tgamsdatasparse self)
{
  DATASTORAGE_tgamsdatasearcher result;

  result = ValueCast(DATASTORAGE_tgamsdatasearcher,
    DATASTORAGE_tgamsdatasparsesearcher_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatasparsesearcher,_P3alloc_object(&
    DATASTORAGE_tgamsdatasparsesearcher_CD)),self));
  return result;
}  /* createsearcher */

Procedure DATASTORAGE_tgamsdatasparse_DOT_insertrecord(
  DATASTORAGE_tgamsdatasparse self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_integer d, d2;
  DATASTORAGE_psparsecell pc, pn;
  DATASTORAGE_psparsecell pdown;

  { register DATASTORAGE_tgamsdatasparsesearcher_OD *_W2=ValueCast(
    DATASTORAGE_tgamsdatasparsesearcher,ads);
    if (self->DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount == 0) { 
      { register SYSTEM_int32 _stop = _W2->
          DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim;
        if ((d = 1) <=  _stop) do {
          (*_W2->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[
            d] = DATASTORAGE_tgamsdatasparse_DOT_getcell(self,d);
          (ValueCast(DATASTORAGE_psparsecell,(*_W2->
            DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d - 1]))->
            pcellright = ValueCast(DATASTORAGE_psparsecell,(*_W2->
            DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d]);
          { register DATASTORAGE_tsparsecell *_W3=ValueCast(
            DATASTORAGE_psparsecell,(*_W2->
            DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d]);
            _W3->cellkey = (*akey)[d - 1];
            _W3->pcelldown = NULL;

          }
        
        } while (d++ !=  _stop);

      }
    } else 
      { register SYSTEM_int32 _stop = self->
          DATASTORAGE_tgamsdatastore_DOT_fdimension;
        if ((d = _W2->
          DATASTORAGE_tgamsdatasparsesearcher_DOT_flastvalid + 1) <=  _stop) do {
          pc = ValueCast(DATASTORAGE_psparsecell,(*_W2->
            DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d]);
          pdown = pc->pcelldown;
          pn = ValueCast(DATASTORAGE_psparsecell,
            DATASTORAGE_tgamsdatasparse_DOT_getcell(self,d));
          pn->cellkey = (*akey)[d - 1];
          if ((*akey)[d - 1] < pc->cellkey) {
            pn->pcelldown = pc;
            (ValueCast(DATASTORAGE_psparsecell,(*_W2->
              DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d - 1]))->
              pcellright = pn;
          } else {
            pn->pcelldown = pdown;
            pc->pcelldown = pn;
          } 
          (*_W2->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[
            d] = pn;
          { register SYSTEM_int32 _stop = self->
              DATASTORAGE_tgamsdatastore_DOT_fdimension;
            if ((d2 = d + 1) <=  _stop) do {
              pn = ValueCast(DATASTORAGE_psparsecell,
                DATASTORAGE_tgamsdatasparse_DOT_getcell(self,d2));
              { register DATASTORAGE_tsparsecell *_W3=pn;
                _W3->cellkey = (*akey)[d2 - 1];
                _W3->pcelldown = NULL;

              }
              (*_W2->
                DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[
                d2] = pn;
              (ValueCast(DATASTORAGE_psparsecell,(*_W2->
                DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[
                d2 - 1]))->pcellright = pn;
            
            } while (d2++ !=  _stop);

          }
          SYSTEM_break(BRK_3);
        
CNT_3:;
        } while (d++ !=  _stop);
BRK_3:;

      }
    if (self->DATASTORAGE_tgamsdatastore_DOT_fdatasize > 0) 
      GMSOBJ_cmove(adata,&(ValueCast(DATASTORAGE_psparsedatacell,(*_W2->
        DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[self->
        DATASTORAGE_tgamsdatastore_DOT_fdimension]))->celldata,self->
        DATASTORAGE_tgamsdatastore_DOT_fdatasize);

  }
}  /* insertrecord */

static Function(DATASTORAGE_psparsecell ) cleanup(
  DATASTORAGE_psparsecell p,
  SYSTEM_integer d,
  DATASTORAGE_tgamsdatasparse *_2self)
{
  DATASTORAGE_psparsecell result;
  DATASTORAGE_psparsecell pr;
  DATASTORAGE_psparsecell pp;
  DATASTORAGE_psparsecell ps;
  DATASTORAGE_psparsecell pdown;

  result = NULL;
  ps = p;
  while (ps != NULL) {
    pdown = ps->pcelldown;
    if (d < (*_2self)->DATASTORAGE_tgamsdatastore_DOT_fdimension) { 
      pr = cleanup(ps->pcellright,d + 1,_2self);
    } else 
      if (DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(ValueCast(
        DATASTORAGE_tgamsdatastore,*_2self),&(ValueCast(
        DATASTORAGE_psparsedatacell,ps))->celldata)) { 
        pr = NULL;
      } else 
        pr = ps;
    if (pr == NULL) { 
      DATASTORAGE_tgamsdatasparse_DOT_freecell(*_2self,ps,d);
    } else {
      if (result == NULL) {
        result = ps;
        pp = result;
      } else {
        pp->pcelldown = ps;
        pp = ps;
      } 
      if (d < (*_2self)->DATASTORAGE_tgamsdatastore_DOT_fdimension) 
        ps->pcellright = pr;
    } 
    ps = pdown;
  
}
  if (result != NULL) 
    pp->pcelldown = NULL;
  return result;
}  /* cleanup */

Procedure DATASTORAGE_tgamsdatasparse_DOT_endassign(
  DATASTORAGE_tgamsdatasparse self)
{
  DATASTORAGE_tgamsdatastore_DOT_endassign(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  if (self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec) 
    (ValueCast(DATASTORAGE_psparsecell,(*self->
      DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0]))->pcellright = 
      cleanup((ValueCast(DATASTORAGE_psparsecell,(*self->
      DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0]))->pcellright,1,&
      self);
}  /* endassign */

Function(SYSTEM_pointer ) DATASTORAGE_tgamsdatasparse_DOT_getcell(
  DATASTORAGE_tgamsdatasparse self,
  SYSTEM_integer d)
{
  SYSTEM_pointer result;

  if (d < self->DATASTORAGE_tgamsdatastore_DOT_fdimension) {
    self->DATASTORAGE_tgamsdatasparse_DOT_fcellcount = self->
      DATASTORAGE_tgamsdatasparse_DOT_fcellcount + 1;
    result = GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,
      DATASTORAGE_szcell);
  } else {
    self->DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount = self->
      DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount + 1;
    result = GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,
      DATASTORAGE_szdatacell + self->
      DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  } 
  return result;
}  /* getcell */

Procedure DATASTORAGE_tgamsdatasparse_DOT_freecell(
  DATASTORAGE_tgamsdatasparse self,
  SYSTEM_pointer p,
  SYSTEM_integer d)
{
  if (d < self->DATASTORAGE_tgamsdatastore_DOT_fdimension) {
    GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,p,
      DATASTORAGE_szcell);
    self->DATASTORAGE_tgamsdatasparse_DOT_fcellcount = self->
      DATASTORAGE_tgamsdatasparse_DOT_fcellcount - 1;
  } else {
    GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,p,
      DATASTORAGE_szdatacell + self->
      DATASTORAGE_tgamsdatastore_DOT_fdatasize);
    self->DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount = self->
      DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount - 1;
  } 
}  /* freecell */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatasparse_DOT_memoryused(
  DATASTORAGE_tgamsdatasparse self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamsdatasparse_DOT_fcellcount * 
    DATASTORAGE_szcell + self->
    DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount * (
    DATASTORAGE_szdatacell + self->
    DATASTORAGE_tgamsdatastore_DOT_fdatasize);
  return result;
}  /* memoryused */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatasparse_DOT_getcount(
  DATASTORAGE_tgamsdatasparse self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount;
  return result;
}  /* getcount */
typedef struct DATASTORAGE_tlinkeddatarec_S *
  DATASTORAGE_plinkeddatarec;
typedef SYSTEM_uint8 _sub_2DATASTORAGE;
typedef SYSTEM_integer _arr_1DATASTORAGE[20];
typedef struct DATASTORAGE_tlinkeddatarec_S {
  DATASTORAGE_plinkeddatarec recnext;
  DATASTORAGE_plinkeddatarec hashnext;
  union{
    struct{
      SYSUTILS_P3_tbytearray recdata;
    } _c1;
    struct{
      _arr_1DATASTORAGE reckeys;
    } _c2;
  } _u;
} DATASTORAGE_tlinkeddatarec;


Constructor(DATASTORAGE_tgamsdatahashed ) 
  DATASTORAGE_tgamsdatahashed_DOT_create(
  DATASTORAGE_tgamsdatahashed self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec)
{
  ValueCast(DATASTORAGE_tgamsdatahashed,
    DATASTORAGE_tgamsdatastore_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatastore,self),adimension,adatasize,adefrec));
  self->DATASTORAGE_tgamsdatahashed_DOT_hl = ValueCast(
    DATASTORAGE_tgamshashlist,DATASTORAGE_tgamshashlist_DOT_create(ValueCast(
    DATASTORAGE_tgamshashlist,_P3alloc_object(&
    DATASTORAGE_tgamshashlist_CD)),adimension,adatasize));
  return self;
}  /* create */

Destructor(DATASTORAGE_tgamsdatahashed ) 
  DATASTORAGE_tgamsdatahashed_DOT_destroy(
  DATASTORAGE_tgamsdatahashed self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    DATASTORAGE_tgamsdatahashed_DOT_hl));
  DATASTORAGE_tgamsdatastore_DOT_destroy(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  return self;
}  /* destroy */

Procedure DATASTORAGE_tgamsdatahashed_DOT_clear(
  DATASTORAGE_tgamsdatahashed self)
{
  DATASTORAGE_tgamshashlist_DOT_clear(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl);
}  /* clear */

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatahashed_DOT_startread(
  DATASTORAGE_tgamsdatahashed self)
{
  SYSTEM_boolean result;

  result = DATASTORAGE_tgamshashlist_DOT_startread(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl);
  return result;
}  /* startread */

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatahashed_DOT_getnextkey(
  DATASTORAGE_tgamsdatahashed self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_P3_pbyte result;

  result = DATASTORAGE_tgamshashlist_DOT_getnextkey(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl,akey);
  return result;
}  /* getnextkey */

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatahashed_DOT_createsearcher(
  DATASTORAGE_tgamsdatahashed self)
{
  DATASTORAGE_tgamsdatasearcher result;

  result = ValueCast(DATASTORAGE_tgamsdatasearcher,
    DATASTORAGE_tgamsdatasearcher_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatasearcher,_P3alloc_object(&
    DATASTORAGE_tgamsdatahashedsearcher_CD)),ValueCast(
    DATASTORAGE_tgamsdatastore,self)));
  return result;
}  /* createsearcher */

Procedure DATASTORAGE_tgamsdatahashed_DOT_loadrecord(
  DATASTORAGE_tgamsdatahashed self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  DATASTORAGE_tgamshashlist_DOT_loadrecord(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl,akey,adata);
}  /* loadrecord */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatahashed_DOT_memoryused(
  DATASTORAGE_tgamsdatahashed self)
{
  SYSTEM_integer result;

  result = DATASTORAGE_tgamshashlist_DOT_memoryused(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl);
  return result;
}  /* memoryused */

Procedure DATASTORAGE_tgamsdatahashed_DOT_insertrecord(
  DATASTORAGE_tgamsdatahashed self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  DATASTORAGE_tgamshashlist_DOT_additem(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl,akey,adata);
}  /* insertrecord */

Procedure DATASTORAGE_tgamsdatahashed_DOT_endassign(
  DATASTORAGE_tgamsdatahashed self)
{
  DATASTORAGE_tgamsdatastore_DOT_endassign(ValueCast(
    DATASTORAGE_tgamsdatastore,self));
  if (self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec) {
    DATASTORAGE_tgamshashlist_DOT_removedefaults(self->
      DATASTORAGE_tgamsdatahashed_DOT_hl,&(*self->
      DATASTORAGE_tgamsdatastore_DOT_pdefrec)[0]);
    self->DATASTORAGE_tgamsdatastore_DOT_fsawdefrec = SYSTEM_false;
  } 
}  /* endassign */

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatahashed_DOT_getcount(
  DATASTORAGE_tgamsdatahashed self)
{
  SYSTEM_integer result;

  result = DATASTORAGE_tgamshashlist_DOT_getcount(self->
    DATASTORAGE_tgamsdatahashed_DOT_hl);
  return result;
}  /* getcount */

Constructor(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatasearcher_DOT_create(
  DATASTORAGE_tgamsdatasearcher self,
  DATASTORAGE_tgamsdatastore ads)
{
  ValueCast(DATASTORAGE_tgamsdatasearcher,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->DATASTORAGE_tgamsdatasearcher_DOT_ds = ads;
  return self;
}  /* create */

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatasearcher_DOT_startsearch(
  DATASTORAGE_tgamsdatasearcher self)
{
  SYSTEM_boolean result;

  result = VirtMethodCall(self->DATASTORAGE_tgamsdatasearcher_DOT_ds, 
    DATASTORAGE_tgamsdatastore_DOT_getcount_T, 1, (self->
    DATASTORAGE_tgamsdatasearcher_DOT_ds)) > 0;
  return result;
}  /* startsearch */

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatatablesearcher_DOT_startsearch(
  DATASTORAGE_tgamsdatatablesearcher self)
{
  SYSTEM_boolean result;

  result = DATASTORAGE_tgamsdatasearcher_DOT_startsearch(ValueCast(
    DATASTORAGE_tgamsdatasearcher,self));
  self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = 0;
  return result;
}  /* startsearch */

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatatablesearcher_DOT_search(
  DATASTORAGE_tgamsdatatablesearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata)
{
  SYSTEM_boolean result;
  cnstdef {close_search = 4};
  SYSTEM_integer l, h, i, c;
  GMSGEN_pbytedataarray prec;

  { register DATASTORAGE_tgamsdatatable_OD *_W2=ValueCast(
    DATASTORAGE_tgamsdatatable,self->
    DATASTORAGE_tgamsdatasearcher_DOT_ds);
    result = SYSTEM_false;
    l = 0;
    if (self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex < 0) 
      self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = 0;
    h = _W2->DATASTORAGE_tgamsdatatable_DOT_xlist->
      DATASTORAGE_treclist_DOT_fcount - 1;
    if (self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex > h) 
      self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = h;
    if (self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex < 0) {
      self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = 0;
      return result;
    } 
    prec = ValueCast(GMSGEN_pbytedataarray,
      DATASTORAGE_treclist_DOT_getitem(_W2->
      DATASTORAGE_tgamsdatatable_DOT_xlist,self->
      DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex));
    c = DATASTORAGE_tgamsdatastore_DOT_comparekeys(ValueCast(
      DATASTORAGE_tgamsdatastore,_W2),akey,ValueCast(
      GMSGEN_pintegerarrayone,prec));
    if (c == 0) 
      goto _Lfound_76;
    if (c > 0) { 
      for (i = 1;i <= (SYSTEM_int32)close_search;++i) {
        l = self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex + 
          i;
        if (l > h) {
          self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = h + 1;
          return result;
        } 
        prec = ValueCast(GMSGEN_pbytedataarray,
          DATASTORAGE_treclist_DOT_getitem(_W2->
          DATASTORAGE_tgamsdatatable_DOT_xlist,l));
        c = DATASTORAGE_tgamsdatastore_DOT_comparekeys(ValueCast(
          DATASTORAGE_tgamsdatastore,_W2),akey,ValueCast(
          GMSGEN_pintegerarrayone,prec));
        if (c == 0) {
          self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = l;
          goto _Lfound_76;
        } 
        if (c < 0) {
          self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = l;
          return result;
        } 
      
      }
    } else 
      for (i = 1;i <= (SYSTEM_int32)close_search;++i) {
        h = self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex - 
          i;
        if (h < 0) {
          self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = 0;
          return result;
        } 
        prec = ValueCast(GMSGEN_pbytedataarray,
          DATASTORAGE_treclist_DOT_getitem(_W2->
          DATASTORAGE_tgamsdatatable_DOT_xlist,h));
        c = DATASTORAGE_tgamsdatastore_DOT_comparekeys(ValueCast(
          DATASTORAGE_tgamsdatastore,_W2),akey,ValueCast(
          GMSGEN_pintegerarrayone,prec));
        if (c == 0) {
          self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = h;
          goto _Lfound_76;
        } 
        if (c > 0) {
          self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = h + 1;
          return result;
        } 
      
      }
    while (l <= h) {
      self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = ValueCast(
        SYSTEM_uint32,l + h) >> 1;
      prec = ValueCast(GMSGEN_pbytedataarray,
        DATASTORAGE_treclist_DOT_getitem(_W2->
        DATASTORAGE_tgamsdatatable_DOT_xlist,self->
        DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex));
      c = DATASTORAGE_tgamsdatastore_DOT_comparekeys(ValueCast(
        DATASTORAGE_tgamsdatastore,_W2),akey,ValueCast(
        GMSGEN_pintegerarrayone,prec));
      if (c > 0) { 
        l = self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex + 1;
      } else {
        if (c == 0) 
          goto _Lfound_76;
        h = self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex - 1;
      } 
    
}
    self->DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex = l;
    return result;
    _Lfound_76:;
    result = SYSTEM_true;
    *apdata = ValueCast(SYSTEM_pointer,&(*prec)[_W2->
      DATASTORAGE_tgamsdatastore_DOT_fkeysize]);

  }
  return result;
}  /* search */

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatafullsearcher_DOT_search(
  DATASTORAGE_tgamsdatafullsearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata)
{
  SYSTEM_boolean result;

  { register DATASTORAGE_tgamsdatafull_OD *_W2=ValueCast(
    DATASTORAGE_tgamsdatafull,self->
    DATASTORAGE_tgamsdatasearcher_DOT_ds);
    *apdata = DATASTORAGE_tgamsdatafull_DOT_getoffset(ValueCast(
      DATASTORAGE_tgamsdatafull,_W2),akey);
    result = *apdata != NULL && !
      DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(ValueCast(
      DATASTORAGE_tgamsdatastore,_W2),ValueCast(SYSTEM_P3_pbyte,*
      apdata));

  }
  return result;
}  /* search */

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_startsearch(
  DATASTORAGE_tgamsdatasparsesearcher self)
{
  SYSTEM_boolean result;
  SYSTEM_integer d;

  result = DATASTORAGE_tgamsdatasearcher_DOT_startsearch(ValueCast(
    DATASTORAGE_tgamsdatasearcher,self));
  (*self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[0] = (*(ValueCast(
    DATASTORAGE_tgamsdatasparse,self->
    DATASTORAGE_tgamsdatasearcher_DOT_ds))->
    DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs)[0];
  if (!result) { 
    (*self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[1] = NULL;
  } else 
    { register SYSTEM_int32 _stop = self->
        DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim;
      if ((d = 1) <=  _stop) do {
        (*self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d] = (ValueCast(
          DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d - 1]))->
          pcellright;
      } while (d++ !=  _stop);

    }
  return result;
}  /* startsearch */

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_search(
  DATASTORAGE_tgamsdatasparsesearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata)
{
  SYSTEM_boolean result;
  SYSTEM_integer d, d2, key;
  DATASTORAGE_psparsecell pc;
  DATASTORAGE_psparsecell pdown;

  result = SYSTEM_false;
  self->DATASTORAGE_tgamsdatasparsesearcher_DOT_flastvalid = 0;
  *apdata = NULL;
  if ((*self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[1] == NULL) 
    return result;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim;
    if ((d = 1) <=  _stop) do {
      key = (*akey)[d - 1];
      pc = ValueCast(DATASTORAGE_psparsecell,(*self->
        DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d]);
      if (pc->cellkey == key) {
        self->DATASTORAGE_tgamsdatasparsesearcher_DOT_flastvalid = d;
        SYSTEM_continue(CNT_4);
      } 
      if (key < pc->cellkey) 
        pc = (ValueCast(DATASTORAGE_psparsecell,(*self->
          DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d - 1]))->
          pcellright;
      pdown = pc->pcelldown;
      while (key > pc->cellkey && pdown != NULL && key >= pdown->
        cellkey) {
        pc = pdown;
        pdown = pc->pcelldown;
      
}
      if ((*self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[
        d] != ValueCast(SYSTEM_pointer,pc)) {
        (*self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d] = 
          pc;
        { register SYSTEM_int32 _stop = self->
            DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim;
          if ((d2 = d + 1) <=  _stop) do {
            (*self->
              DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d2] = (ValueCast(
              DATASTORAGE_psparsecell,(*self->
              DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[d2 - 1]))->
              pcellright;
          } while (d2++ !=  _stop);

        }
      } 
      if (pc->cellkey == key) {
        self->DATASTORAGE_tgamsdatasparsesearcher_DOT_flastvalid = d;
        SYSTEM_continue(CNT_4);
      } 
      return result;
    
CNT_4:;
    } while (d++ !=  _stop);
BRK_4:;

  }
  *apdata = ValueCast(SYSTEM_pointer,&(ValueCast(
    DATASTORAGE_psparsedatacell,(*self->
    DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs)[self->
    DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim]))->celldata);
  result = SYSTEM_true;
  return result;
}  /* search */

Constructor(DATASTORAGE_tgamsdatasparsesearcher ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_create(
  DATASTORAGE_tgamsdatasparsesearcher self,
  DATASTORAGE_tgamsdatasparse ads)
{
  ValueCast(DATASTORAGE_tgamsdatasparsesearcher,
    DATASTORAGE_tgamsdatasearcher_DOT_create(ValueCast(
    DATASTORAGE_tgamsdatasearcher,self),ValueCast(
    DATASTORAGE_tgamsdatastore,ads)));
  self->DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs = 
    DATASTORAGE_tgamsdatastore_DOT_allocptrs(ValueCast(
    DATASTORAGE_tgamsdatastore,ads));
  self->DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim = ads->
    DATASTORAGE_tgamsdatastore_DOT_fdimension;
  return self;
}  /* create */

Destructor(DATASTORAGE_tgamsdatasparsesearcher ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_destroy(
  DATASTORAGE_tgamsdatasparsesearcher self)
{
  DATASTORAGE_tgamsdatastore_DOT_freeptrs(self->
    DATASTORAGE_tgamsdatasearcher_DOT_ds,self->
    DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatahashedsearcher_DOT_search(
  DATASTORAGE_tgamsdatahashedsearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata)
{
  SYSTEM_boolean result;

  { register DATASTORAGE_tgamshashlist_OD *_W2=(ValueCast(
    DATASTORAGE_tgamsdatahashed,self->
    DATASTORAGE_tgamsdatasearcher_DOT_ds))->
    DATASTORAGE_tgamsdatahashed_DOT_hl;
    *apdata = DATASTORAGE_tgamshashlist_DOT_indexof(ValueCast(
      DATASTORAGE_tgamshashlist,_W2),akey);
    result = *apdata != NULL;

  }
  return result;
}  /* search */

Constructor(DATASTORAGE_tgamshashlist ) 
  DATASTORAGE_tgamshashlist_DOT_create(
  DATASTORAGE_tgamshashlist self,
  SYSTEM_integer adim,
  SYSTEM_integer adatasize)
{
  ValueCast(DATASTORAGE_tgamshashlist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->DATASTORAGE_tgamshashlist_DOT_fdimension = adim;
  self->DATASTORAGE_tgamshashlist_DOT_fdatasize = adatasize;
  self->DATASTORAGE_tgamshashlist_DOT_fkeysize = adim * sizeof(
    SYSTEM_longint);
  self->DATASTORAGE_tgamshashlist_DOT_phashtable = NULL;
  self->DATASTORAGE_tgamshashlist_DOT_hashsize = 0;
  self->DATASTORAGE_tgamshashlist_DOT_lnkdata = ValueCast(
    DATASTORAGE_tlinkeddata,DATASTORAGE_tlinkeddata_DOT_create(ValueCast(
    DATASTORAGE_tlinkeddata,_P3alloc_object(&
    DATASTORAGE_tlinkeddata_CD)),adim,adatasize));
  return self;
}  /* create */

Function(SYSTEM_boolean ) DATASTORAGE_tgamshashlist_DOT_additem(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_boolean result;
  SYSTEM_integer hv;
  DATASTORAGE_plinkeddatarec prec;

  if (self->DATASTORAGE_tgamshashlist_DOT_phashtable == NULL || self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata->
    DATASTORAGE_tlinkeddata_DOT_fcount > self->
    DATASTORAGE_tgamshashlist_DOT_rehashcnt) 
    DATASTORAGE_tgamshashlist_DOT_hashall(self);
  hv = DATASTORAGE_tgamshashlist_DOT_hash(self,akey);
  prec = ValueCast(DATASTORAGE_plinkeddatarec,(*self->
    DATASTORAGE_tgamshashlist_DOT_phashtable)[hv]);
  while (prec != NULL) {

    if (!DATASTORAGE_tgamshashlist_DOT_equalkeys(self,akey,ValueCast(
      GMSGEN_pintegerarrayone,&prec->_u._c2.reckeys[0]))) { 
      prec = prec->hashnext;
    } else {
      GMSOBJ_cmove(adata,&prec->_u._c1.recdata[self->
        DATASTORAGE_tgamshashlist_DOT_fkeysize],self->
        DATASTORAGE_tgamshashlist_DOT_fdatasize);
      result = SYSTEM_false;
      return result;
    } 
}
  prec = ValueCast(DATASTORAGE_plinkeddatarec,
    DATASTORAGE_tlinkeddata_DOT_additem(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata,akey,adata));
  prec->hashnext = ValueCast(DATASTORAGE_plinkeddatarec,(*self->
    DATASTORAGE_tgamshashlist_DOT_phashtable)[hv]);
  (*self->DATASTORAGE_tgamshashlist_DOT_phashtable)[hv] = prec;
  result = SYSTEM_true;
  return result;
}  /* additem */

Procedure DATASTORAGE_tgamshashlist_DOT_loadrecord(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  if (self->DATASTORAGE_tgamshashlist_DOT_phashtable != NULL) 
    DATASTORAGE_tgamshashlist_DOT_clearhashlist(self);
  DATASTORAGE_tlinkeddata_DOT_additem(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata,akey,adata);
}  /* loadrecord */

Destructor(DATASTORAGE_tgamshashlist ) 
  DATASTORAGE_tgamshashlist_DOT_destroy(
  DATASTORAGE_tgamshashlist self)
{
  DATASTORAGE_tgamshashlist_DOT_clear(self);
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure DATASTORAGE_tgamshashlist_DOT_clear(
  DATASTORAGE_tgamshashlist self)
{
  DATASTORAGE_tgamshashlist_DOT_clearhashlist(self);
  DATASTORAGE_tlinkeddata_DOT_clear(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata);
}  /* clear */

Procedure DATASTORAGE_tgamshashlist_DOT_clearhashlist(
  DATASTORAGE_tgamshashlist self)
{
  if (self->DATASTORAGE_tgamshashlist_DOT_phashtable != NULL) {
    GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,self->
      DATASTORAGE_tgamshashlist_DOT_phashtable,self->
      DATASTORAGE_tgamshashlist_DOT_hashsize * sizeof(SYSTEM_pointer));
    self->DATASTORAGE_tgamshashlist_DOT_phashtable = NULL;
    self->DATASTORAGE_tgamshashlist_DOT_hashsize = 0;
  } 
}  /* clearhashlist */

Function(SYSTEM_boolean ) DATASTORAGE_tgamshashlist_DOT_equalkeys(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey1,
  GMSGEN_pintegerarrayone akey2)
{
  SYSTEM_boolean result;
  SYSTEM_integer d;

  result = (*akey1)[0] == (*akey2)[0];
  if (!result) 
    return result;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamshashlist_DOT_fdimension;
    if ((d = 2) <=  _stop) do {
      result = (*akey1)[d - 1] == (*akey2)[d - 1];
      if (!result) 
        return result;
    
    } while (d++ !=  _stop);

  }
  return result;
}  /* equalkeys */

Function(SYSTEM_integer ) DATASTORAGE_tgamshashlist_DOT_hash(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  SYSTEM_cardinal v;

  v = (*akey)[0];
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamshashlist_DOT_fdimension;
    if ((d = 2) <=  _stop) do {
      v = 1234593 * v + (*akey)[d - 1] & 2147483647;
    } while (d++ !=  _stop);

  }
  result = v % self->DATASTORAGE_tgamshashlist_DOT_hashsize;
  return result;
}  /* hash */

Procedure DATASTORAGE_tgamshashlist_DOT_hashall(
  DATASTORAGE_tgamshashlist self)
{
  SYSTEM_integer hv;
  DATASTORAGE_plinkeddatarec prec;

  DATASTORAGE_tgamshashlist_DOT_clearhashlist(self);
  DATASTORAGE_tgamshashlist_DOT_hashtablereset(self,self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata->
    DATASTORAGE_tlinkeddata_DOT_fcount);
  { register DATASTORAGE_tlinkeddata_OD *_W2=self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata;
    prec = ValueCast(DATASTORAGE_plinkeddatarec,_W2->
      DATASTORAGE_tlinkeddata_DOT_fhead);
    while (prec != NULL) {
      hv = DATASTORAGE_tgamshashlist_DOT_hash(self,ValueCast(
        GMSGEN_pintegerarrayone,&prec->_u._c2.reckeys[0]));
      prec->hashnext = ValueCast(DATASTORAGE_plinkeddatarec,(*self->
        DATASTORAGE_tgamshashlist_DOT_phashtable)[hv]);
      (*self->DATASTORAGE_tgamshashlist_DOT_phashtable)[hv] = prec;
      prec = prec->recnext;
    
}

  }
}  /* hashall */

static Function(SYSTEM_integer ) DATASTORAGE_calcnexthashsize(
  SYSTEM_integer cnt,
  SYSTEM_integer *nxt)
{
  SYSTEM_integer result;
  cnstdef {hashsize_1 = 997};
  cnstdef {next_1 = 1500};
  cnstdef {hashsize_2 = 9973};
  cnstdef {next_2 = 15000};
  cnstdef {hashsize_3 = 99991};
  cnstdef {next_3 = 150000};
  cnstdef {hashsize_4 = 999979};
  cnstdef {next_4 = 1500000};
  cnstdef {hashsize_5 = 9999991};
  cnstdef {next_5 = SYSTEM_maxint};

  if (cnt >= next_4) {
    result = hashsize_5;
    *nxt = next_5;
  } else 
    if (cnt >= next_3) {
      result = hashsize_4;
      *nxt = next_4;
    } else 
      if (cnt >= 15000) {
        result = hashsize_3;
        *nxt = next_3;
      } else 
        if (cnt >= 1500) {
          result = hashsize_2;
          *nxt = next_2;
        } else {
          result = hashsize_1;
          *nxt = next_1;
        } 
  return result;
}  /* calcnexthashsize */

Procedure DATASTORAGE_tgamshashlist_DOT_hashtablereset(
  DATASTORAGE_tgamshashlist self,
  SYSTEM_integer acnt)
{
  SYSTEM_integer n;

  self->DATASTORAGE_tgamshashlist_DOT_hashsize = 
    DATASTORAGE_calcnexthashsize(acnt,&self->
    DATASTORAGE_tgamshashlist_DOT_rehashcnt);
  self->DATASTORAGE_tgamshashlist_DOT_phashtable = ValueCast(
    SYSTEM_P3_ppointerarray,GMSHEAPNEW_theapmgr_DOT_xgetmem(
    GMSHEAPNEW_gheap,self->DATASTORAGE_tgamshashlist_DOT_hashsize * sizeof(
    SYSTEM_pointer)));
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tgamshashlist_DOT_hashsize - 1;
    if ((n = 0) <=  _stop) do {
      (*self->DATASTORAGE_tgamshashlist_DOT_phashtable)[n] = NULL;
    } while (n++ !=  _stop);

  }
}  /* hashtablereset */

Function(SYSTEM_pointer ) DATASTORAGE_tgamshashlist_DOT_indexof(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_pointer result;
  DATASTORAGE_plinkeddatarec prec;

  if (self->DATASTORAGE_tgamshashlist_DOT_phashtable == NULL) 
    DATASTORAGE_tgamshashlist_DOT_hashall(self);
  prec = ValueCast(DATASTORAGE_plinkeddatarec,(*self->
    DATASTORAGE_tgamshashlist_DOT_phashtable)[
    DATASTORAGE_tgamshashlist_DOT_hash(self,akey)]);
  while (prec != NULL) {

    if (!DATASTORAGE_tgamshashlist_DOT_equalkeys(self,akey,ValueCast(
      GMSGEN_pintegerarrayone,&prec->_u._c2.reckeys[0]))) { 
      prec = prec->hashnext;
    } else {
      result = ValueCast(SYSTEM_pointer,&prec->_u._c1.recdata[self->
        DATASTORAGE_tgamshashlist_DOT_fkeysize]);
      return result;
    } 
}
  result = NULL;
  return result;
}  /* indexof */

Function(SYSTEM_integer ) DATASTORAGE_tgamshashlist_DOT_memoryused(
  DATASTORAGE_tgamshashlist self)
{
  SYSTEM_integer result;

  result = DATASTORAGE_tlinkeddata_DOT_memoryused(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata);
  if (self->DATASTORAGE_tgamshashlist_DOT_hashsize > 0) 
    result = result + self->DATASTORAGE_tgamshashlist_DOT_hashsize * sizeof(
      SYSTEM_pointer);
  return result;
}  /* memoryused */

Function(SYSTEM_integer ) DATASTORAGE_tgamshashlist_DOT_getcount(
  DATASTORAGE_tgamshashlist self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tgamshashlist_DOT_lnkdata->
    DATASTORAGE_tlinkeddata_DOT_fcount;
  return result;
}  /* getcount */

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamshashlist_DOT_getnextkey(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_P3_pbyte result;

  result = DATASTORAGE_tlinkeddata_DOT_getnextkey(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata,&self->
    DATASTORAGE_tgamshashlist_DOT_fcurrreadp,akey);
  return result;
}  /* getnextkey */

Function(SYSTEM_boolean ) DATASTORAGE_tgamshashlist_DOT_startread(
  DATASTORAGE_tgamshashlist self)
{
  SYSTEM_boolean result;

  result = DATASTORAGE_tlinkeddata_DOT_startread(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata,&self->
    DATASTORAGE_tgamshashlist_DOT_fcurrreadp,ValueCast(
    GMSGEN_pintegerarrayone,NULL));
  return result;
}  /* startread */

Procedure DATASTORAGE_tgamshashlist_DOT_removedefaults(
  DATASTORAGE_tgamshashlist self,
  const SYSTEM_untyped *defdata)
{
  if (DATASTORAGE_tlinkeddata_DOT_removedefaults(self->
    DATASTORAGE_tgamshashlist_DOT_lnkdata,defdata)) 
    DATASTORAGE_tgamshashlist_DOT_clearhashlist(self);
}  /* removedefaults */

Constructor(DATASTORAGE_treclist ) DATASTORAGE_treclist_DOT_create(
  DATASTORAGE_treclist self,
  SYSTEM_integer arecsize)
{
  ValueCast(DATASTORAGE_treclist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->DATASTORAGE_treclist_DOT_frecsize = arecsize;
  self->DATASTORAGE_treclist_DOT_flist = NULL;
  self->DATASTORAGE_treclist_DOT_fcapacity = 0;
  self->DATASTORAGE_treclist_DOT_fcount = 0;
  return self;
}  /* create */

Destructor(DATASTORAGE_treclist ) DATASTORAGE_treclist_DOT_destroy(
  DATASTORAGE_treclist self)
{
  DATASTORAGE_treclist_DOT_clear(self);
  SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
    DATASTORAGE_treclist_DOT_flist),0);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure DATASTORAGE_treclist_DOT_setcapacity(
  DATASTORAGE_treclist self,
  SYSTEM_integer newcapacity)
{
  if (newcapacity != self->DATASTORAGE_treclist_DOT_fcapacity) {
    if (newcapacity < self->DATASTORAGE_treclist_DOT_fcount) 
      newcapacity = self->DATASTORAGE_treclist_DOT_fcount;
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      DATASTORAGE_treclist_DOT_flist),newcapacity * sizeof(
      SYSTEM_pointer));
    self->DATASTORAGE_treclist_DOT_fcapacity = newcapacity;
  } 
}  /* setcapacity */

Procedure DATASTORAGE_treclist_DOT_grow(
  DATASTORAGE_treclist self)
{
  SYSTEM_integer delta;

  if (self->DATASTORAGE_treclist_DOT_fcapacity >= 1048576) { 
    delta = self->DATASTORAGE_treclist_DOT_fcapacity /  4;
  } else 
    if (self->DATASTORAGE_treclist_DOT_fcapacity == 0) { 
      delta = 16;
    } else 
      delta = 3 * self->DATASTORAGE_treclist_DOT_fcapacity;
  DATASTORAGE_treclist_DOT_setcapacity(self,self->
    DATASTORAGE_treclist_DOT_fcapacity + delta);
}  /* grow */

Function(SYSTEM_pointer ) DATASTORAGE_treclist_DOT_additem(
  DATASTORAGE_treclist self)
{
  SYSTEM_pointer result;

  if (self->DATASTORAGE_treclist_DOT_fcount == self->
    DATASTORAGE_treclist_DOT_fcapacity) 
    VirtMethodCall(self, DATASTORAGE_treclist_DOT_grow_T, 1, (self));
  result = GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,self->
    DATASTORAGE_treclist_DOT_frecsize);
  (*self->DATASTORAGE_treclist_DOT_flist)[self->
    DATASTORAGE_treclist_DOT_fcount] = result;
  _P3inc0(self->DATASTORAGE_treclist_DOT_fcount);
  return result;
}  /* additem */

Procedure DATASTORAGE_treclist_DOT_clear(
  DATASTORAGE_treclist self)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_treclist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,(*self->
        DATASTORAGE_treclist_DOT_flist)[n],self->
        DATASTORAGE_treclist_DOT_frecsize);
    } while (n++ !=  _stop);

  }
  self->DATASTORAGE_treclist_DOT_fcount = 0;
  DATASTORAGE_treclist_DOT_setcapacity(self,0);
}  /* clear */

Procedure DATASTORAGE_treclist_DOT_exchange(
  DATASTORAGE_treclist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_pointer t;

  t = (*self->DATASTORAGE_treclist_DOT_flist)[index1];
  (*self->DATASTORAGE_treclist_DOT_flist)[index1] = (*self->
    DATASTORAGE_treclist_DOT_flist)[index2];
  (*self->DATASTORAGE_treclist_DOT_flist)[index2] = t;
}  /* exchange */

Procedure DATASTORAGE_treclist_DOT_freeitem(
  DATASTORAGE_treclist self,
  SYSTEM_integer index)
{
  GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,(*self->
    DATASTORAGE_treclist_DOT_flist)[index],self->
    DATASTORAGE_treclist_DOT_frecsize);
}  /* freeitem */

Function(SYSTEM_pointer ) DATASTORAGE_treclist_DOT_getitem(
  DATASTORAGE_treclist self,
  SYSTEM_integer index)
{
  SYSTEM_pointer result;

  result = (*self->DATASTORAGE_treclist_DOT_flist)[index];
  return result;
}  /* getitem */

Function(SYSTEM_pointer ) DATASTORAGE_treclist_DOT_insert(
  DATASTORAGE_treclist self,
  SYSTEM_integer index)
{
  SYSTEM_pointer result;

  if (self->DATASTORAGE_treclist_DOT_fcount == self->
    DATASTORAGE_treclist_DOT_fcapacity) 
    VirtMethodCall(self, DATASTORAGE_treclist_DOT_grow_T, 1, (self));
  if (index < self->DATASTORAGE_treclist_DOT_fcount) 
    SYSTEM_move(&(*self->DATASTORAGE_treclist_DOT_flist)[index],&(*
      self->DATASTORAGE_treclist_DOT_flist)[index + 1],(self->
      DATASTORAGE_treclist_DOT_fcount - index) * sizeof(SYSTEM_pointer));
  result = GMSHEAPNEW_theapmgr_DOT_xgetmem(GMSHEAPNEW_gheap,self->
    DATASTORAGE_treclist_DOT_frecsize);
  (*self->DATASTORAGE_treclist_DOT_flist)[index] = result;
  _P3inc0(self->DATASTORAGE_treclist_DOT_fcount);
  return result;
}  /* insert */

Procedure DATASTORAGE_treclist_DOT_remove(
  DATASTORAGE_treclist self,
  SYSTEM_integer index)
{
  GMSHEAPNEW_theapmgr_DOT_xfreemem(GMSHEAPNEW_gheap,(*self->
    DATASTORAGE_treclist_DOT_flist)[index],self->
    DATASTORAGE_treclist_DOT_frecsize);
  (*self->DATASTORAGE_treclist_DOT_flist)[index] = NULL;
}  /* remove */

Procedure DATASTORAGE_treclist_DOT_cleanup(
  DATASTORAGE_treclist self)
{
  SYSTEM_integer n;
  SYSTEM_integer w;

  w = 0;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_treclist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if ((*self->DATASTORAGE_treclist_DOT_flist)[n] != NULL) {
        if (w != n) 
          (*self->DATASTORAGE_treclist_DOT_flist)[w] = (*self->
            DATASTORAGE_treclist_DOT_flist)[n];
        w = w + 1;
      } 
    } while (n++ !=  _stop);

  }
  self->DATASTORAGE_treclist_DOT_fcount = w;
}  /* cleanup */

Constructor(DATASTORAGE_tlinkeddata ) 
  DATASTORAGE_tlinkeddata_DOT_create(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize)
{
  ValueCast(DATASTORAGE_tlinkeddata,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->DATASTORAGE_tlinkeddata_DOT_myheap = ValueCast(
    GMSHEAPNEW_theapmgr,GMSHEAPNEW_theapmgr_DOT_create(ValueCast(
    GMSHEAPNEW_theapmgr,_P3alloc_object(&GMSHEAPNEW_theapmgr_CD)),
    GMSHEAPNEW_bbmgr,_P3str1("\013TLinkedData")));
  self->DATASTORAGE_tlinkeddata_DOT_fdimension = adimension;
  self->DATASTORAGE_tlinkeddata_DOT_fkeysize = adimension * sizeof(
    SYSTEM_longint);
  self->DATASTORAGE_tlinkeddata_DOT_fdatasize = adatasize;
  self->DATASTORAGE_tlinkeddata_DOT_ftotalsize = sizeof(SYSTEM_pointer) + sizeof(
    SYSTEM_pointer) + self->DATASTORAGE_tlinkeddata_DOT_fkeysize + 
    self->DATASTORAGE_tlinkeddata_DOT_fdatasize;
  self->DATASTORAGE_tlinkeddata_DOT_fhead = NULL;
  self->DATASTORAGE_tlinkeddata_DOT_ftail = NULL;
  self->DATASTORAGE_tlinkeddata_DOT_fcount = 0;
  self->DATASTORAGE_tlinkeddata_DOT_fmaxkey = 0;
  self->DATASTORAGE_tlinkeddata_DOT_fminkey = SYSTEM_maxint;
  return self;
}  /* create */

Destructor(DATASTORAGE_tlinkeddata ) 
  DATASTORAGE_tlinkeddata_DOT_destroy(
  DATASTORAGE_tlinkeddata self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    DATASTORAGE_tlinkeddata_DOT_myheap));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure DATASTORAGE_tlinkeddata_DOT_clear(
  DATASTORAGE_tlinkeddata self)
{
  DATASTORAGE_plinkeddatarec p, pn;

  p = ValueCast(DATASTORAGE_plinkeddatarec,self->
    DATASTORAGE_tlinkeddata_DOT_fhead);
  while (p != NULL) {
    pn = p->recnext;
    GMSHEAPNEW_theapmgr_DOT_xfreemem(self->
      DATASTORAGE_tlinkeddata_DOT_myheap,p,self->
      DATASTORAGE_tlinkeddata_DOT_ftotalsize);
    p = pn;
  
}
  self->DATASTORAGE_tlinkeddata_DOT_fcount = 0;
  self->DATASTORAGE_tlinkeddata_DOT_fhead = NULL;
  self->DATASTORAGE_tlinkeddata_DOT_ftail = NULL;
  self->DATASTORAGE_tlinkeddata_DOT_fmaxkey = 0;
  self->DATASTORAGE_tlinkeddata_DOT_fminkey = SYSTEM_maxint;
}  /* clear */

Function(SYSTEM_integer ) DATASTORAGE_tlinkeddata_DOT_memoryused(
  DATASTORAGE_tlinkeddata self)
{
  SYSTEM_integer result;

  result = self->DATASTORAGE_tlinkeddata_DOT_fcount * self->
    DATASTORAGE_tlinkeddata_DOT_ftotalsize;
  return result;
}  /* memoryused */

Function(SYSTEM_pointer ) DATASTORAGE_tlinkeddata_DOT_additem(
  DATASTORAGE_tlinkeddata self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata)
{
  SYSTEM_pointer result;
  SYSTEM_integer d;
  SYSTEM_integer key;

  result = GMSHEAPNEW_theapmgr_DOT_xgetmem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,self->
    DATASTORAGE_tlinkeddata_DOT_ftotalsize);
  if (self->DATASTORAGE_tlinkeddata_DOT_fhead == NULL) { 
    self->DATASTORAGE_tlinkeddata_DOT_fhead = result;
  } else 
    (ValueCast(DATASTORAGE_plinkeddatarec,self->
      DATASTORAGE_tlinkeddata_DOT_ftail))->recnext = ValueCast(
      DATASTORAGE_plinkeddatarec,result);
  self->DATASTORAGE_tlinkeddata_DOT_ftail = result;
  { register DATASTORAGE_tlinkeddatarec *_W2=ValueCast(
    DATASTORAGE_plinkeddatarec,result);
    _W2->recnext = NULL;
    GMSOBJ_cmove(&(*akey)[0],&_W2->_u._c1.recdata[0],self->
      DATASTORAGE_tlinkeddata_DOT_fkeysize);
    GMSOBJ_cmove(adata,&_W2->_u._c1.recdata[self->
      DATASTORAGE_tlinkeddata_DOT_fkeysize],self->
      DATASTORAGE_tlinkeddata_DOT_fdatasize);

  }
  self->DATASTORAGE_tlinkeddata_DOT_fcount = self->
    DATASTORAGE_tlinkeddata_DOT_fcount + 1;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tlinkeddata_DOT_fdimension;
    if ((d = 1) <=  _stop) do {
      key = (*akey)[d - 1];
      if (key > self->DATASTORAGE_tlinkeddata_DOT_fmaxkey) 
        self->DATASTORAGE_tlinkeddata_DOT_fmaxkey = key;
      if (key < self->DATASTORAGE_tlinkeddata_DOT_fminkey) 
        self->DATASTORAGE_tlinkeddata_DOT_fminkey = key;
    
    } while (d++ !=  _stop);

  }
  return result;
}  /* additem */

static Function(SYSTEM_boolean ) issorted(
  DATASTORAGE_tlinkeddata *_2self)
{
  SYSTEM_boolean result;
  GMSGEN_pintegerarrayone prevkey;
  SYSTEM_integer d;
  SYSTEM_integer kd;
  DATASTORAGE_plinkeddatarec r;

  r = ValueCast(DATASTORAGE_plinkeddatarec,(*_2self)->
    DATASTORAGE_tlinkeddata_DOT_fhead);
  prevkey = ValueCast(GMSGEN_pintegerarrayone,r->_u._c2.reckeys);
  r = r->recnext;
  result = SYSTEM_true;
  while (r != NULL) {
    { register SYSTEM_int32 _stop = (*_2self)->
        DATASTORAGE_tlinkeddata_DOT_fdimension;
      if ((d = 1) <=  _stop) do {
        kd = r->_u._c2.reckeys[d - 1] - (*prevkey)[d - 1];
        if (kd != 0) 
          SYSTEM_break(BRK_5);
      
CNT_5:;
      } while (d++ !=  _stop);
BRK_5:;

    }
    if (kd < 0) {
      result = SYSTEM_false;
      SYSTEM_break(BRK_6);
    } 
    prevkey = ValueCast(GMSGEN_pintegerarrayone,r->_u._c2.reckeys);
    r = r->recnext;
  
CNT_6:;
  }
BRK_6:;
  return result;
}  /* issorted */

Procedure DATASTORAGE_tlinkeddata_DOT_sort(
  DATASTORAGE_tlinkeddata self,
  GMSGEN_pintegerarrayone amap)
{
  SYSTEM_P3_ppointerarray head, tail;
  SYSTEM_integer key;
  SYSTEM_integer allocsize;
  SYSTEM_integer keybase;
  SYSTEM_integer d;
  DATASTORAGE_plinkeddatarec r;

  if (self->DATASTORAGE_tlinkeddata_DOT_fhead == NULL) 
    return;
  if (issorted(&self)) 
    return;
  allocsize = (self->DATASTORAGE_tlinkeddata_DOT_fmaxkey - self->
    DATASTORAGE_tlinkeddata_DOT_fminkey + 1) * sizeof(
    SYSTEM_pointer);
  head = ValueCast(SYSTEM_P3_ppointerarray,
    GMSHEAPNEW_theapmgr_DOT_xgetmem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,allocsize));
  tail = ValueCast(SYSTEM_P3_ppointerarray,
    GMSHEAPNEW_theapmgr_DOT_xgetmem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,allocsize));
  keybase = self->DATASTORAGE_tlinkeddata_DOT_fminkey;
  { register SYSTEM_int32 _stop = self->
      DATASTORAGE_tlinkeddata_DOT_fmaxkey - keybase;
    if ((key = 0) <=  _stop) do {
      (*head)[key] = NULL;
    } while (key++ !=  _stop);

  }
  for (d = self->DATASTORAGE_tlinkeddata_DOT_fdimension;d >= (
    SYSTEM_int32)1;--d) {
    r = ValueCast(DATASTORAGE_plinkeddatarec,self->
      DATASTORAGE_tlinkeddata_DOT_fhead);
    while (r != NULL) {
      if (amap == NULL) { 
        key = r->_u._c2.reckeys[d - 1] - keybase;
      } else 
        key = r->_u._c2.reckeys[(*amap)[d - 1] - 1] - 
          keybase;
      if ((*head)[key] == NULL) { 
        (*head)[key] = r;
      } else 
        (ValueCast(DATASTORAGE_plinkeddatarec,(*tail)[key]))->recnext = 
          r;
      (*tail)[key] = r;
      r = r->recnext;
    
}
    r = NULL;
    for (key = self->DATASTORAGE_tlinkeddata_DOT_fmaxkey - keybase;key >= (
      SYSTEM_int32)0;--key) {
      if ((*head)[key] != NULL) {
        (ValueCast(DATASTORAGE_plinkeddatarec,(*tail)[key]))->recnext = 
          r;
        r = ValueCast(DATASTORAGE_plinkeddatarec,(*head)[key]);
        (*head)[key] = NULL;
      } 
    }
    self->DATASTORAGE_tlinkeddata_DOT_fhead = r;
  
  }
  self->DATASTORAGE_tlinkeddata_DOT_ftail = NULL;
  GMSHEAPNEW_theapmgr_DOT_xfreemem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,head,allocsize);
  GMSHEAPNEW_theapmgr_DOT_xfreemem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,tail,allocsize);
}  /* sort */

Function(SYSTEM_boolean ) DATASTORAGE_tlinkeddata_DOT_startread(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_pointer *p,
  GMSGEN_pintegerarrayone amap)
{
  SYSTEM_boolean result;

  result = self->DATASTORAGE_tlinkeddata_DOT_fcount > 0;
  if (!result) { 
    *p = NULL;
  } else {
    DATASTORAGE_tlinkeddata_DOT_sort(self,amap);
    *p = self->DATASTORAGE_tlinkeddata_DOT_fhead;
  } 
  return result;
}  /* startread */

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tlinkeddata_DOT_getnextkey(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_pointer *p,
  GMSGEN_pintegerarrayone akey)
{
  SYSTEM_P3_pbyte result;

  if (*p == NULL) { 
    result = NULL;
  } else 
    { register DATASTORAGE_tlinkeddatarec *_W2=ValueCast(
      DATASTORAGE_plinkeddatarec,*p);
      GMSOBJ_cmove(&_W2->_u._c1.recdata[0],&(*akey)[0],self->
        DATASTORAGE_tlinkeddata_DOT_fkeysize);
      result = ValueCast(SYSTEM_P3_pbyte,&_W2->_u._c1.recdata[self->
        DATASTORAGE_tlinkeddata_DOT_fkeysize]);
      *p = _W2->recnext;

    }
  return result;
}  /* getnextkey */

Function(SYSTEM_boolean ) DATASTORAGE_tlinkeddata_DOT_getnextrecord(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_pointer *p,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_untyped *data)
{
  SYSTEM_boolean result;

  result = *p != NULL;
  if (result) 
    { register DATASTORAGE_tlinkeddatarec *_W2=ValueCast(
      DATASTORAGE_plinkeddatarec,*p);
      GMSOBJ_cmove(&_W2->_u._c1.recdata[0],&(*akey)[0],self->
        DATASTORAGE_tlinkeddata_DOT_fkeysize);
      GMSOBJ_cmove(&_W2->_u._c1.recdata[self->
        DATASTORAGE_tlinkeddata_DOT_fkeysize],data,self->
        DATASTORAGE_tlinkeddata_DOT_fdatasize);
      *p = _W2->recnext;

    }
  return result;
}  /* getnextrecord */

Function(GMSGEN_pintegerarrayone ) 
  DATASTORAGE_tlinkeddata_DOT_allocindex(
  DATASTORAGE_tlinkeddata self)
{
  GMSGEN_pintegerarrayone result;

  result = ValueCast(GMSGEN_pintegerarrayone,
    GMSHEAPNEW_theapmgr_DOT_xgetmem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,self->
    DATASTORAGE_tlinkeddata_DOT_fkeysize));
  return result;
}  /* allocindex */

Procedure DATASTORAGE_tlinkeddata_DOT_freeindex(
  DATASTORAGE_tlinkeddata self,
  GMSGEN_pintegerarrayone p)
{
  GMSHEAPNEW_theapmgr_DOT_xfreemem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,p,self->
    DATASTORAGE_tlinkeddata_DOT_fkeysize);
}  /* freeindex */

Procedure DATASTORAGE_tlinkeddata_DOT_createlist(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_P3_ppointerarray *alist)
{
  DATASTORAGE_plinkeddatarec prec;
  SYSTEM_integer n;

  *alist = ValueCast(SYSTEM_P3_ppointerarray,
    GMSHEAPNEW_theapmgr_DOT_xgetmem(self->
    DATASTORAGE_tlinkeddata_DOT_myheap,(self->
    DATASTORAGE_tlinkeddata_DOT_fcount + 1) * sizeof(
    SYSTEM_pointer)));
  prec = ValueCast(DATASTORAGE_plinkeddatarec,self->
    DATASTORAGE_tlinkeddata_DOT_fhead);
  n = 0;
  while (prec != NULL) {
    n = n + 1;
    (**alist)[n] = ValueCast(SYSTEM_pointer,&prec->_u._c2.reckeys[0]);
    prec = prec->recnext;
  
}
}  /* createlist */

Function(SYSTEM_boolean ) DATASTORAGE_tlinkeddata_DOT_removedefaults(
  DATASTORAGE_tlinkeddata self,
  const SYSTEM_untyped *defdata)
{
  SYSTEM_boolean result;
  DATASTORAGE_plinkeddatarec prec;
  DATASTORAGE_plinkeddatarec pnext;
  DATASTORAGE_plinkeddatarec prevrec;

  result = SYSTEM_false;
  if (self->DATASTORAGE_tlinkeddata_DOT_fdatasize == 0) 
    return result;
  prec = ValueCast(DATASTORAGE_plinkeddatarec,self->
    DATASTORAGE_tlinkeddata_DOT_fhead);
  self->DATASTORAGE_tlinkeddata_DOT_fhead = NULL;
  prevrec = NULL;
  while (prec != NULL) {

    if (!DATASTORAGE_dataequal(&prec->_u._c1.recdata[self->
      DATASTORAGE_tlinkeddata_DOT_fkeysize],defdata,self->
      DATASTORAGE_tlinkeddata_DOT_fdatasize)) {
      if (prevrec == NULL) { 
        self->DATASTORAGE_tlinkeddata_DOT_fhead = prec;
      } else 
        prevrec->recnext = prec;
      prevrec = prec;
      prec = prec->recnext;
    } else {
      result = SYSTEM_true;
      pnext = prec->recnext;
      GMSHEAPNEW_theapmgr_DOT_xfreemem(self->
        DATASTORAGE_tlinkeddata_DOT_myheap,prec,self->
        DATASTORAGE_tlinkeddata_DOT_ftotalsize);
      self->DATASTORAGE_tlinkeddata_DOT_fcount = self->
        DATASTORAGE_tlinkeddata_DOT_fcount - 1;
      prec = pnext;
    } 
}
  if (prevrec != NULL) 
    prevrec->recnext = NULL;
  self->DATASTORAGE_tlinkeddata_DOT_ftail = prevrec;
  return result;
}  /* removedefaults */

/* unit datastorage */
void _Init_Module_datastorage(void)
{
} /* _Init_Module_datastorage */

void _Final_Module_datastorage(void)
{
} /* _Final_Module_datastorage */

