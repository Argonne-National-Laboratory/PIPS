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
#include "gmsdata.h"


void * const GMSDATA_tgrowarray_VT[] = {(void*)&
  GMSDATA_tgrowarray_DOT_destroy};

/* Class descriptor for 'tgrowarray' */
const SYSTEM_classdescriptor_t GMSDATA_tgrowarray_CD = {
  _P3str1("\012tgrowarray"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSDATA_tgrowarray_OD), GMSDATA_tgrowarray_VT, NULL};


void * const GMSDATA_tgrowarrayfxd_VT[] = {(void*)&
  GMSDATA_tgrowarray_DOT_destroy};

/* Class descriptor for 'tgrowarrayfxd' */
const SYSTEM_classdescriptor_t GMSDATA_tgrowarrayfxd_CD = {
  _P3str1("\015tgrowarrayfxd"), 
  &GMSDATA_tgrowarray_CD, NULL, 0, 
  sizeof(GMSDATA_tgrowarrayfxd_OD), GMSDATA_tgrowarrayfxd_VT, NULL};


void * const GMSDATA_txintlist_VT[] = {(void*)&
  GMSDATA_txintlist_DOT_destroy};

/* Class descriptor for 'txintlist' */
const SYSTEM_classdescriptor_t GMSDATA_txintlist_CD = {
  _P3str1("\011txintlist"), 
  &GMSDATA_tgrowarrayfxd_CD, NULL, 0, 
  sizeof(GMSDATA_txintlist_OD), GMSDATA_txintlist_VT, NULL};


void * const GMSDATA_ttblgamsdata_VT[] = {(void*)&
  GMSDATA_ttblgamsdata_DOT_destroy};

/* Class descriptor for 'ttblgamsdata' */
const SYSTEM_classdescriptor_t GMSDATA_ttblgamsdata_CD = {
  _P3str1("\014ttblgamsdata"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSDATA_ttblgamsdata_OD), GMSDATA_ttblgamsdata_VT, NULL};


void * const GMSDATA_trorcmapper_VT[] = {(void*)&
  GMSDATA_ttblgamsdata_DOT_destroy};

/* Class descriptor for 'trorcmapper' */
const SYSTEM_classdescriptor_t GMSDATA_trorcmapper_CD = {
  _P3str1("\013trorcmapper"), 
  &GMSDATA_ttblgamsdata_CD, NULL, 0, 
  sizeof(GMSDATA_trorcmapper_OD), GMSDATA_trorcmapper_VT, NULL};


Constructor(GMSDATA_tgrowarray ) GMSDATA_tgrowarray_DOT_create(
  GMSDATA_tgrowarray self)
{
  ValueCast(GMSDATA_tgrowarray,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSDATA_tgrowarray_DOT_baseallocated = 0;
  self->GMSDATA_tgrowarray_DOT_pbase = NULL;
  self->GMSDATA_tgrowarray_DOT_pcurrentbuf = NULL;
  self->GMSDATA_tgrowarray_DOT_baseused =  -1;
  return self;
}  /* create */

Destructor(GMSDATA_tgrowarray ) GMSDATA_tgrowarray_DOT_destroy(
  GMSDATA_tgrowarray self)
{
  GMSDATA_tgrowarray_DOT_clear(self);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSDATA_tgrowarray_DOT_clear(
  GMSDATA_tgrowarray self)
{
  while (self->GMSDATA_tgrowarray_DOT_baseused >= 0) {
    _P3freemem((*self->GMSDATA_tgrowarray_DOT_pbase)[self->
      GMSDATA_tgrowarray_DOT_baseused]);
    self->GMSDATA_tgrowarray_DOT_baseused = self->
      GMSDATA_tgrowarray_DOT_baseused - 1;
  
}
  SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
    GMSDATA_tgrowarray_DOT_pbase),0);
  self->GMSDATA_tgrowarray_DOT_baseallocated = 0;
  self->GMSDATA_tgrowarray_DOT_pcurrentbuf = NULL;
}  /* clear */

Function(SYSTEM_int64 ) GMSDATA_tgrowarray_DOT_memoryused(
  GMSDATA_tgrowarray self)
{
  SYSTEM_int64 result;

  if (self->GMSDATA_tgrowarray_DOT_pcurrentbuf == NULL) { 
    result = 0;
  } else 
    result = self->GMSDATA_tgrowarray_DOT_baseallocated * sizeof(
      SYSTEM_pointer) + self->GMSDATA_tgrowarray_DOT_baseused * 
      GMSDATA_bufsize + self->GMSDATA_tgrowarray_DOT_pcurrentbuf->
      bytesused;
  return result;
}  /* memoryused */

Function(SYSTEM_pointer ) GMSDATA_tgrowarray_DOT_reservemem(
  GMSDATA_tgrowarray self,
  SYSTEM_integer l)
{
  SYSTEM_pointer result;

  if (self->GMSDATA_tgrowarray_DOT_pcurrentbuf == NULL || self->
    GMSDATA_tgrowarray_DOT_pcurrentbuf->bytesused + l > 
    GMSDATA_bufsize) {
    self->GMSDATA_tgrowarray_DOT_baseused = self->
      GMSDATA_tgrowarray_DOT_baseused + 1;
    if (self->GMSDATA_tgrowarray_DOT_baseused >= self->
      GMSDATA_tgrowarray_DOT_baseallocated) {
      if (self->GMSDATA_tgrowarray_DOT_baseallocated == 0) { 
        self->GMSDATA_tgrowarray_DOT_baseallocated = 32;
      } else 
        self->GMSDATA_tgrowarray_DOT_baseallocated = 2 * self->
          GMSDATA_tgrowarray_DOT_baseallocated;
      SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
        GMSDATA_tgrowarray_DOT_pbase),self->
        GMSDATA_tgrowarray_DOT_baseallocated * sizeof(SYSTEM_pointer));
    } 
    _P3getmem(self->GMSDATA_tgrowarray_DOT_pcurrentbuf,sizeof(
      GMSDATA_tgadatabuffer));
    (*self->GMSDATA_tgrowarray_DOT_pbase)[self->
      GMSDATA_tgrowarray_DOT_baseused] = self->
      GMSDATA_tgrowarray_DOT_pcurrentbuf;
    self->GMSDATA_tgrowarray_DOT_pcurrentbuf->bytesused = 0;
  } 
  { register GMSDATA_tgadatabuffer *_W2=self->
    GMSDATA_tgrowarray_DOT_pcurrentbuf;
    result = ValueCast(SYSTEM_pointer,&_W2->buffer[_W2->bytesused]);
    _W2->bytesused = _W2->bytesused + l;

  }
  return result;
}  /* reservemem */

Function(SYSTEM_pointer ) GMSDATA_tgrowarray_DOT_reserveandclear(
  GMSDATA_tgrowarray self,
  SYSTEM_integer l)
{
  SYSTEM_pointer result;

  result = GMSDATA_tgrowarray_DOT_reservemem(self,l);
  SYSTEM_P3_fillchar(result,l,0);
  return result;
}  /* reserveandclear */

Constructor(GMSDATA_tgrowarrayfxd ) GMSDATA_tgrowarrayfxd_DOT_create(
  GMSDATA_tgrowarrayfxd self,
  SYSTEM_integer asize)
{
  ValueCast(GMSDATA_tgrowarrayfxd,GMSDATA_tgrowarray_DOT_create(ValueCast(
    GMSDATA_tgrowarray,self)));
  self->GMSDATA_tgrowarrayfxd_DOT_fsize = asize;
  self->GMSDATA_tgrowarrayfxd_DOT_fstorefact = GMSDATA_bufsize /  self->
    GMSDATA_tgrowarrayfxd_DOT_fsize;
  return self;
}  /* create */

Function(SYSTEM_pointer ) GMSDATA_tgrowarrayfxd_DOT_additem(
  GMSDATA_tgrowarrayfxd self,
  const SYSTEM_untyped *r)
{
  SYSTEM_pointer result;

  result = GMSDATA_tgrowarray_DOT_reservemem(ValueCast(
    GMSDATA_tgrowarray,self),self->GMSDATA_tgrowarrayfxd_DOT_fsize);
  GMSOBJ_cmove(r,ValueCast(SYSTEM_P3_pbyte,result),self->
    GMSDATA_tgrowarrayfxd_DOT_fsize);
  self->GMSDATA_tgrowarrayfxd_DOT_fcount = self->
    GMSDATA_tgrowarrayfxd_DOT_fcount + 1;
  return result;
}  /* additem */

Function(GMSGEN_pbytedataarray ) 
  GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(
  GMSDATA_tgrowarrayfxd self,
  SYSTEM_integer n)
{
  GMSGEN_pbytedataarray result;

  result = ValueCast(GMSGEN_pbytedataarray,&(*self->
    GMSDATA_tgrowarray_DOT_pbase)[n /  self->
    GMSDATA_tgrowarrayfxd_DOT_fstorefact]->buffer[n % self->
    GMSDATA_tgrowarrayfxd_DOT_fstorefact * self->
    GMSDATA_tgrowarrayfxd_DOT_fsize]);
  return result;
}  /* getitemptrindx */

Procedure GMSDATA_tgrowarrayfxd_DOT_getitem(
  GMSDATA_tgrowarrayfxd self,
  SYSTEM_integer n,
  SYSTEM_untyped *r)
{
  GMSGEN_pbytedataarray pb;

  pb = GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self,n);
  GMSOBJ_cmove(&(*pb)[0],r,self->GMSDATA_tgrowarrayfxd_DOT_fsize);
}  /* getitem */

Procedure GMSDATA_tgrowarrayfxd_DOT_clear(
  GMSDATA_tgrowarrayfxd self)
{
  GMSDATA_tgrowarray_DOT_clear(ValueCast(GMSDATA_tgrowarray,self));
  self->GMSDATA_tgrowarrayfxd_DOT_fcount = 0;
}  /* clear */

Constructor(GMSDATA_ttblgamsdata ) GMSDATA_ttblgamsdata_DOT_create(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer adim,
  SYSTEM_integer adatasize)
{
  ValueCast(GMSDATA_ttblgamsdata,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSDATA_ttblgamsdata_DOT_ds = ValueCast(GMSDATA_tgrowarrayfxd,
    GMSDATA_tgrowarrayfxd_DOT_create(ValueCast(GMSDATA_tgrowarrayfxd,
    _P3alloc_object(&GMSDATA_tgrowarrayfxd_CD)),adim * sizeof(
    SYSTEM_longint) + adatasize));
  self->GMSDATA_ttblgamsdata_DOT__fdim = adim;
  self->GMSDATA_ttblgamsdata_DOT_findexsize = self->
    GMSDATA_ttblgamsdata_DOT__fdim * sizeof(SYSTEM_longint);
  self->GMSDATA_ttblgamsdata_DOT_fdatasize = adatasize;
  self->GMSDATA_ttblgamsdata_DOT_flist = ValueCast(GMSOBJ_txlist,
    SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,_P3alloc_object(&
    GMSOBJ_txlist_CD))));
  self->GMSDATA_ttblgamsdata_DOT_fissorted = SYSTEM_true;
  self->GMSDATA_ttblgamsdata_DOT_flastindex =  -1;
  return self;
}  /* create */

Destructor(GMSDATA_ttblgamsdata ) GMSDATA_ttblgamsdata_DOT_destroy(
  GMSDATA_ttblgamsdata self)
{
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GMSDATA_ttblgamsdata_DOT_ds));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GMSDATA_ttblgamsdata_DOT_flist));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSDATA_ttblgamsdata_DOT_addrecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  const SYSTEM_untyped *buffer)
{
  GMSDATA_ttblgamsdata_DOT_insertrecord(self,self->
    GMSDATA_ttblgamsdata_DOT_flist->GMSOBJ_txlist_DOT_fcount,inx,
    buffer);
}  /* addrecord */

Procedure GMSDATA_ttblgamsdata_DOT_insertrecord(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  const SYSTEM_integer *inx,
  const SYSTEM_untyped *buffer)
{
  GMSGEN_pbytedataarray p;

  p = ValueCast(GMSGEN_pbytedataarray,
    GMSDATA_tgrowarray_DOT_reservemem(ValueCast(GMSDATA_tgrowarray,
    self->GMSDATA_ttblgamsdata_DOT_ds),self->
    GMSDATA_ttblgamsdata_DOT_ds->GMSDATA_tgrowarrayfxd_DOT_fsize));
  GMSOBJ_cmove(inx,&(*p)[0],self->
    GMSDATA_ttblgamsdata_DOT_findexsize);
  GMSOBJ_cmove(buffer,&(*p)[self->GMSDATA_ttblgamsdata_DOT_findexsize],
    self->GMSDATA_ttblgamsdata_DOT_fdatasize);
  GMSOBJ_txlist_DOT_insert(self->GMSDATA_ttblgamsdata_DOT_flist,n,p);
  self->GMSDATA_ttblgamsdata_DOT_fissorted = SYSTEM_false;
}  /* insertrecord */

Procedure GMSDATA_ttblgamsdata_DOT_sort(
  GMSDATA_ttblgamsdata self)
{
  SYSTEM_integer n;
  SYSTEM_boolean sortneeded;

  if (!self->GMSDATA_ttblgamsdata_DOT_fissorted) {
    sortneeded = SYSTEM_false;
    { register SYSTEM_int32 _stop = self->
        GMSDATA_ttblgamsdata_DOT_flist->GMSOBJ_txlist_DOT_fcount - 2;
      if ((n = 0) <=  _stop) do {
        if (GMSDATA_ttblgamsdata_DOT_compare(self,n,n + 1) > 0) {
          sortneeded = SYSTEM_true;
          SYSTEM_break(BRK_1);
        } 
CNT_1:;
      } while (n++ !=  _stop);
BRK_1:;

    }
    if (sortneeded) 
      GMSDATA_ttblgamsdata_DOT_quicksort(self,0,self->
        GMSDATA_ttblgamsdata_DOT_flist->GMSOBJ_txlist_DOT_fcount - 1);
    self->GMSDATA_ttblgamsdata_DOT_fissorted = SYSTEM_true;
  } 
}  /* sort */

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_compare(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GMSSPECS_ptindex p1, p2;

  result = 0;
  p1 = ValueCast(GMSSPECS_ptindex,GMSOBJ_txlist_DOT_get(self->
    GMSDATA_ttblgamsdata_DOT_flist,index1));
  p2 = ValueCast(GMSSPECS_ptindex,GMSOBJ_txlist_DOT_get(self->
    GMSDATA_ttblgamsdata_DOT_flist,index2));
  { register SYSTEM_int32 _stop = self->GMSDATA_ttblgamsdata_DOT__fdim;
    if ((d = 1) <=  _stop) do {
      result = (*p1)[d - 1] - (*p2)[d - 1];
      if (result != 0) 
        return result;
    
    } while (d++ !=  _stop);

  }
  return result;
}  /* compare */

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_comparewithrecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  SYSTEM_integer n)
{
  SYSTEM_integer result;
  SYSTEM_integer d;
  GMSSPECS_ptindex p1;

  result = 0;
  p1 = ValueCast(GMSSPECS_ptindex,GMSOBJ_txlist_DOT_get(self->
    GMSDATA_ttblgamsdata_DOT_flist,n));
  { register SYSTEM_int32 _stop = self->GMSDATA_ttblgamsdata_DOT__fdim;
    if ((d = 1) <=  _stop) do {
      result = inx[d - 1] - (*p1)[d - 1];
      if (result != 0) 
        return result;
    
    } while (d++ !=  _stop);

  }
  return result;
}  /* comparewithrecord */

Procedure GMSDATA_ttblgamsdata_DOT_exchange(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_pointer p;

  p = GMSOBJ_txlist_DOT_get(self->GMSDATA_ttblgamsdata_DOT_flist,
    index1);
  GMSOBJ_txlist_DOT_put(self->GMSDATA_ttblgamsdata_DOT_flist,index1,
    GMSOBJ_txlist_DOT_get(self->GMSDATA_ttblgamsdata_DOT_flist,index2));
  GMSOBJ_txlist_DOT_put(self->GMSDATA_ttblgamsdata_DOT_flist,index2,p);
}  /* exchange */

Procedure GMSDATA_ttblgamsdata_DOT_quicksort(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer l,
  SYSTEM_integer r)
{
  SYSTEM_integer i, j, p;

  do {
    i = l;
    j = r;
    p = ValueCast(SYSTEM_uint32,l + r) >> 1;
    do {
      while (GMSDATA_ttblgamsdata_DOT_compare(self,i,p) < 0) {

        _P3inc0(i);
}
      while (GMSDATA_ttblgamsdata_DOT_compare(self,j,p) > 0) {

        _P3dec0(j);
}
      if (i <= j) {
        GMSDATA_ttblgamsdata_DOT_exchange(self,i,j);
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
      GMSDATA_ttblgamsdata_DOT_quicksort(self,l,j);
    l = i;
  } while (!(i >= r));
}  /* quicksort */

Procedure GMSDATA_ttblgamsdata_DOT_getrecord(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  SYSTEM_integer *inx,
  SYSTEM_untyped *buffer)
{
  GMSGEN_pbytedataarray p;

  p = ValueCast(GMSGEN_pbytedataarray,GMSOBJ_txlist_DOT_get(self->
    GMSDATA_ttblgamsdata_DOT_flist,n));
  GMSOBJ_cmove(&(*p)[0],inx,self->
    GMSDATA_ttblgamsdata_DOT_findexsize);
  GMSOBJ_cmove(&(*p)[self->GMSDATA_ttblgamsdata_DOT_findexsize],buffer,
    self->GMSDATA_ttblgamsdata_DOT_fdatasize);
}  /* getrecord */

Procedure GMSDATA_ttblgamsdata_DOT_getkeys(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  SYSTEM_integer *inx)
{
  GMSOBJ_cmove(&(*ValueCast(GMSGEN_pbytedataarray,
    GMSOBJ_txlist_DOT_get(self->GMSDATA_ttblgamsdata_DOT_flist,n)))[0],
    inx,self->GMSDATA_ttblgamsdata_DOT_findexsize);
}  /* getkeys */

Procedure GMSDATA_ttblgamsdata_DOT_getdata(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  SYSTEM_untyped *buffer)
{
  GMSOBJ_cmove(&(*ValueCast(GMSGEN_pbytedataarray,
    GMSOBJ_txlist_DOT_get(self->GMSDATA_ttblgamsdata_DOT_flist,n)))[
    self->GMSDATA_ttblgamsdata_DOT_findexsize],buffer,self->
    GMSDATA_ttblgamsdata_DOT_fdatasize);
}  /* getdata */

Function(SYSTEM_pointer ) GMSDATA_ttblgamsdata_DOT_getdataptr(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n)
{
  SYSTEM_pointer result;

  result = ValueCast(SYSTEM_pointer,&(*ValueCast(
    SYSTEM_P3_pintegerarray,GMSOBJ_txlist_DOT_get(self->
    GMSDATA_ttblgamsdata_DOT_flist,n)))[self->
    GMSDATA_ttblgamsdata_DOT__fdim]);
  return result;
}  /* getdataptr */

Function(SYSTEM_boolean ) GMSDATA_ttblgamsdata_DOT_searchrecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  SYSTEM_integer *recnr)
{
  SYSTEM_boolean result;
  SYSTEM_integer l, h, i, c;

  result = SYSTEM_false;
  h = GMSDATA_ttblgamsdata_DOT_getcount(self) - 1;
  if (h < 0) {
    *recnr = 0;
    self->GMSDATA_ttblgamsdata_DOT_flastindex = 0;
    return result;
  } 
  l = 0;
  self->GMSDATA_ttblgamsdata_DOT_flastindex = self->
    GMSDATA_ttblgamsdata_DOT_flastindex + 1;
  if (self->GMSDATA_ttblgamsdata_DOT_flastindex >= 0 && self->
    GMSDATA_ttblgamsdata_DOT_flastindex <= h) {
    c = GMSDATA_ttblgamsdata_DOT_comparewithrecord(self,inx,self->
      GMSDATA_ttblgamsdata_DOT_flastindex);
    if (c == 0) {
      result = SYSTEM_true;
      *recnr = self->GMSDATA_ttblgamsdata_DOT_flastindex;
      return result;
    } 
    if (c < 0) { 
      h = self->GMSDATA_ttblgamsdata_DOT_flastindex - 1;
    } else 
      l = self->GMSDATA_ttblgamsdata_DOT_flastindex + 1;
  } 
  while (l <= h) {
    i = ValueCast(SYSTEM_uint32,l + h) >> 1;
    c = GMSDATA_ttblgamsdata_DOT_comparewithrecord(self,inx,i);
    if (c > 0) { 
      l = i + 1;
    } else 
      if (c != 0) { 
        h = i - 1;
      } else {
        result = SYSTEM_true;
        l = i;
        SYSTEM_break(BRK_2);
      } 
  
CNT_2:;
  }
BRK_2:;
  *recnr = l;
  self->GMSDATA_ttblgamsdata_DOT_flastindex = l;
  return result;
}  /* searchrecord */

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_getcount(
  GMSDATA_ttblgamsdata self)
{
  SYSTEM_integer result;

  result = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  return result;
}  /* getcount */

Function(SYSTEM_boolean ) GMSDATA_ttblgamsdata_DOT_adduniquerecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  const SYSTEM_untyped *buffer)
{
  SYSTEM_boolean result;
  SYSTEM_integer n;

  result = !GMSDATA_ttblgamsdata_DOT_searchrecord(self,inx,&n);
  if (result) 
    GMSDATA_ttblgamsdata_DOT_insertrecord(self,n,inx,buffer);
  return result;
}  /* adduniquerecord */

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_getcapacity(
  GMSDATA_ttblgamsdata self)
{
  SYSTEM_integer result;

  result = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcapacity;
  return result;
}  /* getcapacity */

Procedure GMSDATA_ttblgamsdata_DOT_setcapacity(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n)
{
  GMSOBJ_txlist_DOT_setcapacity(self->GMSDATA_ttblgamsdata_DOT_flist,n);
}  /* setcapacity */

Procedure GMSDATA_ttblgamsdata_DOT_clear(
  GMSDATA_ttblgamsdata self)
{
  GMSDATA_tgrowarrayfxd_DOT_clear(self->GMSDATA_ttblgamsdata_DOT_ds);
  GMSOBJ_txlist_DOT_clear(self->GMSDATA_ttblgamsdata_DOT_flist);
}  /* clear */

Function(SYSTEM_int64 ) GMSDATA_ttblgamsdata_DOT_memoryused(
  GMSDATA_ttblgamsdata self)
{
  SYSTEM_int64 result;

  result = GMSDATA_tgrowarray_DOT_memoryused(ValueCast(
    GMSDATA_tgrowarray,self->GMSDATA_ttblgamsdata_DOT_ds)) + 
    GMSOBJ_txlist_DOT_memoryused(self->GMSDATA_ttblgamsdata_DOT_flist);
  return result;
}  /* memoryused */
typedef struct GMSDATA_trorcrecord_S *GMSDATA_prorcrecord;
typedef struct GMSDATA_trorcrecord_S {
  SYSTEM_integer orgindx,newindx;
} GMSDATA_trorcrecord;


Constructor(GMSDATA_trorcmapper ) GMSDATA_trorcmapper_DOT_create(
  GMSDATA_trorcmapper self,
  SYSTEM_integer adim)
{
  ValueCast(GMSDATA_trorcmapper,GMSDATA_ttblgamsdata_DOT_create(ValueCast(
    GMSDATA_ttblgamsdata,self),adim,sizeof(GMSDATA_trorcrecord)));
  return self;
}  /* create */
cnstdef {GMSDATA_ndummyrecord = 2137483647};

Procedure GMSDATA_trorcmapper_DOT_adddummyrecord(
  GMSDATA_trorcmapper self)
{
  GMSSPECS_tindex index;
  SYSTEM_integer d;
  GMSDATA_trorcrecord r;

  r.orgindx = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  r.newindx = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  { register SYSTEM_int32 _stop = self->GMSDATA_ttblgamsdata_DOT__fdim - 1;
    if ((d = 1) <=  _stop) do {
      index[d - 1] = 0;
    } while (d++ !=  _stop);

  }
  index[self->GMSDATA_ttblgamsdata_DOT__fdim - 1] = 
    GMSDATA_ndummyrecord + self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  GMSDATA_ttblgamsdata_DOT_adduniquerecord(ValueCast(
    GMSDATA_ttblgamsdata,self),index,&r);
}  /* adddummyrecord */

Procedure GMSDATA_trorcmapper_DOT_addrecord(
  GMSDATA_trorcmapper self,
  const SYSTEM_integer *inx)
{
  GMSDATA_trorcrecord r;

  r.orgindx = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  r.newindx = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  GMSDATA_ttblgamsdata_DOT_adduniquerecord(ValueCast(
    GMSDATA_ttblgamsdata,self),inx,&r);
}  /* addrecord */

Function(SYSTEM_boolean ) GMSDATA_trorcmapper_DOT_addanyrecord(
  GMSDATA_trorcmapper self,
  const SYSTEM_integer *inx)
{
  SYSTEM_boolean result;
  GMSDATA_trorcrecord r;

  r.orgindx = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  r.newindx = self->GMSDATA_ttblgamsdata_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  result = GMSDATA_ttblgamsdata_DOT_adduniquerecord(ValueCast(
    GMSDATA_ttblgamsdata,self),inx,&r);
  if (!result) 
    GMSDATA_trorcmapper_DOT_adddummyrecord(self);
  return result;
}  /* addanyrecord */

Procedure GMSDATA_trorcmapper_DOT_finishindex(
  GMSDATA_trorcmapper self,
  SYSTEM_boolean keeppos)
{
  SYSTEM_integer n;

  if (keeppos) 
    GMSDATA_trorcmapper_DOT_sort2(self);
  { register SYSTEM_int32 _stop = self->GMSDATA_ttblgamsdata_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      (ValueCast(GMSDATA_prorcrecord,
        GMSDATA_ttblgamsdata_DOT_getdataptr(ValueCast(
        GMSDATA_ttblgamsdata,self),n)))->newindx = n;
    } while (n++ !=  _stop);

  }
  if (keeppos) 
    GMSDATA_ttblgamsdata_DOT_sort(ValueCast(GMSDATA_ttblgamsdata,self));
}  /* finishindex */

Function(SYSTEM_integer ) GMSDATA_trorcmapper_DOT_getrorc(
  GMSDATA_trorcmapper self,
  const SYSTEM_integer *inx)
{
  SYSTEM_integer result;

  if (!GMSDATA_ttblgamsdata_DOT_searchrecord(ValueCast(
    GMSDATA_ttblgamsdata,self),inx,&result)) { 
    result =  -1;
  } else 
    if (inx[self->GMSDATA_ttblgamsdata_DOT__fdim - 1] >= 
      GMSDATA_ndummyrecord) { 
      result =  -1;
    } else 
      result = (ValueCast(GMSDATA_prorcrecord,
        GMSDATA_ttblgamsdata_DOT_getdataptr(ValueCast(
        GMSDATA_ttblgamsdata,self),result)))->newindx;
  return result;
}  /* getrorc */

Function(SYSTEM_integer ) GMSDATA_trorcmapper_DOT_compare2(
  GMSDATA_trorcmapper self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_integer result;

  result = (ValueCast(GMSDATA_prorcrecord,
    GMSDATA_ttblgamsdata_DOT_getdataptr(ValueCast(GMSDATA_ttblgamsdata,
    self),index1)))->orgindx - (ValueCast(GMSDATA_prorcrecord,
    GMSDATA_ttblgamsdata_DOT_getdataptr(ValueCast(GMSDATA_ttblgamsdata,
    self),index2)))->orgindx;
  return result;
}  /* compare2 */

Procedure GMSDATA_trorcmapper_DOT_quicksort2(
  GMSDATA_trorcmapper self,
  SYSTEM_integer l,
  SYSTEM_integer r)
{
  SYSTEM_integer i, j, p;

  do {
    i = l;
    j = r;
    p = ValueCast(SYSTEM_uint32,l + r) >> 1;
    do {
      while (GMSDATA_trorcmapper_DOT_compare2(self,i,p) < 0) {

        _P3inc0(i);
}
      while (GMSDATA_trorcmapper_DOT_compare2(self,j,p) > 0) {

        _P3dec0(j);
}
      if (i <= j) {
        GMSDATA_ttblgamsdata_DOT_exchange(ValueCast(
          GMSDATA_ttblgamsdata,self),i,j);
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
      GMSDATA_trorcmapper_DOT_quicksort2(self,l,j);
    l = i;
  } while (!(i >= r));
}  /* quicksort2 */

Procedure GMSDATA_trorcmapper_DOT_sort2(
  GMSDATA_trorcmapper self)
{
  if (self->GMSDATA_ttblgamsdata_DOT_flist->GMSOBJ_txlist_DOT_fcount > 1) 
    GMSDATA_trorcmapper_DOT_quicksort2(self,0,self->
      GMSDATA_ttblgamsdata_DOT_flist->GMSOBJ_txlist_DOT_fcount - 1);
  self->GMSDATA_ttblgamsdata_DOT_fissorted = SYSTEM_false;
}  /* sort2 */

Constructor(GMSDATA_txintlist ) GMSDATA_txintlist_DOT_create(
  GMSDATA_txintlist self)
{
  ValueCast(GMSDATA_txintlist,GMSDATA_tgrowarrayfxd_DOT_create(ValueCast(
    GMSDATA_tgrowarrayfxd,self),sizeof(SYSTEM_longint)));
  self->GMSDATA_tgrowarrayfxd_DOT_fcount = 0;
  return self;
}  /* create */

Function(SYSTEM_integer ) GMSDATA_txintlist_DOT_add(
  GMSDATA_txintlist self,
  SYSTEM_integer item)
{
  SYSTEM_integer result;

  result = self->GMSDATA_tgrowarrayfxd_DOT_fcount;
  GMSDATA_tgrowarrayfxd_DOT_additem(ValueCast(GMSDATA_tgrowarrayfxd,
    self),&item);
  return result;
}  /* add */

Destructor(GMSDATA_txintlist ) GMSDATA_txintlist_DOT_destroy(
  GMSDATA_txintlist self)
{
  GMSDATA_tgrowarray_DOT_destroy(ValueCast(GMSDATA_tgrowarray,self));
  return self;
}  /* destroy */

Procedure GMSDATA_txintlist_DOT_exchange(
  GMSDATA_txintlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_P3_pinteger p1, p2;
  SYSTEM_integer t;

  p1 = ValueCast(SYSTEM_P3_pinteger,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(ValueCast(
    GMSDATA_tgrowarrayfxd,self),index1));
  p2 = ValueCast(SYSTEM_P3_pinteger,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(ValueCast(
    GMSDATA_tgrowarrayfxd,self),index2));
  t = *p1;
  *p1 = *p2;
  *p2 = t;
}  /* exchange */

Function(SYSTEM_integer ) GMSDATA_txintlist_DOT_getitems(
  GMSDATA_txintlist self,
  SYSTEM_integer index)
{
  SYSTEM_integer result;

  result = *ValueCast(SYSTEM_P3_pinteger,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(ValueCast(
    GMSDATA_tgrowarrayfxd,self),index));
  return result;
}  /* getitems */

Procedure GMSDATA_txintlist_DOT_setitems(
  GMSDATA_txintlist self,
  SYSTEM_integer index,
  SYSTEM_integer v)
{
  while (index >= self->GMSDATA_tgrowarrayfxd_DOT_fcount) {
    GMSDATA_tgrowarray_DOT_reserveandclear(ValueCast(
      GMSDATA_tgrowarray,self),self->GMSDATA_tgrowarrayfxd_DOT_fsize);
    _P3inc0(self->GMSDATA_tgrowarrayfxd_DOT_fcount);
  
}
  *ValueCast(SYSTEM_P3_pinteger,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(ValueCast(
    GMSDATA_tgrowarrayfxd,self),index)) = v;
}  /* setitems */

/* unit gmsdata */
void _Init_Module_gmsdata(void)
{
} /* _Init_Module_gmsdata */

void _Final_Module_gmsdata(void)
{
} /* _Final_Module_gmsdata */

