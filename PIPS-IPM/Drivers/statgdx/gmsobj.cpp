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


void * const GMSOBJ_tquicksortclass_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy, (void*)&_P3_abstract_call1, (void*)&
  _P3_abstract_call1};

/* Class descriptor for 'tquicksortclass' */
const SYSTEM_classdescriptor_t GMSOBJ_tquicksortclass_CD = {
  _P3str1("\017tquicksortclass"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSOBJ_tquicksortclass_OD), GMSOBJ_tquicksortclass_VT, NULL};


void * const GMSOBJ_txlist_VT[] = {(void*)&GMSOBJ_txlist_DOT_destroy, (void*)&
  GMSOBJ_txlist_DOT_exchange, (void*)&_P3_abstract_call1, (void*)&
  GMSOBJ_txlist_DOT_grow, (void*)&GMSOBJ_txlist_DOT_freeitem};

/* Class descriptor for 'txlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txlist_CD = {
  _P3str1("\006txlist"), 
  &GMSOBJ_tquicksortclass_CD, NULL, 0, 
  sizeof(GMSOBJ_txlist_OD), GMSOBJ_txlist_VT, NULL};


void * const GMSOBJ_txstrings_VT[] = {(void*)&
  GMSOBJ_txlist_DOT_destroy, (void*)&GMSOBJ_txlist_DOT_exchange, (void*)&
  GMSOBJ_txstrings_DOT_compare, (void*)&GMSOBJ_txlist_DOT_grow, (void*)&
  GMSOBJ_txstrings_DOT_freeitem};

/* Class descriptor for 'txstrings' */
const SYSTEM_classdescriptor_t GMSOBJ_txstrings_CD = {
  _P3str1("\011txstrings"), 
  &GMSOBJ_txlist_CD, NULL, 0, 
  sizeof(GMSOBJ_txstrings_OD), GMSOBJ_txstrings_VT, NULL};


void * const GMSOBJ_txpcharlist_VT[] = {(void*)&
  GMSOBJ_txpcharlist_DOT_destroy};

/* Class descriptor for 'txpcharlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txpcharlist_CD = {
  _P3str1("\013txpcharlist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSOBJ_txpcharlist_OD), GMSOBJ_txpcharlist_VT, NULL};


void * const GMSOBJ_txcustomstringlist_VT[] = {(void*)&
  GMSOBJ_txcustomstringlist_DOT_destroy, (void*)&
  GMSOBJ_txcustomstringlist_DOT_exchange, (void*)&
  GMSOBJ_txcustomstringlist_DOT_compare, (void*)&
  GMSOBJ_txcustomstringlist_DOT_grow, (void*)&
  GMSOBJ_txcustomstringlist_DOT_freeobject};

/* Class descriptor for 'txcustomstringlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txcustomstringlist_CD = {
  _P3str1("\022txcustomstringlist"), 
  &GMSOBJ_tquicksortclass_CD, NULL, 0, 
  sizeof(GMSOBJ_txcustomstringlist_OD), GMSOBJ_txcustomstringlist_VT, NULL};


void * const GMSOBJ_txstringlist_VT[] = {(void*)&
  GMSOBJ_txcustomstringlist_DOT_destroy, (void*)&
  GMSOBJ_txcustomstringlist_DOT_exchange, (void*)&
  GMSOBJ_txcustomstringlist_DOT_compare, (void*)&
  GMSOBJ_txcustomstringlist_DOT_grow, (void*)&
  GMSOBJ_txcustomstringlist_DOT_freeobject};

/* Class descriptor for 'txstringlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txstringlist_CD = {
  _P3str1("\014txstringlist"), 
  &GMSOBJ_txcustomstringlist_CD, NULL, 0, 
  sizeof(GMSOBJ_txstringlist_OD), GMSOBJ_txstringlist_VT, NULL};


void * const GMSOBJ_txsortedstringlist_VT[] = {(void*)&
  GMSOBJ_txcustomstringlist_DOT_destroy, (void*)&
  GMSOBJ_txcustomstringlist_DOT_exchange, (void*)&
  GMSOBJ_txcustomstringlist_DOT_compare, (void*)&
  GMSOBJ_txcustomstringlist_DOT_grow, (void*)&
  GMSOBJ_txcustomstringlist_DOT_freeobject};

/* Class descriptor for 'txsortedstringlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txsortedstringlist_CD = {
  _P3str1("\022txsortedstringlist"), 
  &GMSOBJ_txcustomstringlist_CD, NULL, 0, 
  sizeof(GMSOBJ_txsortedstringlist_OD), GMSOBJ_txsortedstringlist_VT, NULL};


void * const GMSOBJ_txstrstrlist_VT[] = {(void*)&
  GMSOBJ_txcustomstringlist_DOT_destroy, (void*)&
  GMSOBJ_txcustomstringlist_DOT_exchange, (void*)&
  GMSOBJ_txcustomstringlist_DOT_compare, (void*)&
  GMSOBJ_txcustomstringlist_DOT_grow, (void*)&
  GMSOBJ_txstrstrlist_DOT_freeobject};

/* Class descriptor for 'txstrstrlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txstrstrlist_CD = {
  _P3str1("\014txstrstrlist"), 
  &GMSOBJ_txsortedstringlist_CD, NULL, 0, 
  sizeof(GMSOBJ_txstrstrlist_OD), GMSOBJ_txstrstrlist_VT, NULL};


void * const GMSOBJ_txhashedstringlist_VT[] = {(void*)&
  GMSOBJ_txhashedstringlist_DOT_destroy, (void*)&
  GMSOBJ_txcustomstringlist_DOT_exchange, (void*)&
  GMSOBJ_txhashedstringlist_DOT_compare, (void*)&
  GMSOBJ_txcustomstringlist_DOT_grow, (void*)&
  GMSOBJ_txcustomstringlist_DOT_freeobject, (void*)&
  GMSOBJ_txhashedstringlist_DOT_equaltoentry, (void*)&
  GMSOBJ_txhashedstringlist_DOT_hashvalue};

/* Class descriptor for 'txhashedstringlist' */
const SYSTEM_classdescriptor_t GMSOBJ_txhashedstringlist_CD = {
  _P3str1("\022txhashedstringlist"), 
  &GMSOBJ_txcustomstringlist_CD, NULL, 0, 
  sizeof(GMSOBJ_txhashedstringlist_OD), GMSOBJ_txhashedstringlist_VT, NULL};


void * const GMSOBJ_txstrpool_VT[] = {(void*)&
  GMSOBJ_txhashedstringlist_DOT_destroy, (void*)&
  GMSOBJ_txcustomstringlist_DOT_exchange, (void*)&
  GMSOBJ_txstrpool_DOT_compare, (void*)&
  GMSOBJ_txcustomstringlist_DOT_grow, (void*)&
  GMSOBJ_txcustomstringlist_DOT_freeobject, (void*)&
  GMSOBJ_txstrpool_DOT_equaltoentry, (void*)&
  GMSOBJ_txhashedstringlist_DOT_hashvalue};

/* Class descriptor for 'txstrpool' */
const SYSTEM_classdescriptor_t GMSOBJ_txstrpool_CD = {
  _P3str1("\011txstrpool"), 
  &GMSOBJ_txhashedstringlist_CD, NULL, 0, 
  sizeof(GMSOBJ_txstrpool_OD), GMSOBJ_txstrpool_VT, NULL};


void * const GMSOBJ_tbooleanbitarray_VT[] = {(void*)&
  GMSOBJ_tbooleanbitarray_DOT_destroy};

/* Class descriptor for 'tbooleanbitarray' */
const SYSTEM_classdescriptor_t GMSOBJ_tbooleanbitarray_CD = {
  _P3str1("\020tbooleanbitarray"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(GMSOBJ_tbooleanbitarray_OD), GMSOBJ_tbooleanbitarray_VT, NULL};


Procedure GMSOBJ_cmove(
  const SYSTEM_untyped *src,
  SYSTEM_untyped *dest,
  SYSTEM_integer len)
{
  GMSGEN_pbytedataarray psrc, pdest;
  SYSTEM_integer k;

  if (len > 32) { 
    SYSTEM_move(src,dest,len);
  } else {
    psrc = ValueCast(GMSGEN_pbytedataarray,src);
    pdest = ValueCast(GMSGEN_pbytedataarray,dest);
    switch (len) {
      case 0: 
        break;
      case 1: 
        (*pdest)[0] = (*psrc)[0];
        break;
      case 2: 
        (*pdest)[0] = (*psrc)[0];
        (*pdest)[1] = (*psrc)[1];
        break;
      case 3: 
        (*pdest)[0] = (*psrc)[0];
        (*pdest)[1] = (*psrc)[1];
        (*pdest)[2] = (*psrc)[2];
        break;
      case 4: 
        (*pdest)[0] = (*psrc)[0];
        (*pdest)[1] = (*psrc)[1];
        (*pdest)[2] = (*psrc)[2];
        (*pdest)[3] = (*psrc)[3];
        break;
      default:
        { register SYSTEM_int32 _stop = len - 1;
          if ((k = 0) <=  _stop) do {
            (*pdest)[k] = (*psrc)[k];
          } while (k++ !=  _stop);

        }
    }
  } 
}  /* cmove */
typedef struct GMSOBJ_tiprec_S {
  union{
    struct{
      SYSTEM_nativeuint n;
    } _c1;
    struct{
      SYSTEM_pointer p;
    } _c2;
  } _u;
} GMSOBJ_tiprec;


Function(SYSTEM_pointer ) GMSOBJ_copyint2ptr(
  SYSTEM_integer i)
{
  SYSTEM_pointer result;
  SYSTEM_uint32 u;
  GMSOBJ_tiprec r;

  u = ValueCast(SYSTEM_uint32,i);
  r._u._c1.n = u;
  result = r._u._c2.p;
  return result;
}  /* copyint2ptr */

Function(SYSTEM_integer ) GMSOBJ_copyptr2int(
  SYSTEM_pointer p)
{
  SYSTEM_integer result;
  SYSTEM_uint32 u;
  GMSOBJ_tiprec r;

  r._u._c2.p = p;
  u = ValueCast(SYSTEM_uint32,r._u._c1.n);
  result = ValueCast(SYSTEM_int32,u);
  return result;
}  /* copyptr2int */

Procedure GMSOBJ_tquicksortclass_DOT_quicksort(
  GMSOBJ_tquicksortclass self,
  SYSTEM_integer l,
  SYSTEM_integer r)
{
  SYSTEM_integer i, j, p;

  do {
    i = l;
    j = r;
    p = ValueCast(SYSTEM_uint32,l + r) >> 1;
    do {
      while (VirtMethodCall(self, GMSOBJ_tquicksortclass_DOT_compare_T, 2, (
        self,i,p)) < 0) {

        _P3inc0(i);
}
      while (VirtMethodCall(self, GMSOBJ_tquicksortclass_DOT_compare_T, 2, (
        self,j,p)) > 0) {

        _P3dec0(j);
}
      if (i <= j) {
        VirtMethodCall(self, GMSOBJ_tquicksortclass_DOT_exchange_T, 1, (
          self,i,j));
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
      GMSOBJ_tquicksortclass_DOT_quicksort(self,l,j);
    l = i;
  } while (!(i >= r));
}  /* quicksort */

Procedure GMSOBJ_tquicksortclass_DOT_sortn(
  GMSOBJ_tquicksortclass self,
  SYSTEM_integer n)
{
  if (n > 1) 
    GMSOBJ_tquicksortclass_DOT_quicksort(self,SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased),n - 1 + SYSTEM_ord(
      self->GMSOBJ_tquicksortclass_DOT_onebased));
}  /* sortn */

Destructor(GMSOBJ_txlist ) GMSOBJ_txlist_DOT_destroy(
  GMSOBJ_txlist self)
{
  GMSOBJ_txlist_DOT_clear(self);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSOBJ_txlist_DOT_freeitem(
  GMSOBJ_txlist self,
  SYSTEM_integer index)
{
}  /* freeitem */

Function(SYSTEM_integer ) GMSOBJ_txlist_DOT_add(
  GMSOBJ_txlist self,
  SYSTEM_pointer item)
{
  SYSTEM_integer result;

  result = self->GMSOBJ_txlist_DOT_fcount;
  if (result == self->GMSOBJ_txlist_DOT_fcapacity) 
    VirtMethodCall(self, GMSOBJ_txlist_DOT_grow_T, 3, (self));
  (*self->GMSOBJ_txlist_DOT_flist)[result] = item;
  _P3inc0(self->GMSOBJ_txlist_DOT_fcount);
  if (self->GMSOBJ_tquicksortclass_DOT_onebased) 
    _P3inc0(result);
  return result;
}  /* add */

Function(SYSTEM_pointer ) GMSOBJ_txlist_DOT_getlast(
  GMSOBJ_txlist self)
{
  SYSTEM_pointer result;

  if (self->GMSOBJ_txlist_DOT_fcount <= 0) { 
    result = NULL;
  } else 
    result = (*self->GMSOBJ_txlist_DOT_flist)[self->
      GMSOBJ_txlist_DOT_fcount - 1];
  return result;
}  /* getlast */

Procedure GMSOBJ_txlist_DOT_clear(
  GMSOBJ_txlist self)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased);
    if ((n = self->GMSOBJ_txlist_DOT_fcount - 1 + SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased)) >=  _stop) do {
      VirtMethodCall(self, GMSOBJ_txlist_DOT_freeitem_T, 4, (self,n));
    } while (n-- !=  _stop);

  }
  self->GMSOBJ_txlist_DOT_fcount = 0;
  GMSOBJ_txlist_DOT_setcapacity(self,0);
}  /* clear */

Procedure GMSOBJ_txlist_DOT_delete(
  GMSOBJ_txlist self,
  SYSTEM_integer index)
{
  VirtMethodCall(self, GMSOBJ_txlist_DOT_freeitem_T, 4, (self,index));
  _P3dec0(self->GMSOBJ_txlist_DOT_fcount);
  if (index < self->GMSOBJ_txlist_DOT_fcount) {
    if (self->GMSOBJ_tquicksortclass_DOT_onebased) 
      _P3dec0(index);
    SYSTEM_move(&(*self->GMSOBJ_txlist_DOT_flist)[index + 1],&(*
      self->GMSOBJ_txlist_DOT_flist)[index],(self->
      GMSOBJ_txlist_DOT_fcount - index) * sizeof(SYSTEM_pointer));
  } 
}  /* delete */

Procedure GMSOBJ_txlist_DOT_exchange(
  GMSOBJ_txlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_pointer item;

  if (self->GMSOBJ_tquicksortclass_DOT_onebased) {
    _P3dec0(index1);
    _P3dec0(index2);
  } 
  item = (*self->GMSOBJ_txlist_DOT_flist)[index1];
  (*self->GMSOBJ_txlist_DOT_flist)[index1] = (*self->
    GMSOBJ_txlist_DOT_flist)[index2];
  (*self->GMSOBJ_txlist_DOT_flist)[index2] = item;
}  /* exchange */

Function(SYSTEM_pointer ) GMSOBJ_txlist_DOT_get(
  GMSOBJ_txlist self,
  SYSTEM_integer index)
{
  SYSTEM_pointer result;

  result = (*self->GMSOBJ_txlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)];
  return result;
}  /* get */

Procedure GMSOBJ_txlist_DOT_grow(
  GMSOBJ_txlist self)
{
  SYSTEM_integer delta;

  if (self->GMSOBJ_txlist_DOT_fcapacity >= 1048576) { 
    delta = self->GMSOBJ_txlist_DOT_fcapacity /  4;
  } else 
    if (self->GMSOBJ_txlist_DOT_fcapacity == 0) { 
      delta = 16;
    } else 
      delta = 7 * self->GMSOBJ_txlist_DOT_fcapacity;
  GMSOBJ_txlist_DOT_setcapacity(self,self->GMSOBJ_txlist_DOT_fcapacity + 
    delta);
}  /* grow */

Function(SYSTEM_integer ) GMSOBJ_txlist_DOT_indexof(
  GMSOBJ_txlist self,
  SYSTEM_pointer item)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if ((*self->GMSOBJ_txlist_DOT_flist)[n] == item) {
        result = n + SYSTEM_ord(self->
          GMSOBJ_tquicksortclass_DOT_onebased);
        return result;
      } 
    } while (n++ !=  _stop);

  }
  result =  -1;
  return result;
}  /* indexof */

Procedure GMSOBJ_txlist_DOT_insert(
  GMSOBJ_txlist self,
  SYSTEM_integer index,
  SYSTEM_pointer item)
{
  if (self->GMSOBJ_txlist_DOT_fcount == self->
    GMSOBJ_txlist_DOT_fcapacity) 
    VirtMethodCall(self, GMSOBJ_txlist_DOT_grow_T, 3, (self));
  if (self->GMSOBJ_tquicksortclass_DOT_onebased) 
    _P3dec0(index);
  if (index < self->GMSOBJ_txlist_DOT_fcount) 
    SYSTEM_move(&(*self->GMSOBJ_txlist_DOT_flist)[index],&(*self->
      GMSOBJ_txlist_DOT_flist)[index + 1],(self->
      GMSOBJ_txlist_DOT_fcount - index) * sizeof(SYSTEM_pointer));
  (*self->GMSOBJ_txlist_DOT_flist)[index] = item;
  _P3inc0(self->GMSOBJ_txlist_DOT_fcount);
}  /* insert */

Procedure GMSOBJ_txlist_DOT_put(
  GMSOBJ_txlist self,
  SYSTEM_integer index,
  SYSTEM_pointer item)
{
  VirtMethodCall(self, GMSOBJ_txlist_DOT_freeitem_T, 4, (self,index));
  (*self->GMSOBJ_txlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)] = item;
}  /* put */

Function(SYSTEM_integer ) GMSOBJ_txlist_DOT_remove(
  GMSOBJ_txlist self,
  SYSTEM_pointer item)
{
  SYSTEM_integer result;

  result = self->GMSOBJ_txlist_DOT_fcount - 1;
  while (result >= 0 && (*self->GMSOBJ_txlist_DOT_flist)[result] != 
    item) {

    _P3dec0(result);
}
  if (result >= SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)) 
    GMSOBJ_txlist_DOT_delete(self,result);
  return result;
}  /* remove */

Procedure GMSOBJ_txlist_DOT_setcapacity(
  GMSOBJ_txlist self,
  SYSTEM_integer newcapacity)
{
  if (newcapacity != self->GMSOBJ_txlist_DOT_fcapacity) {
    if (newcapacity < self->GMSOBJ_txlist_DOT_fcount) 
      newcapacity = self->GMSOBJ_txlist_DOT_fcount;
    SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
      GMSOBJ_txlist_DOT_flist),newcapacity * sizeof(SYSTEM_pointer));
    self->GMSOBJ_txlist_DOT_fcapacity = newcapacity;
  } 
}  /* setcapacity */

Procedure GMSOBJ_txlist_DOT_setcount(
  GMSOBJ_txlist self,
  SYSTEM_integer newcount)
{
  SYSTEM_integer i;

  if (newcount != self->GMSOBJ_txlist_DOT_fcount) {
    if (newcount > self->GMSOBJ_txlist_DOT_fcapacity) 
      GMSOBJ_txlist_DOT_setcapacity(self,newcount);
    if (newcount > self->GMSOBJ_txlist_DOT_fcount) { 
      SYSTEM_P3_fillchar(&(*self->GMSOBJ_txlist_DOT_flist)[self->
        GMSOBJ_txlist_DOT_fcount],(newcount - self->
        GMSOBJ_txlist_DOT_fcount) * sizeof(SYSTEM_pointer),0);
    } else 
      { register SYSTEM_int32 _stop = newcount;
        if ((i = self->GMSOBJ_txlist_DOT_fcount - 1) >=  _stop) do {
          VirtMethodCall(self, GMSOBJ_txlist_DOT_freeitem_T, 4, (self,i));
        } while (i-- !=  _stop);

      }
    self->GMSOBJ_txlist_DOT_fcount = newcount;
  } 
}  /* setcount */

Function(SYSTEM_pointer ) GMSOBJ_txlist_DOT_extract(
  GMSOBJ_txlist self,
  SYSTEM_pointer item)
{
  SYSTEM_pointer result;
  SYSTEM_integer i;

  result = NULL;
  i = GMSOBJ_txlist_DOT_indexof(self,item);
  if (self->GMSOBJ_tquicksortclass_DOT_onebased) 
    _P3dec0(i);
  if (i >= 0) {
    result = item;
    _P3dec0(self->GMSOBJ_txlist_DOT_fcount);
    if (i < self->GMSOBJ_txlist_DOT_fcount) 
      SYSTEM_move(&(*self->GMSOBJ_txlist_DOT_flist)[i + 1],&(*
        self->GMSOBJ_txlist_DOT_flist)[i],(self->
        GMSOBJ_txlist_DOT_fcount - i) * sizeof(SYSTEM_pointer));
  } 
  return result;
}  /* extract */

Function(SYSTEM_int64 ) GMSOBJ_txlist_DOT_memoryused(
  GMSOBJ_txlist self)
{
  SYSTEM_int64 result;

  result = self->GMSOBJ_txlist_DOT_fcapacity * sizeof(SYSTEM_pointer);
  return result;
}  /* memoryused */

Function(SYSTEM_integer ) GMSOBJ_txstrings_DOT_add(
  GMSOBJ_txstrings self,
  const SYSTEM_ansichar *item)
{
  SYSTEM_integer result;

  result = GMSOBJ_txlist_DOT_add(ValueCast(GMSOBJ_txlist,self),
    STRUTILX_newstringm(item,&self->GMSOBJ_txstrings_DOT_fstrmemory));
  return result;
}  /* add */

Function(SYSTEM_integer ) GMSOBJ_txstrings_DOT_compare(
  GMSOBJ_txstrings self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_integer result;

  result = STRUTILX_pstrucmp(ValueCast(SYSTEM_P3_pshortstring,(*self->
    GMSOBJ_txlist_DOT_flist)[index1 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)]),ValueCast(
    SYSTEM_P3_pshortstring,(*self->GMSOBJ_txlist_DOT_flist)[index2 - 
    SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)]));
  return result;
}  /* compare */

Function(SYSTEM_ansichar *) GMSOBJ_txstrings_DOT_extract(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrings self,
  const SYSTEM_ansichar *item)
{
  SYSTEM_integer n;

  n = GMSOBJ_txstrings_DOT_indexof(self,item);
  if (n < 0) { 
    _P3strclr(result);
  } else {
    _P3strcpy(result,_len_ret,item);
    GMSOBJ_txlist_DOT_delete(ValueCast(GMSOBJ_txlist,self),n);
  } 
  return result;
}  /* extract */

Procedure GMSOBJ_txstrings_DOT_freeitem(
  GMSOBJ_txstrings self,
  SYSTEM_integer index)
{
  STRUTILX_disposestringm(ValueCast(SYSTEM_P3_pshortstring,(*self->
    GMSOBJ_txlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)]),&self->
    GMSOBJ_txstrings_DOT_fstrmemory);
}  /* freeitem */

Function(SYSTEM_ansichar *) GMSOBJ_txstrings_DOT_get(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrings self,
  SYSTEM_integer index)
{
  STRUTILX_getstring(result,_len_ret,ValueCast(SYSTEM_P3_pshortstring,(*
    self->GMSOBJ_txlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)]));
  return result;
}  /* get */

Function(SYSTEM_integer ) GMSOBJ_txstrings_DOT_indexof(
  GMSOBJ_txstrings self,
  const SYSTEM_ansichar *item)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GMSOBJ_txlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if (STRUTILX_pstruequal(ValueCast(SYSTEM_P3_pshortstring,(*self->
        GMSOBJ_txlist_DOT_flist)[n]),ValueCast(SYSTEM_P3_pshortstring,
        item))) {
        result = n + SYSTEM_ord(self->
          GMSOBJ_tquicksortclass_DOT_onebased);
        return result;
      } 
    } while (n++ !=  _stop);

  }
  result =  -1;
  return result;
}  /* indexof */

Procedure GMSOBJ_txstrings_DOT_insert(
  GMSOBJ_txstrings self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *item)
{
  GMSOBJ_txlist_DOT_insert(ValueCast(GMSOBJ_txlist,self),index,
    STRUTILX_newstringm(item,&self->GMSOBJ_txstrings_DOT_fstrmemory));
}  /* insert */

Function(SYSTEM_int64 ) GMSOBJ_txstrings_DOT_memoryused(
  GMSOBJ_txstrings self)
{
  SYSTEM_int64 result;

  result = self->GMSOBJ_txstrings_DOT_fstrmemory + 
    GMSOBJ_txlist_DOT_memoryused(ValueCast(GMSOBJ_txlist,self));
  return result;
}  /* memoryused */

Procedure GMSOBJ_txstrings_DOT_put(
  GMSOBJ_txstrings self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring item;

  _P3strcpy(item,255,_ftmp1);
  VirtMethodCall(ValueCast(GMSOBJ_txlist,self), 
    GMSOBJ_txlist_DOT_freeitem_T, 4, (ValueCast(GMSOBJ_txlist,self),
    index));
  (*self->GMSOBJ_txlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)] = STRUTILX_newstringm(item,&
    self->GMSOBJ_txstrings_DOT_fstrmemory);
}  /* put */

Procedure GMSOBJ_txstrings_DOT_sort(
  GMSOBJ_txstrings self)
{
  GMSOBJ_tquicksortclass_DOT_sortn(ValueCast(GMSOBJ_tquicksortclass,
    self),self->GMSOBJ_txlist_DOT_fcount);
}  /* sort */

Constructor(GMSOBJ_txcustomstringlist ) 
  GMSOBJ_txcustomstringlist_DOT_create(
  GMSOBJ_txcustomstringlist self)
{
  ValueCast(GMSOBJ_txcustomstringlist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSOBJ_txcustomstringlist_DOT_flist = NULL;
  self->GMSOBJ_txcustomstringlist_DOT_fcount = 0;
  self->GMSOBJ_txcustomstringlist_DOT_fcapacity = 0;
  self->GMSOBJ_tquicksortclass_DOT_onebased = SYSTEM_false;
  self->GMSOBJ_txcustomstringlist_DOT_fstrmemory = 0;
  return self;
}  /* create */

Destructor(GMSOBJ_txcustomstringlist ) 
  GMSOBJ_txcustomstringlist_DOT_destroy(
  GMSOBJ_txcustomstringlist self)
{
  GMSOBJ_txcustomstringlist_DOT_clear(self);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSOBJ_txcustomstringlist_DOT_clear(
  GMSOBJ_txcustomstringlist self)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased);
    if ((n = self->GMSOBJ_txcustomstringlist_DOT_fcount - 1 + 
      SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)) >=  _stop) do {
      GMSOBJ_txcustomstringlist_DOT_freeitem(self,n);
    } while (n-- !=  _stop);

  }
  self->GMSOBJ_txcustomstringlist_DOT_fcount = 0;
  GMSOBJ_txcustomstringlist_DOT_setcapacity(self,0);
}  /* clear */

Procedure GMSOBJ_txcustomstringlist_DOT_delete(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index)
{
  GMSOBJ_txcustomstringlist_DOT_freeitem(self,index);
  if (self->GMSOBJ_tquicksortclass_DOT_onebased) 
    _P3dec0(index);
  _P3dec0(self->GMSOBJ_txcustomstringlist_DOT_fcount);
  if (index < self->GMSOBJ_txcustomstringlist_DOT_fcount) 
    SYSTEM_move(&(*self->GMSOBJ_txcustomstringlist_DOT_flist)[index + 1],&(*
      self->GMSOBJ_txcustomstringlist_DOT_flist)[index],(self->
      GMSOBJ_txcustomstringlist_DOT_fcount - index) * sizeof(
      GMSOBJ_tstringitem));
}  /* delete */

Procedure GMSOBJ_txcustomstringlist_DOT_freeitem(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index)
{
  STRUTILX_disposestringm((*self->GMSOBJ_txcustomstringlist_DOT_flist)[
    index - SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)].
    fstring,&self->GMSOBJ_txcustomstringlist_DOT_fstrmemory);
  VirtMethodCall(self, GMSOBJ_txcustomstringlist_DOT_freeobject_T, 4, (
    self,index));
}  /* freeitem */

Procedure GMSOBJ_txcustomstringlist_DOT_freeobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index)
{
}  /* freeobject */

Function(SYSTEM_ansichar *) GMSOBJ_txcustomstringlist_DOT_getname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index)
{
  STRUTILX_getstring(result,_len_ret,(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring);
  return result;
}  /* getname */

Procedure GMSOBJ_txcustomstringlist_DOT_setname(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *v)
{
  STRUTILX_strassignm(&(*self->GMSOBJ_txcustomstringlist_DOT_flist)[
    index - SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)].
    fstring,v,&self->GMSOBJ_txcustomstringlist_DOT_fstrmemory);
}  /* setname */

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_add(
  GMSOBJ_txcustomstringlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = GMSOBJ_txcustomstringlist_DOT_addobject(self,s,NULL);
  return result;
}  /* add */

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_addobject(
  GMSOBJ_txcustomstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer)
{
  SYSTEM_integer result;

  result = self->GMSOBJ_txcustomstringlist_DOT_fcount + SYSTEM_ord(
    self->GMSOBJ_tquicksortclass_DOT_onebased);
  GMSOBJ_txcustomstringlist_DOT_insertitem(self,result,s,apointer);
  return result;
}  /* addobject */

Procedure GMSOBJ_txcustomstringlist_DOT_exchange(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  GMSOBJ_tstringitem item;

  if (self->GMSOBJ_tquicksortclass_DOT_onebased) {
    _P3dec0(index1);
    _P3dec0(index2);
  } 
  item = (*self->GMSOBJ_txcustomstringlist_DOT_flist)[index1];
  (*self->GMSOBJ_txcustomstringlist_DOT_flist)[index1] = (*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index2];
  (*self->GMSOBJ_txcustomstringlist_DOT_flist)[index2] = item;
}  /* exchange */

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_compare(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_integer result;

  result = STRUTILX_pstrucmp((*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index1 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring,(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index2 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring);
  return result;
}  /* compare */

Function(SYSTEM_tobject ) GMSOBJ_txcustomstringlist_DOT_getobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index)
{
  SYSTEM_tobject result;

  result = (*self->GMSOBJ_txcustomstringlist_DOT_flist)[index - 
    SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)].fobject;
  return result;
}  /* getobject */

Procedure GMSOBJ_txcustomstringlist_DOT_grow(
  GMSOBJ_txcustomstringlist self)
{
  SYSTEM_integer delta;

  if (self->GMSOBJ_txcustomstringlist_DOT_fcapacity >= 1048576) { 
    delta = self->GMSOBJ_txcustomstringlist_DOT_fcapacity /  4;
  } else 
    if (self->GMSOBJ_txcustomstringlist_DOT_fcapacity == 0) { 
      delta = 16;
    } else 
      delta = 7 * self->GMSOBJ_txcustomstringlist_DOT_fcapacity;
  GMSOBJ_txcustomstringlist_DOT_setcapacity(self,self->
    GMSOBJ_txcustomstringlist_DOT_fcapacity + delta);
}  /* grow */

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_indexof(
  GMSOBJ_txcustomstringlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->
      GMSOBJ_txcustomstringlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if (STRUTILX_pstruequal((*self->
        GMSOBJ_txcustomstringlist_DOT_flist)[n].fstring,ValueCast(
        SYSTEM_P3_pshortstring,s))) {
        result = n + SYSTEM_ord(self->
          GMSOBJ_tquicksortclass_DOT_onebased);
        return result;
      } 
    } while (n++ !=  _stop);

  }
  result =  -1;
  return result;
}  /* indexof */

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_indexofobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_tobject aobject)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->
      GMSOBJ_txcustomstringlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      if ((*self->GMSOBJ_txcustomstringlist_DOT_flist)[n].fobject == 
        aobject) {
        result = n + SYSTEM_ord(self->
          GMSOBJ_tquicksortclass_DOT_onebased);
        return result;
      } 
    } while (n++ !=  _stop);

  }
  result =  -1;
  return result;
}  /* indexofobject */

Procedure GMSOBJ_txcustomstringlist_DOT_insertitem(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer)
{
  if (self->GMSOBJ_txcustomstringlist_DOT_fcount == self->
    GMSOBJ_txcustomstringlist_DOT_fcapacity) 
    VirtMethodCall(self, GMSOBJ_txcustomstringlist_DOT_grow_T, 3, (self));
  if (self->GMSOBJ_tquicksortclass_DOT_onebased) 
    _P3dec0(index);
  if (index < self->GMSOBJ_txcustomstringlist_DOT_fcount) 
    SYSTEM_move(&(*self->GMSOBJ_txcustomstringlist_DOT_flist)[index],&(*
      self->GMSOBJ_txcustomstringlist_DOT_flist)[index + 1],(self->
      GMSOBJ_txcustomstringlist_DOT_fcount - index) * sizeof(
      GMSOBJ_tstringitem));
  { register GMSOBJ_tstringitem *_W2= &(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index];
    _W2->fstring = STRUTILX_newstringm(s,&self->
      GMSOBJ_txcustomstringlist_DOT_fstrmemory);
    _W2->fobject = ValueCast(SYSTEM_tobject,apointer);

  }
  _P3inc0(self->GMSOBJ_txcustomstringlist_DOT_fcount);
}  /* insertitem */

Procedure GMSOBJ_txcustomstringlist_DOT_putobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index,
  SYSTEM_tobject aobject)
{
  (*self->GMSOBJ_txcustomstringlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fobject = aobject;
}  /* putobject */

Procedure GMSOBJ_txcustomstringlist_DOT_setcapacity(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer newcapacity)
{
  if (newcapacity < self->GMSOBJ_txcustomstringlist_DOT_fcount) 
    newcapacity = self->GMSOBJ_txcustomstringlist_DOT_fcount;
  SYSTEM_reallocmem(&PointerCast(SYSTEM_pointer,&self->
    GMSOBJ_txcustomstringlist_DOT_flist),newcapacity * sizeof(
    GMSOBJ_tstringitem));
  self->GMSOBJ_txcustomstringlist_DOT_fcapacity = newcapacity;
}  /* setcapacity */

Function(SYSTEM_int64 ) GMSOBJ_txcustomstringlist_DOT_memoryused(
  GMSOBJ_txcustomstringlist self)
{
  SYSTEM_int64 result;

  result = self->GMSOBJ_txcustomstringlist_DOT_fcapacity * sizeof(
    GMSOBJ_tstringitem) + self->
    GMSOBJ_txcustomstringlist_DOT_fstrmemory;
  return result;
}  /* memoryused */

Procedure GMSOBJ_txstringlist_DOT_insert(
  GMSOBJ_txstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *s)
{
  GMSOBJ_txcustomstringlist_DOT_insertitem(ValueCast(
    GMSOBJ_txcustomstringlist,self),index,s,NULL);
}  /* insert */

Procedure GMSOBJ_txstringlist_DOT_insertobject(
  GMSOBJ_txstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer)
{
  GMSOBJ_txcustomstringlist_DOT_insertitem(ValueCast(
    GMSOBJ_txcustomstringlist,self),index,s,apointer);
}  /* insertobject */

Constructor(GMSOBJ_txsortedstringlist ) 
  GMSOBJ_txsortedstringlist_DOT_create(
  GMSOBJ_txsortedstringlist self)
{
  ValueCast(GMSOBJ_txsortedstringlist,
    GMSOBJ_txcustomstringlist_DOT_create(ValueCast(
    GMSOBJ_txcustomstringlist,self)));
  self->GMSOBJ_txsortedstringlist_DOT_fupdatecount = 0;
  self->GMSOBJ_txsortedstringlist_DOT_fsorted = SYSTEM_true;
  return self;
}  /* create */

Function(SYSTEM_integer ) GMSOBJ_txsortedstringlist_DOT_addobject(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer)
{
  SYSTEM_integer result;

  if (self->GMSOBJ_txsortedstringlist_DOT_fupdatecount != 0 || 
    self->GMSOBJ_txcustomstringlist_DOT_fcount == 0) {
    result = self->GMSOBJ_txcustomstringlist_DOT_fcount + SYSTEM_ord(
      self->GMSOBJ_tquicksortclass_DOT_onebased);
    self->GMSOBJ_txsortedstringlist_DOT_fsorted = SYSTEM_false;
  } else 
    GMSOBJ_txsortedstringlist_DOT_find(self,s,&result);
  GMSOBJ_txcustomstringlist_DOT_insertitem(ValueCast(
    GMSOBJ_txcustomstringlist,self),result,s,apointer);
  return result;
}  /* addobject */

Function(SYSTEM_integer ) GMSOBJ_txsortedstringlist_DOT_add(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = GMSOBJ_txsortedstringlist_DOT_addobject(self,s,NULL);
  return result;
}  /* add */

Procedure GMSOBJ_txsortedstringlist_DOT_beginupdate(
  GMSOBJ_txsortedstringlist self)
{
  _P3inc0(self->GMSOBJ_txsortedstringlist_DOT_fupdatecount);
}  /* beginupdate */

Procedure GMSOBJ_txsortedstringlist_DOT_endupdate(
  GMSOBJ_txsortedstringlist self)
{
  _P3dec0(self->GMSOBJ_txsortedstringlist_DOT_fupdatecount);
  if (self->GMSOBJ_txsortedstringlist_DOT_fupdatecount == 0) 
    GMSOBJ_txsortedstringlist_DOT_setsorted(self,SYSTEM_true);
}  /* endupdate */

Procedure GMSOBJ_txsortedstringlist_DOT_setsorted(
  GMSOBJ_txsortedstringlist self,
  SYSTEM_boolean value)
{
  if (self->GMSOBJ_txsortedstringlist_DOT_fsorted != value) {
    if (value) 
      GMSOBJ_tquicksortclass_DOT_sortn(ValueCast(
        GMSOBJ_tquicksortclass,self),self->
        GMSOBJ_txcustomstringlist_DOT_fcount);
    self->GMSOBJ_txsortedstringlist_DOT_fsorted = value;
  } 
}  /* setsorted */

Function(SYSTEM_boolean ) GMSOBJ_txsortedstringlist_DOT_find(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer *index)
{
  SYSTEM_boolean result;
  SYSTEM_integer l, h, i, c;

  result = SYSTEM_false;
  l = 0;
  h = self->GMSOBJ_txcustomstringlist_DOT_fcount - 1;
  while (l <= h) {
    i = ValueCast(SYSTEM_uint32,l + h) >> 1;
    c = STRUTILX_pstrucmp(ValueCast(SYSTEM_P3_pshortstring,s),(*self->
      GMSOBJ_txcustomstringlist_DOT_flist)[i].fstring);
    if (c > 0) { 
      l = i + 1;
    } else {
      h = i - 1;
      if (c == 0) {
        result = SYSTEM_true;
        l = i;
      } 
    } 
  
}
  *index = l + SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased);
  return result;
}  /* find */

Function(SYSTEM_integer ) GMSOBJ_txsortedstringlist_DOT_indexof(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  if (!self->GMSOBJ_txsortedstringlist_DOT_fsorted) { 
    result = GMSOBJ_txcustomstringlist_DOT_indexof(ValueCast(
      GMSOBJ_txcustomstringlist,self),s);
  } else 
    if (!GMSOBJ_txsortedstringlist_DOT_find(self,s,&result)) 
      result =  -1;
  return result;
}  /* indexof */

static Procedure GMSOBJ_xhashstats(
  SYSTEM_P3_ppointerarray ph,
  SYSTEM_integer hashcount,
  SYSTEM_integer *amin,
  SYSTEM_integer *amax,
  SYSTEM_integer *aavg)
{
  SYSTEM_integer n;
  GMSOBJ_phashrecord p;
  SYSTEM_integer c;

  *amin = SYSTEM_maxint;
  *amax = 0;
  *aavg = 0;
  { register SYSTEM_int32 _stop = hashcount - 1;
    if ((n = 0) <=  _stop) do {
      c = 0;
      p = ValueCast(GMSOBJ_phashrecord,(*ph)[n]);
      while (p != NULL) {
        _P3inc0(c);
        p = p->pnext;
      
}
      if (c < *amin) 
        *amin = c;
      if (c > *amax) 
        *amax = c;
      _P3inc1(*aavg,c);
    
    } while (n++ !=  _stop);

  }
  *aavg = SYSTEM_round(ValueCast(SYSTEM_double,*aavg) /  (hashcount + 1));
}  /* xhashstats */

Constructor(GMSOBJ_txhashedstringlist ) 
  GMSOBJ_txhashedstringlist_DOT_create(
  GMSOBJ_txhashedstringlist self)
{
  ValueCast(GMSOBJ_txhashedstringlist,
    GMSOBJ_txcustomstringlist_DOT_create(ValueCast(
    GMSOBJ_txcustomstringlist,self)));
  self->GMSOBJ_txhashedstringlist_DOT_phash = NULL;
  self->GMSOBJ_txhashedstringlist_DOT_hashcount = 0;
  return self;
}  /* create */

Destructor(GMSOBJ_txhashedstringlist ) 
  GMSOBJ_txhashedstringlist_DOT_destroy(
  GMSOBJ_txhashedstringlist self)
{
  GMSOBJ_txhashedstringlist_DOT_clear(self);
  GMSOBJ_txcustomstringlist_DOT_destroy(ValueCast(
    GMSOBJ_txcustomstringlist,self));
  return self;
}  /* destroy */

Procedure GMSOBJ_txhashedstringlist_DOT_clear(
  GMSOBJ_txhashedstringlist self)
{
  GMSOBJ_txhashedstringlist_DOT_clearhashlist(self);
  GMSOBJ_txcustomstringlist_DOT_clear(ValueCast(
    GMSOBJ_txcustomstringlist,self));
}  /* clear */

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_indexof(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  GMSOBJ_phashrecord ph;

  if (self->GMSOBJ_txhashedstringlist_DOT_phash == NULL) 
    GMSOBJ_txhashedstringlist_DOT_sethashsize(self,self->
      GMSOBJ_txcustomstringlist_DOT_fcount);
  ph = ValueCast(GMSOBJ_phashrecord,(*self->
    GMSOBJ_txhashedstringlist_DOT_phash)[VirtMethodCall(self, 
    GMSOBJ_txhashedstringlist_DOT_hashvalue_T, 6, (self,s))]);
  while (ph != NULL) {
    if (VirtMethodCall(self, 
      GMSOBJ_txhashedstringlist_DOT_equaltoentry_T, 5, (self,s,ph->
      refnr))) 
      SYSTEM_break(BRK_1);
    ph = ph->pnext;
  
CNT_1:;
  }
BRK_1:;
  if (ph == NULL) { 
    result =  -1;
  } else 
    result = ph->refnr + SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased);
  return result;
}  /* indexof */

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_addobject(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer)
{
  SYSTEM_integer result;
  GMSOBJ_phashrecord ph;
  SYSTEM_integer hv;

  if (self->GMSOBJ_txhashedstringlist_DOT_phash == NULL || self->
    GMSOBJ_txcustomstringlist_DOT_fcount > 4 * self->
    GMSOBJ_txhashedstringlist_DOT_hashcount) 
    GMSOBJ_txhashedstringlist_DOT_sethashsize(self,self->
      GMSOBJ_txcustomstringlist_DOT_fcount);
  hv = VirtMethodCall(self, GMSOBJ_txhashedstringlist_DOT_hashvalue_T, 6, (
    self,s));
  ph = ValueCast(GMSOBJ_phashrecord,(*self->
    GMSOBJ_txhashedstringlist_DOT_phash)[hv]);
  while (ph != NULL) {
    if (VirtMethodCall(self, 
      GMSOBJ_txhashedstringlist_DOT_equaltoentry_T, 5, (self,s,ph->
      refnr))) 
      SYSTEM_break(BRK_2);
    ph = ph->pnext;
  
CNT_2:;
  }
BRK_2:;
  if (ph != NULL) { 
    result = ph->refnr + SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased);
  } else {
    result = self->GMSOBJ_txcustomstringlist_DOT_fcount + SYSTEM_ord(
      self->GMSOBJ_tquicksortclass_DOT_onebased);
    GMSOBJ_txcustomstringlist_DOT_insertitem(ValueCast(
      GMSOBJ_txcustomstringlist,self),result,s,apointer);
    _P3getmem(ph,sizeof(GMSOBJ_thashrecord));
    ph->pnext = ValueCast(GMSOBJ_phashrecord,(*self->
      GMSOBJ_txhashedstringlist_DOT_phash)[hv]);
    ph->refnr = result - SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased);
    (*self->GMSOBJ_txhashedstringlist_DOT_phash)[hv] = ph;
  } 
  return result;
}  /* addobject */

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_add(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = GMSOBJ_txhashedstringlist_DOT_addobject(self,s,NULL);
  return result;
}  /* add */

Function(SYSTEM_cardinal ) GMSOBJ_txhashedstringlist_DOT_hashvalue(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *name)
{
  SYSTEM_cardinal result;
  SYSTEM_integer i;

  result = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(name);
    if ((i = 1) <=  _stop) do {
      result = 211 * result + SYSTEM_ord(SYSTEM_upcase(name[i]));
    } while (i++ !=  _stop);

  }
  result = (result & 2147483647) % self->
    GMSOBJ_txhashedstringlist_DOT_hashcount;
  return result;
}  /* hashvalue */

Procedure GMSOBJ_txhashedstringlist_DOT_clearhashlist(
  GMSOBJ_txhashedstringlist self)
{
  SYSTEM_integer n;
  GMSOBJ_phashrecord p1, p2;

  if (self->GMSOBJ_txhashedstringlist_DOT_phash != NULL) {
    { register SYSTEM_int32 _stop = self->
        GMSOBJ_txhashedstringlist_DOT_hashcount - 1;
      if ((n = 0) <=  _stop) do {
        p1 = ValueCast(GMSOBJ_phashrecord,(*self->
          GMSOBJ_txhashedstringlist_DOT_phash)[n]);
        (*self->GMSOBJ_txhashedstringlist_DOT_phash)[n] = NULL;
        while (p1 != NULL) {
          p2 = p1->pnext;
          _P3freemem2(p1,sizeof(GMSOBJ_thashrecord));
          p1 = p2;
        
}
      
      } while (n++ !=  _stop);

    }
    _P3freemem(self->GMSOBJ_txhashedstringlist_DOT_phash);
    self->GMSOBJ_txhashedstringlist_DOT_phash = NULL;
    self->GMSOBJ_txhashedstringlist_DOT_hashcount = 0;
  } 
}  /* clearhashlist */

Procedure GMSOBJ_txhashedstringlist_DOT_sethashsize(
  GMSOBJ_txhashedstringlist self,
  SYSTEM_integer v)
{
  SYSTEM_integer n;
  SYSTEM_integer hv;
  GMSOBJ_phashrecord ph;

  if (v <= 9973) { 
    v = 9973;
  } else 
    v = 99991;
  if (v > 0 && v == self->GMSOBJ_txhashedstringlist_DOT_hashcount) 
    return;
  GMSOBJ_txhashedstringlist_DOT_clearhashlist(self);
  self->GMSOBJ_txhashedstringlist_DOT_hashcount = v;
  _P3getmem(self->GMSOBJ_txhashedstringlist_DOT_phash,sizeof(
    SYSTEM_pointer) * self->GMSOBJ_txhashedstringlist_DOT_hashcount);
  { register SYSTEM_int32 _stop = self->
      GMSOBJ_txhashedstringlist_DOT_hashcount - 1;
    if ((n = 0) <=  _stop) do {
      (*self->GMSOBJ_txhashedstringlist_DOT_phash)[n] = NULL;
    } while (n++ !=  _stop);

  }
  { register SYSTEM_int32 _stop = self->
      GMSOBJ_txcustomstringlist_DOT_fcount - 1 + SYSTEM_ord(self->
      GMSOBJ_tquicksortclass_DOT_onebased);
    if ((n = SYSTEM_ord(self->GMSOBJ_tquicksortclass_DOT_onebased)) <=  _stop) do {
      {
        SYSTEM_shortstring _t1;

        hv = VirtMethodCall(self, 
          GMSOBJ_txhashedstringlist_DOT_hashvalue_T, 6, (self,
          GMSOBJ_txcustomstringlist_DOT_getname(_t1,255,ValueCast(
          GMSOBJ_txcustomstringlist,self),n)));
      }
      _P3getmem(ph,sizeof(GMSOBJ_thashrecord));
      ph->pnext = ValueCast(GMSOBJ_phashrecord,(*self->
        GMSOBJ_txhashedstringlist_DOT_phash)[hv]);
      ph->refnr = n - SYSTEM_ord(self->
        GMSOBJ_tquicksortclass_DOT_onebased);
      (*self->GMSOBJ_txhashedstringlist_DOT_phash)[hv] = ph;
    
    } while (n++ !=  _stop);

  }
}  /* sethashsize */

Procedure GMSOBJ_txhashedstringlist_DOT_sort(
  GMSOBJ_txhashedstringlist self)
{
  GMSOBJ_txhashedstringlist_DOT_clearhashlist(self);
  GMSOBJ_tquicksortclass_DOT_sortn(ValueCast(GMSOBJ_tquicksortclass,
    self),self->GMSOBJ_txcustomstringlist_DOT_fcount);
}  /* sort */

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_compare(
  GMSOBJ_txhashedstringlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_integer result;

  result = STRUTILX_pstrucmp((*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index1 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring,(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index2 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring);
  return result;
}  /* compare */

Procedure GMSOBJ_txhashedstringlist_DOT_hashstats(
  GMSOBJ_txhashedstringlist self,
  SYSTEM_integer *amin,
  SYSTEM_integer *amax,
  SYSTEM_integer *aavg)
{
  GMSOBJ_xhashstats(self->GMSOBJ_txhashedstringlist_DOT_phash,self->
    GMSOBJ_txhashedstringlist_DOT_hashcount,amin,amax,aavg);
}  /* hashstats */

Function(SYSTEM_boolean ) GMSOBJ_txhashedstringlist_DOT_equaltoentry(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer en)
{
  SYSTEM_boolean result;

  result = STRUTILX_pstruequal(ValueCast(SYSTEM_P3_pshortstring,s),(*
    self->GMSOBJ_txcustomstringlist_DOT_flist)[en].fstring);
  return result;
}  /* equaltoentry */

Function(SYSTEM_integer ) GMSOBJ_txstrpool_DOT_compare(
  GMSOBJ_txstrpool self,
  SYSTEM_integer index1,
  SYSTEM_integer index2)
{
  SYSTEM_integer result;

  result = STRUTILX_pstrcmp((*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index1 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring,(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index2 - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fstring);
  return result;
}  /* compare */

Function(SYSTEM_boolean ) GMSOBJ_txstrpool_DOT_equaltoentry(
  GMSOBJ_txstrpool self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer en)
{
  SYSTEM_boolean result;

  result = STRUTILX_pstrequal(ValueCast(SYSTEM_P3_pshortstring,s),(*
    self->GMSOBJ_txcustomstringlist_DOT_flist)[en].fstring);
  return result;
}  /* equaltoentry */

Function(SYSTEM_integer ) GMSOBJ_txstrstrlist_DOT_addobject(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2)
{
  SYSTEM_integer result;

  result = GMSOBJ_txsortedstringlist_DOT_addobject(ValueCast(
    GMSOBJ_txsortedstringlist,self),s1,STRUTILX_newstringm(s2,&self->
    GMSOBJ_txcustomstringlist_DOT_fstrmemory));
  return result;
}  /* addobject */

Procedure GMSOBJ_txstrstrlist_DOT_freeobject(
  GMSOBJ_txstrstrlist self,
  SYSTEM_integer index)
{
  STRUTILX_disposestringm(ValueCast(SYSTEM_P3_pshortstring,(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fobject),&self->
    GMSOBJ_txcustomstringlist_DOT_fstrmemory);
}  /* freeobject */

Function(SYSTEM_boolean ) GMSOBJ_txstrstrlist_DOT_getasboolean(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s1;
  SYSTEM_boolean result;
  SYSTEM_shortstring s2;

  _P3strcpy(s1,255,_ftmp1);
  GMSOBJ_txstrstrlist_DOT_getasstring(s2,255,self,s1);
  if (_P3strcmpE(s2,_P3str1("\000"))) { 
    result = SYSTEM_false;
  } else 
    result = _P3SET_in_1(s2[1],_P3char('Y'),_P3SET_in_1(s2[1],_P3char('y'),
      _P3SET_in_1(s2[1],_P3char('1'),_P3SET_in_1(s2[1],_P3char('t'),
      _P3SET_equal(s2[1],_P3char('T'))))));
  return result;
}  /* getasboolean */

Function(SYSTEM_double ) GMSOBJ_txstrstrlist_DOT_getasdouble(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s1;
  SYSTEM_double result;
  SYSTEM_integer n;
  SYSTEM_shortstring s2;

  _P3strcpy(s1,255,_ftmp1);
  GMSOBJ_txstrstrlist_DOT_getasstring(s2,255,self,s1);
  if (_P3strcmpE(s2,_P3str1("\000"))) { 
    result = 0.0;
  } else {
    _P3val_d(s2,result,&n);
    if (n != 0) 
      result = 0.0;
  } 
  return result;
}  /* getasdouble */

Function(SYSTEM_integer ) GMSOBJ_txstrstrlist_DOT_getasinteger(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s1;
  SYSTEM_integer result;
  SYSTEM_integer n;
  SYSTEM_shortstring s2;

  _P3strcpy(s1,255,_ftmp1);
  GMSOBJ_txstrstrlist_DOT_getasstring(s2,255,self,s1);
  if (_P3strcmpE(s2,_P3str1("\000"))) { 
    result = 0;
  } else {
    _P3val_i(s2,result,&n);
    if (n != 0) 
      result = 0;
  } 
  return result;
}  /* getasinteger */

Function(SYSTEM_ansichar *) GMSOBJ_txstrstrlist_DOT_getasstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring s1;
  SYSTEM_integer n;

  _P3strcpy(s1,255,_ftmp1);
  n = GMSOBJ_txsortedstringlist_DOT_indexof(ValueCast(
    GMSOBJ_txsortedstringlist,self),s1);
  if (n < 0) { 
    _P3strclr(result);
  } else {
    GMSOBJ_txstrstrlist_DOT_getobject(result,_len_ret,self,n);
    if (_P3stccmpE(result,GMSOBJ_non_empty)) 
      _P3strclr(result);
  } 
  return result;
}  /* getasstring */

Function(SYSTEM_ansichar *) GMSOBJ_txstrstrlist_DOT_getobject(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrstrlist self,
  SYSTEM_integer index)
{
  STRUTILX_getstring(result,_len_ret,ValueCast(SYSTEM_P3_pshortstring,(*
    self->GMSOBJ_txcustomstringlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fobject));
  return result;
}  /* getobject */

Procedure GMSOBJ_txstrstrlist_DOT_putobject(
  GMSOBJ_txstrstrlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *_ftmp1)
{
  SYSTEM_shortstring astr;

  _P3strcpy(astr,255,_ftmp1);
  VirtMethodCall(ValueCast(GMSOBJ_txcustomstringlist,self), 
    GMSOBJ_txcustomstringlist_DOT_freeobject_T, 4, (ValueCast(
    GMSOBJ_txcustomstringlist,self),index));
  PointerCast(SYSTEM_P3_pshortstring,&(*self->
    GMSOBJ_txcustomstringlist_DOT_flist)[index - SYSTEM_ord(self->
    GMSOBJ_tquicksortclass_DOT_onebased)].fobject) = 
    STRUTILX_newstringm(astr,&self->
    GMSOBJ_txcustomstringlist_DOT_fstrmemory);
}  /* putobject */

Procedure GMSOBJ_txstrstrlist_DOT_setasstring(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  const SYSTEM_ansichar *_ftmp2)
{
  SYSTEM_shortstring s1, s2;
  SYSTEM_integer n;

  _P3strcpy(s1,255,_ftmp1);
  _P3strcpy(s2,255,_ftmp2);
  n = GMSOBJ_txsortedstringlist_DOT_indexof(ValueCast(
    GMSOBJ_txsortedstringlist,self),s1);
  if (_P3strcmpE(s2,_P3str1("\000"))) { 
    if (n >= 0) 
      GMSOBJ_txcustomstringlist_DOT_delete(ValueCast(
        GMSOBJ_txcustomstringlist,self),n);
  } else 
    if (n >= 0) { 
      GMSOBJ_txstrstrlist_DOT_putobject(self,n,s2);
    } else 
      GMSOBJ_txstrstrlist_DOT_addobject(self,s1,s2);
}  /* setasstring */

Procedure GMSOBJ_txstrstrlist_DOT_setasinteger(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_integer v)
{
  SYSTEM_shortstring s1;

  _P3strcpy(s1,255,_ftmp1);
  if (v == 0) { 
    GMSOBJ_txstrstrlist_DOT_setasstring(self,s1,_P3str1("\000"));
  } else 
    {
      SYSTEM_shortstring _t1;

      GMSOBJ_txstrstrlist_DOT_setasstring(self,s1,STRUTILX_inttostr(
        _t1,255,v));
    }
}  /* setasinteger */

Procedure GMSOBJ_txstrstrlist_DOT_setasdouble(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_double v)
{
  SYSTEM_shortstring s1;
  SYSTEM_shortstring s2;

  _P3strcpy(s1,255,_ftmp1);
  if (v == 0.0) { 
    GMSOBJ_txstrstrlist_DOT_setasstring(self,s1,_P3str1("\000"));
  } else {
    _P3str_d0(v,s2,255);
    GMSOBJ_txstrstrlist_DOT_setasstring(self,s1,s2);
  } 
}  /* setasdouble */

Procedure GMSOBJ_txstrstrlist_DOT_setasboolean(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_boolean v)
{
  SYSTEM_shortstring s1;

  _P3strcpy(s1,255,_ftmp1);
  if (!v) { 
    GMSOBJ_txstrstrlist_DOT_setasstring(self,s1,_P3str1("\000"));
  } else 
    GMSOBJ_txstrstrlist_DOT_setasstring(self,s1,_P3str1("\001T"));
}  /* setasboolean */

Constructor(GMSOBJ_tbooleanbitarray ) 
  GMSOBJ_tbooleanbitarray_DOT_create(
  GMSOBJ_tbooleanbitarray self)
{
  ValueCast(GMSOBJ_tbooleanbitarray,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSOBJ_tbooleanbitarray_DOT_pdata = NULL;
  self->GMSOBJ_tbooleanbitarray_DOT_fallocated = 0;
  self->GMSOBJ_tbooleanbitarray_DOT_fhighindex =  -1;
  return self;
}  /* create */

Destructor(GMSOBJ_tbooleanbitarray ) 
  GMSOBJ_tbooleanbitarray_DOT_destroy(
  GMSOBJ_tbooleanbitarray self)
{
  if (self->GMSOBJ_tbooleanbitarray_DOT_fallocated > 0) 
    _P3freemem2(self->GMSOBJ_tbooleanbitarray_DOT_pdata,self->
      GMSOBJ_tbooleanbitarray_DOT_fallocated);
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSOBJ_tbooleanbitarray_DOT_clear(
  GMSOBJ_tbooleanbitarray self)
{
  if (self->GMSOBJ_tbooleanbitarray_DOT_fallocated > 0) 
    _P3freemem2(self->GMSOBJ_tbooleanbitarray_DOT_pdata,self->
      GMSOBJ_tbooleanbitarray_DOT_fallocated);
  self->GMSOBJ_tbooleanbitarray_DOT_pdata = NULL;
  self->GMSOBJ_tbooleanbitarray_DOT_fallocated = 0;
  self->GMSOBJ_tbooleanbitarray_DOT_fhighindex =  -1;
}  /* clear */

Procedure GMSOBJ_tbooleanbitarray_DOT_getbitmask(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer v,
  SYSTEM_integer *n,
  SYSTEM_byte *m)
{
  *n = ValueCast(SYSTEM_uint32,v) >> 3;
  *m = 1 << (v & 7);
}  /* getbitmask */

Function(SYSTEM_boolean ) GMSOBJ_tbooleanbitarray_DOT_getbit(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer n)
{
  SYSTEM_boolean result;
  SYSTEM_integer p;
  SYSTEM_byte m;

  if (n < 0 || n > self->GMSOBJ_tbooleanbitarray_DOT_fhighindex) { 
    result = SYSTEM_false;
  } else {
    GMSOBJ_tbooleanbitarray_DOT_getbitmask(self,n,&p,&m);
    result = (ValueCast(SYSTEM_int32,(*ValueCast(GMSGEN_pbytedataarray,
      self->GMSOBJ_tbooleanbitarray_DOT_pdata))[p]) & m) != 0;
  } 
  return result;
}  /* getbit */

Procedure GMSOBJ_tbooleanbitarray_DOT_setbit(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer n,
  SYSTEM_boolean v)
{
  SYSTEM_integer p;
  SYSTEM_byte m;

  if (n >= 0) {
    if (n > self->GMSOBJ_tbooleanbitarray_DOT_fhighindex) {
      if (!v) 
        return;
      GMSOBJ_tbooleanbitarray_DOT_sethighindex(self,n);
    } 
    GMSOBJ_tbooleanbitarray_DOT_getbitmask(self,n,&p,&m);
    if (v) { 
      (*self->GMSOBJ_tbooleanbitarray_DOT_pdata)[p] = ValueCast(
        SYSTEM_int32,(*self->GMSOBJ_tbooleanbitarray_DOT_pdata)[p]) | 
        m;
    } else 
      (*self->GMSOBJ_tbooleanbitarray_DOT_pdata)[p] = ValueCast(
        SYSTEM_int32,(*self->GMSOBJ_tbooleanbitarray_DOT_pdata)[p]) & ~
        m;
  } 
}  /* setbit */

Procedure GMSOBJ_tbooleanbitarray_DOT_sethighindex(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer v)
{
  SYSTEM_pointer newmem;
  SYSTEM_integer newmemsize;
  SYSTEM_integer delta;

  if (v > self->GMSOBJ_tbooleanbitarray_DOT_fhighindex) {
    newmemsize = (v + 8) /  8;
    if (newmemsize > self->GMSOBJ_tbooleanbitarray_DOT_fallocated) {
      delta = 0;
      do {
        if (self->GMSOBJ_tbooleanbitarray_DOT_fallocated == 0) { 
          _P3inc1(delta,256);
        } else 
          if (self->GMSOBJ_tbooleanbitarray_DOT_fallocated < 8192) { 
            _P3inc1(delta,self->GMSOBJ_tbooleanbitarray_DOT_fallocated);
          } else 
            _P3inc1(delta,self->GMSOBJ_tbooleanbitarray_DOT_fallocated /  4);
      } while (!(newmemsize < self->
        GMSOBJ_tbooleanbitarray_DOT_fallocated + delta));
      newmemsize = self->GMSOBJ_tbooleanbitarray_DOT_fallocated + 
        delta;
      _P3getmem(newmem,newmemsize);
      SYSTEM_P3_fillchar(newmem,newmemsize,0);
      if (self->GMSOBJ_tbooleanbitarray_DOT_fallocated != 0) {
        SYSTEM_move(*self->GMSOBJ_tbooleanbitarray_DOT_pdata,newmem,
          self->GMSOBJ_tbooleanbitarray_DOT_fallocated);
        _P3freemem2(self->GMSOBJ_tbooleanbitarray_DOT_pdata,self->
          GMSOBJ_tbooleanbitarray_DOT_fallocated);
      } 
      self->GMSOBJ_tbooleanbitarray_DOT_pdata = ValueCast(
        GMSGEN_pbytedataarray,newmem);
      self->GMSOBJ_tbooleanbitarray_DOT_fallocated = newmemsize;
    } 
    self->GMSOBJ_tbooleanbitarray_DOT_fhighindex = v;
  } 
}  /* sethighindex */

Procedure GMSOBJ_tbooleanbitarray_DOT_iterate(
  GMSOBJ_tbooleanbitarray self,
  GMSOBJ_tbooleanbitarrayiterfunction func)
{
  SYSTEM_integer n;
  SYSTEM_integer p;
  SYSTEM_integer pp;
  SYSTEM_byte m;

  if (self->GMSOBJ_tbooleanbitarray_DOT_fhighindex >= 0) {
    GMSOBJ_tbooleanbitarray_DOT_getbitmask(self,self->
      GMSOBJ_tbooleanbitarray_DOT_fhighindex,&pp,&m);
    { register SYSTEM_int32 _stop = pp;
      if ((p = 0) <=  _stop) do {
        m = (*self->GMSOBJ_tbooleanbitarray_DOT_pdata)[p];
        if (m == 0) 
          SYSTEM_continue(CNT_3);
        n = p * 8;
        do {
          if ((ValueCast(SYSTEM_int32,m) & 1) != 0) 
            if (!(*func)(n)) 
              return;
          n = n + 1;
          m = m >> 1;
        } while (!(m == 0));
      
CNT_3:;
      } while (p++ !=  _stop);
BRK_3:;

    }
  } 
  (*func)( -1);
}  /* iterate */

Procedure GMSOBJ_tbooleanbitarray_DOT_iteratedown(
  GMSOBJ_tbooleanbitarray self,
  GMSOBJ_tbooleanbitarrayiterfunction func)
{
  SYSTEM_integer n;
  SYSTEM_integer p;
  SYSTEM_integer pp;
  SYSTEM_byte m;

  if (self->GMSOBJ_tbooleanbitarray_DOT_fhighindex >= 0) {
    GMSOBJ_tbooleanbitarray_DOT_getbitmask(self,self->
      GMSOBJ_tbooleanbitarray_DOT_fhighindex,&pp,&m);
    for (p = pp;p >= (SYSTEM_int32)0;--p) {
      m = (*self->GMSOBJ_tbooleanbitarray_DOT_pdata)[p];
      if (m == 0) 
        SYSTEM_continue(CNT_4);
      n = (p + 1) * 8 - 1;
      do {
        if ((ValueCast(SYSTEM_int32,m) & 128) != 0) 
          if (!(*func)(n)) 
            return;
        n = n - 1;
        m = m << 1;
      } while (!(m == 0));
    
    CNT_4:;
}
BRK_4:;
  } 
  (*func)( -1);
}  /* iteratedown */

Function(SYSTEM_int64 ) GMSOBJ_tbooleanbitarray_DOT_memoryused(
  GMSOBJ_tbooleanbitarray self)
{
  SYSTEM_int64 result;

  result = self->GMSOBJ_tbooleanbitarray_DOT_fallocated;
  return result;
}  /* memoryused */

Function(SYSTEM_integer ) GMSOBJ_txpcharlist_DOT_add(
  GMSOBJ_txpcharlist self,
  SYSTEM_P3_pansichar p,
  SYSTEM_integer l)
{
  SYSTEM_integer result;
  SYSTEM_P3_pansichar ps, pw;
  SYSTEM_integer size;

  size = sizeof(SYSTEM_longint) + l;
  _P3getmem(ps,size);
  self->GMSOBJ_txpcharlist_DOT_fstrmemory = self->
    GMSOBJ_txpcharlist_DOT_fstrmemory + size;
  pw = ps;
  *ValueCast(SYSTEM_P3_pinteger,pw) = l;
  _P3inc1(pw,sizeof(SYSTEM_longint));
  SYSTEM_move(p,pw,l);
  result = GMSOBJ_txlist_DOT_add(self->GMSOBJ_txpcharlist_DOT_flist,ps);
  return result;
}  /* add */

Procedure GMSOBJ_txpcharlist_DOT_clear(
  GMSOBJ_txpcharlist self)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->GMSOBJ_txpcharlist_DOT_flist->
      GMSOBJ_txlist_DOT_fcount - 1 + SYSTEM_ord(
      GMSOBJ_txpcharlist_DOT_getonebased(self));
    if ((n = SYSTEM_ord(GMSOBJ_txpcharlist_DOT_getonebased(self))) <=  _stop) do {
      _P3freemem0(GMSOBJ_txlist_DOT_get(self->
        GMSOBJ_txpcharlist_DOT_flist,n));
    } while (n++ !=  _stop);

  }
  GMSOBJ_txlist_DOT_clear(self->GMSOBJ_txpcharlist_DOT_flist);
  self->GMSOBJ_txpcharlist_DOT_fstrmemory = 0;
}  /* clear */

Function(SYSTEM_integer ) GMSOBJ_txpcharlist_DOT_count(
  GMSOBJ_txpcharlist self)
{
  SYSTEM_integer result;

  result = self->GMSOBJ_txpcharlist_DOT_flist->
    GMSOBJ_txlist_DOT_fcount;
  return result;
}  /* count */

Constructor(GMSOBJ_txpcharlist ) GMSOBJ_txpcharlist_DOT_create(
  GMSOBJ_txpcharlist self)
{
  ValueCast(GMSOBJ_txpcharlist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->GMSOBJ_txpcharlist_DOT_flist = ValueCast(GMSOBJ_txlist,
    SYSTEM_tobject_DOT_create(ValueCast(SYSTEM_tobject,_P3alloc_object(&
    GMSOBJ_txlist_CD))));
  self->GMSOBJ_txpcharlist_DOT_fstrmemory = 0;
  return self;
}  /* create */

Destructor(GMSOBJ_txpcharlist ) GMSOBJ_txpcharlist_DOT_destroy(
  GMSOBJ_txpcharlist self)
{
  GMSOBJ_txpcharlist_DOT_clear(self);
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    GMSOBJ_txpcharlist_DOT_flist));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure GMSOBJ_txpcharlist_DOT_getitem(
  GMSOBJ_txpcharlist self,
  SYSTEM_integer index,
  SYSTEM_P3_pansichar *p,
  SYSTEM_integer *l)
{
  *p = ValueCast(SYSTEM_P3_pansichar,GMSOBJ_txlist_DOT_get(self->
    GMSOBJ_txpcharlist_DOT_flist,index));
  *l = *ValueCast(SYSTEM_P3_pinteger,*p);
  _P3inc1(*p,sizeof(SYSTEM_longint));
}  /* getitem */

Function(SYSTEM_boolean ) GMSOBJ_txpcharlist_DOT_getonebased(
  GMSOBJ_txpcharlist self)
{
  SYSTEM_boolean result;

  result = self->GMSOBJ_txpcharlist_DOT_flist->
    GMSOBJ_tquicksortclass_DOT_onebased;
  return result;
}  /* getonebased */

Function(SYSTEM_int64 ) GMSOBJ_txpcharlist_DOT_memoryused(
  GMSOBJ_txpcharlist self)
{
  SYSTEM_int64 result;

  result = GMSOBJ_txlist_DOT_memoryused(self->
    GMSOBJ_txpcharlist_DOT_flist) + self->
    GMSOBJ_txpcharlist_DOT_fstrmemory;
  return result;
}  /* memoryused */

Procedure GMSOBJ_txpcharlist_DOT_setonebased(
  GMSOBJ_txpcharlist self,
  SYSTEM_boolean v)
{
  self->GMSOBJ_txpcharlist_DOT_flist->
    GMSOBJ_tquicksortclass_DOT_onebased = v;
}  /* setonebased */

/* unit gmsobj */
void _Init_Module_gmsobj(void)
{
} /* _Init_Module_gmsobj */

void _Final_Module_gmsobj(void)
{
} /* _Final_Module_gmsobj */

