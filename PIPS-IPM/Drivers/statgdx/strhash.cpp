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
#include "gmsdata.h"
#include "strhash.h"


void * const STRHASH_txstrhashlist_VT[] = {(void*)&
  STRHASH_txstrhashlist_DOT_destroy, (void*)&
  STRHASH_txstrhashlist_DOT_hash, (void*)&
  STRHASH_txstrhashlist_DOT_entryequal, (void*)&
  STRHASH_txstrhashlist_DOT_compare, (void*)&
  STRHASH_txstrhashlist_DOT_freeitem};

/* Class descriptor for 'txstrhashlist' */
const SYSTEM_classdescriptor_t STRHASH_txstrhashlist_CD = {
  _P3str1("\015txstrhashlist"), 
  &SYSTEM_tobject_CD, NULL, 0, 
  sizeof(STRHASH_txstrhashlist_OD), STRHASH_txstrhashlist_VT, NULL};


void * const STRHASH_txcsstrhashlist_VT[] = {(void*)&
  STRHASH_txstrhashlist_DOT_destroy, (void*)&
  STRHASH_txcsstrhashlist_DOT_hash, (void*)&
  STRHASH_txcsstrhashlist_DOT_entryequal, (void*)&
  STRHASH_txstrhashlist_DOT_compare, (void*)&
  STRHASH_txstrhashlist_DOT_freeitem};

/* Class descriptor for 'txcsstrhashlist' */
const SYSTEM_classdescriptor_t STRHASH_txcsstrhashlist_CD = {
  _P3str1("\017txcsstrhashlist"), 
  &STRHASH_txstrhashlist_CD, NULL, 0, 
  sizeof(STRHASH_txcsstrhashlist_OD), STRHASH_txcsstrhashlist_VT, NULL};


static Function(SYSTEM_P3_pshortstring ) STRHASH_newstringx(
  const SYSTEM_ansichar *s)
{
  SYSTEM_P3_pshortstring result;

  _P3getmem(result,ValueCast(SYSTEM_int32,SYSTEM_length(s)) + 1);
  _P3strcpy(*result,255,s);
  return result;
}  /* newstringx */

static Procedure STRHASH_disposestringx(
  SYSTEM_P3_pshortstring ps)
{
  _P3freemem2(ps,ValueCast(SYSTEM_int32,SYSTEM_length(*ps)) + 1);
}  /* disposestringx */

static Procedure STRHASH_strassignx(
  SYSTEM_P3_pshortstring *p,
  const SYSTEM_ansichar *s)
{
  STRHASH_disposestringx(*p);
  *p = STRHASH_newstringx(s);
}  /* strassignx */
typedef struct STRHASH_thashbucket_S *STRHASH_phashbucket;
typedef struct STRHASH_thashbucket_S {
  SYSTEM_P3_pshortstring strp;
  STRHASH_phashbucket nxtbuck;
  SYSTEM_integer strnr;
  SYSTEM_tobject obj;
} STRHASH_thashbucket;


Constructor(STRHASH_txstrhashlist ) STRHASH_txstrhashlist_DOT_create(
  STRHASH_txstrhashlist self)
{
  ValueCast(STRHASH_txstrhashlist,SYSTEM_tobject_DOT_create(ValueCast(
    SYSTEM_tobject,self)));
  self->STRHASH_txstrhashlist_DOT_fcount = 0;
  self->STRHASH_txstrhashlist_DOT_buckets = ValueCast(
    GMSDATA_tgrowarrayfxd,GMSDATA_tgrowarrayfxd_DOT_create(ValueCast(
    GMSDATA_tgrowarrayfxd,_P3alloc_object(&GMSDATA_tgrowarrayfxd_CD)),sizeof(
    STRHASH_thashbucket)));
  STRHASH_txstrhashlist_DOT_clearhashtable(self);
  self->STRHASH_txstrhashlist_DOT_onebased = SYSTEM_false;
  self->STRHASH_txstrhashlist_DOT_sortmap = NULL;
  self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_false;
  return self;
}  /* create */

Destructor(STRHASH_txstrhashlist ) STRHASH_txstrhashlist_DOT_destroy(
  STRHASH_txstrhashlist self)
{
  STRHASH_txstrhashlist_DOT_clear(self);
  _P3freemem(self->STRHASH_txstrhashlist_DOT_phashtable);
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    STRHASH_txstrhashlist_DOT_sortmap));
  SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,self->
    STRHASH_txstrhashlist_DOT_buckets));
  SYSTEM_tobject_DOT_destroy(ValueCast(SYSTEM_tobject,self));
  return self;
}  /* destroy */

Procedure STRHASH_txstrhashlist_DOT_clear(
  STRHASH_txstrhashlist self)
{
  SYSTEM_integer n;

  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_fcount - 1 + SYSTEM_ord(self->
      STRHASH_txstrhashlist_DOT_onebased);
    if ((n = SYSTEM_ord(self->STRHASH_txstrhashlist_DOT_onebased)) <=  _stop) do {
      VirtMethodCall(self, STRHASH_txstrhashlist_DOT_freeitem_T, 4, (
        self,n));
    } while (n++ !=  _stop);

  }
  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      STRHASH_disposestringx((ValueCast(STRHASH_phashbucket,
        GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
        STRHASH_txstrhashlist_DOT_buckets,n)))->strp);
    } while (n++ !=  _stop);

  }
  GMSDATA_tgrowarrayfxd_DOT_clear(self->
    STRHASH_txstrhashlist_DOT_buckets);
  self->STRHASH_txstrhashlist_DOT_fcount = 0;
  STRHASH_txstrhashlist_DOT_clearhashtable(self);
  SYSUTILS_P3_freeandnil(&self->STRHASH_txstrhashlist_DOT_sortmap);
  self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_false;
}  /* clear */

Procedure STRHASH_txstrhashlist_DOT_hashtablereset(
  STRHASH_txstrhashlist self,
  SYSTEM_integer acnt)
{
  cnstdef {hashsize_1 = 97};
  cnstdef {next_1 = 150};
  cnstdef {hashsize_2 = 9973};
  cnstdef {next_2 = 10000};
  cnstdef {hashsize_3 = 99991};
  cnstdef {next_3 = 100000};
  cnstdef {hashsize_4 = 999979};
  cnstdef {next_4 = SYSTEM_maxint};
  SYSTEM_integer n;

  if (acnt >= next_3) {
    self->STRHASH_txstrhashlist_DOT_hashsize = hashsize_4;
    self->STRHASH_txstrhashlist_DOT_rehashcnt = next_4;
  } else 
    if (acnt >= 10000) {
      self->STRHASH_txstrhashlist_DOT_hashsize = hashsize_3;
      self->STRHASH_txstrhashlist_DOT_rehashcnt = next_3;
    } else 
      if (acnt >= 150) {
        self->STRHASH_txstrhashlist_DOT_hashsize = hashsize_2;
        self->STRHASH_txstrhashlist_DOT_rehashcnt = next_2;
      } else {
        self->STRHASH_txstrhashlist_DOT_hashsize = hashsize_1;
        self->STRHASH_txstrhashlist_DOT_rehashcnt = next_1;
      } 
  _P3getmem(self->STRHASH_txstrhashlist_DOT_phashtable,self->
    STRHASH_txstrhashlist_DOT_hashsize * sizeof(SYSTEM_pointer));
  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_hashsize - 1;
    if ((n = 0) <=  _stop) do {
      (*self->STRHASH_txstrhashlist_DOT_phashtable)[n] = NULL;
    } while (n++ !=  _stop);

  }
}  /* hashtablereset */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_hash(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer i;

  result = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      result = 211 * result + SYSTEM_ord(SYSTEM_upcase(s[i]));
    } while (i++ !=  _stop);

  }
  result = (result & 2147483647) % self->
    STRHASH_txstrhashlist_DOT_hashsize;
  return result;
}  /* hash */

Function(SYSTEM_boolean ) STRHASH_txstrhashlist_DOT_entryequal(
  STRHASH_txstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2)
{
  SYSTEM_boolean result;

  result = STRUTILX_pstruequal(ps1,ps2);
  return result;
}  /* entryequal */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_compare(
  STRHASH_txstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2)
{
  SYSTEM_integer result;

  result = STRUTILX_pstrucmp(ps1,ps2);
  return result;
}  /* compare */

Procedure STRHASH_txstrhashlist_DOT_hashall(
  STRHASH_txstrhashlist self)
{
  SYSTEM_integer hv;
  STRHASH_phashbucket pbuck;
  SYSTEM_integer n;

  _P3freemem(self->STRHASH_txstrhashlist_DOT_phashtable);
  STRHASH_txstrhashlist_DOT_hashtablereset(self,self->
    STRHASH_txstrhashlist_DOT_fcount);
  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      pbuck = ValueCast(STRHASH_phashbucket,
        GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
        STRHASH_txstrhashlist_DOT_buckets,n));
      { register STRHASH_thashbucket *_W2=pbuck;
        hv = VirtMethodCall(self, STRHASH_txstrhashlist_DOT_hash_T, 1, (
          self,*_W2->strp));
        _W2->nxtbuck = ValueCast(STRHASH_phashbucket,(*self->
          STRHASH_txstrhashlist_DOT_phashtable)[hv]);
        (*self->STRHASH_txstrhashlist_DOT_phashtable)[hv] = pbuck;

      }
    
    } while (n++ !=  _stop);

  }
}  /* hashall */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_addobject(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_tobject aobj)
{
  SYSTEM_integer result;
  SYSTEM_integer hv;
  STRHASH_phashbucket pbuck;

  if (self->STRHASH_txstrhashlist_DOT_fcount >= self->
    STRHASH_txstrhashlist_DOT_rehashcnt) 
    STRHASH_txstrhashlist_DOT_hashall(self);
  hv = VirtMethodCall(self, STRHASH_txstrhashlist_DOT_hash_T, 1, (self,
    s));
  pbuck = ValueCast(STRHASH_phashbucket,(*self->
    STRHASH_txstrhashlist_DOT_phashtable)[hv]);
  while (pbuck != NULL) {

    if (!VirtMethodCall(self, STRHASH_txstrhashlist_DOT_entryequal_T, 2, (
      self,pbuck->strp,ValueCast(SYSTEM_P3_pshortstring,s)))) { 
      pbuck = pbuck->nxtbuck;
    } else {
      result = pbuck->strnr + SYSTEM_ord(self->
        STRHASH_txstrhashlist_DOT_onebased);
      return result;
    } 
}
  pbuck = ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarray_DOT_reservemem(ValueCast(GMSDATA_tgrowarray,
    self->STRHASH_txstrhashlist_DOT_buckets),sizeof(
    STRHASH_thashbucket)));
  { register STRHASH_thashbucket *_W2=pbuck;
    _W2->nxtbuck = ValueCast(STRHASH_phashbucket,(*self->
      STRHASH_txstrhashlist_DOT_phashtable)[hv]);
    (*self->STRHASH_txstrhashlist_DOT_phashtable)[hv] = pbuck;
    _W2->strnr = self->STRHASH_txstrhashlist_DOT_fcount;
    result = self->STRHASH_txstrhashlist_DOT_fcount + SYSTEM_ord(self->
      STRHASH_txstrhashlist_DOT_onebased);
    if (self->STRHASH_txstrhashlist_DOT_sortmap != NULL) {
      GMSDATA_txintlist_DOT_setitems(self->
        STRHASH_txstrhashlist_DOT_sortmap,self->
        STRHASH_txstrhashlist_DOT_fcount,self->
        STRHASH_txstrhashlist_DOT_fcount);
      self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_false;
    } 
    self->STRHASH_txstrhashlist_DOT_fcount = self->
      STRHASH_txstrhashlist_DOT_fcount + 1;
    _W2->strp = STRHASH_newstringx(s);
    _W2->obj = aobj;

  }
  return result;
}  /* addobject */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_storeobject(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_tobject aobj)
{
  SYSTEM_integer result;
  STRHASH_phashbucket pbuck;

  if (self->STRHASH_txstrhashlist_DOT_phashtable != NULL) 
    STRHASH_txstrhashlist_DOT_clearhashtable(self);
  pbuck = ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarray_DOT_reservemem(ValueCast(GMSDATA_tgrowarray,
    self->STRHASH_txstrhashlist_DOT_buckets),sizeof(
    STRHASH_thashbucket)));
  { register STRHASH_thashbucket *_W2=pbuck;
    _W2->nxtbuck = NULL;
    _W2->strnr = self->STRHASH_txstrhashlist_DOT_fcount;
    result = self->STRHASH_txstrhashlist_DOT_fcount + SYSTEM_ord(self->
      STRHASH_txstrhashlist_DOT_onebased);
    if (self->STRHASH_txstrhashlist_DOT_sortmap != NULL) {
      GMSDATA_txintlist_DOT_setitems(self->
        STRHASH_txstrhashlist_DOT_sortmap,self->
        STRHASH_txstrhashlist_DOT_fcount,self->
        STRHASH_txstrhashlist_DOT_fcount);
      self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_false;
    } 
    self->STRHASH_txstrhashlist_DOT_fcount = self->
      STRHASH_txstrhashlist_DOT_fcount + 1;
    _W2->strp = STRHASH_newstringx(s);
    _W2->obj = aobj;

  }
  return result;
}  /* storeobject */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_indexof(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer hv;
  STRHASH_phashbucket pbuck;

  if (self->STRHASH_txstrhashlist_DOT_phashtable == NULL) 
    STRHASH_txstrhashlist_DOT_hashall(self);
  hv = VirtMethodCall(self, STRHASH_txstrhashlist_DOT_hash_T, 1, (self,
    s));
  pbuck = ValueCast(STRHASH_phashbucket,(*self->
    STRHASH_txstrhashlist_DOT_phashtable)[hv]);
  while (pbuck != NULL) {

    if (!VirtMethodCall(self, STRHASH_txstrhashlist_DOT_entryequal_T, 2, (
      self,pbuck->strp,ValueCast(SYSTEM_P3_pshortstring,s)))) { 
      pbuck = pbuck->nxtbuck;
    } else {
      result = pbuck->strnr + SYSTEM_ord(self->
        STRHASH_txstrhashlist_DOT_onebased);
      return result;
    } 
}
  result =  -1;
  return result;
}  /* indexof */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_add(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;

  result = STRHASH_txstrhashlist_DOT_addobject(self,s,ValueCast(
    SYSTEM_tobject,NULL));
  return result;
}  /* add */

Function(SYSTEM_ansichar *) STRHASH_txstrhashlist_DOT_getstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  STRHASH_txstrhashlist self,
  SYSTEM_integer n)
{
  _P3strcpy(result,_len_ret,*(ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased))))->strp);
  return result;
}  /* getstring */

Function(SYSTEM_ansichar *) STRHASH_txstrhashlist_DOT_getsortedstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  STRHASH_txstrhashlist self,
  SYSTEM_integer n)
{
  if (!self->STRHASH_txstrhashlist_DOT_fsorted) 
    STRHASH_txstrhashlist_DOT_sort(self);
  _P3strcpy(result,_len_ret,*(ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,GMSDATA_txintlist_DOT_getitems(
    self->STRHASH_txstrhashlist_DOT_sortmap,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased)))))->strp);
  return result;
}  /* getsortedstring */

Function(SYSTEM_tobject ) STRHASH_txstrhashlist_DOT_getobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n)
{
  SYSTEM_tobject result;

  result = (ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased))))->obj;
  return result;
}  /* getobject */

Function(SYSTEM_tobject ) STRHASH_txstrhashlist_DOT_getsortedobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n)
{
  SYSTEM_tobject result;

  if (!self->STRHASH_txstrhashlist_DOT_fsorted) 
    STRHASH_txstrhashlist_DOT_sort(self);
  result = (ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,GMSDATA_txintlist_DOT_getitems(
    self->STRHASH_txstrhashlist_DOT_sortmap,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased)))))->obj;
  return result;
}  /* getsortedobject */

Procedure STRHASH_txstrhashlist_DOT_setobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n,
  SYSTEM_tobject aobj)
{
  (ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased))))->obj = aobj;
}  /* setobject */

Procedure STRHASH_txstrhashlist_DOT_setsortedobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n,
  SYSTEM_tobject aobj)
{
  if (!self->STRHASH_txstrhashlist_DOT_fsorted) 
    STRHASH_txstrhashlist_DOT_sort(self);
  (ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,GMSDATA_txintlist_DOT_getitems(
    self->STRHASH_txstrhashlist_DOT_sortmap,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased)))))->obj = aobj;
}  /* setsortedobject */

Procedure STRHASH_txstrhashlist_DOT_savetostream(
  STRHASH_txstrhashlist self,
  GMSSTRM_txstream s)
{
  SYSTEM_integer n;

  GMSSTRM_txstream_DOT_writeinteger(s,self->
    STRHASH_txstrhashlist_DOT_fcount);
  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      GMSSTRM_txstream_DOT_writepstring(s,(ValueCast(
        STRHASH_phashbucket,GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(
        self->STRHASH_txstrhashlist_DOT_buckets,n)))->strp);
    } while (n++ !=  _stop);

  }
}  /* savetostream */

Procedure STRHASH_txstrhashlist_DOT_loadfromstream(
  STRHASH_txstrhashlist self,
  GMSSTRM_txstream s)
{
  SYSTEM_integer n;
  SYSTEM_integer cnt;
  SYSTEM_shortstring s2;

  STRHASH_txstrhashlist_DOT_clear(self);
  cnt = VirtMethodCall(s, GMSSTRM_txstream_DOT_readinteger_T, 7, (s));
  { register SYSTEM_int32 _stop = cnt - 1;
    if ((n = 0) <=  _stop) do {
      GMSSTRM_txstream_DOT_readstring(s2,255,s);
      STRHASH_txstrhashlist_DOT_storeobject(self,s2,ValueCast(
        SYSTEM_tobject,NULL));
    
    } while (n++ !=  _stop);

  }
}  /* loadfromstream */

static Procedure quicksort(
  SYSTEM_integer l,
  SYSTEM_integer r,
  STRHASH_txstrhashlist *_2self)
{
  SYSTEM_integer i, j, p;
  SYSTEM_P3_pshortstring rp;

  do {
    i = l;
    j = r;
    p = ValueCast(SYSTEM_uint32,l + r) >> 1;
    do {
      rp = (ValueCast(STRHASH_phashbucket,
        GMSDATA_tgrowarrayfxd_DOT_getitemptrindx((*_2self)->
        STRHASH_txstrhashlist_DOT_buckets,
        GMSDATA_txintlist_DOT_getitems((*_2self)->
        STRHASH_txstrhashlist_DOT_sortmap,p))))->strp;
      while (VirtMethodCall(*_2self, 
        STRHASH_txstrhashlist_DOT_compare_T, 3, (*_2self,(ValueCast(
        STRHASH_phashbucket,GMSDATA_tgrowarrayfxd_DOT_getitemptrindx((*
        _2self)->STRHASH_txstrhashlist_DOT_buckets,
        GMSDATA_txintlist_DOT_getitems((*_2self)->
        STRHASH_txstrhashlist_DOT_sortmap,i))))->strp,rp)) < 0) {

        _P3inc0(i);
}
      while (VirtMethodCall(*_2self, 
        STRHASH_txstrhashlist_DOT_compare_T, 3, (*_2self,(ValueCast(
        STRHASH_phashbucket,GMSDATA_tgrowarrayfxd_DOT_getitemptrindx((*
        _2self)->STRHASH_txstrhashlist_DOT_buckets,
        GMSDATA_txintlist_DOT_getitems((*_2self)->
        STRHASH_txstrhashlist_DOT_sortmap,j))))->strp,rp)) > 0) {

        _P3dec0(j);
}
      if (i <= j) {
        GMSDATA_txintlist_DOT_exchange((*_2self)->
          STRHASH_txstrhashlist_DOT_sortmap,i,j);
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

Procedure STRHASH_txstrhashlist_DOT_sort(
  STRHASH_txstrhashlist self)
{
  SYSTEM_integer n;
  SYSTEM_P3_pshortstring psn;
  SYSTEM_P3_pshortstring psn1;

  if (self->STRHASH_txstrhashlist_DOT_sortmap == NULL) {
    self->STRHASH_txstrhashlist_DOT_sortmap = ValueCast(
      GMSDATA_txintlist,GMSDATA_txintlist_DOT_create(ValueCast(
      GMSDATA_txintlist,_P3alloc_object(&GMSDATA_txintlist_CD))));
    { register SYSTEM_int32 _stop = self->
        STRHASH_txstrhashlist_DOT_fcount - 1;
      if ((n = 0) <=  _stop) do {
        GMSDATA_txintlist_DOT_setitems(self->
          STRHASH_txstrhashlist_DOT_sortmap,n,n);
      } while (n++ !=  _stop);

    }
    self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_false;
  } 
  if (!self->STRHASH_txstrhashlist_DOT_fsorted) {
    if (self->STRHASH_txstrhashlist_DOT_fcount >= 2) {
      psn = (ValueCast(STRHASH_phashbucket,
        GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
        STRHASH_txstrhashlist_DOT_buckets,0)))->strp;
      { register SYSTEM_int32 _stop = self->
          STRHASH_txstrhashlist_DOT_fcount - 2;
        if ((n = 0) <=  _stop) do {
          psn1 = (ValueCast(STRHASH_phashbucket,
            GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
            STRHASH_txstrhashlist_DOT_buckets,n + 1)))->strp;
          if (VirtMethodCall(self, STRHASH_txstrhashlist_DOT_compare_T, 3, (
            self,psn,psn1)) > 0) {
            quicksort(0,self->STRHASH_txstrhashlist_DOT_fcount - 1,&
              self);
            SYSTEM_break(BRK_1);
          } 
          psn = psn1;
        
CNT_1:;
        } while (n++ !=  _stop);
BRK_1:;

      }
    } 
    self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_true;
  } 
}  /* sort */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_getstringlength(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n)
{
  SYSTEM_integer result;

  result = *ValueCast(SYSTEM_P3_pbyte,(ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,n - SYSTEM_ord(self->
    STRHASH_txstrhashlist_DOT_onebased))))->strp);
  return result;
}  /* getstringlength */

Procedure STRHASH_txstrhashlist_DOT_clearhashtable(
  STRHASH_txstrhashlist self)
{
  _P3freemem(self->STRHASH_txstrhashlist_DOT_phashtable);
  self->STRHASH_txstrhashlist_DOT_phashtable = NULL;
  self->STRHASH_txstrhashlist_DOT_hashsize = 0;
  self->STRHASH_txstrhashlist_DOT_rehashcnt = 0;
}  /* clearhashtable */

Procedure STRHASH_txstrhashlist_DOT_freeitem(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n)
{
}  /* freeitem */

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_memoryused(
  STRHASH_txstrhashlist self)
{
  SYSTEM_integer result;
  SYSTEM_integer n;

  result = 0;
  { register SYSTEM_int32 _stop = self->
      STRHASH_txstrhashlist_DOT_fcount - 1;
    if ((n = 0) <=  _stop) do {
      result = result + SYSTEM_length(*(ValueCast(STRHASH_phashbucket,
        GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
        STRHASH_txstrhashlist_DOT_buckets,n)))->strp);
    } while (n++ !=  _stop);

  }
  result = result + GMSDATA_tgrowarray_DOT_memoryused(ValueCast(
    GMSDATA_tgrowarray,self->STRHASH_txstrhashlist_DOT_buckets));
  if (self->STRHASH_txstrhashlist_DOT_phashtable != NULL) 
    result = result + self->STRHASH_txstrhashlist_DOT_hashsize * sizeof(
      SYSTEM_pointer);
  if (self->STRHASH_txstrhashlist_DOT_sortmap != NULL) 
    result = result + GMSDATA_tgrowarray_DOT_memoryused(ValueCast(
      GMSDATA_tgrowarray,self->STRHASH_txstrhashlist_DOT_sortmap));
  return result;
}  /* memoryused */

Procedure STRHASH_txstrhashlist_DOT_renameentry(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer hv0;
  SYSTEM_integer hv1;
  STRHASH_phashbucket pbuck;
  STRHASH_phashbucket prevbuck;

  n = n - SYSTEM_ord(self->STRHASH_txstrhashlist_DOT_onebased);
  if (self->STRHASH_txstrhashlist_DOT_fsorted) {
    SYSUTILS_P3_freeandnil(&self->STRHASH_txstrhashlist_DOT_sortmap);
    self->STRHASH_txstrhashlist_DOT_fsorted = SYSTEM_false;
  } 
  if (self->STRHASH_txstrhashlist_DOT_phashtable != NULL) {
    hv0 = VirtMethodCall(self, STRHASH_txstrhashlist_DOT_hash_T, 1, (
      self,*(ValueCast(STRHASH_phashbucket,
      GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
      STRHASH_txstrhashlist_DOT_buckets,n)))->strp));
    hv1 = VirtMethodCall(self, STRHASH_txstrhashlist_DOT_hash_T, 1, (
      self,s));
    if (hv0 != hv1) {
      prevbuck = NULL;
      pbuck = ValueCast(STRHASH_phashbucket,(*self->
        STRHASH_txstrhashlist_DOT_phashtable)[hv0]);
      while (SYSTEM_true) {
        if (pbuck->strnr == n) 
          SYSTEM_break(BRK_2);
        prevbuck = pbuck;
        pbuck = pbuck->nxtbuck;
      
CNT_2:;
      }
BRK_2:;
      if (prevbuck == NULL) { 
        (*self->STRHASH_txstrhashlist_DOT_phashtable)[hv0] = pbuck->
          nxtbuck;
      } else 
        prevbuck->nxtbuck = pbuck->nxtbuck;
      pbuck->nxtbuck = ValueCast(STRHASH_phashbucket,(*self->
        STRHASH_txstrhashlist_DOT_phashtable)[hv1]);
      (*self->STRHASH_txstrhashlist_DOT_phashtable)[hv1] = pbuck;
    } 
  } 
  STRHASH_strassignx(&(ValueCast(STRHASH_phashbucket,
    GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(self->
    STRHASH_txstrhashlist_DOT_buckets,n)))->strp,s);
}  /* renameentry */

Function(SYSTEM_boolean ) STRHASH_txcsstrhashlist_DOT_entryequal(
  STRHASH_txcsstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2)
{
  SYSTEM_boolean result;

  result = STRUTILX_pstrequal(ps1,ps2);
  return result;
}  /* entryequal */

Function(SYSTEM_integer ) STRHASH_txcsstrhashlist_DOT_hash(
  STRHASH_txcsstrhashlist self,
  const SYSTEM_ansichar *s)
{
  SYSTEM_integer result;
  SYSTEM_integer i;

  result = 0;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      result = 211 * result + SYSTEM_ord(s[i]);
    } while (i++ !=  _stop);

  }
  result = (result & 2147483647) % self->
    STRHASH_txstrhashlist_DOT_hashsize;
  return result;
}  /* hash */

/* unit strhash */
void _Init_Module_strhash(void)
{
} /* _Init_Module_strhash */

void _Final_Module_strhash(void)
{
} /* _Final_Module_strhash */

