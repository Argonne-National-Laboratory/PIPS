#ifndef _P3___strhash___H
#define _P3___strhash___H

typedef struct STRHASH_txstrhashlist_OD_S* STRHASH_txstrhashlist; /* sy_class */
typedef struct STRHASH_txstrhashlist_OD_S {  /* Objects of 'txstrhashlist' */
  SYSTEM_classreference_t CD;  /* = &STRHASH_txstrhashlist_CD */
  GMSDATA_tgrowarrayfxd STRHASH_txstrhashlist_DOT_buckets;
  SYSTEM_P3_ppointerarray STRHASH_txstrhashlist_DOT_phashtable;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_hashsize;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_rehashcnt;
  GMSDATA_txintlist STRHASH_txstrhashlist_DOT_sortmap;
  SYSTEM_boolean STRHASH_txstrhashlist_DOT_fsorted;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_fcount;
  SYSTEM_boolean STRHASH_txstrhashlist_DOT_onebased;
} STRHASH_txstrhashlist_OD;


Procedure STRHASH_txstrhashlist_DOT_clearhashtable(
  STRHASH_txstrhashlist self);

Procedure STRHASH_txstrhashlist_DOT_hashtablereset(
  STRHASH_txstrhashlist self,
  SYSTEM_integer acnt);

Prototype Function(SYSTEM_integer ) (*STRHASH_txstrhashlist_DOT_hash_T)(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_hash(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s);

Prototype Function(SYSTEM_boolean ) (*
  STRHASH_txstrhashlist_DOT_entryequal_T)(
  STRHASH_txstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2);

Function(SYSTEM_boolean ) STRHASH_txstrhashlist_DOT_entryequal(
  STRHASH_txstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2);

Prototype Function(SYSTEM_integer ) (*
  STRHASH_txstrhashlist_DOT_compare_T)(
  STRHASH_txstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_compare(
  STRHASH_txstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2);

Procedure STRHASH_txstrhashlist_DOT_hashall(
  STRHASH_txstrhashlist self);

Function(SYSTEM_tobject ) STRHASH_txstrhashlist_DOT_getobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Function(SYSTEM_tobject ) STRHASH_txstrhashlist_DOT_getsortedobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Procedure STRHASH_txstrhashlist_DOT_setobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n,
  SYSTEM_tobject aobj);

Procedure STRHASH_txstrhashlist_DOT_setsortedobject(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n,
  SYSTEM_tobject aobj);

Function(SYSTEM_ansichar *) STRHASH_txstrhashlist_DOT_getstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) STRHASH_txstrhashlist_DOT_getsortedstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Procedure STRHASH_txstrhashlist_DOT_sort(
  STRHASH_txstrhashlist self);

Constructor(STRHASH_txstrhashlist ) STRHASH_txstrhashlist_DOT_create(
  STRHASH_txstrhashlist self);

Destructor(STRHASH_txstrhashlist ) STRHASH_txstrhashlist_DOT_destroy(
  STRHASH_txstrhashlist self);

Procedure STRHASH_txstrhashlist_DOT_clear(
  STRHASH_txstrhashlist self);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_storeobject(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_tobject aobj);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_addobject(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_tobject aobj);

Prototype Procedure (*STRHASH_txstrhashlist_DOT_freeitem_T)(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Procedure STRHASH_txstrhashlist_DOT_freeitem(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_add(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_indexof(
  STRHASH_txstrhashlist self,
  const SYSTEM_ansichar *s);

Procedure STRHASH_txstrhashlist_DOT_loadfromstream(
  STRHASH_txstrhashlist self,
  GMSSTRM_txstream s);

Procedure STRHASH_txstrhashlist_DOT_savetostream(
  STRHASH_txstrhashlist self,
  GMSSTRM_txstream s);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_getstringlength(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n);

Function(SYSTEM_integer ) STRHASH_txstrhashlist_DOT_memoryused(
  STRHASH_txstrhashlist self);

Procedure STRHASH_txstrhashlist_DOT_renameentry(
  STRHASH_txstrhashlist self,
  SYSTEM_integer n,
  const SYSTEM_ansichar *s);
extern void * const STRHASH_txstrhashlist_VT[];
extern const SYSTEM_classdescriptor_t STRHASH_txstrhashlist_CD;


typedef struct STRHASH_txcsstrhashlist_OD_S* STRHASH_txcsstrhashlist; /* sy_class */
typedef struct STRHASH_txcsstrhashlist_OD_S {  /* Objects of 'txcsstrhashlist' */
  SYSTEM_classreference_t CD;  /* = &STRHASH_txcsstrhashlist_CD */
  GMSDATA_tgrowarrayfxd STRHASH_txstrhashlist_DOT_buckets;
  SYSTEM_P3_ppointerarray STRHASH_txstrhashlist_DOT_phashtable;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_hashsize;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_rehashcnt;
  GMSDATA_txintlist STRHASH_txstrhashlist_DOT_sortmap;
  SYSTEM_boolean STRHASH_txstrhashlist_DOT_fsorted;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_fcount;
  SYSTEM_boolean STRHASH_txstrhashlist_DOT_onebased;
} STRHASH_txcsstrhashlist_OD;


Function(SYSTEM_integer ) STRHASH_txcsstrhashlist_DOT_hash(
  STRHASH_txcsstrhashlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_boolean ) STRHASH_txcsstrhashlist_DOT_entryequal(
  STRHASH_txcsstrhashlist self,
  SYSTEM_P3_pshortstring ps1,
  SYSTEM_P3_pshortstring ps2);
extern void * const STRHASH_txcsstrhashlist_VT[];
extern const SYSTEM_classdescriptor_t STRHASH_txcsstrhashlist_CD;



extern void _Init_Module_strhash(void);
extern void _Final_Module_strhash(void);

#endif /* ! defined _P3___strhash___H */
