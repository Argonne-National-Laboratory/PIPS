#ifndef _P3___gmsdata___H
#define _P3___gmsdata___H

cnstdef {GMSDATA_bufsize = 16384};
typedef SYSTEM_uint16 _sub_1GMSDATA;
typedef SYSTEM_byte _arr_0GMSDATA[16384];
typedef struct GMSDATA_tgadatabuffer_S {
  SYSTEM_integer bytesused;
  SYSTEM_integer filler;
  _arr_0GMSDATA buffer;
} GMSDATA_tgadatabuffer;

typedef GMSDATA_tgadatabuffer *GMSDATA_pgadatabuffer;
typedef SYSTEM_uint32 _sub_2GMSDATA;
typedef GMSDATA_pgadatabuffer GMSDATA_tgadataarray[268435455];
typedef GMSDATA_tgadataarray *GMSDATA_pgadataarray;
typedef struct GMSDATA_tgrowarray_OD_S* GMSDATA_tgrowarray; /* sy_class */
typedef struct GMSDATA_tgrowarray_OD_S {  /* Objects of 'tgrowarray' */
  SYSTEM_classreference_t CD;  /* = &GMSDATA_tgrowarray_CD */
  GMSDATA_pgadataarray GMSDATA_tgrowarray_DOT_pbase;
  SYSTEM_integer GMSDATA_tgrowarray_DOT_baseallocated;
  SYSTEM_integer GMSDATA_tgrowarray_DOT_baseused;
  GMSDATA_pgadatabuffer GMSDATA_tgrowarray_DOT_pcurrentbuf;
} GMSDATA_tgrowarray_OD;


Constructor(GMSDATA_tgrowarray ) GMSDATA_tgrowarray_DOT_create(
  GMSDATA_tgrowarray self);

Destructor(GMSDATA_tgrowarray ) GMSDATA_tgrowarray_DOT_destroy(
  GMSDATA_tgrowarray self);

Procedure GMSDATA_tgrowarray_DOT_clear(
  GMSDATA_tgrowarray self);

Function(SYSTEM_pointer ) GMSDATA_tgrowarray_DOT_reservemem(
  GMSDATA_tgrowarray self,
  SYSTEM_integer l);

Function(SYSTEM_pointer ) GMSDATA_tgrowarray_DOT_reserveandclear(
  GMSDATA_tgrowarray self,
  SYSTEM_integer l);

Function(SYSTEM_int64 ) GMSDATA_tgrowarray_DOT_memoryused(
  GMSDATA_tgrowarray self);
extern void * const GMSDATA_tgrowarray_VT[];
extern const SYSTEM_classdescriptor_t GMSDATA_tgrowarray_CD;


typedef struct GMSDATA_tgrowarrayfxd_OD_S* GMSDATA_tgrowarrayfxd; /* sy_class */
typedef struct GMSDATA_tgrowarrayfxd_OD_S {  /* Objects of 'tgrowarrayfxd' */
  SYSTEM_classreference_t CD;  /* = &GMSDATA_tgrowarrayfxd_CD */
  GMSDATA_pgadataarray GMSDATA_tgrowarray_DOT_pbase;
  SYSTEM_integer GMSDATA_tgrowarray_DOT_baseallocated;
  SYSTEM_integer GMSDATA_tgrowarray_DOT_baseused;
  GMSDATA_pgadatabuffer GMSDATA_tgrowarray_DOT_pcurrentbuf;
  SYSTEM_integer GMSDATA_tgrowarrayfxd_DOT_fsize;
  SYSTEM_integer GMSDATA_tgrowarrayfxd_DOT_fstorefact;
  SYSTEM_integer GMSDATA_tgrowarrayfxd_DOT_fcount;
} GMSDATA_tgrowarrayfxd_OD;


Constructor(GMSDATA_tgrowarrayfxd ) GMSDATA_tgrowarrayfxd_DOT_create(
  GMSDATA_tgrowarrayfxd self,
  SYSTEM_integer asize);

Function(SYSTEM_pointer ) GMSDATA_tgrowarrayfxd_DOT_additem(
  GMSDATA_tgrowarrayfxd self,
  const SYSTEM_untyped *r);

Function(GMSGEN_pbytedataarray ) 
  GMSDATA_tgrowarrayfxd_DOT_getitemptrindx(
  GMSDATA_tgrowarrayfxd self,
  SYSTEM_integer n);

Procedure GMSDATA_tgrowarrayfxd_DOT_getitem(
  GMSDATA_tgrowarrayfxd self,
  SYSTEM_integer n,
  SYSTEM_untyped *r);

Procedure GMSDATA_tgrowarrayfxd_DOT_clear(
  GMSDATA_tgrowarrayfxd self);
extern void * const GMSDATA_tgrowarrayfxd_VT[];
extern const SYSTEM_classdescriptor_t GMSDATA_tgrowarrayfxd_CD;


typedef struct GMSDATA_txintlist_OD_S* GMSDATA_txintlist; /* sy_class */
typedef struct GMSDATA_txintlist_OD_S {  /* Objects of 'txintlist' */
  SYSTEM_classreference_t CD;  /* = &GMSDATA_txintlist_CD */
  GMSDATA_pgadataarray GMSDATA_tgrowarray_DOT_pbase;
  SYSTEM_integer GMSDATA_tgrowarray_DOT_baseallocated;
  SYSTEM_integer GMSDATA_tgrowarray_DOT_baseused;
  GMSDATA_pgadatabuffer GMSDATA_tgrowarray_DOT_pcurrentbuf;
  SYSTEM_integer GMSDATA_tgrowarrayfxd_DOT_fsize;
  SYSTEM_integer GMSDATA_tgrowarrayfxd_DOT_fstorefact;
  SYSTEM_integer GMSDATA_tgrowarrayfxd_DOT_fcount;
} GMSDATA_txintlist_OD;


Function(SYSTEM_integer ) GMSDATA_txintlist_DOT_getitems(
  GMSDATA_txintlist self,
  SYSTEM_integer index);

Procedure GMSDATA_txintlist_DOT_setitems(
  GMSDATA_txintlist self,
  SYSTEM_integer index,
  SYSTEM_integer v);

Constructor(GMSDATA_txintlist ) GMSDATA_txintlist_DOT_create(
  GMSDATA_txintlist self);

Destructor(GMSDATA_txintlist ) GMSDATA_txintlist_DOT_destroy(
  GMSDATA_txintlist self);

Function(SYSTEM_integer ) GMSDATA_txintlist_DOT_add(
  GMSDATA_txintlist self,
  SYSTEM_integer item);

Procedure GMSDATA_txintlist_DOT_exchange(
  GMSDATA_txintlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);
extern void * const GMSDATA_txintlist_VT[];
extern const SYSTEM_classdescriptor_t GMSDATA_txintlist_CD;


typedef struct GMSDATA_ttblgamsdata_OD_S* GMSDATA_ttblgamsdata; /* sy_class */
typedef struct GMSDATA_ttblgamsdata_OD_S {  /* Objects of 'ttblgamsdata' */
  SYSTEM_classreference_t CD;  /* = &GMSDATA_ttblgamsdata_CD */
  GMSDATA_tgrowarrayfxd GMSDATA_ttblgamsdata_DOT_ds;
  GMSOBJ_txlist GMSDATA_ttblgamsdata_DOT_flist;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT__fdim;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT_findexsize;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT_fdatasize;
  SYSTEM_boolean GMSDATA_ttblgamsdata_DOT_fissorted;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT_flastindex;
} GMSDATA_ttblgamsdata_OD;


Procedure GMSDATA_ttblgamsdata_DOT_quicksort(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer l,
  SYSTEM_integer r);

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_compare(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_comparewithrecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  SYSTEM_integer n);

Procedure GMSDATA_ttblgamsdata_DOT_exchange(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_getcount(
  GMSDATA_ttblgamsdata self);

Procedure GMSDATA_ttblgamsdata_DOT_insertrecord(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  const SYSTEM_integer *inx,
  const SYSTEM_untyped *buffer);

Function(SYSTEM_integer ) GMSDATA_ttblgamsdata_DOT_getcapacity(
  GMSDATA_ttblgamsdata self);

Procedure GMSDATA_ttblgamsdata_DOT_setcapacity(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n);

Constructor(GMSDATA_ttblgamsdata ) GMSDATA_ttblgamsdata_DOT_create(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer adim,
  SYSTEM_integer adatasize);

Destructor(GMSDATA_ttblgamsdata ) GMSDATA_ttblgamsdata_DOT_destroy(
  GMSDATA_ttblgamsdata self);

Procedure GMSDATA_ttblgamsdata_DOT_addrecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  const SYSTEM_untyped *buffer);

Function(SYSTEM_boolean ) GMSDATA_ttblgamsdata_DOT_adduniquerecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  const SYSTEM_untyped *buffer);

Procedure GMSDATA_ttblgamsdata_DOT_getrecord(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  SYSTEM_integer *inx,
  SYSTEM_untyped *buffer);

Procedure GMSDATA_ttblgamsdata_DOT_sort(
  GMSDATA_ttblgamsdata self);

Procedure GMSDATA_ttblgamsdata_DOT_getkeys(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  SYSTEM_integer *inx);

Procedure GMSDATA_ttblgamsdata_DOT_getdata(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n,
  SYSTEM_untyped *buffer);

Function(SYSTEM_pointer ) GMSDATA_ttblgamsdata_DOT_getdataptr(
  GMSDATA_ttblgamsdata self,
  SYSTEM_integer n);

Function(SYSTEM_boolean ) GMSDATA_ttblgamsdata_DOT_searchrecord(
  GMSDATA_ttblgamsdata self,
  const SYSTEM_integer *inx,
  SYSTEM_integer *recnr);

Procedure GMSDATA_ttblgamsdata_DOT_clear(
  GMSDATA_ttblgamsdata self);

Function(SYSTEM_int64 ) GMSDATA_ttblgamsdata_DOT_memoryused(
  GMSDATA_ttblgamsdata self);
extern void * const GMSDATA_ttblgamsdata_VT[];
extern const SYSTEM_classdescriptor_t GMSDATA_ttblgamsdata_CD;


typedef struct GMSDATA_trorcmapper_OD_S* GMSDATA_trorcmapper; /* sy_class */
typedef struct GMSDATA_trorcmapper_OD_S {  /* Objects of 'trorcmapper' */
  SYSTEM_classreference_t CD;  /* = &GMSDATA_trorcmapper_CD */
  GMSDATA_tgrowarrayfxd GMSDATA_ttblgamsdata_DOT_ds;
  GMSOBJ_txlist GMSDATA_ttblgamsdata_DOT_flist;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT__fdim;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT_findexsize;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT_fdatasize;
  SYSTEM_boolean GMSDATA_ttblgamsdata_DOT_fissorted;
  SYSTEM_integer GMSDATA_ttblgamsdata_DOT_flastindex;
} GMSDATA_trorcmapper_OD;


Procedure GMSDATA_trorcmapper_DOT_sort2(
  GMSDATA_trorcmapper self);

Function(SYSTEM_integer ) GMSDATA_trorcmapper_DOT_compare2(
  GMSDATA_trorcmapper self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Procedure GMSDATA_trorcmapper_DOT_quicksort2(
  GMSDATA_trorcmapper self,
  SYSTEM_integer l,
  SYSTEM_integer r);

Constructor(GMSDATA_trorcmapper ) GMSDATA_trorcmapper_DOT_create(
  GMSDATA_trorcmapper self,
  SYSTEM_integer adim);

Procedure GMSDATA_trorcmapper_DOT_addrecord(
  GMSDATA_trorcmapper self,
  const SYSTEM_integer *inx);

Procedure GMSDATA_trorcmapper_DOT_adddummyrecord(
  GMSDATA_trorcmapper self);

Function(SYSTEM_boolean ) GMSDATA_trorcmapper_DOT_addanyrecord(
  GMSDATA_trorcmapper self,
  const SYSTEM_integer *inx);

Procedure GMSDATA_trorcmapper_DOT_finishindex(
  GMSDATA_trorcmapper self,
  SYSTEM_boolean keeppos);

Function(SYSTEM_integer ) GMSDATA_trorcmapper_DOT_getrorc(
  GMSDATA_trorcmapper self,
  const SYSTEM_integer *inx);
extern void * const GMSDATA_trorcmapper_VT[];
extern const SYSTEM_classdescriptor_t GMSDATA_trorcmapper_CD;



extern void _Init_Module_gmsdata(void);
extern void _Final_Module_gmsdata(void);

#endif /* ! defined _P3___gmsdata___H */
