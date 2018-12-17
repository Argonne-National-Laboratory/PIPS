#ifndef _P3___datastorage___H
#define _P3___datastorage___H

typedef struct DATASTORAGE_tgamsdatasearcher_OD_S* 
  DATASTORAGE_tgamsdatasearcher; /* sy_class */
extern void * const DATASTORAGE_tgamsdatasearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasearcher_CD;


typedef struct DATASTORAGE_tgamsdatastore_OD_S* 
  DATASTORAGE_tgamsdatastore; /* sy_class */
typedef struct DATASTORAGE_tgamsdatastore_OD_S {  /* Objects of 'tgamsdatastore' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatastore_CD */
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_ftotalsize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_frecnr;
  GMSGEN_pbytedataarray DATASTORAGE_tgamsdatastore_DOT_pdefrec;
  DATASTORAGE_tgamsdatasearcher DATASTORAGE_tgamsdatastore_DOT_fasrch;
  SYSTEM_boolean DATASTORAGE_tgamsdatastore_DOT_fsawdefrec;
} DATASTORAGE_tgamsdatastore_OD;


Constructor(DATASTORAGE_tgamsdatastore ) 
  DATASTORAGE_tgamsdatastore_DOT_create(
  DATASTORAGE_tgamsdatastore self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec);

Function(SYSTEM_P3_ppointerarray ) 
  DATASTORAGE_tgamsdatastore_DOT_allocptrs(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_freeptrs(
  DATASTORAGE_tgamsdatastore self,
  SYSTEM_P3_ppointerarray p);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatastore_DOT_comparekeys(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone key1,
  GMSGEN_pintegerarrayone key2);

Prototype Function(SYSTEM_integer ) (*
  DATASTORAGE_tgamsdatastore_DOT_getcount_T)(
  DATASTORAGE_tgamsdatastore self);

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatastore_DOT_isdefaultdata(
  DATASTORAGE_tgamsdatastore self,
  const SYSTEM_untyped *adata);

Destructor(DATASTORAGE_tgamsdatastore ) 
  DATASTORAGE_tgamsdatastore_DOT_destroy(
  DATASTORAGE_tgamsdatastore self);

Function(GMSGEN_pintegerarrayone ) 
  DATASTORAGE_tgamsdatastore_DOT_allocindex(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_freeindex(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone p);

Prototype Procedure (*DATASTORAGE_tgamsdatastore_DOT_clear_T)(
  DATASTORAGE_tgamsdatastore self);

Prototype Procedure (*DATASTORAGE_tgamsdatastore_DOT_insertrecord_T)(
  DATASTORAGE_tgamsdatastore self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Prototype Procedure (*DATASTORAGE_tgamsdatastore_DOT_loadrecord_T)(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Prototype Function(SYSTEM_boolean ) (*
  DATASTORAGE_tgamsdatastore_DOT_startread_T)(
  DATASTORAGE_tgamsdatastore self);

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatastore_DOT_startread(
  DATASTORAGE_tgamsdatastore self);

Prototype Function(SYSTEM_P3_pbyte ) (*
  DATASTORAGE_tgamsdatastore_DOT_getnextkey_T)(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone akey);

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatastore_DOT_getnextrecord(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_untyped *adata);

Prototype Procedure (*DATASTORAGE_tgamsdatastore_DOT_endread_T)(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_endread(
  DATASTORAGE_tgamsdatastore self);

Prototype Procedure (*DATASTORAGE_tgamsdatastore_DOT_startassign_T)(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_startassign(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_assignrecord(
  DATASTORAGE_tgamsdatastore self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Prototype Procedure (*DATASTORAGE_tgamsdatastore_DOT_endassign_T)(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_endassign(
  DATASTORAGE_tgamsdatastore self);

Prototype Function(SYSTEM_integer ) (*
  DATASTORAGE_tgamsdatastore_DOT_memoryused_T)(
  DATASTORAGE_tgamsdatastore self);

Prototype Function(DATASTORAGE_tgamsdatasearcher ) (*
  DATASTORAGE_tgamsdatastore_DOT_createsearcher_T)(
  DATASTORAGE_tgamsdatastore self);

Procedure DATASTORAGE_tgamsdatastore_DOT_verify(
  DATASTORAGE_tgamsdatastore self,
  SYSTEM_boolean print);
extern void * const DATASTORAGE_tgamsdatastore_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatastore_CD;


typedef struct DATASTORAGE_tgamsdatasearcher_OD_S {  /* Objects of 'tgamsdatasearcher' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatasearcher_CD */
  DATASTORAGE_tgamsdatastore DATASTORAGE_tgamsdatasearcher_DOT_ds;
} DATASTORAGE_tgamsdatasearcher_OD;


Constructor(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatasearcher_DOT_create(
  DATASTORAGE_tgamsdatasearcher self,
  DATASTORAGE_tgamsdatastore ads);

Prototype Function(SYSTEM_boolean ) (*
  DATASTORAGE_tgamsdatasearcher_DOT_startsearch_T)(
  DATASTORAGE_tgamsdatasearcher self);

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatasearcher_DOT_startsearch(
  DATASTORAGE_tgamsdatasearcher self);

Prototype Function(SYSTEM_boolean ) (*
  DATASTORAGE_tgamsdatasearcher_DOT_search_T)(
  DATASTORAGE_tgamsdatasearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata);
extern void * const DATASTORAGE_tgamsdatasearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasearcher_CD;


typedef struct DATASTORAGE_treclist_OD_S* DATASTORAGE_treclist; /* sy_class */
typedef struct DATASTORAGE_treclist_OD_S {  /* Objects of 'treclist' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_treclist_CD */
  SYSTEM_integer DATASTORAGE_treclist_DOT_fcapacity;
  SYSTEM_integer DATASTORAGE_treclist_DOT_frecsize;
  SYSTEM_integer DATASTORAGE_treclist_DOT_fcount;
  SYSTEM_P3_ppointerarray DATASTORAGE_treclist_DOT_flist;
} DATASTORAGE_treclist_OD;


Function(SYSTEM_pointer ) DATASTORAGE_treclist_DOT_getitem(
  DATASTORAGE_treclist self,
  SYSTEM_integer index);

Procedure DATASTORAGE_treclist_DOT_setcapacity(
  DATASTORAGE_treclist self,
  SYSTEM_integer newcapacity);

Prototype Procedure (*DATASTORAGE_treclist_DOT_grow_T)(
  DATASTORAGE_treclist self);

Procedure DATASTORAGE_treclist_DOT_grow(
  DATASTORAGE_treclist self);

Procedure DATASTORAGE_treclist_DOT_freeitem(
  DATASTORAGE_treclist self,
  SYSTEM_integer index);

Constructor(DATASTORAGE_treclist ) DATASTORAGE_treclist_DOT_create(
  DATASTORAGE_treclist self,
  SYSTEM_integer arecsize);

Destructor(DATASTORAGE_treclist ) DATASTORAGE_treclist_DOT_destroy(
  DATASTORAGE_treclist self);

Function(SYSTEM_pointer ) DATASTORAGE_treclist_DOT_additem(
  DATASTORAGE_treclist self);

Procedure DATASTORAGE_treclist_DOT_clear(
  DATASTORAGE_treclist self);

Function(SYSTEM_pointer ) DATASTORAGE_treclist_DOT_insert(
  DATASTORAGE_treclist self,
  SYSTEM_integer index);

Procedure DATASTORAGE_treclist_DOT_exchange(
  DATASTORAGE_treclist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Procedure DATASTORAGE_treclist_DOT_remove(
  DATASTORAGE_treclist self,
  SYSTEM_integer index);

Procedure DATASTORAGE_treclist_DOT_cleanup(
  DATASTORAGE_treclist self);
extern void * const DATASTORAGE_treclist_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_treclist_CD;


typedef struct DATASTORAGE_tgamsdatatablesearcher_OD_S* 
  DATASTORAGE_tgamsdatatablesearcher; /* sy_class */
extern void * const DATASTORAGE_tgamsdatatablesearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatatablesearcher_CD;


typedef struct DATASTORAGE_tgamsdatatable_OD_S* 
  DATASTORAGE_tgamsdatatable; /* sy_class */
typedef struct DATASTORAGE_tgamsdatatable_OD_S {  /* Objects of 'tgamsdatatable' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatatable_CD */
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_ftotalsize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_frecnr;
  GMSGEN_pbytedataarray DATASTORAGE_tgamsdatastore_DOT_pdefrec;
  DATASTORAGE_tgamsdatasearcher DATASTORAGE_tgamsdatastore_DOT_fasrch;
  SYSTEM_boolean DATASTORAGE_tgamsdatastore_DOT_fsawdefrec;
  DATASTORAGE_treclist DATASTORAGE_tgamsdatatable_DOT_xlist;
} DATASTORAGE_tgamsdatatable_OD;


Function(SYSTEM_integer ) DATASTORAGE_tgamsdatatable_DOT_getcount(
  DATASTORAGE_tgamsdatatable self);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatatable_DOT_getcapacity(
  DATASTORAGE_tgamsdatatable self);

Procedure DATASTORAGE_tgamsdatatable_DOT_setcapacity(
  DATASTORAGE_tgamsdatatable self,
  SYSTEM_integer n);

Constructor(DATASTORAGE_tgamsdatatable ) 
  DATASTORAGE_tgamsdatatable_DOT_create(
  DATASTORAGE_tgamsdatatable self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec);

Destructor(DATASTORAGE_tgamsdatatable ) 
  DATASTORAGE_tgamsdatatable_DOT_destroy(
  DATASTORAGE_tgamsdatatable self);

Procedure DATASTORAGE_tgamsdatatable_DOT_clear(
  DATASTORAGE_tgamsdatatable self);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatatable_DOT_memoryused(
  DATASTORAGE_tgamsdatatable self);

Procedure DATASTORAGE_tgamsdatatable_DOT_insertrecord(
  DATASTORAGE_tgamsdatatable self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Procedure DATASTORAGE_tgamsdatatable_DOT_loadrecord(
  DATASTORAGE_tgamsdatatable self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatatable_DOT_getnextkey(
  DATASTORAGE_tgamsdatatable self,
  GMSGEN_pintegerarrayone akey);

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatatable_DOT_createsearcher(
  DATASTORAGE_tgamsdatatable self);

Procedure DATASTORAGE_tgamsdatatable_DOT_endassign(
  DATASTORAGE_tgamsdatatable self);

Procedure DATASTORAGE_tgamsdatatable_DOT_sort(
  DATASTORAGE_tgamsdatatable self);
extern void * const DATASTORAGE_tgamsdatatable_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatatable_CD;


typedef struct DATASTORAGE_tgamsdatatablesearcher_OD_S {  /* Objects of 'tgamsdatatablesearcher' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatatablesearcher_CD */
  DATASTORAGE_tgamsdatastore DATASTORAGE_tgamsdatasearcher_DOT_ds;
  SYSTEM_integer DATASTORAGE_tgamsdatatablesearcher_DOT_flastindex;
} DATASTORAGE_tgamsdatatablesearcher_OD;


Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatatablesearcher_DOT_startsearch(
  DATASTORAGE_tgamsdatatablesearcher self);

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatatablesearcher_DOT_search(
  DATASTORAGE_tgamsdatatablesearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata);
extern void * const DATASTORAGE_tgamsdatatablesearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatatablesearcher_CD;


typedef struct DATASTORAGE_tintegermapping_OD_S* 
  DATASTORAGE_tintegermapping; /* sy_class */
typedef struct DATASTORAGE_tintegermapping_OD_S {  /* Objects of 'tintegermapping' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tintegermapping_CD */
  SYSTEM_integer DATASTORAGE_tintegermapping_DOT_fhighestindex;
  SYSTEM_integer DATASTORAGE_tintegermapping_DOT_fcapacity;
  SYSTEM_P3_pintegerarray DATASTORAGE_tintegermapping_DOT_pmap;
} DATASTORAGE_tintegermapping_OD;


Procedure DATASTORAGE_tintegermapping_DOT_setmapping(
  DATASTORAGE_tintegermapping self,
  SYSTEM_integer f,
  SYSTEM_integer t);

Function(SYSTEM_integer ) DATASTORAGE_tintegermapping_DOT_getmapping(
  DATASTORAGE_tintegermapping self,
  SYSTEM_integer f);

Constructor(DATASTORAGE_tintegermapping ) 
  DATASTORAGE_tintegermapping_DOT_create(
  DATASTORAGE_tintegermapping self);

Destructor(DATASTORAGE_tintegermapping ) 
  DATASTORAGE_tintegermapping_DOT_destroy(
  DATASTORAGE_tintegermapping self);

Function(SYSTEM_integer ) DATASTORAGE_tintegermapping_DOT_memoryused(
  DATASTORAGE_tintegermapping self);
extern void * const DATASTORAGE_tintegermapping_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tintegermapping_CD;


typedef SYSTEM_uint32 _sub_0DATASTORAGE;
typedef DATASTORAGE_tintegermapping DATASTORAGE_tintegermappingarray[10000000];
typedef DATASTORAGE_tintegermappingarray *
  DATASTORAGE_pintegermappingarray;
typedef struct DATASTORAGE_tgamsdatafull_OD_S* 
  DATASTORAGE_tgamsdatafull; /* sy_class */
typedef struct DATASTORAGE_tgamsdatafull_OD_S {  /* Objects of 'tgamsdatafull' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatafull_CD */
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_ftotalsize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_frecnr;
  GMSGEN_pbytedataarray DATASTORAGE_tgamsdatastore_DOT_pdefrec;
  DATASTORAGE_tgamsdatasearcher DATASTORAGE_tgamsdatastore_DOT_fasrch;
  SYSTEM_boolean DATASTORAGE_tgamsdatastore_DOT_fsawdefrec;
  SYSTEM_P3_pbyte DATASTORAGE_tgamsdatafull_DOT_pdata;
  SYSTEM_integer DATASTORAGE_tgamsdatafull_DOT_fcount;
  GMSGEN_pintegerarrayone DATASTORAGE_tgamsdatafull_DOT_xmult;
  GMSGEN_pintegerarrayone DATASTORAGE_tgamsdatafull_DOT_xmin;
  GMSGEN_pintegerarrayone DATASTORAGE_tgamsdatafull_DOT_xcnt;
  SYSTEM_integer DATASTORAGE_tgamsdatafull_DOT_fallocsize;
  DATASTORAGE_pintegermappingarray DATASTORAGE_tgamsdatafull_DOT_finxmap;
} DATASTORAGE_tgamsdatafull_OD;


Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatafull_DOT_getoffset(
  DATASTORAGE_tgamsdatafull self,
  GMSGEN_pintegerarrayone akey);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatafull_DOT_getcount(
  DATASTORAGE_tgamsdatafull self);

Constructor(DATASTORAGE_tgamsdatafull ) 
  DATASTORAGE_tgamsdatafull_DOT_create(
  DATASTORAGE_tgamsdatafull self,
  DATASTORAGE_tgamsdatastore ds);

Destructor(DATASTORAGE_tgamsdatafull ) 
  DATASTORAGE_tgamsdatafull_DOT_destroy(
  DATASTORAGE_tgamsdatafull self);

Procedure DATASTORAGE_tgamsdatafull_DOT_allocatememory(
  DATASTORAGE_tgamsdatafull self,
  DATASTORAGE_tgamsdatastore ds);

Procedure DATASTORAGE_tgamsdatafull_DOT_clear(
  DATASTORAGE_tgamsdatafull self);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatafull_DOT_memoryused(
  DATASTORAGE_tgamsdatafull self);

Procedure DATASTORAGE_tgamsdatafull_DOT_loadrecord(
  DATASTORAGE_tgamsdatafull self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatafull_DOT_createsearcher(
  DATASTORAGE_tgamsdatafull self);
extern void * const DATASTORAGE_tgamsdatafull_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatafull_CD;


typedef struct DATASTORAGE_tgamsdatafullsearcher_OD_S* 
  DATASTORAGE_tgamsdatafullsearcher; /* sy_class */
typedef struct DATASTORAGE_tgamsdatafullsearcher_OD_S {  /* Objects of 'tgamsdatafullsearcher' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatafullsearcher_CD */
  DATASTORAGE_tgamsdatastore DATASTORAGE_tgamsdatasearcher_DOT_ds;
} DATASTORAGE_tgamsdatafullsearcher_OD;


Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatafullsearcher_DOT_search(
  DATASTORAGE_tgamsdatafullsearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata);
extern void * const DATASTORAGE_tgamsdatafullsearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatafullsearcher_CD;


typedef struct DATASTORAGE_tgamsdatasparse_OD_S* 
  DATASTORAGE_tgamsdatasparse; /* sy_class */
typedef struct DATASTORAGE_tgamsdatasparse_OD_S {  /* Objects of 'tgamsdatasparse' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatasparse_CD */
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_ftotalsize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_frecnr;
  GMSGEN_pbytedataarray DATASTORAGE_tgamsdatastore_DOT_pdefrec;
  DATASTORAGE_tgamsdatasearcher DATASTORAGE_tgamsdatastore_DOT_fasrch;
  SYSTEM_boolean DATASTORAGE_tgamsdatastore_DOT_fsawdefrec;
  SYSTEM_P3_ppointerarray DATASTORAGE_tgamsdatasparse_DOT_fwrkptrs;
  SYSTEM_integer DATASTORAGE_tgamsdatasparse_DOT_fcellcount;
  SYSTEM_integer DATASTORAGE_tgamsdatasparse_DOT_fcelldatacount;
} DATASTORAGE_tgamsdatasparse_OD;


Function(SYSTEM_integer ) DATASTORAGE_tgamsdatasparse_DOT_getcount(
  DATASTORAGE_tgamsdatasparse self);

Function(SYSTEM_pointer ) DATASTORAGE_tgamsdatasparse_DOT_getcell(
  DATASTORAGE_tgamsdatasparse self,
  SYSTEM_integer d);

Procedure DATASTORAGE_tgamsdatasparse_DOT_freecell(
  DATASTORAGE_tgamsdatasparse self,
  SYSTEM_pointer p,
  SYSTEM_integer d);

Constructor(DATASTORAGE_tgamsdatasparse ) 
  DATASTORAGE_tgamsdatasparse_DOT_create(
  DATASTORAGE_tgamsdatasparse self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec);

Destructor(DATASTORAGE_tgamsdatasparse ) 
  DATASTORAGE_tgamsdatasparse_DOT_destroy(
  DATASTORAGE_tgamsdatasparse self);

Procedure DATASTORAGE_tgamsdatasparse_DOT_clear(
  DATASTORAGE_tgamsdatasparse self);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatasparse_DOT_memoryused(
  DATASTORAGE_tgamsdatasparse self);

Procedure DATASTORAGE_tgamsdatasparse_DOT_insertrecord(
  DATASTORAGE_tgamsdatasparse self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Procedure DATASTORAGE_tgamsdatasparse_DOT_loadrecord(
  DATASTORAGE_tgamsdatasparse self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatasparse_DOT_startread(
  DATASTORAGE_tgamsdatasparse self);

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatasparse_DOT_getnextkey(
  DATASTORAGE_tgamsdatasparse self,
  GMSGEN_pintegerarrayone akey);

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatasparse_DOT_createsearcher(
  DATASTORAGE_tgamsdatasparse self);

Procedure DATASTORAGE_tgamsdatasparse_DOT_endassign(
  DATASTORAGE_tgamsdatasparse self);
extern void * const DATASTORAGE_tgamsdatasparse_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasparse_CD;


typedef struct DATASTORAGE_tgamsdatasparsesearcher_OD_S* 
  DATASTORAGE_tgamsdatasparsesearcher; /* sy_class */
typedef struct DATASTORAGE_tgamsdatasparsesearcher_OD_S {  /* Objects of 'tgamsdatasparsesearcher' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatasparsesearcher_CD */
  DATASTORAGE_tgamsdatastore DATASTORAGE_tgamsdatasearcher_DOT_ds;
  SYSTEM_P3_ppointerarray DATASTORAGE_tgamsdatasparsesearcher_DOT_fsearchptrs;
  SYSTEM_integer DATASTORAGE_tgamsdatasparsesearcher_DOT__fdim;
  SYSTEM_integer DATASTORAGE_tgamsdatasparsesearcher_DOT_flastvalid;
} DATASTORAGE_tgamsdatasparsesearcher_OD;


Constructor(DATASTORAGE_tgamsdatasparsesearcher ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_create(
  DATASTORAGE_tgamsdatasparsesearcher self,
  DATASTORAGE_tgamsdatasparse ads);

Destructor(DATASTORAGE_tgamsdatasparsesearcher ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_destroy(
  DATASTORAGE_tgamsdatasparsesearcher self);

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_startsearch(
  DATASTORAGE_tgamsdatasparsesearcher self);

Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatasparsesearcher_DOT_search(
  DATASTORAGE_tgamsdatasparsesearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata);
extern void * const DATASTORAGE_tgamsdatasparsesearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatasparsesearcher_CD;


typedef struct DATASTORAGE_tlinkeddata_OD_S* DATASTORAGE_tlinkeddata; /* sy_class */
typedef struct DATASTORAGE_tlinkeddata_OD_S {  /* Objects of 'tlinkeddata' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tlinkeddata_CD */
  GMSHEAPNEW_theapmgr DATASTORAGE_tlinkeddata_DOT_myheap;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_fmaxkey;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_fminkey;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_ftotalsize;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_fdatasize;
  SYSTEM_pointer DATASTORAGE_tlinkeddata_DOT_fhead;
  SYSTEM_pointer DATASTORAGE_tlinkeddata_DOT_ftail;
  SYSTEM_integer DATASTORAGE_tlinkeddata_DOT_fcount;
} DATASTORAGE_tlinkeddata_OD;


Function(SYSTEM_boolean ) DATASTORAGE_tlinkeddata_DOT_removedefaults(
  DATASTORAGE_tlinkeddata self,
  const SYSTEM_untyped *defdata);

Constructor(DATASTORAGE_tlinkeddata ) 
  DATASTORAGE_tlinkeddata_DOT_create(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize);

Destructor(DATASTORAGE_tlinkeddata ) 
  DATASTORAGE_tlinkeddata_DOT_destroy(
  DATASTORAGE_tlinkeddata self);

Procedure DATASTORAGE_tlinkeddata_DOT_clear(
  DATASTORAGE_tlinkeddata self);

Function(SYSTEM_integer ) DATASTORAGE_tlinkeddata_DOT_memoryused(
  DATASTORAGE_tlinkeddata self);

Function(GMSGEN_pintegerarrayone ) 
  DATASTORAGE_tlinkeddata_DOT_allocindex(
  DATASTORAGE_tlinkeddata self);

Procedure DATASTORAGE_tlinkeddata_DOT_freeindex(
  DATASTORAGE_tlinkeddata self,
  GMSGEN_pintegerarrayone p);

Function(SYSTEM_pointer ) DATASTORAGE_tlinkeddata_DOT_additem(
  DATASTORAGE_tlinkeddata self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Function(SYSTEM_boolean ) DATASTORAGE_tlinkeddata_DOT_startread(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_pointer *p,
  GMSGEN_pintegerarrayone amap);

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tlinkeddata_DOT_getnextkey(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_pointer *p,
  GMSGEN_pintegerarrayone akey);

Function(SYSTEM_boolean ) DATASTORAGE_tlinkeddata_DOT_getnextrecord(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_pointer *p,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_untyped *data);

Procedure DATASTORAGE_tlinkeddata_DOT_sort(
  DATASTORAGE_tlinkeddata self,
  GMSGEN_pintegerarrayone amap);

Procedure DATASTORAGE_tlinkeddata_DOT_createlist(
  DATASTORAGE_tlinkeddata self,
  SYSTEM_P3_ppointerarray *alist);
extern void * const DATASTORAGE_tlinkeddata_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tlinkeddata_CD;


typedef struct DATASTORAGE_tgamshashlist_OD_S* 
  DATASTORAGE_tgamshashlist; /* sy_class */
typedef struct DATASTORAGE_tgamshashlist_OD_S {  /* Objects of 'tgamshashlist' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamshashlist_CD */
  SYSTEM_P3_ppointerarray DATASTORAGE_tgamshashlist_DOT_phashtable;
  DATASTORAGE_tlinkeddata DATASTORAGE_tgamshashlist_DOT_lnkdata;
  SYSTEM_integer DATASTORAGE_tgamshashlist_DOT_hashsize;
  SYSTEM_integer DATASTORAGE_tgamshashlist_DOT_rehashcnt;
  SYSTEM_integer DATASTORAGE_tgamshashlist_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tgamshashlist_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tgamshashlist_DOT_fdatasize;
  SYSTEM_pointer DATASTORAGE_tgamshashlist_DOT_fcurrreadp;
} DATASTORAGE_tgamshashlist_OD;


Procedure DATASTORAGE_tgamshashlist_DOT_hashtablereset(
  DATASTORAGE_tgamshashlist self,
  SYSTEM_integer acnt);

Procedure DATASTORAGE_tgamshashlist_DOT_hashall(
  DATASTORAGE_tgamshashlist self);

Function(SYSTEM_integer ) DATASTORAGE_tgamshashlist_DOT_getcount(
  DATASTORAGE_tgamshashlist self);

Function(SYSTEM_integer ) DATASTORAGE_tgamshashlist_DOT_hash(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey);

Function(SYSTEM_boolean ) DATASTORAGE_tgamshashlist_DOT_equalkeys(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey1,
  GMSGEN_pintegerarrayone akey2);

Procedure DATASTORAGE_tgamshashlist_DOT_removedefaults(
  DATASTORAGE_tgamshashlist self,
  const SYSTEM_untyped *defdata);

Constructor(DATASTORAGE_tgamshashlist ) 
  DATASTORAGE_tgamshashlist_DOT_create(
  DATASTORAGE_tgamshashlist self,
  SYSTEM_integer adim,
  SYSTEM_integer adatasize);

Destructor(DATASTORAGE_tgamshashlist ) 
  DATASTORAGE_tgamshashlist_DOT_destroy(
  DATASTORAGE_tgamshashlist self);

Procedure DATASTORAGE_tgamshashlist_DOT_clear(
  DATASTORAGE_tgamshashlist self);

Function(SYSTEM_pointer ) DATASTORAGE_tgamshashlist_DOT_indexof(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey);

Function(SYSTEM_boolean ) DATASTORAGE_tgamshashlist_DOT_additem(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Procedure DATASTORAGE_tgamshashlist_DOT_loadrecord(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Function(SYSTEM_boolean ) DATASTORAGE_tgamshashlist_DOT_startread(
  DATASTORAGE_tgamshashlist self);

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamshashlist_DOT_getnextkey(
  DATASTORAGE_tgamshashlist self,
  GMSGEN_pintegerarrayone akey);

Procedure DATASTORAGE_tgamshashlist_DOT_clearhashlist(
  DATASTORAGE_tgamshashlist self);

Function(SYSTEM_integer ) DATASTORAGE_tgamshashlist_DOT_memoryused(
  DATASTORAGE_tgamshashlist self);
extern void * const DATASTORAGE_tgamshashlist_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamshashlist_CD;


typedef struct DATASTORAGE_tgamsdatahashedsearcher_OD_S* 
  DATASTORAGE_tgamsdatahashedsearcher; /* sy_class */
extern void * const DATASTORAGE_tgamsdatahashedsearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatahashedsearcher_CD;


typedef struct DATASTORAGE_tgamsdatahashed_OD_S* 
  DATASTORAGE_tgamsdatahashed; /* sy_class */
typedef struct DATASTORAGE_tgamsdatahashed_OD_S {  /* Objects of 'tgamsdatahashed' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatahashed_CD */
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdimension;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fkeysize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_ftotalsize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_fdatasize;
  SYSTEM_integer DATASTORAGE_tgamsdatastore_DOT_frecnr;
  GMSGEN_pbytedataarray DATASTORAGE_tgamsdatastore_DOT_pdefrec;
  DATASTORAGE_tgamsdatasearcher DATASTORAGE_tgamsdatastore_DOT_fasrch;
  SYSTEM_boolean DATASTORAGE_tgamsdatastore_DOT_fsawdefrec;
  DATASTORAGE_tgamshashlist DATASTORAGE_tgamsdatahashed_DOT_hl;
} DATASTORAGE_tgamsdatahashed_OD;


Function(SYSTEM_integer ) DATASTORAGE_tgamsdatahashed_DOT_getcount(
  DATASTORAGE_tgamsdatahashed self);

Constructor(DATASTORAGE_tgamsdatahashed ) 
  DATASTORAGE_tgamsdatahashed_DOT_create(
  DATASTORAGE_tgamsdatahashed self,
  SYSTEM_integer adimension,
  SYSTEM_integer adatasize,
  const SYSTEM_untyped *adefrec);

Destructor(DATASTORAGE_tgamsdatahashed ) 
  DATASTORAGE_tgamsdatahashed_DOT_destroy(
  DATASTORAGE_tgamsdatahashed self);

Procedure DATASTORAGE_tgamsdatahashed_DOT_clear(
  DATASTORAGE_tgamsdatahashed self);

Procedure DATASTORAGE_tgamsdatahashed_DOT_insertrecord(
  DATASTORAGE_tgamsdatahashed self,
  DATASTORAGE_tgamsdatasearcher ads,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Procedure DATASTORAGE_tgamsdatahashed_DOT_loadrecord(
  DATASTORAGE_tgamsdatahashed self,
  GMSGEN_pintegerarrayone akey,
  const SYSTEM_untyped *adata);

Function(SYSTEM_boolean ) DATASTORAGE_tgamsdatahashed_DOT_startread(
  DATASTORAGE_tgamsdatahashed self);

Function(SYSTEM_P3_pbyte ) DATASTORAGE_tgamsdatahashed_DOT_getnextkey(
  DATASTORAGE_tgamsdatahashed self,
  GMSGEN_pintegerarrayone akey);

Function(SYSTEM_integer ) DATASTORAGE_tgamsdatahashed_DOT_memoryused(
  DATASTORAGE_tgamsdatahashed self);

Function(DATASTORAGE_tgamsdatasearcher ) 
  DATASTORAGE_tgamsdatahashed_DOT_createsearcher(
  DATASTORAGE_tgamsdatahashed self);

Procedure DATASTORAGE_tgamsdatahashed_DOT_endassign(
  DATASTORAGE_tgamsdatahashed self);
extern void * const DATASTORAGE_tgamsdatahashed_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatahashed_CD;


typedef struct DATASTORAGE_tgamsdatahashedsearcher_OD_S {  /* Objects of 'tgamsdatahashedsearcher' */
  SYSTEM_classreference_t CD;  /* = &DATASTORAGE_tgamsdatahashedsearcher_CD */
  DATASTORAGE_tgamsdatastore DATASTORAGE_tgamsdatasearcher_DOT_ds;
} DATASTORAGE_tgamsdatahashedsearcher_OD;


Function(SYSTEM_boolean ) 
  DATASTORAGE_tgamsdatahashedsearcher_DOT_search(
  DATASTORAGE_tgamsdatahashedsearcher self,
  GMSGEN_pintegerarrayone akey,
  SYSTEM_pointer *apdata);
extern void * const DATASTORAGE_tgamsdatahashedsearcher_VT[];
extern const SYSTEM_classdescriptor_t DATASTORAGE_tgamsdatahashedsearcher_CD;



extern void _Init_Module_datastorage(void);
extern void _Final_Module_datastorage(void);

#endif /* ! defined _P3___datastorage___H */
