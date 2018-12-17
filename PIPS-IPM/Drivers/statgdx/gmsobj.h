#ifndef _P3___gmsobj___H
#define _P3___gmsobj___H

typedef struct GMSOBJ_tquicksortclass_OD_S* GMSOBJ_tquicksortclass; /* sy_class */
typedef struct GMSOBJ_tquicksortclass_OD_S {  /* Objects of 'tquicksortclass' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_tquicksortclass_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
} GMSOBJ_tquicksortclass_OD;


Procedure GMSOBJ_tquicksortclass_DOT_quicksort(
  GMSOBJ_tquicksortclass self,
  SYSTEM_integer l,
  SYSTEM_integer r);

Prototype Procedure (*GMSOBJ_tquicksortclass_DOT_exchange_T)(
  GMSOBJ_tquicksortclass self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Prototype Function(SYSTEM_integer ) (*
  GMSOBJ_tquicksortclass_DOT_compare_T)(
  GMSOBJ_tquicksortclass self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Procedure GMSOBJ_tquicksortclass_DOT_sortn(
  GMSOBJ_tquicksortclass self,
  SYSTEM_integer n);
extern void * const GMSOBJ_tquicksortclass_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_tquicksortclass_CD;


typedef struct GMSOBJ_txlist_OD_S* GMSOBJ_txlist; /* sy_class */
typedef struct GMSOBJ_txlist_OD_S {  /* Objects of 'txlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txlist_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txlist_DOT_fcapacity;
  SYSTEM_integer GMSOBJ_txlist_DOT_fcount;
  SYSTEM_P3_ppointerarray GMSOBJ_txlist_DOT_flist;
} GMSOBJ_txlist_OD;


Function(SYSTEM_pointer ) GMSOBJ_txlist_DOT_get(
  GMSOBJ_txlist self,
  SYSTEM_integer index);

Function(SYSTEM_pointer ) GMSOBJ_txlist_DOT_getlast(
  GMSOBJ_txlist self);

Procedure GMSOBJ_txlist_DOT_put(
  GMSOBJ_txlist self,
  SYSTEM_integer index,
  SYSTEM_pointer item);

Procedure GMSOBJ_txlist_DOT_setcapacity(
  GMSOBJ_txlist self,
  SYSTEM_integer newcapacity);

Procedure GMSOBJ_txlist_DOT_setcount(
  GMSOBJ_txlist self,
  SYSTEM_integer newcount);

Prototype Procedure (*GMSOBJ_txlist_DOT_grow_T)(
  GMSOBJ_txlist self);

Procedure GMSOBJ_txlist_DOT_grow(
  GMSOBJ_txlist self);

Prototype Procedure (*GMSOBJ_txlist_DOT_freeitem_T)(
  GMSOBJ_txlist self,
  SYSTEM_integer index);

Procedure GMSOBJ_txlist_DOT_freeitem(
  GMSOBJ_txlist self,
  SYSTEM_integer index);

Destructor(GMSOBJ_txlist ) GMSOBJ_txlist_DOT_destroy(
  GMSOBJ_txlist self);

Function(SYSTEM_integer ) GMSOBJ_txlist_DOT_add(
  GMSOBJ_txlist self,
  SYSTEM_pointer item);

Procedure GMSOBJ_txlist_DOT_clear(
  GMSOBJ_txlist self);

Procedure GMSOBJ_txlist_DOT_delete(
  GMSOBJ_txlist self,
  SYSTEM_integer index);

Function(SYSTEM_pointer ) GMSOBJ_txlist_DOT_extract(
  GMSOBJ_txlist self,
  SYSTEM_pointer item);

Function(SYSTEM_integer ) GMSOBJ_txlist_DOT_indexof(
  GMSOBJ_txlist self,
  SYSTEM_pointer item);

Procedure GMSOBJ_txlist_DOT_insert(
  GMSOBJ_txlist self,
  SYSTEM_integer index,
  SYSTEM_pointer item);

Function(SYSTEM_integer ) GMSOBJ_txlist_DOT_remove(
  GMSOBJ_txlist self,
  SYSTEM_pointer item);

Function(SYSTEM_int64 ) GMSOBJ_txlist_DOT_memoryused(
  GMSOBJ_txlist self);

Procedure GMSOBJ_txlist_DOT_exchange(
  GMSOBJ_txlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);
extern void * const GMSOBJ_txlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txlist_CD;


typedef struct GMSOBJ_txstrings_OD_S* GMSOBJ_txstrings; /* sy_class */
typedef struct GMSOBJ_txstrings_OD_S {  /* Objects of 'txstrings' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txstrings_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txlist_DOT_fcapacity;
  SYSTEM_integer GMSOBJ_txlist_DOT_fcount;
  SYSTEM_P3_ppointerarray GMSOBJ_txlist_DOT_flist;
  SYSTEM_int64 GMSOBJ_txstrings_DOT_fstrmemory;
} GMSOBJ_txstrings_OD;


Procedure GMSOBJ_txstrings_DOT_put(
  GMSOBJ_txstrings self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *_ftmp1);

Function(SYSTEM_ansichar *) GMSOBJ_txstrings_DOT_get(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrings self,
  SYSTEM_integer index);

Procedure GMSOBJ_txstrings_DOT_freeitem(
  GMSOBJ_txstrings self,
  SYSTEM_integer index);

Function(SYSTEM_integer ) GMSOBJ_txstrings_DOT_add(
  GMSOBJ_txstrings self,
  const SYSTEM_ansichar *item);

Function(SYSTEM_ansichar *) GMSOBJ_txstrings_DOT_extract(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrings self,
  const SYSTEM_ansichar *item);

Function(SYSTEM_integer ) GMSOBJ_txstrings_DOT_indexof(
  GMSOBJ_txstrings self,
  const SYSTEM_ansichar *item);

Procedure GMSOBJ_txstrings_DOT_insert(
  GMSOBJ_txstrings self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *item);

Function(SYSTEM_integer ) GMSOBJ_txstrings_DOT_compare(
  GMSOBJ_txstrings self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Procedure GMSOBJ_txstrings_DOT_sort(
  GMSOBJ_txstrings self);

Function(SYSTEM_int64 ) GMSOBJ_txstrings_DOT_memoryused(
  GMSOBJ_txstrings self);
extern void * const GMSOBJ_txstrings_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txstrings_CD;


typedef struct GMSOBJ_txpcharlist_OD_S* GMSOBJ_txpcharlist; /* sy_class */
typedef struct GMSOBJ_txpcharlist_OD_S {  /* Objects of 'txpcharlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txpcharlist_CD */
  GMSOBJ_txlist GMSOBJ_txpcharlist_DOT_flist;
  SYSTEM_int64 GMSOBJ_txpcharlist_DOT_fstrmemory;
} GMSOBJ_txpcharlist_OD;


Function(SYSTEM_boolean ) GMSOBJ_txpcharlist_DOT_getonebased(
  GMSOBJ_txpcharlist self);

Procedure GMSOBJ_txpcharlist_DOT_setonebased(
  GMSOBJ_txpcharlist self,
  SYSTEM_boolean v);

Constructor(GMSOBJ_txpcharlist ) GMSOBJ_txpcharlist_DOT_create(
  GMSOBJ_txpcharlist self);

Destructor(GMSOBJ_txpcharlist ) GMSOBJ_txpcharlist_DOT_destroy(
  GMSOBJ_txpcharlist self);

Function(SYSTEM_integer ) GMSOBJ_txpcharlist_DOT_add(
  GMSOBJ_txpcharlist self,
  SYSTEM_P3_pansichar p,
  SYSTEM_integer l);

Procedure GMSOBJ_txpcharlist_DOT_getitem(
  GMSOBJ_txpcharlist self,
  SYSTEM_integer index,
  SYSTEM_P3_pansichar *p,
  SYSTEM_integer *l);

Procedure GMSOBJ_txpcharlist_DOT_clear(
  GMSOBJ_txpcharlist self);

Function(SYSTEM_integer ) GMSOBJ_txpcharlist_DOT_count(
  GMSOBJ_txpcharlist self);

Function(SYSTEM_int64 ) GMSOBJ_txpcharlist_DOT_memoryused(
  GMSOBJ_txpcharlist self);
extern void * const GMSOBJ_txpcharlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txpcharlist_CD;


typedef struct GMSOBJ_tstringitem_S {
  SYSTEM_P3_pshortstring fstring;
  SYSTEM_tobject fobject;
} GMSOBJ_tstringitem;

typedef GMSOBJ_tstringitem *GMSOBJ_pstringitem;
typedef GMSOBJ_tstringitem GMSOBJ_tstringitemlist[10000001];
typedef GMSOBJ_tstringitemlist *GMSOBJ_pstringitemlist;
typedef struct GMSOBJ_txcustomstringlist_OD_S* 
  GMSOBJ_txcustomstringlist; /* sy_class */
typedef struct GMSOBJ_txcustomstringlist_OD_S {  /* Objects of 'txcustomstringlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txcustomstringlist_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcount;
  GMSOBJ_pstringitemlist GMSOBJ_txcustomstringlist_DOT_flist;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcapacity;
  SYSTEM_int64 GMSOBJ_txcustomstringlist_DOT_fstrmemory;
} GMSOBJ_txcustomstringlist_OD;


Procedure GMSOBJ_txcustomstringlist_DOT_setname(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *v);

Procedure GMSOBJ_txcustomstringlist_DOT_setcapacity(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer newcapacity);

Function(SYSTEM_ansichar *) GMSOBJ_txcustomstringlist_DOT_getname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index);

Function(SYSTEM_tobject ) GMSOBJ_txcustomstringlist_DOT_getobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index);

Procedure GMSOBJ_txcustomstringlist_DOT_putobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index,
  SYSTEM_tobject aobject);

Procedure GMSOBJ_txcustomstringlist_DOT_insertitem(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer);

Prototype Procedure (*GMSOBJ_txcustomstringlist_DOT_grow_T)(
  GMSOBJ_txcustomstringlist self);

Procedure GMSOBJ_txcustomstringlist_DOT_grow(
  GMSOBJ_txcustomstringlist self);

Prototype Procedure (*GMSOBJ_txcustomstringlist_DOT_freeobject_T)(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index);

Procedure GMSOBJ_txcustomstringlist_DOT_freeobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index);

Constructor(GMSOBJ_txcustomstringlist ) 
  GMSOBJ_txcustomstringlist_DOT_create(
  GMSOBJ_txcustomstringlist self);

Destructor(GMSOBJ_txcustomstringlist ) 
  GMSOBJ_txcustomstringlist_DOT_destroy(
  GMSOBJ_txcustomstringlist self);

Procedure GMSOBJ_txcustomstringlist_DOT_delete(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index);

Procedure GMSOBJ_txcustomstringlist_DOT_freeitem(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index);

Procedure GMSOBJ_txcustomstringlist_DOT_clear(
  GMSOBJ_txcustomstringlist self);

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_add(
  GMSOBJ_txcustomstringlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_addobject(
  GMSOBJ_txcustomstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer);

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_indexof(
  GMSOBJ_txcustomstringlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_indexofobject(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_tobject aobject);

Procedure GMSOBJ_txcustomstringlist_DOT_exchange(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Function(SYSTEM_integer ) GMSOBJ_txcustomstringlist_DOT_compare(
  GMSOBJ_txcustomstringlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Function(SYSTEM_int64 ) GMSOBJ_txcustomstringlist_DOT_memoryused(
  GMSOBJ_txcustomstringlist self);
extern void * const GMSOBJ_txcustomstringlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txcustomstringlist_CD;


typedef struct GMSOBJ_txstringlist_OD_S* GMSOBJ_txstringlist; /* sy_class */
typedef struct GMSOBJ_txstringlist_OD_S {  /* Objects of 'txstringlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txstringlist_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcount;
  GMSOBJ_pstringitemlist GMSOBJ_txcustomstringlist_DOT_flist;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcapacity;
  SYSTEM_int64 GMSOBJ_txcustomstringlist_DOT_fstrmemory;
} GMSOBJ_txstringlist_OD;


Procedure GMSOBJ_txstringlist_DOT_insert(
  GMSOBJ_txstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *s);

Procedure GMSOBJ_txstringlist_DOT_insertobject(
  GMSOBJ_txstringlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer);
extern void * const GMSOBJ_txstringlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txstringlist_CD;


typedef struct GMSOBJ_txsortedstringlist_OD_S* 
  GMSOBJ_txsortedstringlist; /* sy_class */
typedef struct GMSOBJ_txsortedstringlist_OD_S {  /* Objects of 'txsortedstringlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txsortedstringlist_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcount;
  GMSOBJ_pstringitemlist GMSOBJ_txcustomstringlist_DOT_flist;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcapacity;
  SYSTEM_int64 GMSOBJ_txcustomstringlist_DOT_fstrmemory;
  SYSTEM_integer GMSOBJ_txsortedstringlist_DOT_fupdatecount;
  SYSTEM_boolean GMSOBJ_txsortedstringlist_DOT_fsorted;
} GMSOBJ_txsortedstringlist_OD;


Procedure GMSOBJ_txsortedstringlist_DOT_setsorted(
  GMSOBJ_txsortedstringlist self,
  SYSTEM_boolean value);

Constructor(GMSOBJ_txsortedstringlist ) 
  GMSOBJ_txsortedstringlist_DOT_create(
  GMSOBJ_txsortedstringlist self);

Function(SYSTEM_boolean ) GMSOBJ_txsortedstringlist_DOT_find(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer *index);

Function(SYSTEM_integer ) GMSOBJ_txsortedstringlist_DOT_add(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GMSOBJ_txsortedstringlist_DOT_addobject(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer);

Function(SYSTEM_integer ) GMSOBJ_txsortedstringlist_DOT_indexof(
  GMSOBJ_txsortedstringlist self,
  const SYSTEM_ansichar *s);

Procedure GMSOBJ_txsortedstringlist_DOT_beginupdate(
  GMSOBJ_txsortedstringlist self);

Procedure GMSOBJ_txsortedstringlist_DOT_endupdate(
  GMSOBJ_txsortedstringlist self);
extern void * const GMSOBJ_txsortedstringlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txsortedstringlist_CD;


typedef struct GMSOBJ_txstrstrlist_OD_S* GMSOBJ_txstrstrlist; /* sy_class */
typedef struct GMSOBJ_txstrstrlist_OD_S {  /* Objects of 'txstrstrlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txstrstrlist_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcount;
  GMSOBJ_pstringitemlist GMSOBJ_txcustomstringlist_DOT_flist;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcapacity;
  SYSTEM_int64 GMSOBJ_txcustomstringlist_DOT_fstrmemory;
  SYSTEM_integer GMSOBJ_txsortedstringlist_DOT_fupdatecount;
  SYSTEM_boolean GMSOBJ_txsortedstringlist_DOT_fsorted;
} GMSOBJ_txstrstrlist_OD;


Function(SYSTEM_ansichar *) GMSOBJ_txstrstrlist_DOT_getobject(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrstrlist self,
  SYSTEM_integer index);

Procedure GMSOBJ_txstrstrlist_DOT_putobject(
  GMSOBJ_txstrstrlist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *_ftmp1);

Procedure GMSOBJ_txstrstrlist_DOT_setasstring(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  const SYSTEM_ansichar *_ftmp2);

Procedure GMSOBJ_txstrstrlist_DOT_setasinteger(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_integer v);

Procedure GMSOBJ_txstrstrlist_DOT_setasdouble(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_double v);

Procedure GMSOBJ_txstrstrlist_DOT_setasboolean(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_boolean v);

Function(SYSTEM_ansichar *) GMSOBJ_txstrstrlist_DOT_getasstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1);

Function(SYSTEM_integer ) GMSOBJ_txstrstrlist_DOT_getasinteger(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1);

Function(SYSTEM_double ) GMSOBJ_txstrstrlist_DOT_getasdouble(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1);

Function(SYSTEM_boolean ) GMSOBJ_txstrstrlist_DOT_getasboolean(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *_ftmp1);

Procedure GMSOBJ_txstrstrlist_DOT_freeobject(
  GMSOBJ_txstrstrlist self,
  SYSTEM_integer index);

Function(SYSTEM_integer ) GMSOBJ_txstrstrlist_DOT_addobject(
  GMSOBJ_txstrstrlist self,
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);
extern void * const GMSOBJ_txstrstrlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txstrstrlist_CD;


cnstdef {GMSOBJ_non_empty = _P3char('=')};
typedef struct GMSOBJ_thashrecord_S *GMSOBJ_phashrecord;
typedef struct GMSOBJ_thashrecord_S {
  GMSOBJ_phashrecord pnext;
  SYSTEM_integer refnr;
} GMSOBJ_thashrecord;

typedef struct GMSOBJ_txhashedstringlist_OD_S* 
  GMSOBJ_txhashedstringlist; /* sy_class */
typedef struct GMSOBJ_txhashedstringlist_OD_S {  /* Objects of 'txhashedstringlist' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txhashedstringlist_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcount;
  GMSOBJ_pstringitemlist GMSOBJ_txcustomstringlist_DOT_flist;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcapacity;
  SYSTEM_int64 GMSOBJ_txcustomstringlist_DOT_fstrmemory;
  SYSTEM_P3_ppointerarray GMSOBJ_txhashedstringlist_DOT_phash;
  SYSTEM_integer GMSOBJ_txhashedstringlist_DOT_hashcount;
} GMSOBJ_txhashedstringlist_OD;


Prototype Function(SYSTEM_boolean ) (*
  GMSOBJ_txhashedstringlist_DOT_equaltoentry_T)(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer en);

Function(SYSTEM_boolean ) GMSOBJ_txhashedstringlist_DOT_equaltoentry(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer en);

Procedure GMSOBJ_txhashedstringlist_DOT_clearhashlist(
  GMSOBJ_txhashedstringlist self);

Procedure GMSOBJ_txhashedstringlist_DOT_sethashsize(
  GMSOBJ_txhashedstringlist self,
  SYSTEM_integer v);

Prototype Function(SYSTEM_cardinal ) (*
  GMSOBJ_txhashedstringlist_DOT_hashvalue_T)(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *name);

Function(SYSTEM_cardinal ) GMSOBJ_txhashedstringlist_DOT_hashvalue(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *name);

Constructor(GMSOBJ_txhashedstringlist ) 
  GMSOBJ_txhashedstringlist_DOT_create(
  GMSOBJ_txhashedstringlist self);

Destructor(GMSOBJ_txhashedstringlist ) 
  GMSOBJ_txhashedstringlist_DOT_destroy(
  GMSOBJ_txhashedstringlist self);

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_compare(
  GMSOBJ_txhashedstringlist self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);

Procedure GMSOBJ_txhashedstringlist_DOT_clear(
  GMSOBJ_txhashedstringlist self);

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_indexof(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_add(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GMSOBJ_txhashedstringlist_DOT_addobject(
  GMSOBJ_txhashedstringlist self,
  const SYSTEM_ansichar *s,
  SYSTEM_pointer apointer);

Procedure GMSOBJ_txhashedstringlist_DOT_sort(
  GMSOBJ_txhashedstringlist self);

Procedure GMSOBJ_txhashedstringlist_DOT_hashstats(
  GMSOBJ_txhashedstringlist self,
  SYSTEM_integer *amin,
  SYSTEM_integer *amax,
  SYSTEM_integer *aavg);
extern void * const GMSOBJ_txhashedstringlist_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txhashedstringlist_CD;


typedef struct GMSOBJ_txstrpool_OD_S* GMSOBJ_txstrpool; /* sy_class */
typedef struct GMSOBJ_txstrpool_OD_S {  /* Objects of 'txstrpool' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_txstrpool_CD */
  SYSTEM_boolean GMSOBJ_tquicksortclass_DOT_onebased;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcount;
  GMSOBJ_pstringitemlist GMSOBJ_txcustomstringlist_DOT_flist;
  SYSTEM_integer GMSOBJ_txcustomstringlist_DOT_fcapacity;
  SYSTEM_int64 GMSOBJ_txcustomstringlist_DOT_fstrmemory;
  SYSTEM_P3_ppointerarray GMSOBJ_txhashedstringlist_DOT_phash;
  SYSTEM_integer GMSOBJ_txhashedstringlist_DOT_hashcount;
} GMSOBJ_txstrpool_OD;


Function(SYSTEM_boolean ) GMSOBJ_txstrpool_DOT_equaltoentry(
  GMSOBJ_txstrpool self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer en);

Function(SYSTEM_integer ) GMSOBJ_txstrpool_DOT_compare(
  GMSOBJ_txstrpool self,
  SYSTEM_integer index1,
  SYSTEM_integer index2);
extern void * const GMSOBJ_txstrpool_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_txstrpool_CD;



Prototype Function(SYSTEM_boolean ) (*
  GMSOBJ_tbooleanbitarrayiterfunction)(
SYSTEM_integer n);

typedef struct GMSOBJ_tbooleanbitarray_OD_S* GMSOBJ_tbooleanbitarray; /* sy_class */
typedef struct GMSOBJ_tbooleanbitarray_OD_S {  /* Objects of 'tbooleanbitarray' */
  SYSTEM_classreference_t CD;  /* = &GMSOBJ_tbooleanbitarray_CD */
  GMSGEN_pbytedataarray GMSOBJ_tbooleanbitarray_DOT_pdata;
  SYSTEM_integer GMSOBJ_tbooleanbitarray_DOT_fallocated;
  SYSTEM_integer GMSOBJ_tbooleanbitarray_DOT_fhighindex;
} GMSOBJ_tbooleanbitarray_OD;


Procedure GMSOBJ_tbooleanbitarray_DOT_setbit(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer n,
  SYSTEM_boolean v);

Function(SYSTEM_boolean ) GMSOBJ_tbooleanbitarray_DOT_getbit(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer n);

Procedure GMSOBJ_tbooleanbitarray_DOT_sethighindex(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer v);

Procedure GMSOBJ_tbooleanbitarray_DOT_getbitmask(
  GMSOBJ_tbooleanbitarray self,
  SYSTEM_integer v,
  SYSTEM_integer *n,
  SYSTEM_byte *m);

Constructor(GMSOBJ_tbooleanbitarray ) 
  GMSOBJ_tbooleanbitarray_DOT_create(
  GMSOBJ_tbooleanbitarray self);

Destructor(GMSOBJ_tbooleanbitarray ) 
  GMSOBJ_tbooleanbitarray_DOT_destroy(
  GMSOBJ_tbooleanbitarray self);

Procedure GMSOBJ_tbooleanbitarray_DOT_clear(
  GMSOBJ_tbooleanbitarray self);

Procedure GMSOBJ_tbooleanbitarray_DOT_iterate(
  GMSOBJ_tbooleanbitarray self,
  GMSOBJ_tbooleanbitarrayiterfunction func);

Procedure GMSOBJ_tbooleanbitarray_DOT_iteratedown(
  GMSOBJ_tbooleanbitarray self,
  GMSOBJ_tbooleanbitarrayiterfunction func);

Function(SYSTEM_int64 ) GMSOBJ_tbooleanbitarray_DOT_memoryused(
  GMSOBJ_tbooleanbitarray self);
extern void * const GMSOBJ_tbooleanbitarray_VT[];
extern const SYSTEM_classdescriptor_t GMSOBJ_tbooleanbitarray_CD;



Procedure GMSOBJ_cmove(
  const SYSTEM_untyped *src,
  SYSTEM_untyped *dest,
  SYSTEM_integer len);

Function(SYSTEM_pointer ) GMSOBJ_copyint2ptr(
  SYSTEM_integer i);

Function(SYSTEM_integer ) GMSOBJ_copyptr2int(
  SYSTEM_pointer p);

extern void _Init_Module_gmsobj(void);
extern void _Final_Module_gmsobj(void);

#endif /* ! defined _P3___gmsobj___H */
