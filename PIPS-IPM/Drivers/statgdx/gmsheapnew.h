#ifndef _P3___gmsheapnew___H
#define _P3___gmsheapnew___H


Prototype Procedure (*GMSHEAPNEW_tmemoryreportproc)(
SYSTEM_double mem);

cnstdef {GMSHEAPNEW_bigblocksize = 524288};
cnstdef {GMSHEAPNEW_heapgranularity = 8};
cnstdef {GMSHEAPNEW_lastslot = 32};
typedef SYSTEM_uint8 GMSHEAPNEW_theapslotnr;
typedef struct GMSHEAPNEW_tsmallblock_S *GMSHEAPNEW_psmallblock;
typedef struct GMSHEAPNEW_tsmallblock_S {
  GMSHEAPNEW_psmallblock nextsmallblock;
} GMSHEAPNEW_tsmallblock;

typedef struct GMSHEAPNEW_tslotrecord_S {
  GMSHEAPNEW_psmallblock firstfree;
  SYSTEM_int64 getcount;
  SYSTEM_int64 freecount;
  SYSTEM_int64 listcount;
} GMSHEAPNEW_tslotrecord;

typedef struct GMSHEAPNEW_tlargeblock_S *GMSHEAPNEW_plargeblock;
typedef struct GMSHEAPNEW_tlargeblock_S {
  SYSTEM_integer freeslots;
  SYSTEM_P3_pbyte initialptr;
  SYSTEM_P3_pbyte currptr;
} GMSHEAPNEW_tlargeblock;

typedef struct GMSHEAPNEW_theapmgr_OD_S* GMSHEAPNEW_theapmgr; /* sy_class */
extern void * const GMSHEAPNEW_theapmgr_VT[];
extern const SYSTEM_classdescriptor_t GMSHEAPNEW_theapmgr_CD;


typedef struct GMSHEAPNEW_tbigblockmgr_OD_S* GMSHEAPNEW_tbigblockmgr; /* sy_class */
typedef struct GMSHEAPNEW_tbigblockmgr_OD_S {  /* Objects of 'tbigblockmgr' */
  SYSTEM_classreference_t CD;  /* = &GMSHEAPNEW_tbigblockmgr_CD */
  SYSTEM_P3_pshortstring GMSHEAPNEW_tbigblockmgr_DOT_spname;
  SYSTEM_int64 GMSHEAPNEW_tbigblockmgr_DOT_othermemory;
  SYSTEM_int64 GMSHEAPNEW_tbigblockmgr_DOT_highmark;
  GMSOBJ_txlist GMSHEAPNEW_tbigblockmgr_DOT_freelist;
  GMSOBJ_txlist GMSHEAPNEW_tbigblockmgr_DOT_mgrlist;
} GMSHEAPNEW_tbigblockmgr_OD;


Function(SYSTEM_pointer ) GMSHEAPNEW_tbigblockmgr_DOT_getbigblock(
  GMSHEAPNEW_tbigblockmgr self);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_releasebigblock(
  GMSHEAPNEW_tbigblockmgr self,
  SYSTEM_pointer p);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_reducememorysize(
  GMSHEAPNEW_tbigblockmgr self,
  SYSTEM_int64 delta);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_increasememorysize(
  GMSHEAPNEW_tbigblockmgr self,
  SYSTEM_int64 delta);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_registerheapmgr(
  GMSHEAPNEW_tbigblockmgr self,
  GMSHEAPNEW_theapmgr h);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_removeheapmgr(
  GMSHEAPNEW_tbigblockmgr self,
  GMSHEAPNEW_theapmgr h);

Function(GMSHEAPNEW_theapmgr ) GMSHEAPNEW_tbigblockmgr_DOT_getheapmgr(
  GMSHEAPNEW_tbigblockmgr self,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSHEAPNEW_tbigblockmgr_DOT_getname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSHEAPNEW_tbigblockmgr self);

Constructor(GMSHEAPNEW_tbigblockmgr ) 
  GMSHEAPNEW_tbigblockmgr_DOT_create(
  GMSHEAPNEW_tbigblockmgr self,
  const SYSTEM_ansichar *name);

Destructor(GMSHEAPNEW_tbigblockmgr ) 
  GMSHEAPNEW_tbigblockmgr_DOT_destroy(
  GMSHEAPNEW_tbigblockmgr self);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_clear(
  GMSHEAPNEW_tbigblockmgr self);

Function(SYSTEM_integer ) GMSHEAPNEW_tbigblockmgr_DOT_count(
  GMSHEAPNEW_tbigblockmgr self);

Procedure GMSHEAPNEW_tbigblockmgr_DOT_getbigstats(
  GMSHEAPNEW_tbigblockmgr self,
  SYSTEM_int64 *sizeothermemory,
  SYSTEM_int64 *sizehighmark,
  SYSTEM_int64 *cntfree);

Function(SYSTEM_int64 ) GMSHEAPNEW_tbigblockmgr_DOT_getfreeslotspace(
  GMSHEAPNEW_tbigblockmgr self);
extern void * const GMSHEAPNEW_tbigblockmgr_VT[];
extern const SYSTEM_classdescriptor_t GMSHEAPNEW_tbigblockmgr_CD;


typedef GMSHEAPNEW_tslotrecord _arr_0GMSHEAPNEW[32];
typedef struct GMSHEAPNEW_theapmgr_OD_S {  /* Objects of 'theapmgr' */
  SYSTEM_classreference_t CD;  /* = &GMSHEAPNEW_theapmgr_CD */
  GMSHEAPNEW_plargeblock GMSHEAPNEW_theapmgr_DOT_workbuffer;
  GMSHEAPNEW_tbigblockmgr GMSHEAPNEW_theapmgr_DOT_blockmgr;
  _arr_0GMSHEAPNEW GMSHEAPNEW_theapmgr_DOT_slots;
  SYSTEM_int64 GMSHEAPNEW_theapmgr_DOT_highmark,GMSHEAPNEW_theapmgr_DOT_othermemory,GMSHEAPNEW_theapmgr_DOT_otherget,GMSHEAPNEW_theapmgr_DOT_otherfree;
  SYSTEM_int64 GMSHEAPNEW_theapmgr_DOT_otherget64,GMSHEAPNEW_theapmgr_DOT_otherfree64;
  SYSTEM_int64 GMSHEAPNEW_theapmgr_DOT_realloccnt,GMSHEAPNEW_theapmgr_DOT_reallocused,GMSHEAPNEW_theapmgr_DOT_realloccnt64,GMSHEAPNEW_theapmgr_DOT_reallocused64;
  GMSOBJ_txlist GMSHEAPNEW_theapmgr_DOT_wrkbuffs;
  GMSOBJ_txlist GMSHEAPNEW_theapmgr_DOT_active;
  SYSTEM_P3_pshortstring GMSHEAPNEW_theapmgr_DOT_spname;
} GMSHEAPNEW_theapmgr_OD;


Function(GMSHEAPNEW_plargeblock ) 
  GMSHEAPNEW_theapmgr_DOT_getworkbuffer(
  GMSHEAPNEW_theapmgr self);

Procedure GMSHEAPNEW_theapmgr_DOT_releaseworkbuffer(
  GMSHEAPNEW_theapmgr self,
  GMSHEAPNEW_plargeblock p);

Procedure GMSHEAPNEW_theapmgr_DOT_reducememorysize(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_int64 delta);

Procedure GMSHEAPNEW_theapmgr_DOT_increasememorysize(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_int64 delta);

Function(SYSTEM_ansichar *) GMSHEAPNEW_theapmgr_DOT_getname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSHEAPNEW_theapmgr self);

Constructor(GMSHEAPNEW_theapmgr ) GMSHEAPNEW_theapmgr_DOT_create(
  GMSHEAPNEW_theapmgr self,
  GMSHEAPNEW_tbigblockmgr m,
  const SYSTEM_ansichar *name);

Destructor(GMSHEAPNEW_theapmgr ) GMSHEAPNEW_theapmgr_DOT_destroy(
  GMSHEAPNEW_theapmgr self);

Procedure GMSHEAPNEW_theapmgr_DOT_clear(
  GMSHEAPNEW_theapmgr self);

Function(SYSTEM_pointer ) GMSHEAPNEW_theapmgr_DOT_gmsgetmem(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_word slot);

Procedure GMSHEAPNEW_theapmgr_DOT_gmsfreemem(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer p,
  SYSTEM_word slot);

Function(SYSTEM_pointer ) GMSHEAPNEW_theapmgr_DOT_xgetmem(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_integer size);

Function(SYSTEM_pointer ) GMSHEAPNEW_theapmgr_DOT_xgetmemnc(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_integer size);

Function(SYSTEM_pointer ) GMSHEAPNEW_theapmgr_DOT_xallocmem(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_integer size);

Function(SYSTEM_pointer ) GMSHEAPNEW_theapmgr_DOT_xallocmemnc(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_integer size);

Function(SYSTEM_pointer ) GMSHEAPNEW_theapmgr_DOT_xgetmem64(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_int64 size);

Procedure GMSHEAPNEW_theapmgr_DOT_xfreemem(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer p,
  SYSTEM_integer size);

Procedure GMSHEAPNEW_theapmgr_DOT_xfreememnc(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer p,
  SYSTEM_integer size);

Procedure GMSHEAPNEW_theapmgr_DOT_xfreememandnil(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer *p,
  SYSTEM_integer size);

Procedure GMSHEAPNEW_theapmgr_DOT_xfreemem64(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer p,
  SYSTEM_int64 size);

Procedure GMSHEAPNEW_theapmgr_DOT_xfreemem64andnil(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer *p,
  SYSTEM_int64 size);

Procedure GMSHEAPNEW_theapmgr_DOT_xreallocmem(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer *p,
  SYSTEM_integer oldsize,
  SYSTEM_integer newsize);

Procedure GMSHEAPNEW_theapmgr_DOT_xreallocmemnc(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer *p,
  SYSTEM_integer oldsize,
  SYSTEM_integer newsize);

Procedure GMSHEAPNEW_theapmgr_DOT_xreallocmem64(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_pointer *p,
  SYSTEM_int64 oldsize,
  SYSTEM_int64 newsize);

Procedure GMSHEAPNEW_theapmgr_DOT_getslotcnts(
  GMSHEAPNEW_theapmgr self,
  GMSHEAPNEW_theapslotnr slot,
  SYSTEM_int64 *cntget,
  SYSTEM_int64 *cntfree,
  SYSTEM_int64 *cntavail);

Procedure GMSHEAPNEW_theapmgr_DOT_getblockstats(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_int64 *cntwrkbuffs,
  SYSTEM_int64 *cntactive,
  SYSTEM_int64 *sizeothermemory,
  SYSTEM_int64 *sizehighmark);

Procedure GMSHEAPNEW_theapmgr_DOT_getotherstats(
  GMSHEAPNEW_theapmgr self,
  SYSTEM_boolean do64,
  SYSTEM_int64 *cntget,
  SYSTEM_int64 *cntfree,
  SYSTEM_int64 *cntrealloc,
  SYSTEM_int64 *sizerused);

Function(SYSTEM_int64 ) GMSHEAPNEW_theapmgr_DOT_getfreeslotspace(
  GMSHEAPNEW_theapmgr self);
extern void * const GMSHEAPNEW_theapmgr_VT[];
extern const SYSTEM_classdescriptor_t GMSHEAPNEW_theapmgr_CD;



Procedure GMSHEAPNEW_gmscreatedefaultheap(void);

Procedure GMSHEAPNEW_gmsreleasedefaultheap(void);

Function(SYSTEM_double ) GMSHEAPNEW_gmsmemoryused(void);

Function(SYSTEM_double ) GMSHEAPNEW_gmsmemoryfree(void);

Function(SYSTEM_boolean ) GMSHEAPNEW_setmemorylimit(
  SYSTEM_double limit);

Function(SYSTEM_double ) GMSHEAPNEW_getmemorylimit(void);

Procedure GMSHEAPNEW_setmemoryreportproc(
  GMSHEAPNEW_tmemoryreportproc f);
extern GMSHEAPNEW_tbigblockmgr GMSHEAPNEW_bbmgr;
extern GMSHEAPNEW_theapmgr GMSHEAPNEW_gheap;

extern void _Init_Module_gmsheapnew(void);
extern void _Final_Module_gmsheapnew(void);

#endif /* ! defined _P3___gmsheapnew___H */
