#ifndef _P3___gxfile___H
#define _P3___gxfile___H

extern _P3STR_7 GXFILE_baduel_prefix;
extern _P3STR_7 GXFILE_badstr_prefix;
extern _P3STR_15 GXFILE_strgdxcompress;
extern _P3STR_15 GXFILE_strgdxconvert;
typedef GXDEFS_tgdxuelindex *GXFILE_puelindex;
typedef struct GXFILE_tdfilter_OD_S* GXFILE_tdfilter; /* sy_class */
typedef struct GXFILE_tdfilter_OD_S {  /* Objects of 'tdfilter' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tdfilter_CD */
  SYSTEM_integer GXFILE_tdfilter_DOT_filtnumber;
  SYSTEM_integer GXFILE_tdfilter_DOT_filtmaxuel;
  GMSOBJ_tbooleanbitarray GXFILE_tdfilter_DOT_filtmap;
  SYSTEM_boolean GXFILE_tdfilter_DOT_filtsorted;
} GXFILE_tdfilter_OD;


Constructor(GXFILE_tdfilter ) GXFILE_tdfilter_DOT_create(
  GXFILE_tdfilter self,
  SYSTEM_integer nr,
  SYSTEM_integer userhigh);

Destructor(GXFILE_tdfilter ) GXFILE_tdfilter_DOT_destroy(
  GXFILE_tdfilter self);

Function(SYSTEM_boolean ) GXFILE_tdfilter_DOT_infilter(
  GXFILE_tdfilter self,
  SYSTEM_integer v);

Function(SYSTEM_int64 ) GXFILE_tdfilter_DOT_memoryused(
  GXFILE_tdfilter self);
extern void * const GXFILE_tdfilter_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tdfilter_CD;


typedef struct GXFILE_tfilterlist_OD_S* GXFILE_tfilterlist; /* sy_class */
typedef struct GXFILE_tfilterlist_OD_S {  /* Objects of 'tfilterlist' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tfilterlist_CD */
  GMSOBJ_txlist GXFILE_tfilterlist_DOT_flist;
} GXFILE_tfilterlist_OD;


Constructor(GXFILE_tfilterlist ) GXFILE_tfilterlist_DOT_create(
  GXFILE_tfilterlist self);

Destructor(GXFILE_tfilterlist ) GXFILE_tfilterlist_DOT_destroy(
  GXFILE_tfilterlist self);

Procedure GXFILE_tfilterlist_DOT_addfilter(
  GXFILE_tfilterlist self,
  GXFILE_tdfilter f);

Procedure GXFILE_tfilterlist_DOT_deletefilter(
  GXFILE_tfilterlist self,
  SYSTEM_integer ix);

Function(GXFILE_tdfilter ) GXFILE_tfilterlist_DOT_findfilter(
  GXFILE_tfilterlist self,
  SYSTEM_integer nr);

Function(SYSTEM_int64 ) GXFILE_tfilterlist_DOT_memoryused(
  GXFILE_tfilterlist self);
extern void * const GXFILE_tfilterlist_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tfilterlist_CD;


typedef SYSTEM_byte GXFILE_tgdxdaction; /* Anonymous */ enum{GXFILE_dm_unmapped,GXFILE_dm_strict,GXFILE_dm_filter,
  GXFILE_dm_expand};
typedef struct GXFILE_tdomain_S {
  GXFILE_tdfilter dfilter;
  GXFILE_tgdxdaction daction;
} GXFILE_tdomain;

typedef GXFILE_tdomain GXFILE_tdomainlist[20];
typedef struct GXFILE_tgdxsymbrecord_S *GXFILE_pgdxsymbrecord;
typedef struct GXFILE_tgdxsymbrecord_S {
  SYSTEM_integer ssynr;
  SYSTEM_int64 sposition;
  SYSTEM_integer sdim;
  SYSTEM_integer sdatacount;
  SYSTEM_integer serrors;
  GMSSPECS_tgdxdatatype sdatatype;
  SYSTEM_integer suserinfo;
  SYSTEM_boolean ssettext;
  SYSTEM_shortstring sexpltxt;
  SYSTEM_boolean siscompressed;
  SYSTEM_P3_pintegerarray sdomsymbols;
  SYSTEM_P3_pintegerarray sdomstrings;
  GMSOBJ_txstrings scommentslist;
  SYSTEM_boolean sscalarfrst;
  GMSOBJ_tbooleanbitarray ssetbitmap;
} GXFILE_tgdxsymbrecord;

typedef SYSTEM_byte GXFILE_tgdxintvaltyp; /* Anonymous */ enum{GXFILE_vm_valund,GXFILE_vm_valna,GXFILE_vm_valpin,
  GXFILE_vm_valmin,GXFILE_vm_valeps,GXFILE_vm_zero,GXFILE_vm_one,
  GXFILE_vm_mone,GXFILE_vm_half,GXFILE_vm_two,GXFILE_vm_normal};
typedef SYSTEM_byte GXFILE_tgxfilemode; /* Anonymous */ enum{GXFILE_f_not_open,GXFILE_fr_init,GXFILE_fw_init,
  GXFILE_fw_dom_raw,GXFILE_fw_dom_map,GXFILE_fw_dom_str,
  GXFILE_fw_raw_data,GXFILE_fw_map_data,GXFILE_fw_str_data,
  GXFILE_f_raw_elem,GXFILE_f_map_elem,GXFILE_f_str_elem,
  GXFILE_fr_raw_data,GXFILE_fr_map_data,GXFILE_fr_mapr_data,
  GXFILE_fr_str_data,GXFILE_fr_filter,GXFILE_fr_slice};
typedef _P3SET_31 GXFILE_tgxmodeset;
extern GXFILE_tgxmodeset GXFILE_anywritemode;
extern GXFILE_tgxmodeset GXFILE_anyreadmode;
typedef SYSTEM_byte GXFILE_tgdxelemsize; /* Anonymous */ enum{GXFILE_sz_byte,GXFILE_sz_word,GXFILE_sz_integer};
typedef struct GXFILE_tintegermapping_OD_S* GXFILE_tintegermapping; /* sy_class */
typedef struct GXFILE_tintegermapping_OD_S {  /* Objects of 'tintegermapping' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tintegermapping_CD */
  SYSTEM_integer GXFILE_tintegermapping_DOT_fhighestindex;
  SYSTEM_integer GXFILE_tintegermapping_DOT_fcapacity;
  SYSTEM_P3_pintegerarray GXFILE_tintegermapping_DOT_pmap;
} GXFILE_tintegermapping_OD;


Procedure GXFILE_tintegermapping_DOT_setmapping(
  GXFILE_tintegermapping self,
  SYSTEM_integer f,
  SYSTEM_integer t);

Function(SYSTEM_integer ) GXFILE_tintegermapping_DOT_getmapping(
  GXFILE_tintegermapping self,
  SYSTEM_integer f);

Constructor(GXFILE_tintegermapping ) GXFILE_tintegermapping_DOT_create(
  GXFILE_tintegermapping self);

Destructor(GXFILE_tintegermapping ) GXFILE_tintegermapping_DOT_destroy(
  GXFILE_tintegermapping self);

Function(SYSTEM_int64 ) GXFILE_tintegermapping_DOT_memoryused(
  GXFILE_tintegermapping self);
extern void * const GXFILE_tintegermapping_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tintegermapping_CD;


typedef SYSTEM_byte GXFILE_tuelusermapstatus; /* Anonymous */ enum{GXFILE_map_unknown,GXFILE_map_unsorted,GXFILE_map_sorted,
  GXFILE_map_sortgrow,GXFILE_map_sortfull};
typedef struct GXFILE_tueltable_OD_S* GXFILE_tueltable; /* sy_class */
typedef struct GXFILE_tueltable_OD_S {  /* Objects of 'tueltable' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tueltable_CD */
  GMSDATA_tgrowarrayfxd STRHASH_txstrhashlist_DOT_buckets;
  SYSTEM_P3_ppointerarray STRHASH_txstrhashlist_DOT_phashtable;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_hashsize;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_rehashcnt;
  GMSDATA_txintlist STRHASH_txstrhashlist_DOT_sortmap;
  SYSTEM_boolean STRHASH_txstrhashlist_DOT_fsorted;
  SYSTEM_integer STRHASH_txstrhashlist_DOT_fcount;
  SYSTEM_boolean STRHASH_txstrhashlist_DOT_onebased;
  GXFILE_tintegermapping GXFILE_tueltable_DOT_usruel2ent;
  GXFILE_tuelusermapstatus GXFILE_tueltable_DOT_fmaptouserstatus;
} GXFILE_tueltable_OD;


Function(GXFILE_tuelusermapstatus ) 
  GXFILE_tueltable_DOT_getmaptouserstatus(
  GXFILE_tueltable self);

Procedure GXFILE_tueltable_DOT_resetmaptouserstatus(
  GXFILE_tueltable self);

Constructor(GXFILE_tueltable ) GXFILE_tueltable_DOT_create(
  GXFILE_tueltable self);

Destructor(GXFILE_tueltable ) GXFILE_tueltable_DOT_destroy(
  GXFILE_tueltable self);

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_addusrnew(
  GXFILE_tueltable self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_addusrindxnew(
  GXFILE_tueltable self,
  const SYSTEM_ansichar *s,
  SYSTEM_integer uelnr);

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_newusruel(
  GXFILE_tueltable self,
  SYSTEM_integer en);

Procedure GXFILE_tueltable_DOT_loadfromstream(
  GXFILE_tueltable self,
  GMSSTRM_txstream s);

Function(SYSTEM_integer ) GXFILE_tueltable_DOT_getusermap(
  GXFILE_tueltable self,
  SYSTEM_integer n);

Function(SYSTEM_int64 ) GXFILE_tueltable_DOT_memoryused(
  GXFILE_tueltable self);
extern void * const GXFILE_tueltable_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tueltable_CD;


typedef struct GXFILE_tacronym_OD_S* GXFILE_tacronym; /* sy_class */
typedef struct GXFILE_tacronym_OD_S {  /* Objects of 'tacronym' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tacronym_CD */
  SYSTEM_P3_pshortstring GXFILE_tacronym_DOT_acrname;
  SYSTEM_P3_pshortstring GXFILE_tacronym_DOT_acrtext;
  SYSTEM_integer GXFILE_tacronym_DOT_acrmap;
  SYSTEM_integer GXFILE_tacronym_DOT_acrreadmap;
  SYSTEM_boolean GXFILE_tacronym_DOT_acrautogen;
  SYSTEM_int64 GXFILE_tacronym_DOT_fstrmemory;
} GXFILE_tacronym_OD;


Constructor(GXFILE_tacronym ) GXFILE_tacronym_DOT_create(
  GXFILE_tacronym self,
  const SYSTEM_ansichar *name,
  const SYSTEM_ansichar *text,
  SYSTEM_integer map);

Constructor(GXFILE_tacronym ) GXFILE_tacronym_DOT_createfromstream(
  GXFILE_tacronym self,
  GMSSTRM_txstream s);

Destructor(GXFILE_tacronym ) GXFILE_tacronym_DOT_destroy(
  GXFILE_tacronym self);

Procedure GXFILE_tacronym_DOT_savetostream(
  GXFILE_tacronym self,
  GMSSTRM_txstream s);

Function(SYSTEM_int64 ) GXFILE_tacronym_DOT_memoryused(
  GXFILE_tacronym self);
extern void * const GXFILE_tacronym_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tacronym_CD;


typedef struct GXFILE_tacronymlist_OD_S* GXFILE_tacronymlist; /* sy_class */
typedef struct GXFILE_tacronymlist_OD_S {  /* Objects of 'tacronymlist' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tacronymlist_CD */
  GMSOBJ_txlist GXFILE_tacronymlist_DOT_flist;
} GXFILE_tacronymlist_OD;


Constructor(GXFILE_tacronymlist ) GXFILE_tacronymlist_DOT_create(
  GXFILE_tacronymlist self);

Destructor(GXFILE_tacronymlist ) GXFILE_tacronymlist_DOT_destroy(
  GXFILE_tacronymlist self);

Function(SYSTEM_integer ) GXFILE_tacronymlist_DOT_findentry(
  GXFILE_tacronymlist self,
  SYSTEM_integer map);

Function(SYSTEM_integer ) GXFILE_tacronymlist_DOT_findname(
  GXFILE_tacronymlist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GXFILE_tacronymlist_DOT_addentry(
  GXFILE_tacronymlist self,
  const SYSTEM_ansichar *name,
  const SYSTEM_ansichar *text,
  SYSTEM_integer map);

Procedure GXFILE_tacronymlist_DOT_checkentry(
  GXFILE_tacronymlist self,
  SYSTEM_integer map);

Procedure GXFILE_tacronymlist_DOT_savetostream(
  GXFILE_tacronymlist self,
  GMSSTRM_txstream s);

Procedure GXFILE_tacronymlist_DOT_loadfromstream(
  GXFILE_tacronymlist self,
  GMSSTRM_txstream s);

Function(SYSTEM_int64 ) GXFILE_tacronymlist_DOT_memoryused(
  GXFILE_tacronymlist self);
extern void * const GXFILE_tacronymlist_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tacronymlist_CD;


typedef SYSTEM_byte _enm_0GXFILE; /* Anonymous */ enum{GXFILE_stat_notopen,GXFILE_stat_read,GXFILE_stat_write};
typedef SYSTEM_double _arr_1GXFILE[11];
typedef SYSTEM_byte _enm_2GXFILE; /* Anonymous */ enum{GXFILE_trl_none,GXFILE_trl_errors,GXFILE_trl_some,GXFILE_trl_all};
typedef GXFILE_tgdxelemsize _arr_3GXFILE[20];
typedef GXFILE_tintegermapping _arr_4GXFILE[20];
typedef GXFILE_tintegermapping _arr_5GXFILE[20];
typedef GMSOBJ_tbooleanbitarray _arr_6GXFILE[20];
typedef struct GXFILE_tgxfileobj_OD_S* GXFILE_tgxfileobj; /* sy_class */
typedef struct GXFILE_tgxfileobj_OD_S {  /* Objects of 'tgxfileobj' */
  SYSTEM_classreference_t CD;  /* = &GXFILE_tgxfileobj_CD */
  GMSSTRM_tmibufferedstream GXFILE_tgxfileobj_DOT_ffile;
  GXFILE_tgxfilemode GXFILE_tgxfileobj_DOT_fmode;
  GXFILE_tgxfilemode GXFILE_tgxfileobj_DOT_fmode_aftreg;
  _enm_0GXFILE GXFILE_tgxfileobj_DOT_fstatus;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_fcomprlev;
  GXFILE_tueltable GXFILE_tgxfileobj_DOT_ueltable;
  GMSOBJ_txstrpool GXFILE_tgxfileobj_DOT_settextlist;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_fcurrentdim;
  GMSSPECS_tindex GXFILE_tgxfileobj_DOT_lastelem;
  GMSSPECS_tindex GXFILE_tgxfileobj_DOT_prevelem;
  GMSSPECS_tindex GXFILE_tgxfileobj_DOT_minelem;
  GMSSPECS_tindex GXFILE_tgxfileobj_DOT_maxelem;
  GXDEFS_tgdxstrindex GXFILE_tgxfileobj_DOT_laststrelem;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_datasize;
  GMSSPECS_tvarvaltype GXFILE_tgxfileobj_DOT_lastdatafield;
  STRHASH_txstrhashlist GXFILE_tgxfileobj_DOT_namelist;
  STRHASH_txstrhashlist GXFILE_tgxfileobj_DOT_domainstrlist;
  DATASTORAGE_tlinkeddata GXFILE_tgxfileobj_DOT_sortlist;
  GMSDATA_ttblgamsdata GXFILE_tgxfileobj_DOT_errorlist;
  GXFILE_pgdxsymbrecord GXFILE_tgxfileobj_DOT_cursyptr;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_errcnt;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_errcnttotal;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_lasterror;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_lastreperror;
  GXFILE_tfilterlist GXFILE_tgxfileobj_DOT_filterlist;
  GXFILE_tdfilter GXFILE_tgxfileobj_DOT_curfilter;
  GXFILE_tdomainlist GXFILE_tgxfileobj_DOT_domainlist;
  _arr_1GXFILE GXFILE_tgxfileobj_DOT_intvaluemap,GXFILE_tgxfileobj_DOT_readintvaluemap;
  _enm_2GXFILE GXFILE_tgxfileobj_DOT_tracelevel;
  SYSTEM_shortstring GXFILE_tgxfileobj_DOT_tracestr;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_versionread;
  SYSTEM_shortstring GXFILE_tgxfileobj_DOT_fproducer;
  SYSTEM_shortstring GXFILE_tgxfileobj_DOT_fproducer2;
  SYSTEM_shortstring GXFILE_tgxfileobj_DOT_filesystemid;
  SYSTEM_int64 GXFILE_tgxfileobj_DOT_majorindexposition;
  SYSTEM_int64 GXFILE_tgxfileobj_DOT_nextwriteposition;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_datacount;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_nrmappedadded;
  _arr_3GXFILE GXFILE_tgxfileobj_DOT_elemtype;
  _P3STR_63 GXFILE_tgxfileobj_DOT_majcontext;
  _arr_4GXFILE GXFILE_tgxfileobj_DOT_sliceindxs;
  _arr_5GXFILE GXFILE_tgxfileobj_DOT_slicerevmap;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_slicesynr;
  GXDEFS_tgdxstrindex GXFILE_tgxfileobj_DOT_sliceelems;
  SYSTEM_pointer GXFILE_tgxfileobj_DOT_readptr;
  SYSTEM_boolean GXFILE_tgxfileobj_DOT_douncompress;
  SYSTEM_boolean GXFILE_tgxfileobj_DOT_compressout;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_deltaforwrite;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_deltaforread;
  SYSTEM_double GXFILE_tgxfileobj_DOT_zvalacr;
  GXFILE_tacronymlist GXFILE_tgxfileobj_DOT_acronymlist;
  _arr_6GXFILE GXFILE_tgxfileobj_DOT_wrbitmaps;
  SYSTEM_boolean GXFILE_tgxfileobj_DOT_readuniverse;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_universenr;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_uelcntorig;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_autoconvert;
  SYSTEM_integer GXFILE_tgxfileobj_DOT_nextautoacronym;
  SYSTEM_boolean GXFILE_tgxfileobj_DOT_appendactive;
  GXDEFS_tdomainindexproc GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp;
  GXDEFS_tdatastorefiltproc GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp;
  SYSTEM_boolean GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_callbyref;
  SYSTEM_boolean GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_callbyref;
} GXFILE_tgxfileobj_OD;


Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_preparesymbolwrite(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *caller,
  const SYSTEM_ansichar *aname,
  const SYSTEM_ansichar *atext,
  SYSTEM_integer adim,
  SYSTEM_integer atype,
  SYSTEM_integer auserinfo);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_preparesymbolread(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *caller,
  SYSTEM_integer synr,
  const SYSTEM_integer *adomainnrs,
  GXFILE_tgxfilemode newmode);

Procedure GXFILE_tgxfileobj_DOT_initerrors(
  GXFILE_tgxfileobj self);

Procedure GXFILE_tgxfileobj_DOT_seterror(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n);

Procedure GXFILE_tgxfileobj_DOT_reporterror(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_errorcondition(
  GXFILE_tgxfileobj self,
  SYSTEM_boolean cnd,
  SYSTEM_integer n);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_majorcheckmode(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *routine,
  const _P3set_elem *ms);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_checkmode(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *routine,
  const _P3set_elem *ms);

Procedure GXFILE_tgxfileobj_DOT_writetrace(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *s);

Procedure GXFILE_tgxfileobj_DOT_initdowrite(
  GXFILE_tgxfileobj self,
  SYSTEM_integer nrrecs);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_dowrite(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *aelements,
  const SYSTEM_double *avals);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_doread(
  GXFILE_tgxfileobj self,
  SYSTEM_double *avals,
  SYSTEM_integer *afdim);

Procedure GXFILE_tgxfileobj_DOT_addtoerrorlist(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *aelements,
  const SYSTEM_double *avals);

Procedure GXFILE_tgxfileobj_DOT_getdefaultrecord(
  GXFILE_tgxfileobj self,
  SYSTEM_double *avals);

Function(SYSTEM_double ) GXFILE_tgxfileobj_DOT_acronymremap(
  GXFILE_tgxfileobj self,
  SYSTEM_double v);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_isgoodnewsymbol(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_boolean ) GXFILE_tgxfileobj_DOT_resultwillbesorted(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *adomainnrs);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenreadxx(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *afn,
  SYSTEM_integer filemode,
  SYSTEM_integer *errnr);

Procedure  STDCALL GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_fc(
  GXFILE_tgxfileobj self,
  SYSTEM_integer rawindex,
  SYSTEM_integer mappedindex,
  SYSTEM_pointer uptr);

Function(SYSTEM_integer )  STDCALL 
  GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_fc(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *indx,
  const SYSTEM_double *vals,
  SYSTEM_pointer uptr);

Constructor(GXFILE_tgxfileobj ) GXFILE_tgxfileobj_DOT_create(
  GXFILE_tgxfileobj self,
  SYSTEM_ansichar *errmsg);

Destructor(GXFILE_tgxfileobj ) GXFILE_tgxfileobj_DOT_destroy(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenwrite(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *producer,
  SYSTEM_integer *errnr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenwriteex(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *producer,
  SYSTEM_integer compr,
  SYSTEM_integer *errnr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenread(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  SYSTEM_integer *errnr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxopenappend(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *producer,
  SYSTEM_integer *errnr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsysteminfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *sycnt,
  SYSTEM_integer *uelcnt);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymboldim(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_ansichar *syid,
  SYSTEM_integer *dimen,
  SYSTEM_integer *typ);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolinfox(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *reccnt,
  SYSTEM_integer *userinfo,
  SYSTEM_ansichar *expltxt);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfindsymbol(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  SYSTEM_integer *synr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetuel(
  GXFILE_tgxfileobj self,
  SYSTEM_integer uelnr,
  SYSTEM_ansichar *uel);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetlasterror(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxerrorcount(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterrawstart(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregistermapstart(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterstrstart(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterraw(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *uel);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregistermap(
  GXFILE_tgxfileobj self,
  SYSTEM_integer umap,
  const SYSTEM_ansichar *uel);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterstr(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *uel,
  SYSTEM_integer *uelnr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelregisterdone(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawriterawstart(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  const SYSTEM_ansichar *expltxt,
  SYSTEM_integer dimen,
  SYSTEM_integer typ,
  SYSTEM_integer userinfo);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritemapstart(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  const SYSTEM_ansichar *expltxt,
  SYSTEM_integer dimen,
  SYSTEM_integer typ,
  SYSTEM_integer userinfo);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritestrstart(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *syid,
  const SYSTEM_ansichar *expltxt,
  SYSTEM_integer dimen,
  SYSTEM_integer typ,
  SYSTEM_integer userinfo);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawriteraw(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *keyint,
  const SYSTEM_double *values);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritemap(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *keyint,
  const SYSTEM_double *values);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritestr(
  GXFILE_tgxfileobj self,
  const SYSTEM_shortstring *keystr,
  const SYSTEM_double *values);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxaddsettext(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *txt,
  SYSTEM_integer *txtnr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetelemtext(
  GXFILE_tgxfileobj self,
  SYSTEM_integer txtnr,
  SYSTEM_ansichar *txt,
  SYSTEM_integer *node);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsettextnodenr(
  GXFILE_tgxfileobj self,
  SYSTEM_integer txtnr,
  SYSTEM_integer node);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatawritedone(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdataerrorcount(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdataerrorrecord(
  GXFILE_tgxfileobj self,
  SYSTEM_integer recnr,
  SYSTEM_integer *keyint,
  SYSTEM_double *values);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadrawstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *nrrecs);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadmapstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *nrrecs);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadstrstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *nrrecs);

Function(SYSTEM_integer ) 
  GXFILE_tgxfileobj_DOT_gdxdatareadfilteredstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_integer *filteraction,
  SYSTEM_integer *nrrecs);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadraw(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *keyint,
  SYSTEM_double *values,
  SYSTEM_integer *dimfrst);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadmap(
  GXFILE_tgxfileobj self,
  SYSTEM_integer recnr,
  SYSTEM_integer *keyint,
  SYSTEM_double *values,
  SYSTEM_integer *dimfrst);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadstr(
  GXFILE_tgxfileobj self,
  SYSTEM_shortstring *keystr,
  SYSTEM_double *values,
  SYSTEM_integer *dimfrst);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadslicestart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *elemcounts);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadslice(
  GXFILE_tgxfileobj self,
  const SYSTEM_shortstring *uelfilterstr,
  SYSTEM_integer *dimen,
  GXDEFS_tdatastoreproc dp);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatasliceuels(
  GXFILE_tgxfileobj self,
  const SYSTEM_integer *slicekeyint,
  SYSTEM_shortstring *keystr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareaddone(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterregisterstart(
  GXFILE_tgxfileobj self,
  SYSTEM_integer filternr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterregister(
  GXFILE_tgxfileobj self,
  SYSTEM_integer uelmap);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterregisterdone(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfilterexists(
  GXFILE_tgxfileobj self,
  SYSTEM_integer filternr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxumuelinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *uelcnt,
  SYSTEM_integer *highmap);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxumuelget(
  GXFILE_tgxfileobj self,
  SYSTEM_integer uelnr,
  SYSTEM_ansichar *uel,
  SYSTEM_integer *uelmap);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxumfinduel(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *uel,
  SYSTEM_integer *uelnr,
  SYSTEM_integer *uelmap);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsethastext(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxresetspecialvalues(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsetspecialvalues(
  GXFILE_tgxfileobj self,
  const SYSTEM_double *avals);

Function(SYSTEM_integer ) 
  GXFILE_tgxfileobj_DOT_gdxsetreadspecialvalues(
  GXFILE_tgxfileobj self,
  const SYSTEM_double *avals);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetspecialvalues(
  GXFILE_tgxfileobj self,
  SYSTEM_double *avals);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxmapvalue(
  GXFILE_tgxfileobj self,
  SYSTEM_double d,
  SYSTEM_integer *sv);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsettracelevel(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfileversion(
  GXFILE_tgxfileobj self,
  SYSTEM_ansichar *filestr,
  SYSTEM_ansichar *producestr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxuelmaxlength(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbmaxlength(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbindxmaxlength(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *lengthinfo);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymindex(
  GXFILE_tgxfileobj self,
  SYSTEM_double v);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymnextnr(
  GXFILE_tgxfileobj self,
  SYSTEM_integer nv);

Function(SYSTEM_double ) GXFILE_tgxfileobj_DOT_gdxacronymvalue(
  GXFILE_tgxfileobj self,
  SYSTEM_integer aindx);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymcount(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymname(
  GXFILE_tgxfileobj self,
  SYSTEM_double v,
  SYSTEM_ansichar *aname);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymgetinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  SYSTEM_ansichar *aname,
  SYSTEM_ansichar *txt,
  SYSTEM_integer *aindx);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymadd(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *aname,
  const SYSTEM_ansichar *txt,
  SYSTEM_integer aindx);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymsetinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  const SYSTEM_ansichar *aname,
  const SYSTEM_ansichar *txt,
  SYSTEM_integer aindx);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxacronymgetmapping(
  GXFILE_tgxfileobj self,
  SYSTEM_integer n,
  SYSTEM_integer *orgindx,
  SYSTEM_integer *newindx,
  SYSTEM_integer *autoindex);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymboladdcomment(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_ansichar *txt);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolgetcomment(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer n,
  SYSTEM_ansichar *txt);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolsetdomain(
  GXFILE_tgxfileobj self,
  const SYSTEM_shortstring *domainids);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolsetdomainx(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_shortstring *domainids);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolgetdomain(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer *domainsynrs);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxsymbolgetdomainx(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_shortstring *domainids);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxaddalias(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *id1,
  const SYSTEM_ansichar *id2);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxerrorstr(
  GXFILE_tgxfileobj self,
  SYSTEM_integer errnr,
  SYSTEM_ansichar *errmsg);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetdllversion(
  GXFILE_tgxfileobj self,
  SYSTEM_ansichar *v);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxfileinfo(
  GXFILE_tgxfileobj self,
  SYSTEM_integer *filever,
  SYSTEM_integer *comprlev);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxclose(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxcurrentdim(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxautoconvert(
  GXFILE_tgxfileobj self,
  SYSTEM_integer nv);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadrawfast(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  GXDEFS_tdatastoreproc dp,
  SYSTEM_integer *nrrecs);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  const SYSTEM_shortstring *uelfilterstr,
  GXDEFS_tdatastorefiltproc dp);

Function(SYSTEM_int64 ) GXFILE_tgxfileobj_DOT_gdxgetmemoryused(
  GXFILE_tgxfileobj self);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxgetdomainelements(
  GXFILE_tgxfileobj self,
  SYSTEM_integer synr,
  SYSTEM_integer dimpos,
  SYSTEM_integer filternr,
  GXDEFS_tdomainindexproc dp,
  SYSTEM_integer *nrelem,
  SYSTEM_pointer uptr);

Function(SYSTEM_integer ) GXFILE_tgxfileobj_DOT_gdxrenameuel(
  GXFILE_tgxfileobj self,
  const SYSTEM_ansichar *oldname,
  const SYSTEM_ansichar *newname);
extern void * const GXFILE_tgxfileobj_VT[];
extern const SYSTEM_classdescriptor_t GXFILE_tgxfileobj_CD;


extern SYSTEM_shortstring GXFILE_dllloadpath;

extern void _Init_Module_gxfile(void);
extern void _Final_Module_gxfile(void);

#endif /* ! defined _P3___gxfile___H */
