#ifndef _P3___gxdefs___H
#define _P3___gxdefs___H

cnstdef {GXDEFS_domc_unmapped =  -2};
cnstdef {GXDEFS_domc_expand =  -1};
cnstdef {GXDEFS_domc_strict = 0};
typedef SYSTEM_pointer GXDEFS_pgxfile;
typedef GMSSPECS_tindex GXDEFS_tgdxuelindex;
typedef GMSSPECS_tstrindex GXDEFS_tgdxstrindex;
typedef GMSSPECS_tvarreca GXDEFS_tgdxvalues;
typedef SYSTEM_double GXDEFS_tgdxsvals[7];

Prototype Procedure ( STDCALL *GXDEFS_tdatastoreproc)(
const SYSTEM_integer *indx,
const SYSTEM_double *vals);


Prototype Function(SYSTEM_integer ) ( STDCALL *
  GXDEFS_tdatastorefiltproc)(
const SYSTEM_integer *indx,
const SYSTEM_double *vals,
SYSTEM_pointer uptr);


Prototype Procedure ( STDCALL *GXDEFS_tdomainindexproc)(
SYSTEM_integer rawindex,
SYSTEM_integer mappedindex,
SYSTEM_pointer uptr);


Prototype Function(SYSTEM_integer ) ( STDCALL *
  GXDEFS_tdatastorefiltproc_f)(
const SYSTEM_integer *indx,
const SYSTEM_double *vals,
SYSTEM_int64 *uptr);


Prototype Procedure ( STDCALL *GXDEFS_tdomainindexproc_f)(
SYSTEM_integer *rawindex,
SYSTEM_integer *mappedindex,
SYSTEM_int64 *uptr);

typedef _P3STR_7 _arr_0GXDEFS[5];
extern _arr_0GXDEFS GXDEFS_gdxdatatypstr;
typedef _P3STR_15 _arr_1GXDEFS[5];
extern _arr_1GXDEFS GXDEFS_gdxdatatypstrl;
typedef SYSTEM_integer _arr_2GXDEFS[5];
extern _arr_2GXDEFS GXDEFS_datatypsize;
typedef _P3STR_7 _arr_3GXDEFS[7];
extern _arr_3GXDEFS GXDEFS_gdxspecialvaluesstr;

Function(SYSTEM_boolean ) GXDEFS_canbequoted(
  const SYSTEM_ansichar *s);

Function(SYSTEM_boolean ) GXDEFS_gooduelstring(
  const SYSTEM_ansichar *s);

extern void _Init_Module_gxdefs(void);
extern void _Final_Module_gxdefs(void);

#endif /* ! defined _P3___gxdefs___H */
