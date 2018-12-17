#ifndef _P3___gmsglob___H
#define _P3___gmsglob___H

typedef SYSTEM_byte GMSGLOB_tssymbol; /* Anonymous */ enum{GMSGLOB_ssyeq,GMSGLOB_ssygt,GMSGLOB_ssyge,GMSGLOB_ssylt,
  GMSGLOB_ssyle,GMSGLOB_ssyne,GMSGLOB_ssyplus,GMSGLOB_ssysubtr,
  GMSGLOB_ssymult,GMSGLOB_ssydiv,GMSGLOB_ssylagpp,GMSGLOB_ssylagmm,
  GMSGLOB_ssyasspar,GMSGLOB_ssyassequ,GMSGLOB_ssyassdol,GMSGLOB_ssyor,
  GMSGLOB_ssyxor,GMSGLOB_ssyno,GMSGLOB_ssyyes,GMSGLOB_ssyna,
  GMSGLOB_ssyinf,GMSGLOB_ssyeps,GMSGLOB_ssysum,GMSGLOB_ssyprod,
  GMSGLOB_ssysmin,GMSGLOB_ssysmax,GMSGLOB_ssysca,GMSGLOB_ssyacr,
  GMSGLOB_ssymod,GMSGLOB_ssyset,GMSGLOB_ssypar,GMSGLOB_ssyvar,
  GMSGLOB_ssyequ,GMSGLOB_ssyfile,GMSGLOB_ssypro,GMSGLOB_ssypre,
  GMSGLOB_ssymac,GMSGLOB_ssyfunc,GMSGLOB_ssyendloop,GMSGLOB_ssyendif,
  GMSGLOB_ssyendwhile,GMSGLOB_ssyendfor,GMSGLOB_ssyfre,GMSGLOB_ssybin,
  GMSGLOB_ssypos,GMSGLOB_ssyneg,GMSGLOB_ssyint,GMSGLOB_ssysos1,
  GMSGLOB_ssysos2,GMSGLOB_ssysemi,GMSGLOB_ssysemiint,GMSGLOB_ssymin,
  GMSGLOB_ssymax,GMSGLOB_ssyeque,GMSGLOB_ssyequg,GMSGLOB_ssyequl,
  GMSGLOB_ssyequn,GMSGLOB_ssyequx,GMSGLOB_ssyequc,GMSGLOB_ssyequb,
  GMSGLOB_ssysetm,GMSGLOB_ssysets,GMSGLOB_ssydisp,GMSGLOB_ssyabort,
  GMSGLOB_ssyexec,GMSGLOB_ssyload,GMSGLOB_ssyunload,
  GMSGLOB_ssyloadpoint,GMSGLOB_ssyloadhandle,GMSGLOB_ssyloaddc,
  GMSGLOB_ssyunloaddi,GMSGLOB_ssyunloadidx,GMSGLOB_ssyput,
  GMSGLOB_ssyptl,GMSGLOB_ssyphd,GMSGLOB_ssypclear,GMSGLOB_ssyppg,
  GMSGLOB_ssypcl,GMSGLOB_ssyround,GMSGLOB_ssysquare,GMSGLOB_ssycurly,
  GMSGLOB_ssyimp,GMSGLOB_ssyeqv,GMSGLOB_ssypbruce,GMSGLOB_ssyundf,
  GMSGLOB_ssyother};
typedef _P3STR_15 _arr_0GMSGLOB[86];
extern _arr_0GMSGLOB GMSGLOB_ssymboltext;
typedef GMSSPECS_tvarreca _arr_1GMSGLOB[10];
extern _arr_1GMSGLOB GMSGLOB_defrecvar;
typedef GMSGLOB_tssymbol _sub_3GMSGLOB;
typedef GMSSPECS_tvarreca _arr_2GMSGLOB[7];
extern _arr_2GMSGLOB GMSGLOB_defrecequ;

extern void _Init_Module_gmsglob(void);
extern void _Final_Module_gmsglob(void);

#endif /* ! defined _P3___gmsglob___H */
