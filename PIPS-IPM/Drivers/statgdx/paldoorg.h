#ifndef _P3___paldoorg___H
#define _P3___paldoorg___H

cnstdef {PALDOORG_maxlicense = 65};
cnstdef {PALDOORG_licensemaxsub = 28};
extern _P3STR_15 PALDOORG_cmexliccodes;
cnstdef {PALDOORG_graceeval = 30};
cnstdef {PALDOORG_gracemaint = 30};
cnstdef {PALDOORG_gracebeta = 60};
cnstdef {PALDOORG_demomaxrow = 300};
cnstdef {PALDOORG_demomaxcol = 300};
cnstdef {PALDOORG_demomaxnz = 2000};
cnstdef {PALDOORG_demomaxnlnz = 1000};
cnstdef {PALDOORG_demomaxdisc = 50};
cnstdef {PALDOORG_demoglobalmaxrow = 10};
cnstdef {PALDOORG_demoglobalmaxcol = 10};
typedef _P3STR_95 PALDOORG_tlicstring;
typedef struct PALDOORG_tpalobject_OD_S* PALDOORG_tpalobject; /* sy_class */
typedef struct PALDOORG_tpalobject_OD_S {  /* Objects of 'tpalobject' */
  SYSTEM_classreference_t CD;  /* = &PALDOORG_tpalobject_CD */
  _P3STR_95 PALDOORG_tpalobject_DOT_gdlrelcpr;
  _P3STR_31 PALDOORG_tpalobject_DOT_gdlreldat;
  _P3STR_3 PALDOORG_tpalobject_DOT_gdlrelmaj;
  _P3STR_3 PALDOORG_tpalobject_DOT_gdlrelmin;
  _P3STR_3 PALDOORG_tpalobject_DOT_gdlrelgold;
  _P3STR_3 PALDOORG_tpalobject_DOT_gdlrelplc;
  _P3STR_31 PALDOORG_tpalobject_DOT_gdlrelplt;
  _P3STR_3 PALDOORG_tpalobject_DOT_gdlbldcod;
  _P3STR_15 PALDOORG_tpalobject_DOT_gdlrevision;
  _P3STR_31 PALDOORG_tpalobject_DOT_gdlsysnam;
  _P3STR_3 PALDOORG_tpalobject_DOT_gdlsysver;
  _P3STR_15 PALDOORG_tpalobject_DOT_gdllicdat;
  SYSTEM_integer PALDOORG_tpalobject_DOT_gdllicjul;
  _P3STR_95 PALDOORG_tpalobject_DOT_gdlauditline;
  SYSTEM_integer PALDOORG_tpalobject_DOT_juliantoday;
  SYSTEM_integer PALDOORG_tpalobject_DOT_currentsearchhelper,PALDOORG_tpalobject_DOT_optionsearchhelper,PALDOORG_tpalobject_DOT_licenseactsub,PALDOORG_tpalobject_DOT_licensestatus,PALDOORG_tpalobject_DOT_licenselevel;
  SYSTEM_integer PALDOORG_tpalobject_DOT_licenseversion;
  PALDOORG_tlicstring PALDOORG_tpalobject_DOT_license1;
  PALDOORG_tlicstring PALDOORG_tpalobject_DOT_license2;
  PALDOORG_tlicstring PALDOORG_tpalobject_DOT_license3;
  PALDOORG_tlicstring PALDOORG_tpalobject_DOT_license4;
  PALDOORG_tlicstring PALDOORG_tpalobject_DOT_license5;
  PALDOORG_tlicstring PALDOORG_tpalobject_DOT_license6;
  SYSTEM_integer PALDOORG_tpalobject_DOT_subsyssecondary;
  _P3STR_31 PALDOORG_tpalobject_DOT_subsyscode;
  SYSTEM_integer PALDOORG_tpalobject_DOT_checksum,PALDOORG_tpalobject_DOT_isglobal;
  SYSTEM_integer PALDOORG_tpalobject_DOT_rowcnt,PALDOORG_tpalobject_DOT_colcnt,PALDOORG_tpalobject_DOT_nzcnt,PALDOORG_tpalobject_DOT_nlnzcnt,PALDOORG_tpalobject_DOT_disccnt;
  SYSTEM_pointer PALDOORG_tpalobject_DOT_ml;
} PALDOORG_tpalobject_OD;


Procedure PALDOORG_tpalobject_DOT_makedemolicense(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_licensegetmaxsubsys(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_msgadd(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *msg);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_tampercheck(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_licensecheckv1to3(
  PALDOORG_tpalobject self,
  SYSTEM_integer v1,
  SYSTEM_integer v2,
  SYSTEM_integer v3);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_licensecheckinternal(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_licensechecksubinternal(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *msg,
  SYSTEM_integer numcodes,
  const SYSTEM_ansichar *codes);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_licensegetsubstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_licensegetsubeval(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_licensegetsubmaint(
  PALDOORG_tpalobject self);

Constructor(PALDOORG_tpalobject ) PALDOORG_tpalobject_DOT_create(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *msg);

Constructor(PALDOORG_tpalobject ) PALDOORG_tpalobject_DOT_createx(
  PALDOORG_tpalobject self);

Destructor(PALDOORG_tpalobject ) PALDOORG_tpalobject_DOT_destroy(
  PALDOORG_tpalobject self);

Procedure PALDOORG_tpalobject_DOT_gutsofcreate(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensereadw(
  PALDOORG_tpalobject self,
  SYSTEM_text *f);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensescipw(
  PALDOORG_tpalobject self,
  SYSTEM_text *f);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensedisplay(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self,
  SYSTEM_integer i);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensegetdc(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicensegetacademic(
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensemessage(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensereadu(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *filename,
  SYSTEM_ansichar *msg,
  SYSTEM_integer *rc);

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicensevalidateforplatform(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *pf);

Procedure PALDOORG_tpalobject_DOT_pallicensesetsubsearch(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensegetsubnext(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetsubeval(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetsubmaint(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetevaldate(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetmaintdate(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) 
  PALDOORG_tpalobject_DOT_pallicensegetsubmaintdate(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetjulbase(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetjullice(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensechecksize(
  PALDOORG_tpalobject self,
  SYSTEM_integer m,
  SYSTEM_integer n,
  SYSTEM_integer nz,
  SYSTEM_integer nlnz,
  SYSTEM_integer ndisc);

Procedure PALDOORG_tpalobject_DOT_pallicenseclear(
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palolderlicense3(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicenseexist(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensewritec1(
  PALDOORG_tpalobject self,
  SYSTEM_text *f);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensewritec2(
  PALDOORG_tpalobject self,
  SYSTEM_text *f);

Procedure PALDOORG_tpalobject_DOT_pallicensedemo(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_paltampercheck(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetlevel(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_pallicensegetversion(
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetinstitution(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetinstdc(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_longint ) PALDOORG_tpalobject_DOT_pallicensegetkey(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palsecuritycheck(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensegetdates(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *lcode,
  SYSTEM_longint *eval,
  SYSTEM_longint *maint);

Function(SYSTEM_longint ) 
  PALDOORG_tpalobject_DOT_pallicensegetclipcomponent(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *lcode);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_pallicensegetid(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetlicensee(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetleveltext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetvendor(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetplatformtext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetmudtext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensestatusmessage(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensechecksubx(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *sname,
  const SYSTEM_ansichar *_ftmp1,
  SYSTEM_integer *daysleft);

Procedure PALDOORG_tpalobject_DOT_palsetauditline(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *auditline);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palauditrun(
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetcpr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgetver(
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetrel(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetgold(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetcod(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgethdr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgetjul(
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetlicdat(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetbldcod(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetreldat(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_ansichar *) PALDOORG_tpalobject_DOT_palgetrevision(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palisbeta(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_palisalfa(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgettoday(
  PALDOORG_tpalobject self);

Function(SYSTEM_integer ) PALDOORG_tpalobject_DOT_palgetjuliandays(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *s);

Procedure PALDOORG_tpalobject_DOT_palauditfields(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *auditline,
  SYSTEM_ansichar *v1,
  SYSTEM_ansichar *v2,
  SYSTEM_ansichar *v3);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_palgetshortauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);

Procedure PALDOORG_tpalobject_DOT_pallicenseregistergams(
  PALDOORG_tpalobject self,
  SYSTEM_integer linenr,
  const SYSTEM_ansichar *liceline);

Procedure PALDOORG_tpalobject_DOT_pallicenseregistergamsdone(
  PALDOORG_tpalobject self);

Procedure PALDOORG_tpalobject_DOT_pallicenseregistersystem(
  PALDOORG_tpalobject self,
  SYSTEM_integer numcodes,
  const SYSTEM_ansichar *codes,
  SYSTEM_integer magicnum,
  SYSTEM_integer globalflag);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensevalidation(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensecheck(
  PALDOORG_tpalobject self,
  SYSTEM_integer m,
  SYSTEM_integer n,
  SYSTEM_integer nz,
  SYSTEM_integer nlnz,
  SYSTEM_integer ndisc);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicensegetmessage(
  PALDOORG_tpalobject self,
  SYSTEM_ansichar *msg);

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicenseisdemocheckout(
  PALDOORG_tpalobject self);

Procedure PALDOORG_tpalobject_DOT_pallicensedemolimits(
  PALDOORG_tpalobject self,
  SYSTEM_integer *m,
  SYSTEM_integer *n,
  SYSTEM_integer *nz,
  SYSTEM_integer *nlnz,
  SYSTEM_integer *ndisc);

Procedure PALDOORG_tpalobject_DOT_pallicenseglobaldemolimits(
  PALDOORG_tpalobject self,
  SYSTEM_integer *m,
  SYSTEM_integer *n);

Function(SYSTEM_boolean ) PALDOORG_tpalobject_DOT_pallicenseisacademic(
  PALDOORG_tpalobject self);

Function(SYSTEM_boolean ) 
  PALDOORG_tpalobject_DOT_pallicensechecksubsys(
  PALDOORG_tpalobject self,
  const SYSTEM_ansichar *codes);

Function(SYSTEM_ansichar *) 
  PALDOORG_tpalobject_DOT_pallicensegetplatform(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  PALDOORG_tpalobject self);
extern void * const PALDOORG_tpalobject_VT[];
extern const SYSTEM_classdescriptor_t PALDOORG_tpalobject_CD;



extern void _Init_Module_paldoorg(void);
extern void _Final_Module_paldoorg(void);

#endif /* ! defined _P3___paldoorg___H */
