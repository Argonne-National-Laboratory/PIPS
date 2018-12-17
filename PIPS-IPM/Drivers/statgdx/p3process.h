#ifndef _P3___p3process___H
#define _P3___p3process___H

typedef SYSTEM_byte P3PROCESS_tkillhow; /* Anonymous */ enum{P3PROCESS_soft,P3PROCESS_hard};
typedef struct P3PROCESS_tprocinfo_S {
  SYSTEM_cardinal pid;
  SYSTEM_cardinal tid;
  SYSTEM_nativeuint hprocess;
} P3PROCESS_tprocinfo;

typedef P3PROCESS_tprocinfo *P3PROCESS_pprocinfo;
typedef struct P3PROCESS_texecarglist_OD_S* P3PROCESS_texecarglist; /* sy_class */
typedef struct P3PROCESS_texecarglist_OD_S {  /* Objects of 'texecarglist' */
  SYSTEM_classreference_t CD;  /* = &P3PROCESS_texecarglist_CD */
  SYSTEM_integer P3PROCESS_texecarglist_DOT_fcapacity;
  SYSTEM_boolean P3PROCESS_texecarglist_DOT_finheritedhandles;
  SYSTEM_integer P3PROCESS_texecarglist_DOT_fcount;
  SYSTEM_P3_ppointerarray P3PROCESS_texecarglist_DOT_flist;
} P3PROCESS_texecarglist_OD;


Function(SYSTEM_ansichar *) P3PROCESS_texecarglist_DOT_get(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  P3PROCESS_texecarglist self,
  SYSTEM_integer index);

Function(SYSTEM_ansichar *) P3PROCESS_texecarglist_DOT_getlast(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  P3PROCESS_texecarglist self);

Procedure P3PROCESS_texecarglist_DOT_put(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *_ftmp1);

Procedure P3PROCESS_texecarglist_DOT_setcapacity(
  P3PROCESS_texecarglist self,
  SYSTEM_integer newcapacity);

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_split(
  P3PROCESS_texecarglist self,
  SYSTEM_integer append,
  const SYSTEM_ansichar *s);

Procedure P3PROCESS_texecarglist_DOT_grow(
  P3PROCESS_texecarglist self);

Procedure P3PROCESS_texecarglist_DOT_freeitem(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index);

Destructor(P3PROCESS_texecarglist ) P3PROCESS_texecarglist_DOT_destroy(
  P3PROCESS_texecarglist self);

Constructor(P3PROCESS_texecarglist ) P3PROCESS_texecarglist_DOT_create(
  P3PROCESS_texecarglist self);

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_add(
  P3PROCESS_texecarglist self,
  const SYSTEM_ansichar *item);

Procedure P3PROCESS_texecarglist_DOT_clear(
  P3PROCESS_texecarglist self);

Procedure P3PROCESS_texecarglist_DOT_delete(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index);

Procedure P3PROCESS_texecarglist_DOT_insert(
  P3PROCESS_texecarglist self,
  SYSTEM_integer index,
  const SYSTEM_ansichar *item);

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_splitappend(
  P3PROCESS_texecarglist self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) P3PROCESS_texecarglist_DOT_splitprepend(
  P3PROCESS_texecarglist self,
  const SYSTEM_ansichar *s);
extern void * const P3PROCESS_texecarglist_VT[];
extern const SYSTEM_classdescriptor_t P3PROCESS_texecarglist_CD;



Function(SYSTEM_integer ) P3PROCESS_p3systemp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3system(
  const SYSTEM_ansichar *cmd,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3system2(
  const SYSTEM_ansichar *progname,
  const SYSTEM_ansichar *progparams,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3systeml(
  const SYSTEM_ansichar *progname,
  P3PROCESS_texecarglist progparams,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3execp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3exec(
  const SYSTEM_ansichar *cmd,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3exec2(
  const SYSTEM_ansichar *progname,
  const SYSTEM_ansichar *progparams,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3execl(
  const SYSTEM_ansichar *progname,
  P3PROCESS_texecarglist progparams,
  SYSTEM_integer *progrc);

Function(SYSTEM_integer ) P3PROCESS_p3asyncexecp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean newconsole,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg);

Function(SYSTEM_integer ) P3PROCESS_p3asyncsystemp(
  SYSTEM_P3_pansichar cmdptr,
  SYSTEM_boolean newconsole,
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_ansichar *msg);

Function(SYSTEM_integer ) P3PROCESS_p3asyncstatus(
  P3PROCESS_tprocinfo *procinfo,
  SYSTEM_integer *progrc,
  SYSTEM_ansichar *msg);

Function(SYSTEM_cardinal ) P3PROCESS_p3getpid(void);

Function(SYSTEM_boolean ) P3PROCESS_p3ispidvalid(
  SYSTEM_cardinal pid);

Function(SYSTEM_boolean ) P3PROCESS_p3ispidrunning(
  SYSTEM_cardinal pid);

Function(SYSTEM_boolean ) P3PROCESS_p3killprocess(
  const P3PROCESS_tprocinfo *procinfo,
  P3PROCESS_tkillhow how);

Function(SYSTEM_boolean ) P3PROCESS_p3killprocgrouptp(
  const P3PROCESS_tprocinfo *procinfo,
  P3PROCESS_tkillhow how);

Function(SYSTEM_boolean ) P3PROCESS_p3killprocgrouptk(
  const P3PROCESS_tprocinfo *procinfo,
  P3PROCESS_tkillhow how);

Function(SYSTEM_P3_pansichar ) P3PROCESS_p3getcommandline(void);
cnstdef {P3PROCESS_p3ctrlhandlerok = 0};
cnstdef {P3PROCESS_p3ctrlhandlerwasempty = 1};
cnstdef {P3PROCESS_p3ctrlhandlersysfail = 2};

Prototype Procedure (*P3PROCESS_tctrlhandler)(void);


Function(SYSTEM_integer ) P3PROCESS_p3installctrlhandler(
  P3PROCESS_tctrlhandler newhandler);

Function(SYSTEM_integer ) P3PROCESS_p3uninstallctrlhandler(void);

Function(P3PROCESS_tctrlhandler ) P3PROCESS_p3getctrlhandler(void);

Procedure P3PROCESS_p3defaultshowwindow(void);

Procedure P3PROCESS_p3setshowwindow(
  SYSTEM_integer showwindow);

Function(SYSTEM_integer ) P3PROCESS_p3getshowwindow(void);

Function(SYSTEM_integer ) P3PROCESS_p3getnumberofprocessors(void);

extern void _Init_Module_p3process(void);
extern void _Final_Module_p3process(void);

#endif /* ! defined _P3___p3process___H */
