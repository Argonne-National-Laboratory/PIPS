#ifndef _P3___runner___H
#define _P3___runner___H

cnstdef {RUNNER_ec_cannot_modify = 1};
cnstdef {RUNNER_ec_process_active = 2};
cnstdef {RUNNER_ec_empty_cmd_line = 3};
typedef SYSTEM_byte RUNNER_tvisible; /* Anonymous */ enum{RUNNER_vis_hide,RUNNER_vis_minimized,RUNNER_vis_normal};
typedef struct RUNNER_tmsghandler_OD_S* RUNNER_tmsghandler; /* sy_class */
typedef struct RUNNER_tmsghandler_OD_S {  /* Objects of 'tmsghandler' */
  SYSTEM_classreference_t CD;  /* = &RUNNER_tmsghandler_CD */
  SYSTEM_integer RUNNER_tmsghandler_DOT_fverbose;
  SYSTEM_shortstring RUNNER_tmsghandler_DOT_fmsgpfx;
} RUNNER_tmsghandler_OD;


Constructor(RUNNER_tmsghandler ) RUNNER_tmsghandler_DOT_create(
  RUNNER_tmsghandler self,
  const SYSTEM_ansichar *msgpfx);

Procedure RUNNER_tmsghandler_DOT_errormessage(
  RUNNER_tmsghandler self,
  SYSTEM_integer ec,
  const SYSTEM_ansichar *s);

Procedure RUNNER_tmsghandler_DOT_logmessage(
  RUNNER_tmsghandler self,
  const SYSTEM_ansichar *s);

Procedure RUNNER_tmsghandler_DOT_debugmessage(
  RUNNER_tmsghandler self,
  const SYSTEM_ansichar *s);
extern void * const RUNNER_tmsghandler_VT[];
extern const SYSTEM_classdescriptor_t RUNNER_tmsghandler_CD;


typedef struct RUNNER_trunner_OD_S* RUNNER_trunner; /* sy_class */
typedef struct RUNNER_trunner_OD_S {  /* Objects of 'trunner' */
  SYSTEM_classreference_t CD;  /* = &RUNNER_trunner_CD */
  RUNNER_tmsghandler RUNNER_trunner_DOT_fmsghandler;
  SYSTEM_shortstring RUNNER_trunner_DOT_fexecutable;
  GMSOBJ_txstrings RUNNER_trunner_DOT_fparams;
  SYSTEM_shortstring RUNNER_trunner_DOT_fworkdir;
  SYSTEM_P3_pansichar RUNNER_trunner_DOT_fcommandline;
  SYSTEM_boolean RUNNER_trunner_DOT_fisrunning;
  SYSTEM_boolean RUNNER_trunner_DOT_finherithandles;
  SYSTEM_boolean RUNNER_trunner_DOT_fuseshell;
  RUNNER_tvisible RUNNER_trunner_DOT_fvisible;
  SYSTEM_integer RUNNER_trunner_DOT_fprogrc;
} RUNNER_trunner_OD;


Function(SYSTEM_boolean ) RUNNER_trunner_DOT_errorwhenrunning(
  RUNNER_trunner self,
  const SYSTEM_ansichar *s);

Procedure RUNNER_trunner_DOT_setexecutable(
  RUNNER_trunner self,
  const SYSTEM_ansichar *v);

Procedure RUNNER_trunner_DOT_setworkdir(
  RUNNER_trunner self,
  const SYSTEM_ansichar *v);

Procedure RUNNER_trunner_DOT_setinherithandles(
  RUNNER_trunner self,
  SYSTEM_boolean v);

Procedure RUNNER_trunner_DOT_setuseshell(
  RUNNER_trunner self,
  SYSTEM_boolean v);

Procedure RUNNER_trunner_DOT_commandlinechanged(
  RUNNER_trunner self);

Function(SYSTEM_integer ) RUNNER_trunner_DOT_getverbose(
  RUNNER_trunner self);

Procedure RUNNER_trunner_DOT_setverbose(
  RUNNER_trunner self,
  SYSTEM_integer v);

Procedure RUNNER_trunner_DOT_setvisible(
  RUNNER_trunner self,
  RUNNER_tvisible v);

Constructor(RUNNER_trunner ) RUNNER_trunner_DOT_create(
  RUNNER_trunner self);

Destructor(RUNNER_trunner ) RUNNER_trunner_DOT_destroy(
  RUNNER_trunner self);

Procedure RUNNER_trunner_DOT_paramsadd(
  RUNNER_trunner self,
  const SYSTEM_ansichar *v);

Procedure RUNNER_trunner_DOT_paramsclear(
  RUNNER_trunner self);

Function(SYSTEM_integer ) RUNNER_trunner_DOT_paramscount(
  RUNNER_trunner self);

Function(SYSTEM_P3_pansichar ) RUNNER_trunner_DOT_commandline(
  RUNNER_trunner self);

Function(SYSTEM_integer ) RUNNER_trunner_DOT_startandwait(
  RUNNER_trunner self);
extern void * const RUNNER_trunner_VT[];
extern const SYSTEM_classdescriptor_t RUNNER_trunner_CD;



extern void _Init_Module_runner(void);
extern void _Final_Module_runner(void);

#endif /* ! defined _P3___runner___H */
