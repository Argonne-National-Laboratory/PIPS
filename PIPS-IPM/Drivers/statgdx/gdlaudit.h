#ifndef _P3___gdlaudit___H
#define _P3___gdlaudit___H


Procedure GDLAUDIT_gdlsetauditline(
  const SYSTEM_ansichar *auditline);

Procedure GDLAUDIT_gdlsetauditlinelib(
  const SYSTEM_ansichar *auditline);

Function(SYSTEM_boolean ) GDLAUDIT_gdlauditrun(void);

Function(SYSTEM_ansichar *) GDLAUDIT_gdlgetauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_ansichar *) GDLAUDIT_gdlgetshortauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Procedure GDLAUDIT_gdlauditfields(
  const SYSTEM_ansichar *auditline,
  SYSTEM_ansichar *v1,
  SYSTEM_ansichar *v2,
  SYSTEM_ansichar *v3);

extern void _Init_Module_gdlaudit(void);
extern void _Final_Module_gdlaudit(void);

#endif /* ! defined _P3___gdlaudit___H */
