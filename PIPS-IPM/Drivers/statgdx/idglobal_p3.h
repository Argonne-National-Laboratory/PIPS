#ifndef _P3___idglobal_p3___H
#define _P3___idglobal_p3___H


Function(SYSTEM_cardinal ) IDGLOBAL_P3_gettickcount(void);

Function(SYSTEM_cardinal ) IDGLOBAL_P3_gettickdiff(
  SYSTEM_cardinal aoldtickcount,
  SYSTEM_cardinal anewtickcount);

extern void _Init_Module_idglobal_p3(void);
extern void _Final_Module_idglobal_p3(void);

#endif /* ! defined _P3___idglobal_p3___H */
