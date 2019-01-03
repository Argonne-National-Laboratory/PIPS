#include "p3io.h"
#include "p3platform.h"
#include "p3utils.h"
#include "system_p3.h"
#include "p3process.h"
#include "p3library.h"
#include "math_p3.h"
#include "p3ieeefp.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "p3threads.h"
#include "idglobal_p3.h"
#include "paldoorg.h"
#include "gdlaudit.h"

static PALDOORG_tpalobject GDLAUDIT_palobj;
static SYSTEM_shortstring GDLAUDIT_msg;

Procedure GDLAUDIT_gdlsetauditline(
  const SYSTEM_ansichar *auditline)
{
  PALDOORG_tpalobject_DOT_palsetauditline(GDLAUDIT_palobj,auditline);
}  /* gdlsetauditline */

Procedure GDLAUDIT_gdlsetauditlinelib(
  const SYSTEM_ansichar *auditline)
{
  PALDOORG_tpalobject_DOT_palsetauditline(GDLAUDIT_palobj,auditline);
}  /* gdlsetauditlinelib */

Function(SYSTEM_boolean ) GDLAUDIT_gdlauditrun(void)
{
  SYSTEM_boolean result;

  result = PALDOORG_tpalobject_DOT_palauditrun(GDLAUDIT_palobj);
  return result;
}  /* gdlauditrun */

Function(SYSTEM_ansichar *) GDLAUDIT_gdlgetauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  PALDOORG_tpalobject_DOT_palgetauditline(result,_len_ret,
    GDLAUDIT_palobj);
  return result;
}  /* gdlgetauditline */

Function(SYSTEM_ansichar *) GDLAUDIT_gdlgetshortauditline(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret)
{
  PALDOORG_tpalobject_DOT_palgetshortauditline(result,_len_ret,
    GDLAUDIT_palobj);
  return result;
}  /* gdlgetshortauditline */

Procedure GDLAUDIT_gdlauditfields(
  const SYSTEM_ansichar *auditline,
  SYSTEM_ansichar *v1,
  SYSTEM_ansichar *v2,
  SYSTEM_ansichar *v3)
{
  PALDOORG_tpalobject_DOT_palauditfields(GDLAUDIT_palobj,auditline,v1,
    v2,v3);
}  /* gdlauditfields */

/* unit gdlaudit */
void _Init_Module_gdlaudit(void)
{
  GDLAUDIT_palobj = ValueCast(PALDOORG_tpalobject,
    PALDOORG_tpalobject_DOT_create(ValueCast(PALDOORG_tpalobject,
    _P3alloc_object(&PALDOORG_tpalobject_CD)),GDLAUDIT_msg));
} /* _Init_Module_gdlaudit */

void _Final_Module_gdlaudit(void)
{
  SYSUTILS_P3_freeandnil(&GDLAUDIT_palobj);
} /* _Final_Module_gdlaudit */

