#include "p3io.h"
#include "gmsspecs.h"
#include "gxdefs.h"

_arr_0GXDEFS GXDEFS_gdxdatatypstr = {{3,'S','e','t'}, {3,'P','a','r'}, {3,'V','a','r'}, {3,'E','q','u'}, {5,'A','l','i','a','s'}};
_arr_1GXDEFS GXDEFS_gdxdatatypstrl = {{3,'S','e','t'}, {9,'P','a','r','a','m','e','t','e','r'}, {8,'V','a','r','i','a','b','l','e'}, {8,'E','q','u','a','t','i','o','n'}, {5,'A','l','i','a','s'}};
_arr_2GXDEFS GXDEFS_datatypsize = {1, 1, 5, 5, 0};
_arr_3GXDEFS GXDEFS_gdxspecialvaluesstr = {{4,'U','n','d','f'}, {2,'N','A'}, {4,'+','I','n','f'}, {4,'-','I','n','f'}, {3,'E','p','s'}, {1,'0'}, {5,'A','c','r','o','N'}};

Function(SYSTEM_boolean ) GXDEFS_canbequoted(
  const SYSTEM_ansichar *s)
{
  SYSTEM_boolean result;
  SYSTEM_boolean saw_single, saw_double;
  SYSTEM_integer i;
  SYSTEM_ansichar ch;

  result = SYSTEM_false;
  saw_single = SYSTEM_false;
  saw_double = SYSTEM_false;
  { register SYSTEM_int32 _stop = SYSTEM_length(s);
    if ((i = 1) <=  _stop) do {
      ch = s[i];
      if (ch == _P3char('\'')) {
        if (saw_double) 
          return result;
        saw_single = SYSTEM_true;
      } else 
        if (ch == _P3char('\"')) {
          if (saw_single) 
            return result;
          saw_double = SYSTEM_true;
        } else 
          if (ch < _P3char(' ')) 
            return result;
    
    } while (i++ !=  _stop);

  }
  result = SYSTEM_true;
  return result;
}  /* canbequoted */

Function(SYSTEM_boolean ) GXDEFS_gooduelstring(
  const SYSTEM_ansichar *s)
{
  SYSTEM_boolean result;

  result = SYSTEM_length(s) <= GMSSPECS_maxnamelen;
  if (result) 
    result = GXDEFS_canbequoted(s);
  return result;
}  /* gooduelstring */

/* unit gxdefs */
void _Init_Module_gxdefs(void)
{
} /* _Init_Module_gxdefs */

void _Final_Module_gxdefs(void)
{
} /* _Final_Module_gxdefs */

