#include "p3io.h"
#include "gmsspecs.h"

SYSTEM_double GMSSPECS_old_valund = 1.0e30;
SYSTEM_double GMSSPECS_old_valna = 2.0e30;
SYSTEM_double GMSSPECS_old_valpin = 3.0e30;
SYSTEM_double GMSSPECS_old_valmin = 4.0e30;
SYSTEM_double GMSSPECS_old_valeps = 5.0e30;
SYSTEM_double GMSSPECS_old_valacr = 10.0e30;
SYSTEM_double GMSSPECS_valund = 1.0e300;
SYSTEM_double GMSSPECS_valna = 2.0e300;
SYSTEM_double GMSSPECS_valpin = 3.0e300;
SYSTEM_double GMSSPECS_valmin = 4.0e300;
SYSTEM_double GMSSPECS_valeps = 5.0e300;
SYSTEM_double GMSSPECS_valacr = 10.0e300;
SYSTEM_integer GMSSPECS_valnaint = 2100000000;
SYSTEM_double GMSSPECS_valiup = 3.0e300;
SYSTEM_double GMSSPECS_valbig = 1.0e299;
SYSTEM_double GMSSPECS_valsmall = -1.0e299;
SYSTEM_double GMSSPECS_valtiny = 1.0e-250;
SYSTEM_double GMSSPECS_defiterlim = 2.0e9;
_P3SET_15 GMSSPECS_varstypx = {56,0};
_P3SET_15 GMSSPECS_varstypi = {198,3};
_arr_0GMSSPECS GMSSPECS_varstyptxt = {{8,'u','n','k','n','o','w','n',' '}, {8,'b','i','n','a','r','y',' ',' '}, {8,'i','n','t','e','g','e','r',' '}, {8,'p','o','s','i','t','i','v','e'}, {8,'n','e','g','a','t','i','v','e'}, {8,'f','r','e','e',' ',' ',' ',' '}, {8,'s','o','s','1',' ',' ',' ',' '}, {8,'s','o','s','2',' ',' ',' ',' '}, {8,'s','e','m','i','c','o','n','t'}, {8,'s','e','m','i','i','n','t',' '}};
_arr_1GMSSPECS GMSSPECS_equstypinfo;
_arr_2GMSSPECS GMSSPECS_sufftxt = {{1,'L'}, {1,'M'}, {2,'L','O'}, {2,'U','P'}, {5,'S','C','A','L','E'}};
_arr_3GMSSPECS GMSSPECS_equstyp = {{5,' ','=','E','=',' '}, {5,' ','=','G','=',' '}, {5,' ','=','L','=',' '}, {5,' ','=','N','=',' '}, {5,' ','=','X','=',' '}, {5,' ','=','C','=',' '}, {5,' ','=','B','=',' '}};
_arr_5GMSSPECS GMSSPECS_equctyp = {_P3char('E'), _P3char('G'), _P3char('L'), _P3char('N'), _P3char('X'), _P3char('C'), _P3char('B')};
_arr_7GMSSPECS GMSSPECS_varstyp = {{1,'x'}, {1,'b'}, {1,'i'}, {3,'s','1','s'}, {3,'s','2','s'}, {2,'s','c'}, {2,'s','i'}};
_arr_9GMSSPECS GMSSPECS_solprinttxt = {{9,'0',' ','S','u','m','m','a','r','y'}, {8,'1',' ','R','e','p','o','r','t'}, {7,'2',' ','Q','u','i','e','t'}};
_arr_11GMSSPECS GMSSPECS_handlestattxt = {{9,'0',' ','U','n','k','n','o','w','n'}, {9,'1',' ','R','u','n','n','i','n','g'}, {7,'2',' ','R','e','a','d','y'}, {9,'3',' ','F','a','i','l','u','r','e'}};
_arr_13GMSSPECS GMSSPECS_solvelinktxt = {{14,'0',' ','C','h','a','i','n',' ','S','c','r','i','p','t'}, {13,'1',' ','C','a','l','l',' ','S','c','r','i','p','t'}, {13,'2',' ','C','a','l','l',' ','M','o','d','u','l','e'}, {12,'3',' ','A','s','y','n','c',' ','G','r','i','d'}, {16,'4',' ','A','s','y','n','c',' ','S','i','m','u','l','a','t','e'}, {14,'5',' ','L','o','a','d',' ','L','i','b','r','a','r','y'}};
_arr_15GMSSPECS GMSSPECS_solvestatustxt = {{28,'1',' ','N','o','r','m','a','l',' ','C','o','m','p','l','e','t','i','o','n',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'2',' ','I','t','e','r','a','t','i','o','n',' ','I','n','t','e','r','r','u','p','t',' ',' ',' ',' ',' ',' ',' '}, {28,'3',' ','R','e','s','o','u','r','c','e',' ','I','n','t','e','r','r','u','p','t',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'4',' ','T','e','r','m','i','n','a','t','e','d',' ','B','y',' ','S','o','l','v','e','r',' ',' ',' ',' ',' ',' '}, {28,'5',' ','E','v','a','l','u','a','t','i','o','n',' ','I','n','t','e','r','r','u','p','t',' ',' ',' ',' ',' ',' '}, {28,'6',' ','C','a','p','a','b','i','l','i','t','y',' ','P','r','o','b','l','e','m','s',' ',' ',' ',' ',' ',' ',' '}, {28,'7',' ','L','i','c','e','n','s','i','n','g',' ','P','r','o','b','l','e','m','s',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'8',' ','U','s','e','r',' ','I','n','t','e','r','r','u','p','t',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'9',' ','S','e','t','u','p',' ','F','a','i','l','u','r','e',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','0',' ','S','o','l','v','e','r',' ','F','a','i','l','u','r','e',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','1',' ','I','n','t','e','r','n','a','l',' ','S','o','l','v','e','r',' ','F','a','i','l','u','r','e',' ',' '}, {28,'1','2',' ','S','o','l','v','e',' ','P','r','o','c','e','s','s','i','n','g',' ','S','k','i','p','p','e','d',' '}, {28,'1','3',' ','S','y','s','t','e','m',' ','F','a','i','l','u','r','e',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}};
_arr_17GMSSPECS GMSSPECS_modelstatustxt = {{28,'1',' ','O','p','t','i','m','a','l',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'2',' ','L','o','c','a','l','l','y',' ','O','p','t','i','m','a','l',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'3',' ','U','n','b','o','u','n','d','e','d',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'4',' ','I','n','f','e','a','s','i','b','l','e',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'5',' ','L','o','c','a','l','l','y',' ','I','n','f','e','a','s','i','b','l','e',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'6',' ','I','n','t','e','r','m','e','d','i','a','t','e',' ','I','n','f','e','a','s','i','b','l','e',' ',' ',' '}, {28,'7',' ','F','e','a','s','i','b','l','e',' ','S','o','l','u','t','i','o','n',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'8',' ','I','n','t','e','g','e','r',' ','S','o','l','u','t','i','o','n',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'9',' ','I','n','t','e','r','m','e','d','i','a','t','e',' ','N','o','n','-','I','n','t','e','g','e','r',' ',' '}, {28,'1','0',' ','I','n','t','e','g','e','r',' ','I','n','f','e','a','s','i','b','l','e',' ',' ',' ',' ',' ',' ',' '}, {28,'1','1',' ','L','i','c','e','n','s','i','n','g',' ','P','r','o','b','l','e','m',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','2',' ','E','r','r','o','r',' ','U','n','k','n','o','w','n',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','3',' ','E','r','r','o','r',' ','N','o',' ','S','o','l','u','t','i','o','n',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','4',' ','N','o',' ','S','o','l','u','t','i','o','n',' ','R','e','t','u','r','n','e','d',' ',' ',' ',' ',' '}, {28,'1','5',' ','S','o','l','v','e','d',' ','U','n','i','q','u','e',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','6',' ','S','o','l','v','e','d',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','7',' ','S','o','l','v','e','d',' ','S','i','n','g','u','l','a','r',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '}, {28,'1','8',' ','U','n','b','o','u','n','d','e','d',' ','-',' ','N','o',' ','S','o','l','u','t','i','o','n',' ',' '}, {28,'1','9',' ','I','n','f','e','a','s','i','b','l','e',' ','-',' ','N','o',' ','S','o','l','u','t','i','o','n',' '}};
_arr_19GMSSPECS GMSSPECS_solvetrigger = {0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 1};
_arr_21GMSSPECS GMSSPECS_modeltrigger = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0};
_arr_23GMSSPECS GMSSPECS_modelsolution = {1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0};

Function(GMSSPECS_tgmsvalue ) GMSSPECS_mapval(
  SYSTEM_double x)
{
  GMSSPECS_tgmsvalue result;
  SYSTEM_longint k;

  if (x < GMSSPECS_valund) {
    result = GMSSPECS_xvreal;
    return result;
  } 
  if (x >= GMSSPECS_valacr) { 
    result = GMSSPECS_xvacr;
  } else {
    x = x /  GMSSPECS_valund;
    k = SYSTEM_round(x);
    if (SYSTEM_abs_r(k - x) > 1.0e-5) { 
      result = GMSSPECS_xvund;
    } else 
      switch (k) {
        case 1: 
          result = GMSSPECS_xvund;
          break;
        case 2: 
          result = GMSSPECS_xvna;
          break;
        case 3: 
          result = GMSSPECS_xvpin;
          break;
        case 4: 
          result = GMSSPECS_xvmin;
          break;
        case 5: 
          result = GMSSPECS_xveps;
          break;
        default:
          result = GMSSPECS_xvacr;
      }
  } 
  return result;
}  /* mapval */

Function(GMSSPECS_txgmsvalue ) GMSSPECS_xmapval(
  SYSTEM_double x)
{
  GMSSPECS_txgmsvalue result;

  if (x < GMSSPECS_valund) { 
    if (x < 0) { 
      result = GMSSPECS_vneg;
    } else 
      if (x == 0) { 
        result = GMSSPECS_vzero;
      } else 
        result = GMSSPECS_vpos;
  } else 
    result = ValueCast(GMSSPECS_txgmsvalue,SYSTEM_ord(GMSSPECS_mapval(
      x)) + 2);
  return result;
}  /* xmapval */

Function(SYSTEM_double ) GMSSPECS_old_new_val(
  SYSTEM_double x)
{
  SYSTEM_double result;
  SYSTEM_longint k;

  if (x < GMSSPECS_old_valund) {
    result = x;
    return result;
  } 
  if (x >= GMSSPECS_old_valacr) {
    k = SYSTEM_round(x /  GMSSPECS_old_valacr);
    result = GMSSPECS_valacr * k;
    return result;
  } 
  x = x /  GMSSPECS_old_valund;
  k = SYSTEM_round(x);
  if (SYSTEM_abs_r(k - x) > 1.0e-5) { 
    result = GMSSPECS_valund;
  } else 
    switch (k) {
      case 1: 
        result = GMSSPECS_valund;
        break;
      case 2: 
        result = GMSSPECS_valna;
        break;
      case 3: 
        result = GMSSPECS_valpin;
        break;
      case 4: 
        result = GMSSPECS_valmin;
        break;
      case 5: 
        result = GMSSPECS_valeps;
        break;
      default:
        result = GMSSPECS_valund;
    }
  return result;
}  /* old_new_val */

/* unit gmsspecs */
void _Init_Module_gmsspecs(void)
{
} /* _Init_Module_gmsspecs */

void _Final_Module_gmsspecs(void)
{
} /* _Final_Module_gmsspecs */

