#include "p3io.h"
#include "gmsspecs.h"
#include "gmsglob.h"

_arr_0GMSGLOB GMSGLOB_ssymboltext = {{2,'e','q'}, {2,'g','t'}, {2,'g','e'}, {2,'l','t'}, {2,'l','e'}, {2,'n','e'}, {4,'p','l','u','s'}, {5,'s','u','b','t','r'}, {4,'m','u','l','t'}, {3,'d','i','v'}, {5,'l','a','g','p','p'}, {5,'l','a','g','m','m'}, {6,'a','s','s','p','a','r'}, {6,'a','s','s','e','q','u'}, {6,'a','s','s','d','o','l'}, {2,'o','r'}, {3,'x','o','r'}, {2,'n','o'}, {3,'y','e','s'}, {2,'n','a'}, {3,'i','n','f'}, {3,'e','p','s'}, {3,'s','u','m'}, {4,'p','r','o','d'}, {4,'s','m','i','n'}, {4,'s','m','a','x'}, {3,'s','c','a'}, {3,'a','c','r'}, {3,'m','o','d'}, {3,'s','e','t'}, {3,'p','a','r'}, {3,'v','a','r'}, {3,'e','q','u'}, {4,'f','i','l','e'}, {3,'p','r','o'}, {3,'p','r','e'}, {3,'m','a','c'}, {4,'f','u','n','c'}, {7,'e','n','d','l','o','o','p'}, {5,'e','n','d','i','f'}, {8,'e','n','d','w','h','i','l','e'}, {6,'e','n','d','f','o','r'}, {3,'f','r','e'}, {3,'b','i','n'}, {3,'p','o','s'}, {3,'n','e','g'}, {3,'i','n','t'}, {4,'s','o','s','1'}, {4,'s','o','s','2'}, {4,'s','e','m','i'}, {7,'s','e','m','i','i','n','t'}, {3,'m','i','n'}, {3,'m','a','x'}, {4,'e','q','u','e'}, {4,'e','q','u','g'}, {4,'e','q','u','l'}, {4,'e','q','u','n'}, {4,'e','q','u','x'}, {4,'e','q','u','c'}, {4,'e','q','u','b'}, {3,'s','e','t'}, {9,'s','i','n','g','l','e','t','o','n'}, {4,'d','i','s','p'}, {5,'a','b','o','r','t'}, {4,'e','x','e','c'}, {4,'l','o','a','d'}, {6,'u','n','l','o','a','d'}, {9,'l','o','a','d','p','o','i','n','t'}, {10,'l','o','a','d','h','a','n','d','l','e'}, {6,'l','o','a','d','d','c'}, {8,'u','n','l','o','a','d','d','i'}, {9,'u','n','l','o','a','d','i','d','x'}, {3,'p','u','t'}, {3,'p','t','l'}, {3,'p','h','d'}, {6,'p','c','l','e','a','r'}, {3,'p','p','g'}, {3,'p','c','l'}, {5,'r','o','u','n','d'}, {6,'s','q','u','a','r','e'}, {5,'c','u','r','l','y'}, {3,'i','m','p'}, {3,'e','q','v'}, {6,'p','b','r','u','c','e'}, {4,'u','n','d','f'}, {5,'o','t','h','e','r'}};
_arr_1GMSGLOB GMSGLOB_defrecvar;
_arr_2GMSGLOB GMSGLOB_defrecequ;

static Procedure GMSGLOB_initdefaultrecords(void)
{
  GMSSPECS_tvarvaltype f;
  GMSSPECS_tvarstyp v;
  typedef GMSGLOB_tssymbol _sub_0INITDEFAULTRECORDS;
  _sub_0INITDEFAULTRECORDS e;

  if ((v = GMSSPECS_stypunknwn) <= (GMSSPECS_stypsemiint)) do {
    for (f = GMSSPECS_vallevel;f <= GMSSPECS_valupper;++f) {
      GMSGLOB_defrecvar[v][f] = 0.0;
    }
    GMSGLOB_defrecvar[v][GMSSPECS_valscale] = 1.0;
  
  } while (v++ != (GMSSPECS_stypsemiint));

  GMSGLOB_defrecvar[GMSSPECS_stypbin][GMSSPECS_valupper] = 1.0;
  GMSGLOB_defrecvar[GMSSPECS_stypint][GMSSPECS_valupper] = 
    GMSSPECS_valiup;
  GMSGLOB_defrecvar[GMSSPECS_styppos][GMSSPECS_valupper] = 
    GMSSPECS_valpin;
  GMSGLOB_defrecvar[GMSSPECS_stypneg][GMSSPECS_vallower] = 
    GMSSPECS_valmin;
  GMSGLOB_defrecvar[GMSSPECS_stypfre][GMSSPECS_vallower] = 
    GMSSPECS_valmin;
  GMSGLOB_defrecvar[GMSSPECS_stypfre][GMSSPECS_valupper] = 
    GMSSPECS_valpin;
  GMSGLOB_defrecvar[GMSSPECS_stypsos1][GMSSPECS_valupper] = 
    GMSSPECS_valpin;
  GMSGLOB_defrecvar[GMSSPECS_stypsos2][GMSSPECS_valupper] = 
    GMSSPECS_valpin;
  GMSGLOB_defrecvar[GMSSPECS_stypsemi][GMSSPECS_vallower] = 1.0;
  GMSGLOB_defrecvar[GMSSPECS_stypsemi][GMSSPECS_valupper] = 
    GMSSPECS_valpin;
  GMSGLOB_defrecvar[GMSSPECS_stypsemiint][GMSSPECS_vallower] = 1.0;
  GMSGLOB_defrecvar[GMSSPECS_stypsemiint][GMSSPECS_valupper] = 
    GMSSPECS_valiup;
  for (e = GMSGLOB_ssyeque;e <= GMSGLOB_ssyequb;++e) {
    for (f = GMSSPECS_vallevel;f <= GMSSPECS_valupper;++f) {
      GMSGLOB_defrecequ[e - 53][f] = 0.0;
    }
    GMSGLOB_defrecequ[e - 53][GMSSPECS_valscale] = 1.0;
  
  }
  GMSGLOB_defrecequ[1][GMSSPECS_valupper] = GMSSPECS_valpin;
  GMSGLOB_defrecequ[2][GMSSPECS_vallower] = GMSSPECS_valmin;
  GMSGLOB_defrecequ[3][GMSSPECS_vallower] = GMSSPECS_valmin;
  GMSGLOB_defrecequ[3][GMSSPECS_valupper] = GMSSPECS_valpin;
  GMSGLOB_defrecequ[5][GMSSPECS_valupper] = GMSSPECS_valpin;
}  /* initdefaultrecords */
static GMSGLOB_tssymbol GMSGLOB_e;

/* unit gmsglob */
void _Init_Module_gmsglob(void)
{
  for (GMSGLOB_e = GMSGLOB_ssyeque;GMSGLOB_e <= GMSGLOB_ssyequb;++
    GMSGLOB_e) {
    GMSSPECS_equstypinfo[ValueCast(GMSSPECS_tequstyp,SYSTEM_ord(
      GMSGLOB_e) - 53)] = SYSTEM_ord(GMSGLOB_e);
  }
  GMSGLOB_initdefaultrecords();
} /* _Init_Module_gmsglob */

void _Final_Module_gmsglob(void)
{
} /* _Final_Module_gmsglob */

