#ifndef _P3___gmsglobx___H
#define _P3___gmsglobx___H

typedef SYSTEM_byte GMSGLOBX_tsycsetconstants; /* Anonymous */ enum{GMSGLOBX_sycmodeltypes,GMSGLOBX_sycgamsparameters,
  GMSGLOBX_sycgamsparametersynonyms,
  GMSGLOBX_sycgamsparametersynonymmap,GMSGLOBX_sycdollaroptions,
  GMSGLOBX_sycgamsfunctions,GMSGLOBX_sycsystemsuffixes,
  GMSGLOBX_sycempty,GMSGLOBX_sycpredefinedsymbols,
  GMSGLOBX_sycgussmodelattributes,GMSGLOBX_sycsetconstants,
  GMSGLOBX_sycsolvernames,GMSGLOBX_sycplatforms,GMSGLOBX_sycvendors,
  GMSGLOBX_syccomponents,GMSGLOBX_sycclipcodes,
  GMSGLOBX_sycgamslicenses,GMSGLOBX_sycgamslicensetypes,
  GMSGLOBX_syccomponentsolvermap,GMSGLOBX_sycclipcomponentmap,
  GMSGLOBX_sycsolverplatformmap,GMSGLOBX_sycsolvertypeplatformmap};
cnstdef {GMSGLOBX_maxsetconstants = 22};
cnstdef {GMSGLOBX_maxsetconstantslength = 23};

Function(SYSTEM_ansichar *) GMSGLOBX_setconstantskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_setconstantstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_setconstantslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxsolvernames = 145};
cnstdef {GMSGLOBX_maxsolvernameslength = 12};

Function(SYSTEM_ansichar *) GMSGLOBX_solvernameskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_solvernamestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_solvernameslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxplatforms = 13};
cnstdef {GMSGLOBX_maxplatformslength = 3};

Function(SYSTEM_ansichar *) GMSGLOBX_platformskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_platformstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_platformstext2(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_platformslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxvendors = 21};
cnstdef {GMSGLOBX_maxvendorslength = 1};

Function(SYSTEM_ansichar *) GMSGLOBX_vendorskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_vendorstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_vendorslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxcomponents = 47};
cnstdef {GMSGLOBX_maxcomponentslength = 14};

Function(SYSTEM_ansichar *) GMSGLOBX_componentskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_componentstext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_componentslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxclipcodes = 62};
cnstdef {GMSGLOBX_maxclipcodeslength = 2};

Function(SYSTEM_ansichar *) GMSGLOBX_clipcodeskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_clipcodestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_clipcodeslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxgamslicenses = 5};
cnstdef {GMSGLOBX_maxgamslicenseslength = 2};

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicenseskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicensestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_gamslicenseslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxgamslicensetypes = 20};
cnstdef {GMSGLOBX_maxgamslicensetypeslength = 1};

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicensetypeskey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_gamslicensetypestext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_gamslicensetypeslookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxmodeltypesx = 16};
cnstdef {GMSGLOBX_maxmodeltypesxlength = 6};

Function(SYSTEM_ansichar *) GMSGLOBX_modeltypesxkey(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_ansichar *) GMSGLOBX_modeltypesxtext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer n);

Function(SYSTEM_integer ) GMSGLOBX_modeltypesxlookup(
  const SYSTEM_ansichar *_ftmp1);
cnstdef {GMSGLOBX_maxcomponentsolvermap = 162};

Function(SYSTEM_integer ) GMSGLOBX_componentsolvermapmap(
  SYSTEM_integer i,
  SYSTEM_integer j);
cnstdef {GMSGLOBX_maxclipcomponentmap = 51};

Function(SYSTEM_integer ) GMSGLOBX_clipcomponentmapmap(
  SYSTEM_integer i,
  SYSTEM_integer j);
cnstdef {GMSGLOBX_maxsolverplatformmap = 815};

Function(SYSTEM_integer ) GMSGLOBX_solverplatformmapmap(
  SYSTEM_integer i,
  SYSTEM_integer j);
cnstdef {GMSGLOBX_maxsolvertypeplatformmap = 2831};

Function(SYSTEM_integer ) GMSGLOBX_solvertypeplatformmapmap(
  SYSTEM_integer i,
  SYSTEM_integer j);

Function(SYSTEM_ansichar *) GMSGLOBX_hostplatform(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

extern void _Init_Module_gmsglobx(void);
extern void _Final_Module_gmsglobx(void);

#endif /* ! defined _P3___gmsglobx___H */
