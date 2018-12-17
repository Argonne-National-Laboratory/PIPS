#ifndef _P3___p3platform___H
#define _P3___p3platform___H

typedef SYSTEM_byte P3PLATFORM_tosfiletype; /* Anonymous */ enum{P3PLATFORM_osfilewin,P3PLATFORM_osfileunix,
  P3PLATFORM_osfilemissing};
typedef _P3STR_31 _arr_0P3PLATFORM[3];
extern _arr_0P3PLATFORM P3PLATFORM_osfiletypetext;
typedef SYSTEM_byte P3PLATFORM_tosplatform; /* Anonymous */ enum{P3PLATFORM_oswindowsdos,P3PLATFORM_oswindows95,
  P3PLATFORM_oswindowsnt,P3PLATFORM_oswindows64emt,P3PLATFORM_osaix,
  P3PLATFORM_oshpux,P3PLATFORM_osirix,P3PLATFORM_oslinux,
  P3PLATFORM_ossunos_sparc32,P3PLATFORM_ososf1,P3PLATFORM_osdarwin,
  P3PLATFORM_oslinux86_64,P3PLATFORM_ossunos_i86pc,
  P3PLATFORM_ossunos_sparc64,P3PLATFORM_osdarwin_i386,
  P3PLATFORM_osdarwin_x64,P3PLATFORM_osbluegene,P3PLATFORM_osmissing};
typedef _P3STR_63 _arr_1P3PLATFORM[18];
extern _arr_1P3PLATFORM P3PLATFORM_osplatformtext;
typedef _P3STR_7 _arr_2P3PLATFORM[18];
extern _arr_2P3PLATFORM P3PLATFORM_osdllextension;
typedef _P3STR_3 _arr_3P3PLATFORM[18];
extern _arr_3P3PLATFORM P3PLATFORM_osdllprefix;

Function(P3PLATFORM_tosfiletype ) P3PLATFORM_osfiletype(void);

Function(P3PLATFORM_tosplatform ) P3PLATFORM_osplatform(void);

Function(SYSTEM_ansichar *) P3PLATFORM_osnullfilename(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_ansichar *) P3PLATFORM_osconsolename(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_ansichar *) P3PLATFORM_oslanguagepascal(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_ansichar *) P3PLATFORM_oslanguagec(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_boolean ) P3PLATFORM_nativeislittleendian(void);

extern void _Init_Module_p3platform(void);
extern void _Final_Module_p3platform(void);

#endif /* ! defined _P3___p3platform___H */
