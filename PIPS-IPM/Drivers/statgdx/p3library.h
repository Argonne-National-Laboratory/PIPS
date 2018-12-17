#ifndef _P3___p3library___H
#define _P3___p3library___H

typedef SYSTEM_pointer P3LIBRARY_tlibhandle;

Function(P3LIBRARY_tlibhandle ) P3LIBRARY_p3loadlibrary(
  const SYSTEM_ansichar *lib,
  SYSTEM_ansichar *loadmsg);

Function(SYSTEM_pointer ) P3LIBRARY_p3getprocaddress(
  P3LIBRARY_tlibhandle handle,
  const SYSTEM_ansichar *name);

Function(SYSTEM_boolean ) P3LIBRARY_p3freelibrary(
  P3LIBRARY_tlibhandle handle);

Function(SYSTEM_boolean ) P3LIBRARY_p3libhandleisnil(
  P3LIBRARY_tlibhandle handle);

Function(SYSTEM_ansichar *) P3LIBRARY_p3makelibname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *path,
  const SYSTEM_ansichar *base);

Function(P3LIBRARY_tlibhandle ) P3LIBRARY_p3nillibhandle(void);

Function(SYSTEM_ansichar *) P3LIBRARY_p3libraryext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_ansichar *) P3LIBRARY_p3libraryprefix(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

extern void _Init_Module_p3library(void);
extern void _Final_Module_p3library(void);

#endif /* ! defined _P3___p3library___H */
