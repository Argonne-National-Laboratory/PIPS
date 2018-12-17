#ifndef _P3___xcompress___H
#define _P3___xcompress___H

typedef SYSTEM_pointer XCOMPRESS_pgzfile;

Function(SYSTEM_boolean ) XCOMPRESS_loadzliblibrary(
  const SYSTEM_ansichar *fn,
  SYSTEM_ansichar *loadmsg);

Procedure XCOMPRESS_unloadzliblibrary(void);

Function(SYSTEM_boolean ) XCOMPRESS_zlibdllloaded(void);

Function(SYSTEM_integer ) XCOMPRESS_compress(
  SYSTEM_pointer pdest,
  CLIBTYPES_clib_ulong *ldest,
  SYSTEM_pointer psrc,
  CLIBTYPES_clib_ulong lsrc);

Function(SYSTEM_integer ) XCOMPRESS_uncompress(
  SYSTEM_pointer pdest,
  CLIBTYPES_clib_ulong *ldest,
  SYSTEM_pointer psrc,
  CLIBTYPES_clib_ulong lsrc);

Function(XCOMPRESS_pgzfile ) XCOMPRESS_gzreadopen(
  const SYSTEM_ansichar *fn);

Function(SYSTEM_integer ) XCOMPRESS_gzread(
  XCOMPRESS_pgzfile pgz,
  SYSTEM_untyped *buf,
  SYSTEM_longword ldest);

Function(SYSTEM_integer ) XCOMPRESS_gzreadclose(
  XCOMPRESS_pgzfile *pgz);

extern void _Init_Module_xcompress(void);
extern void _Final_Module_xcompress(void);

#endif /* ! defined _P3___xcompress___H */
