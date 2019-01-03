#ifndef _P3___p3utils___H
#define _P3___p3utils___H


Function(SYSTEM_integer ) P3UTILS_p3chmod(
  const SYSTEM_ansichar *path,
  SYSTEM_integer mode);

Function(SYSTEM_double ) P3UTILS_realtrunc(
  SYSTEM_double x);

Function(SYSTEM_double ) P3UTILS_realround(
  SYSTEM_double x);

Function(SYSTEM_ansichar *) P3UTILS_floattoe(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double y,
  SYSTEM_integer decimals);

Function(SYSTEM_ansichar *) P3UTILS_replacefileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *extension);

Function(SYSTEM_ansichar *) P3UTILS_completefileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *extension);

Function(SYSTEM_ansichar *) P3UTILS_paramstrzero(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_boolean ) P3UTILS_prefixpath(
  const SYSTEM_ansichar *s);

Function(SYSTEM_boolean ) P3UTILS_prefixloadpath(
  const SYSTEM_ansichar *dir);

Function(SYSTEM_boolean ) P3UTILS_p3setenv(
  const SYSTEM_ansichar *name,
  const SYSTEM_ansichar *val);

Procedure P3UTILS_p3unsetenv(
  const SYSTEM_ansichar *name);

Function(SYSTEM_boolean ) P3UTILS_p3issetenv(
  const SYSTEM_ansichar *name);

Procedure P3UTILS_p3setconsoletitle(
  const SYSTEM_ansichar *s);

Procedure P3UTILS_p3nopopups(void);
typedef SYSTEM_pointer P3UTILS_tp3filehandle;
typedef SYSTEM_byte P3UTILS_tp3fileopenaction; /* Anonymous */ enum{P3UTILS_p3openread,P3UTILS_p3openwrite,P3UTILS_p3openreadwrite};

Function(SYSTEM_integer ) P3UTILS_p3fileopen(
  const SYSTEM_ansichar *fname,
  P3UTILS_tp3fileopenaction mode,
  P3UTILS_tp3filehandle *h);

Function(SYSTEM_integer ) P3UTILS_p3fileclose(
  P3UTILS_tp3filehandle *h);

Function(SYSTEM_integer ) P3UTILS_p3fileread(
  P3UTILS_tp3filehandle h,
  SYSTEM_untyped *buffer,
  SYSTEM_longword buflen,
  SYSTEM_longword *numread);

Function(SYSTEM_integer ) P3UTILS_p3filewrite(
  P3UTILS_tp3filehandle h,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword buflen,
  SYSTEM_longword *numwritten);

Function(SYSTEM_integer ) P3UTILS_p3filegetsize(
  P3UTILS_tp3filehandle h,
  SYSTEM_int64 *filesize);

Function(SYSTEM_integer ) P3UTILS_p3filesetpointer(
  P3UTILS_tp3filehandle h,
  SYSTEM_int64 distance,
  SYSTEM_int64 *newpointer,
  SYSTEM_longword whence);

Function(SYSTEM_integer ) P3UTILS_p3filegetpointer(
  P3UTILS_tp3filehandle h,
  SYSTEM_int64 *filepointer);
cnstdef {P3UTILS_p3_file_begin = 0};
cnstdef {P3UTILS_p3_file_current = 1};
cnstdef {P3UTILS_p3_file_end = 2};

Function(SYSTEM_int64 ) P3UTILS_p3allocmemsize64(void);

Function(SYSTEM_pointer ) P3UTILS_p3allocmem64(
  SYSTEM_int64 size);

Procedure P3UTILS_p3fillchar64(
  SYSTEM_untyped *p,
  SYSTEM_int64 size,
  SYSTEM_byte fillvalue);

Procedure P3UTILS_p3freemem64(
  SYSTEM_pointer *p,
  SYSTEM_int64 size);

Procedure P3UTILS_p3getmem64(
  SYSTEM_pointer *p,
  SYSTEM_int64 size);

Procedure P3UTILS_p3move64(
  const SYSTEM_untyped *source,
  SYSTEM_untyped *dest,
  SYSTEM_int64 sz);

Procedure P3UTILS_p3reallocmem64(
  SYSTEM_pointer *p,
  SYSTEM_int64 size);

Prototype Function(SYSTEM_boolean ) (*P3UTILS_thavedatacb)(
SYSTEM_untyped *buf,
SYSTEM_integer len,
SYSTEM_pointer usermem);


Procedure P3UTILS_p3getfromurl(
  const SYSTEM_ansichar *servername,
  const SYSTEM_ansichar *filename,
  SYSTEM_word port,
  P3UTILS_thavedatacb havedata,
  SYSTEM_pointer usermem,
  SYSTEM_ansichar *msg);
cnstdef {P3UTILS_os_unknown = 0};
cnstdef {P3UTILS_os_win9x = 1};
cnstdef {P3UTILS_os_winnt = 2};
cnstdef {P3UTILS_os_wine = 3};
cnstdef {P3UTILS_os_hx = 4};

Function(SYSTEM_integer ) P3UTILS_p3getwindowsversion(void);

Function(SYSTEM_ansichar *) P3UTILS_p3getcomputername(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_ansichar *) P3UTILS_p3getusername(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_boolean ) P3UTILS_p3senddatamessage(
  SYSTEM_boolean broadcast,
  const SYSTEM_ansichar *_ftmp1,
  const SYSTEM_ansichar *_ftmp2);

Function(SYSTEM_ansichar *) P3UTILS_p3pushdeflocale(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Procedure P3UTILS_p3popdeflocale(
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) P3UTILS_p3getexecname(
  SYSTEM_ansichar *execname,
  SYSTEM_ansichar *msg);

Function(SYSTEM_integer ) P3UTILS_p3getlibname(
  SYSTEM_ansichar *libname,
  SYSTEM_ansichar *msg);

Function(SYSTEM_integer ) P3UTILS_p3someioresult(void);

extern void _Init_Module_p3utils(void);
extern void _Final_Module_p3utils(void);

#endif /* ! defined _P3___p3utils___H */
