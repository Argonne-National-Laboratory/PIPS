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
#include "clibtypes.h"
#include "gmslibname.h"
#include "xcompress.h"

static P3LIBRARY_tlibhandle XCOMPRESS_zlibhandle;

Prototype Function(SYSTEM_integer ) ( CDECL *_Func_0XCOMPRESS)(
SYSTEM_pointer pdest,
CLIBTYPES_clib_ulong *ldest,
SYSTEM_pointer psrc,
CLIBTYPES_clib_ulong lsrc);

static _Func_0XCOMPRESS XCOMPRESS_pcompress;

Prototype Function(SYSTEM_integer ) ( CDECL *_Func_1XCOMPRESS)(
SYSTEM_pointer pdest,
CLIBTYPES_clib_ulong *ldest,
SYSTEM_pointer psrc,
CLIBTYPES_clib_ulong lsrc);

static _Func_1XCOMPRESS XCOMPRESS_puncompress;

Prototype Function(XCOMPRESS_pgzfile ) ( CDECL *_Func_2XCOMPRESS)(
SYSTEM_P3_pchar fn,
SYSTEM_P3_pchar mode);

static _Func_2XCOMPRESS XCOMPRESS_pgzreadopen;

Prototype Function(SYSTEM_integer ) ( CDECL *_Func_3XCOMPRESS)(
XCOMPRESS_pgzfile pgz,
SYSTEM_pointer pdest,
SYSTEM_cardinal ldest);

static _Func_3XCOMPRESS XCOMPRESS_pgzread;

Prototype Function(SYSTEM_integer ) ( CDECL *_Func_4XCOMPRESS)(
XCOMPRESS_pgzfile pgz);

static _Func_4XCOMPRESS XCOMPRESS_pgzreadclose;

Function(SYSTEM_boolean ) XCOMPRESS_zlibdllloaded(void)
{
  SYSTEM_boolean result;

  result = XCOMPRESS_zlibhandle != ValueCast(SYSTEM_pointer,0);
  return result;
}  /* zlibdllloaded */

static Function(SYSTEM_pointer ) loadentry(
  const SYSTEM_ansichar *n,
  SYSTEM_ansichar *_2wfn,
  SYSTEM_ansichar *_2loadmsg)
{
  SYSTEM_pointer result;

  if (_P3strcmpN(_2loadmsg,_P3str1("\000"))) { 
    result = NULL;
  } else {
    {
      SYSTEM_shortstring _t1;

      result = P3LIBRARY_p3getprocaddress(XCOMPRESS_zlibhandle,
        SYSUTILS_P3_lowercase(_t1,255,n));
    }
    if (result == NULL) 
      {
        _P3STR_255 _t1;
        _P3STR_255 _t2;

        _P3strcat(_2loadmsg,255,_P3strcat(_t2,255,_P3strcat(_t1,255,_P3str1("\021Entry not found: "),
          n),_P3str1("\004 in ")),_2wfn);
      }
  } 
  return result;
}  /* loadentry */

Function(SYSTEM_boolean ) XCOMPRESS_loadzliblibrary(
  const SYSTEM_ansichar *fn,
  SYSTEM_ansichar *loadmsg)
{
  SYSTEM_boolean result;
  SYSTEM_shortstring wfn, basename;
  SYSTEM_shortstring path;

  _P3strclr(loadmsg);
  if (XCOMPRESS_zlibhandle == ValueCast(SYSTEM_pointer,0)) {
    SYSUTILS_P3_extractfilepath(path,255,fn);
    SYSUTILS_P3_extractfilename(basename,255,fn);
    if (_P3strcmpE(basename,_P3str1("\000"))) 
      _P3strcpy(basename,255,_P3str1("\010gmszlib1"));
    GMSLIBNAME_gamslibnamep3(wfn,255,basename);
    _P3strcat(wfn,255,path,wfn);
    XCOMPRESS_zlibhandle = P3LIBRARY_p3loadlibrary(wfn,loadmsg);
    if (XCOMPRESS_zlibhandle != ValueCast(SYSTEM_pointer,0) && 
      _P3strcmpE(loadmsg,_P3str1("\000"))) {
      PointerCast(SYSTEM_pointer,&XCOMPRESS_pcompress) = loadentry(_P3str1("\010compress"),
        wfn,loadmsg);
      PointerCast(SYSTEM_pointer,&XCOMPRESS_puncompress) = loadentry(_P3str1("\012uncompress"),
        wfn,loadmsg);
      PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzreadopen) = loadentry(_P3str1("\006gzopen"),
        wfn,loadmsg);
      PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzread) = loadentry(_P3str1("\006gzread"),
        wfn,loadmsg);
      PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzreadclose) = loadentry(_P3str1("\007gzclose"),
        wfn,loadmsg);
    } 
  } 
  if (_P3strcmpN(loadmsg,_P3str1("\000"))) {
    PointerCast(SYSTEM_pointer,&XCOMPRESS_pcompress) = NULL;
    PointerCast(SYSTEM_pointer,&XCOMPRESS_puncompress) = NULL;
    PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzreadopen) = NULL;
    PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzread) = NULL;
    PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzreadclose) = NULL;
  } 
  result = ValueCast(SYSTEM_pointer,XCOMPRESS_pcompress) != NULL;
  return result;
}  /* loadzliblibrary */

Procedure XCOMPRESS_unloadzliblibrary(void)
{
  if (XCOMPRESS_zlibhandle != ValueCast(SYSTEM_pointer,0)) {
    P3LIBRARY_p3freelibrary(XCOMPRESS_zlibhandle);
    XCOMPRESS_zlibhandle = ValueCast(SYSTEM_pointer,0);
  } 
  PointerCast(SYSTEM_pointer,&XCOMPRESS_pcompress) = NULL;
  PointerCast(SYSTEM_pointer,&XCOMPRESS_puncompress) = NULL;
  PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzreadopen) = NULL;
  PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzread) = NULL;
  PointerCast(SYSTEM_pointer,&XCOMPRESS_pgzreadclose) = NULL;
}  /* unloadzliblibrary */

Function(SYSTEM_integer ) XCOMPRESS_compress(
  SYSTEM_pointer pdest,
  CLIBTYPES_clib_ulong *ldest,
  SYSTEM_pointer psrc,
  CLIBTYPES_clib_ulong lsrc)
{
  SYSTEM_integer result;

  result = (*XCOMPRESS_pcompress)(pdest,ldest,psrc,lsrc);
  return result;
}  /* compress */

Function(SYSTEM_integer ) XCOMPRESS_uncompress(
  SYSTEM_pointer pdest,
  CLIBTYPES_clib_ulong *ldest,
  SYSTEM_pointer psrc,
  CLIBTYPES_clib_ulong lsrc)
{
  SYSTEM_integer result;

  result = (*XCOMPRESS_puncompress)(pdest,ldest,psrc,lsrc);
  return result;
}  /* uncompress */

Function(XCOMPRESS_pgzfile ) XCOMPRESS_gzreadopen(
  const SYSTEM_ansichar *fn)
{
  XCOMPRESS_pgzfile result;
  SYSTEM_shortstring sfn;
  SYSTEM_shortstring smode;

  _P3strcat(sfn,255,fn,_P3str1("\001\000"));
  _P3strcpy(smode,255,_P3str1("\003rb\000"));
  result = (*XCOMPRESS_pgzreadopen)(ValueCast(SYSTEM_P3_pansichar,&sfn[1]),ValueCast(
    SYSTEM_P3_pansichar,&smode[1]));
  return result;
}  /* gzreadopen */

Function(SYSTEM_integer ) XCOMPRESS_gzread(
  XCOMPRESS_pgzfile pgz,
  SYSTEM_untyped *buf,
  SYSTEM_longword ldest)
{
  SYSTEM_integer result;

  result = (*XCOMPRESS_pgzread)(pgz,ValueCast(SYSTEM_pointer,buf),
    ldest);
  return result;
}  /* gzread */

Function(SYSTEM_integer ) XCOMPRESS_gzreadclose(
  XCOMPRESS_pgzfile *pgz)
{
  SYSTEM_integer result;

  result = (*XCOMPRESS_pgzreadclose)(*pgz);
  *pgz = NULL;
  return result;
}  /* gzreadclose */

/* unit xcompress */
void _Init_Module_xcompress(void)
{
  XCOMPRESS_zlibhandle = ValueCast(SYSTEM_pointer,0);
  XCOMPRESS_unloadzliblibrary();
} /* _Init_Module_xcompress */

void _Final_Module_xcompress(void)
{
} /* _Final_Module_xcompress */

