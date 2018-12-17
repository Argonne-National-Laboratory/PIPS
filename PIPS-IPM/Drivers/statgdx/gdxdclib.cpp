#define _P3_LIBRARY /* dll/shared object file*/

#include "p3io.h"
#include "p3platform.h"
#include "system_p3.h"
#include "p3private.h"
#include "exceptions.h"
#include "sysutils_p3.h"
#include "math_p3.h"
#include "p3library.h"
#include "p3utils.h"
#include "p3process.h"
#include "p3ieeefp.h"
#include "p3threads.h"
#include "idglobal_p3.h"
#include "gmsspecs.h"
#include "gmsgen.h"
#include "pchutil.h"
#include "gxdefs.h"
#include "gmsglob.h"
#include "strutilx.h"
#include "gmsobj.h"
#include "gmsdata.h"
#include "clibtypes.h"
#include "gmslibname.h"
#include "xcompress.h"
#include "gmsstrm.h"
#include "strhash.h"
#include "gmsheapnew.h"
#include "datastorage.h"
#include "gmsglobx.h"
#include "paldoorg.h"
#include "gdlaudit.h"
#include "runner.h"
#include "gxfile.h"


Extern_C _P3_DllExport Procedure  STDCALL xcreate(
    SYSTEM_pointer *pgdx)
{
    SYSTEM_shortstring msg;

    *pgdx = ValueCast(GXFILE_tgxfileobj,GXFILE_tgxfileobj_DOT_create(ValueCast(
      GXFILE_tgxfileobj,_P3alloc_object(&GXFILE_tgxfileobj_CD)),msg));
}  /* xcreate */

static Function(SYSTEM_integer )  STDCALL GDXDCLIB_checkp(
    SYSTEM_pointer p,
    SYSTEM_P3_pansichar msgbuf,
    SYSTEM_integer msgbuflen)
{
    SYSTEM_integer result;

    if (NULL == p) {
      PCHUTIL_strpcopyn(msgbuf,_P3str1("\033Error while creating object"),
        msgbuflen);
      result = 0;
    } else {
      (*ValueCast(GMSGEN_pansichararray,msgbuf))[0] = _P3char('\000');
      result = 1;
    } 
    return result;
}  /* checkp */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxcreate(
    SYSTEM_pointer *pgdx,
    SYSTEM_P3_pansichar msgbuf,
    SYSTEM_integer msgbuflen)
{
    SYSTEM_integer result;

    xcreate(pgdx);
    result = GDXDCLIB_checkp(*pgdx,msgbuf,msgbuflen);
    return result;
}  /* gdxcreate */

Extern_C _P3_DllExport Procedure  STDCALL gdxxcreate(
    SYSTEM_pointer *pgdx)
{
    xcreate(pgdx);
}  /* gdxxcreate */

Extern_C _P3_DllExport Procedure  STDCALL x2free(
    SYSTEM_pointer *pgdx)
{
    if (*pgdx != NULL) {
      SYSTEM_tobject_DOT_free(ValueCast(SYSTEM_tobject,*pgdx));
      *pgdx = NULL;
    } 
}  /* x2free */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxfree(
    SYSTEM_pointer *pgdx)
{
    SYSTEM_integer result;

    x2free(pgdx);
    if (NULL == *pgdx) { 
      result = 1;
    } else 
      result = 0;
    return result;
}  /* gdxfree */

Extern_C _P3_DllExport Procedure  STDCALL gdxxfree(
    SYSTEM_pointer *pgdx)
{
    x2free(pgdx);
}  /* gdxxfree */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL xapiversion(
    SYSTEM_integer api,
    SYSTEM_ansichar *msg,
    SYSTEM_integer *comp)
{
    SYSTEM_integer result;

    result = 0;
    *comp = 0;
    {
      SYSTEM_shortstring _t1;
      _P3STR_255 _t2;

      _P3strcat(msg,255,_P3strcat(_t2,255,_P3str1("\100gdxdclib: The API is too old for the used library, API version: "),
        SYSUTILS_P3_inttostr(_t1,255,api)),_P3str1("\024, library version: 7"));
    }
    if (api >= 7) {
      result = 1;
      if (api == 7) {
        *comp = 1;
        _P3strcpy(msg,255,_P3str1("\067gdxdclib: API version and library version are the same."));
      } else {
        *comp = 3;
        _P3strcpy(msg,255,_P3str1("\061gdxdclib: API version is newer than this library."));
      } 
      return result;
    } 
    if (api == 7) {
      result = 1;
      *comp = 2;
      _P3strcpy(msg,255,_P3str1("\106gdxdclib: Client version is compatible to this version of the library."));
      return result;
    } 
    return result;
}  /* xapiversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cxapiversion(
    SYSTEM_integer api,
    SYSTEM_P3_pansichar msg,
    SYSTEM_integer *comp)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_msg;

    result = xapiversion(api,local_msg,comp);
    SYSUTILS_P3_strpcopy(msg,local_msg);
    return result;
}  /* cxapiversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxxapiversion(
    SYSTEM_integer api,
    SYSTEM_ansichar *msg,
    SYSTEM_integer *comp)
{
    SYSTEM_integer result;

    result = xapiversion(api,msg,comp);
    return result;
}  /* gdxxapiversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL xcheck(
    const SYSTEM_ansichar *funcn,
    SYSTEM_integer clnrarg,
    SYSTEM_P3_pintegerarray clsign,
    SYSTEM_ansichar *msg)
{
    SYSTEM_integer result;
    SYSTEM_integer i, entryindex;
    cnstdef {maxarg = 6};
    cnstdef {maxentry = 87};
    cnstdef {maxlen = 24};
    typedef SYSTEM_uint8 _sub_1XCHECK;
    typedef _P3STR_31 _arr_0XCHECK[87];
    static _arr_0XCHECK entryname = {{13,'g','d','x','A','c','r','o','n','y','m','A','d','d'}, {15,'g','d','x','A','c','r','o','n','y','m','C','o','u','n','t'}, {17,'g','d','x','A','c','r','o','n','y','m','G','e','t','I','n','f','o'}, {20,'g','d','x','A','c','r','o','n','y','m','G','e','t','M','a','p','p','i','n','g'}, {15,'g','d','x','A','c','r','o','n','y','m','I','n','d','e','x'}, {14,'g','d','x','A','c','r','o','n','y','m','N','a','m','e'}, {16,'g','d','x','A','c','r','o','n','y','m','N','e','x','t','N','r'}, {17,'g','d','x','A','c','r','o','n','y','m','S','e','t','I','n','f','o'}, {15,'g','d','x','A','c','r','o','n','y','m','V','a','l','u','e'}, {11,'g','d','x','A','d','d','A','l','i','a','s'}, {13,'g','d','x','A','d','d','S','e','t','T','e','x','t'}, {14,'g','d','x','A','u','t','o','C','o','n','v','e','r','t'}, {8,'g','d','x','C','l','o','s','e'}, {17,'g','d','x','D','a','t','a','E','r','r','o','r','C','o','u','n','t'}, {18,'g','d','x','D','a','t','a','E','r','r','o','r','R','e','c','o','r','d'}, {15,'g','d','x','D','a','t','a','R','e','a','d','D','o','n','e'}, {24,'g','d','x','D','a','t','a','R','e','a','d','F','i','l','t','e','r','e','d','S','t','a','r','t'}, {14,'g','d','x','D','a','t','a','R','e','a','d','M','a','p'}, {19,'g','d','x','D','a','t','a','R','e','a','d','M','a','p','S','t','a','r','t'}, {14,'g','d','x','D','a','t','a','R','e','a','d','R','a','w'}, {18,'g','d','x','D','a','t','a','R','e','a','d','R','a','w','F','a','s','t'}, {22,'g','d','x','D','a','t','a','R','e','a','d','R','a','w','F','a','s','t','F','i','l','t'}, {19,'g','d','x','D','a','t','a','R','e','a','d','R','a','w','S','t','a','r','t'}, {16,'g','d','x','D','a','t','a','R','e','a','d','S','l','i','c','e'}, {21,'g','d','x','D','a','t','a','R','e','a','d','S','l','i','c','e','S','t','a','r','t'}, {14,'g','d','x','D','a','t','a','R','e','a','d','S','t','r'}, {19,'g','d','x','D','a','t','a','R','e','a','d','S','t','r','S','t','a','r','t'}, {16,'g','d','x','D','a','t','a','S','l','i','c','e','U','E','L','S'}, {16,'g','d','x','D','a','t','a','W','r','i','t','e','D','o','n','e'}, {15,'g','d','x','D','a','t','a','W','r','i','t','e','M','a','p'}, {20,'g','d','x','D','a','t','a','W','r','i','t','e','M','a','p','S','t','a','r','t'}, {15,'g','d','x','D','a','t','a','W','r','i','t','e','R','a','w'}, {20,'g','d','x','D','a','t','a','W','r','i','t','e','R','a','w','S','t','a','r','t'}, {15,'g','d','x','D','a','t','a','W','r','i','t','e','S','t','r'}, {20,'g','d','x','D','a','t','a','W','r','i','t','e','S','t','r','S','t','a','r','t'}, {16,'g','d','x','G','e','t','D','L','L','V','e','r','s','i','o','n'}, {13,'g','d','x','E','r','r','o','r','C','o','u','n','t'}, {11,'g','d','x','E','r','r','o','r','S','t','r'}, {11,'g','d','x','F','i','l','e','I','n','f','o'}, {14,'g','d','x','F','i','l','e','V','e','r','s','i','o','n'}, {15,'g','d','x','F','i','l','t','e','r','E','x','i','s','t','s'}, {17,'g','d','x','F','i','l','t','e','r','R','e','g','i','s','t','e','r'}, {21,'g','d','x','F','i','l','t','e','r','R','e','g','i','s','t','e','r','D','o','n','e'}, {22,'g','d','x','F','i','l','t','e','r','R','e','g','i','s','t','e','r','S','t','a','r','t'}, {13,'g','d','x','F','i','n','d','S','y','m','b','o','l'}, {14,'g','d','x','G','e','t','E','l','e','m','T','e','x','t'}, {15,'g','d','x','G','e','t','L','a','s','t','E','r','r','o','r'}, {16,'g','d','x','G','e','t','M','e','m','o','r','y','U','s','e','d'}, {19,'g','d','x','G','e','t','S','p','e','c','i','a','l','V','a','l','u','e','s'}, {9,'g','d','x','G','e','t','U','E','L'}, {11,'g','d','x','M','a','p','V','a','l','u','e'}, {13,'g','d','x','O','p','e','n','A','p','p','e','n','d'}, {11,'g','d','x','O','p','e','n','R','e','a','d'}, {12,'g','d','x','O','p','e','n','W','r','i','t','e'}, {14,'g','d','x','O','p','e','n','W','r','i','t','e','E','x'}, {21,'g','d','x','R','e','s','e','t','S','p','e','c','i','a','l','V','a','l','u','e','s'}, {13,'g','d','x','S','e','t','H','a','s','T','e','x','t'}, {23,'g','d','x','S','e','t','R','e','a','d','S','p','e','c','i','a','l','V','a','l','u','e','s'}, {19,'g','d','x','S','e','t','S','p','e','c','i','a','l','V','a','l','u','e','s'}, {16,'g','d','x','S','e','t','T','e','x','t','N','o','d','e','N','r'}, {16,'g','d','x','S','e','t','T','r','a','c','e','L','e','v','e','l'}, {20,'g','d','x','S','y','m','b','I','n','d','x','M','a','x','L','e','n','g','t','h'}, {16,'g','d','x','S','y','m','b','M','a','x','L','e','n','g','t','h'}, {19,'g','d','x','S','y','m','b','o','l','A','d','d','C','o','m','m','e','n','t'}, {19,'g','d','x','S','y','m','b','o','l','G','e','t','C','o','m','m','e','n','t'}, {18,'g','d','x','S','y','m','b','o','l','G','e','t','D','o','m','a','i','n'}, {19,'g','d','x','S','y','m','b','o','l','G','e','t','D','o','m','a','i','n','X'}, {12,'g','d','x','S','y','m','b','o','l','D','i','m'}, {13,'g','d','x','S','y','m','b','o','l','I','n','f','o'}, {14,'g','d','x','S','y','m','b','o','l','I','n','f','o','X'}, {18,'g','d','x','S','y','m','b','o','l','S','e','t','D','o','m','a','i','n'}, {19,'g','d','x','S','y','m','b','o','l','S','e','t','D','o','m','a','i','n','X'}, {13,'g','d','x','S','y','s','t','e','m','I','n','f','o'}, {15,'g','d','x','U','E','L','M','a','x','L','e','n','g','t','h'}, {18,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','D','o','n','e'}, {17,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','M','a','p'}, {22,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','M','a','p','S','t','a','r','t'}, {17,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','R','a','w'}, {22,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','R','a','w','S','t','a','r','t'}, {17,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','S','t','r'}, {22,'g','d','x','U','E','L','R','e','g','i','s','t','e','r','S','t','r','S','t','a','r','t'}, {12,'g','d','x','U','M','F','i','n','d','U','E','L'}, {11,'g','d','x','U','M','U','e','l','G','e','t'}, {12,'g','d','x','U','M','U','e','l','I','n','f','o'}, {20,'g','d','x','G','e','t','D','o','m','a','i','n','E','l','e','m','e','n','t','s'}, {13,'g','d','x','C','u','r','r','e','n','t','D','i','m'}, {12,'g','d','x','R','e','n','a','m','e','U','E','L'}};
    typedef SYSTEM_uint8 _sub_3XCHECK;
    typedef SYSTEM_integer _arr_2XCHECK[87];
    static _arr_2XCHECK dllnrarg = {3, 0, 4, 4, 1, 2, 1, 4, 1, 2, 2, 1, 0, 0, 3, 0, 3, 4, 2, 3, 3, 3, 2, 3, 2, 3, 2, 2, 0, 2, 5, 2, 5, 2, 5, 1, 0, 2, 2, 2, 1, 1, 0, 1, 2, 3, 0, 0, 1, 2, 2, 3, 2, 3, 4, 0, 1, 1, 1, 2, 2, 2, 0, 2, 3, 2, 2, 1, 4, 4, 1, 2, 2, 0, 0, 2, 0, 1, 0, 2, 0, 3, 3, 2, 6, 0, 2};
    typedef SYSTEM_uint8 _sub_5XCHECK;
    typedef SYSTEM_uint8 _sub_7XCHECK;
    typedef SYSTEM_integer _arr_6XCHECK[7];
    typedef _arr_6XCHECK _arr_4XCHECK[87];
    static _arr_4XCHECK dllsign = {{3, 11, 11, 3,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 12, 12, 4,  -1,  -1}, {3, 3, 4, 4, 4,  -1,  -1}, 
      {3, 13,  -1,  -1,  -1,  -1,  -1}, {3, 13, 12,  -1,  -1,  -1,  -1}, 
      {3, 3,  -1,  -1,  -1,  -1,  -1}, {3, 3, 11, 11, 3,  -1,  -1}, 
      {13, 3,  -1,  -1,  -1,  -1,  -1}, {3, 11, 11,  -1,  -1,  -1,  -1}, 
      {3, 11, 4,  -1,  -1,  -1,  -1}, {3, 3,  -1,  -1,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 52, 54,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 51, 4,  -1,  -1,  -1}, {3, 3, 52, 54, 4,  -1,  -1}, 
      {3, 3, 4,  -1,  -1,  -1,  -1}, {3, 52, 54, 4,  -1,  -1,  -1}, 
      {3, 3, 59, 4,  -1,  -1,  -1}, {3, 3, 55, 59,  -1,  -1,  -1}, 
      {3, 3, 4,  -1,  -1,  -1,  -1}, {3, 55, 4, 59,  -1,  -1,  -1}, 
      {3, 3, 52,  -1,  -1,  -1,  -1}, {3, 56, 54, 4,  -1,  -1,  -1}, 
      {3, 3, 4,  -1,  -1,  -1,  -1}, {3, 51, 56,  -1,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 51, 53,  -1,  -1,  -1,  -1}, {3, 11, 11, 3, 3, 3,  -1}, 
      {3, 51, 53,  -1,  -1,  -1,  -1}, {3, 11, 11, 3, 3, 3,  -1}, 
      {3, 55, 53,  -1,  -1,  -1,  -1}, {3, 11, 11, 3, 3, 3,  -1}, 
      {3, 12,  -1,  -1,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 12,  -1,  -1,  -1,  -1}, {3, 4, 4,  -1,  -1,  -1,  -1}, 
      {3, 12, 12,  -1,  -1,  -1,  -1}, {3, 3,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3,  -1,  -1,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3,  -1,  -1,  -1,  -1,  -1}, {3, 11, 4,  -1,  -1,  -1,  -1}, 
      {3, 3, 12, 4,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {23,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 58,  -1,  -1,  -1,  -1,  -1}, {3, 3, 12,  -1,  -1,  -1,  -1}, 
      {3, 13, 4,  -1,  -1,  -1,  -1}, {3, 11, 11, 4,  -1,  -1,  -1}, 
      {3, 11, 4,  -1,  -1,  -1,  -1}, {3, 11, 11, 4,  -1,  -1,  -1}, 
      {3, 11, 11, 3, 4,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3,  -1,  -1,  -1,  -1,  -1}, {3, 57,  -1,  -1,  -1,  -1,  -1}, 
      {3, 57,  -1,  -1,  -1,  -1,  -1}, {3, 3, 3,  -1,  -1,  -1,  -1}, 
      {3, 3, 11,  -1,  -1,  -1,  -1}, {3, 3, 52,  -1,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 11,  -1,  -1,  -1,  -1}, {3, 3, 3, 12,  -1,  -1,  -1}, 
      {3, 3, 52,  -1,  -1,  -1,  -1}, {3, 3, 56,  -1,  -1,  -1,  -1}, 
      {3, 3,  -1,  -1,  -1,  -1,  -1}, {3, 3, 12, 4, 4,  -1,  -1}, 
      {3, 3, 4, 4, 12,  -1,  -1}, {3, 55,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 55,  -1,  -1,  -1,  -1}, {3, 4, 4,  -1,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 3, 11,  -1,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 11,  -1,  -1,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 11, 4,  -1,  -1,  -1,  -1}, {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 11, 4, 4,  -1,  -1,  -1}, {3, 3, 12, 4,  -1,  -1,  -1}, 
      {3, 4, 4,  -1,  -1,  -1,  -1}, {3, 3, 3, 3, 59, 4, 1}, 
      {3,  -1,  -1,  -1,  -1,  -1,  -1}, 
      {3, 11, 11,  -1,  -1,  -1,  -1}};

    entryindex = 0;
    for (i = 1;i <= (SYSTEM_int32)maxentry;++i) {
      if (_P3strcmpE(entryname[i - 1],funcn)) {
        entryindex = i;
        SYSTEM_break(BRK_1);
      } 
    CNT_1:;
}
BRK_1:;
    if (entryindex == 0) {
      result = 0;
      {
        _P3STR_255 _t1;

        _P3strcat(msg,255,_P3strcat(_t1,255,_P3str1("\012gdxdclib: "),
          funcn),_P3str1("\040 cannot be found in the library."));
      }
      return result;
    } 
    result = 1;
    _P3strclr(msg);
    if (dllnrarg[entryindex - 1] != clnrarg) {
      result = 0;
      {
        _P3STR_255 _t1;
        _P3STR_255 _t2;
        SYSTEM_shortstring _t3;
        _P3STR_255 _t4;
        _P3STR_255 _t5;
        SYSTEM_shortstring _t6;
        _P3STR_255 _t7;

        _P3strcat(msg,255,_P3strcat(_t7,255,_P3strcat(_t5,255,
          _P3strcat(_t4,255,_P3strcat(_t2,255,_P3strcat(_t1,255,_P3str1("\012gdxdclib: "),
          funcn),_P3str1("\060 has wrong number of arguments, the API expects ")),
          SYSUTILS_P3_inttostr(_t3,255,clnrarg)),_P3str1("\014 but it has ")),
          SYSUTILS_P3_inttostr(_t6,255,dllnrarg[entryindex - 1])),_P3str1("\020 in the library."));
      }
    } else 
      { register SYSTEM_int32 _stop = dllnrarg[entryindex - 1];
        if ((i = 0) <=  _stop) do {
          if (dllsign[entryindex - 1][i] != (*clsign)[i]) {
            result = 0;
            if (_P3strcmpE(msg,_P3str1("\000"))) { 
              {
                _P3STR_255 _t1;
                _P3STR_255 _t2;
                SYSTEM_shortstring _t3;

                _P3strcat(msg,255,_P3strcat(_t2,255,_P3strcat(
                  _t1,255,_P3str1("\012gdxdclib: "),funcn),_P3str1("\046 has wrong argument type for argument ")),
                  SYSUTILS_P3_inttostr(_t3,255,i));
              }
            } else 
              {
                _P3STR_255 _t1;
                SYSTEM_shortstring _t2;

                _P3strcat(msg,255,_P3strcat(_t1,255,msg,_P3str1("\001,")),
                  SYSUTILS_P3_inttostr(_t2,255,i));
              }
          } 
        } while (i++ !=  _stop);

      }
    return result;
}  /* xcheck */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cxcheck(
    SYSTEM_P3_pansichar funcn,
    SYSTEM_integer clnrarg,
    SYSTEM_P3_pintegerarray clsign,
    SYSTEM_P3_pansichar msg)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_msg;

    {
      SYSTEM_shortstring _t1;

      result = xcheck(PCHUTIL_pchartostr(_t1,255,funcn),clnrarg,
        clsign,local_msg);
    }
    SYSUTILS_P3_strpcopy(msg,local_msg);
    return result;
}  /* cxcheck */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxxcheck(
    const SYSTEM_ansichar *funcn,
    SYSTEM_integer clnrarg,
    SYSTEM_P3_pintegerarray clsign,
    SYSTEM_ansichar *msg)
{
    SYSTEM_integer result;

    result = xcheck(funcn,clnrarg,clsign,msg);
    return result;
}  /* gdxxcheck */

Extern_C _P3_DllExport Procedure  STDCALL gdxsetloadpath(
    const SYSTEM_ansichar *s)
{
    _P3strcpy(GXFILE_dllloadpath,255,s);
}  /* gdxsetloadpath */

Extern_C _P3_DllExport Procedure  STDCALL cgdxsetloadpath(
    SYSTEM_P3_pansichar ps)
{
    {
      SYSTEM_shortstring _t1;

      gdxsetloadpath(PCHUTIL_pchartostr(_t1,255,ps));
    }
}  /* cgdxsetloadpath */

Extern_C _P3_DllExport Procedure  STDCALL gdxgetloadpath(
    SYSTEM_ansichar *s)
{
    _P3strcpy(s,255,GXFILE_dllloadpath);
}  /* gdxgetloadpath */

Extern_C _P3_DllExport Procedure  STDCALL cgdxgetloadpath(
    SYSTEM_P3_pansichar ps)
{
    SYSTEM_shortstring s;

    gdxgetloadpath(s);
    SYSUTILS_P3_strpcopy(ps,s);
}  /* cgdxgetloadpath */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymadd(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *aname,
    const SYSTEM_ansichar *txt,
    SYSTEM_integer aindx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymadd(ValueCast(
      GXFILE_tgxfileobj,pgdx),aname,txt,aindx);
    return result;
}  /* gdxacronymadd */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymcount(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymcount(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxacronymcount */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymgetinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    SYSTEM_ansichar *aname,
    SYSTEM_ansichar *txt,
    SYSTEM_integer *aindx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymgetinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),n,aname,txt,aindx);
    return result;
}  /* gdxacronymgetinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymgetmapping(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    SYSTEM_integer *orgindx,
    SYSTEM_integer *newindx,
    SYSTEM_integer *autoindex)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymgetmapping(ValueCast(
      GXFILE_tgxfileobj,pgdx),n,orgindx,newindx,autoindex);
    return result;
}  /* gdxacronymgetmapping */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymindex(
    SYSTEM_pointer pgdx,
    SYSTEM_double v)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymindex(ValueCast(
      GXFILE_tgxfileobj,pgdx),v);
    return result;
}  /* gdxacronymindex */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymname(
    SYSTEM_pointer pgdx,
    SYSTEM_double v,
    SYSTEM_ansichar *aname)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymname(ValueCast(
      GXFILE_tgxfileobj,pgdx),v,aname);
    return result;
}  /* gdxacronymname */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymnextnr(
    SYSTEM_pointer pgdx,
    SYSTEM_integer nv)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymnextnr(ValueCast(
      GXFILE_tgxfileobj,pgdx),nv);
    return result;
}  /* gdxacronymnextnr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxacronymsetinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    const SYSTEM_ansichar *aname,
    const SYSTEM_ansichar *txt,
    SYSTEM_integer aindx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymsetinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),n,aname,txt,aindx);
    return result;
}  /* gdxacronymsetinfo */

Extern_C _P3_DllExport Function(SYSTEM_double )  STDCALL 
    gdxacronymvalue(
    SYSTEM_pointer pgdx,
    SYSTEM_integer aindx)
{
    SYSTEM_double result;

    result = GXFILE_tgxfileobj_DOT_gdxacronymvalue(ValueCast(
      GXFILE_tgxfileobj,pgdx),aindx);
    return result;
}  /* gdxacronymvalue */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxaddalias(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *id1,
    const SYSTEM_ansichar *id2)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxaddalias(ValueCast(
      GXFILE_tgxfileobj,pgdx),id1,id2);
    return result;
}  /* gdxaddalias */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxaddsettext(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *txt,
    SYSTEM_integer *txtnr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxaddsettext(ValueCast(
      GXFILE_tgxfileobj,pgdx),txt,txtnr);
    return result;
}  /* gdxaddsettext */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxautoconvert(
    SYSTEM_pointer pgdx,
    SYSTEM_integer nv)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxautoconvert(ValueCast(
      GXFILE_tgxfileobj,pgdx),nv);
    return result;
}  /* gdxautoconvert */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxclose(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxclose(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxclose */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdataerrorcount(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdataerrorcount(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxdataerrorcount */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdataerrorrecord(
    SYSTEM_pointer pgdx,
    SYSTEM_integer recnr,
    SYSTEM_integer *keyint,
    SYSTEM_double *values)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdataerrorrecord(ValueCast(
      GXFILE_tgxfileobj,pgdx),recnr,keyint,values);
    return result;
}  /* gdxdataerrorrecord */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareaddone(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareaddone(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxdatareaddone */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadfilteredstart(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    const SYSTEM_integer *filteraction,
    SYSTEM_integer *nrrecs)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadfilteredstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,filteraction,nrrecs);
    return result;
}  /* gdxdatareadfilteredstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadmap(
    SYSTEM_pointer pgdx,
    SYSTEM_integer recnr,
    SYSTEM_integer *keyint,
    SYSTEM_double *values,
    SYSTEM_integer *dimfrst)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadmap(ValueCast(
      GXFILE_tgxfileobj,pgdx),recnr,keyint,values,dimfrst);
    return result;
}  /* gdxdatareadmap */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadmapstart(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *nrrecs)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadmapstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,nrrecs);
    return result;
}  /* gdxdatareadmapstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadraw(
    SYSTEM_pointer pgdx,
    SYSTEM_integer *keyint,
    SYSTEM_double *values,
    SYSTEM_integer *dimfrst)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadraw(ValueCast(
      GXFILE_tgxfileobj,pgdx),keyint,values,dimfrst);
    return result;
}  /* gdxdatareadraw */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadrawfast(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    GXDEFS_tdatastoreproc dp,
    SYSTEM_integer *nrrecs)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadrawfast(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,dp,nrrecs);
    return result;
}  /* gdxdatareadrawfast */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadrawfastfilt(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    const SYSTEM_shortstring *uelfilterstr,
    GXDEFS_tdatastorefiltproc dp)
{
    SYSTEM_integer result;

    (ValueCast(GXFILE_tgxfileobj,pgdx))->
      GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_callbyref = 
      SYSTEM_false;
    result = GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,uelfilterstr,dp);
    return result;
}  /* gdxdatareadrawfastfilt */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadrawstart(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *nrrecs)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadrawstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,nrrecs);
    return result;
}  /* gdxdatareadrawstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadslice(
    SYSTEM_pointer pgdx,
    const SYSTEM_shortstring *uelfilterstr,
    SYSTEM_integer *dimen,
    GXDEFS_tdatastoreproc dp)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadslice(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelfilterstr,dimen,dp);
    return result;
}  /* gdxdatareadslice */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadslicestart(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *elemcounts)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadslicestart(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,elemcounts);
    return result;
}  /* gdxdatareadslicestart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadstr(
    SYSTEM_pointer pgdx,
    SYSTEM_shortstring *keystr,
    SYSTEM_double *values,
    SYSTEM_integer *dimfrst)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadstr(ValueCast(
      GXFILE_tgxfileobj,pgdx),keystr,values,dimfrst);
    return result;
}  /* gdxdatareadstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatareadstrstart(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *nrrecs)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadstrstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,nrrecs);
    return result;
}  /* gdxdatareadstrstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatasliceuels(
    SYSTEM_pointer pgdx,
    const SYSTEM_integer *slicekeyint,
    SYSTEM_shortstring *keystr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatasliceuels(ValueCast(
      GXFILE_tgxfileobj,pgdx),slicekeyint,keystr);
    return result;
}  /* gdxdatasliceuels */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawritedone(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawritedone(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxdatawritedone */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawritemap(
    SYSTEM_pointer pgdx,
    const SYSTEM_integer *keyint,
    const SYSTEM_double *values)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawritemap(ValueCast(
      GXFILE_tgxfileobj,pgdx),keyint,values);
    return result;
}  /* gdxdatawritemap */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawritemapstart(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *syid,
    const SYSTEM_ansichar *expltxt,
    SYSTEM_integer dimen,
    SYSTEM_integer typ,
    SYSTEM_integer userinfo)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawritemapstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),syid,expltxt,dimen,typ,userinfo);
    return result;
}  /* gdxdatawritemapstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawriteraw(
    SYSTEM_pointer pgdx,
    const SYSTEM_integer *keyint,
    const SYSTEM_double *values)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawriteraw(ValueCast(
      GXFILE_tgxfileobj,pgdx),keyint,values);
    return result;
}  /* gdxdatawriteraw */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawriterawstart(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *syid,
    const SYSTEM_ansichar *expltxt,
    SYSTEM_integer dimen,
    SYSTEM_integer typ,
    SYSTEM_integer userinfo)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawriterawstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),syid,expltxt,dimen,typ,userinfo);
    return result;
}  /* gdxdatawriterawstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawritestr(
    SYSTEM_pointer pgdx,
    const SYSTEM_shortstring *keystr,
    const SYSTEM_double *values)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawritestr(ValueCast(
      GXFILE_tgxfileobj,pgdx),keystr,values);
    return result;
}  /* gdxdatawritestr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxdatawritestrstart(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *syid,
    const SYSTEM_ansichar *expltxt,
    SYSTEM_integer dimen,
    SYSTEM_integer typ,
    SYSTEM_integer userinfo)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxdatawritestrstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),syid,expltxt,dimen,typ,userinfo);
    return result;
}  /* gdxdatawritestrstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxgetdllversion(
    SYSTEM_pointer pgdx,
    SYSTEM_ansichar *v)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxgetdllversion(ValueCast(
      GXFILE_tgxfileobj,pgdx),v);
    return result;
}  /* gdxgetdllversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxerrorcount(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxerrorcount(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxerrorcount */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxerrorstr(
    SYSTEM_pointer pgdx,
    SYSTEM_integer errnr,
    SYSTEM_ansichar *errmsg)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxerrorstr(ValueCast(
      GXFILE_tgxfileobj,pgdx),errnr,errmsg);
    return result;
}  /* gdxerrorstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxfileinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer *filever,
    SYSTEM_integer *comprlev)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfileinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),filever,comprlev);
    return result;
}  /* gdxfileinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxfileversion(
    SYSTEM_pointer pgdx,
    SYSTEM_ansichar *filestr,
    SYSTEM_ansichar *producestr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfileversion(ValueCast(
      GXFILE_tgxfileobj,pgdx),filestr,producestr);
    return result;
}  /* gdxfileversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxfilterexists(
    SYSTEM_pointer pgdx,
    SYSTEM_integer filternr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfilterexists(ValueCast(
      GXFILE_tgxfileobj,pgdx),filternr);
    return result;
}  /* gdxfilterexists */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxfilterregister(
    SYSTEM_pointer pgdx,
    SYSTEM_integer uelmap)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfilterregister(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelmap);
    return result;
}  /* gdxfilterregister */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxfilterregisterdone(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfilterregisterdone(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxfilterregisterdone */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxfilterregisterstart(
    SYSTEM_pointer pgdx,
    SYSTEM_integer filternr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfilterregisterstart(ValueCast(
      GXFILE_tgxfileobj,pgdx),filternr);
    return result;
}  /* gdxfilterregisterstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxfindsymbol(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *syid,
    SYSTEM_integer *synr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxfindsymbol(ValueCast(
      GXFILE_tgxfileobj,pgdx),syid,synr);
    return result;
}  /* gdxfindsymbol */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxgetelemtext(
    SYSTEM_pointer pgdx,
    SYSTEM_integer txtnr,
    SYSTEM_ansichar *txt,
    SYSTEM_integer *node)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxgetelemtext(ValueCast(
      GXFILE_tgxfileobj,pgdx),txtnr,txt,node);
    return result;
}  /* gdxgetelemtext */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxgetlasterror(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxgetlasterror(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxgetlasterror */

Extern_C _P3_DllExport Function(SYSTEM_int64 )  STDCALL 
    gdxgetmemoryused(
    SYSTEM_pointer pgdx)
{
    SYSTEM_int64 result;

    result = GXFILE_tgxfileobj_DOT_gdxgetmemoryused(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxgetmemoryused */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxgetspecialvalues(
    SYSTEM_pointer pgdx,
    SYSTEM_double *avals)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxgetspecialvalues(ValueCast(
      GXFILE_tgxfileobj,pgdx),avals);
    return result;
}  /* gdxgetspecialvalues */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxgetuel(
    SYSTEM_pointer pgdx,
    SYSTEM_integer uelnr,
    SYSTEM_ansichar *uel)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxgetuel(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelnr,uel);
    return result;
}  /* gdxgetuel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxmapvalue(
    SYSTEM_pointer pgdx,
    SYSTEM_double d,
    SYSTEM_integer *sv)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxmapvalue(ValueCast(
      GXFILE_tgxfileobj,pgdx),d,sv);
    return result;
}  /* gdxmapvalue */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxopenappend(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *filename,
    const SYSTEM_ansichar *producer,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxopenappend(ValueCast(
      GXFILE_tgxfileobj,pgdx),filename,producer,errnr);
    return result;
}  /* gdxopenappend */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxopenread(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *filename,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxopenread(ValueCast(
      GXFILE_tgxfileobj,pgdx),filename,errnr);
    return result;
}  /* gdxopenread */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxopenwrite(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *filename,
    const SYSTEM_ansichar *producer,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxopenwrite(ValueCast(
      GXFILE_tgxfileobj,pgdx),filename,producer,errnr);
    return result;
}  /* gdxopenwrite */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxopenwriteex(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *filename,
    const SYSTEM_ansichar *producer,
    SYSTEM_integer compr,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxopenwriteex(ValueCast(
      GXFILE_tgxfileobj,pgdx),filename,producer,compr,errnr);
    return result;
}  /* gdxopenwriteex */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxresetspecialvalues(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxresetspecialvalues(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxresetspecialvalues */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsethastext(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsethastext(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr);
    return result;
}  /* gdxsethastext */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsetreadspecialvalues(
    SYSTEM_pointer pgdx,
    const SYSTEM_double *avals)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsetreadspecialvalues(ValueCast(
      GXFILE_tgxfileobj,pgdx),avals);
    return result;
}  /* gdxsetreadspecialvalues */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsetspecialvalues(
    SYSTEM_pointer pgdx,
    const SYSTEM_double *avals)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsetspecialvalues(ValueCast(
      GXFILE_tgxfileobj,pgdx),avals);
    return result;
}  /* gdxsetspecialvalues */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsettextnodenr(
    SYSTEM_pointer pgdx,
    SYSTEM_integer txtnr,
    SYSTEM_integer node)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsettextnodenr(ValueCast(
      GXFILE_tgxfileobj,pgdx),txtnr,node);
    return result;
}  /* gdxsettextnodenr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsettracelevel(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    const SYSTEM_ansichar *s)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsettracelevel(ValueCast(
      GXFILE_tgxfileobj,pgdx),n,s);
    return result;
}  /* gdxsettracelevel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbindxmaxlength(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *lengthinfo)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbindxmaxlength(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,lengthinfo);
    return result;
}  /* gdxsymbindxmaxlength */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbmaxlength(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbmaxlength(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxsymbmaxlength */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymboladdcomment(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    const SYSTEM_ansichar *txt)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymboladdcomment(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,txt);
    return result;
}  /* gdxsymboladdcomment */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolgetcomment(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer n,
    SYSTEM_ansichar *txt)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolgetcomment(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,n,txt);
    return result;
}  /* gdxsymbolgetcomment */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolgetdomain(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *domainsynrs)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolgetdomain(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,domainsynrs);
    return result;
}  /* gdxsymbolgetdomain */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolgetdomainx(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_shortstring *domainids)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolgetdomainx(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,domainids);
    return result;
}  /* gdxsymbolgetdomainx */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxsymboldim(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymboldim(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr);
    return result;
}  /* gdxsymboldim */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_ansichar *syid,
    SYSTEM_integer *dimen,
    SYSTEM_integer *typ)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,syid,dimen,typ);
    return result;
}  /* gdxsymbolinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolinfox(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *reccnt,
    SYSTEM_integer *userinfo,
    SYSTEM_ansichar *expltxt)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolinfox(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,reccnt,userinfo,expltxt);
    return result;
}  /* gdxsymbolinfox */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolsetdomain(
    SYSTEM_pointer pgdx,
    const SYSTEM_shortstring *domainids)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolsetdomain(ValueCast(
      GXFILE_tgxfileobj,pgdx),domainids);
    return result;
}  /* gdxsymbolsetdomain */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsymbolsetdomainx(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    const SYSTEM_shortstring *domainids)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolsetdomainx(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,domainids);
    return result;
}  /* gdxsymbolsetdomainx */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxsysteminfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer *sycnt,
    SYSTEM_integer *uelcnt)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxsysteminfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),sycnt,uelcnt);
    return result;
}  /* gdxsysteminfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelmaxlength(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelmaxlength(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxuelmaxlength */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregisterdone(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregisterdone(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxuelregisterdone */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregistermap(
    SYSTEM_pointer pgdx,
    SYSTEM_integer umap,
    const SYSTEM_ansichar *uel)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregistermap(ValueCast(
      GXFILE_tgxfileobj,pgdx),umap,uel);
    return result;
}  /* gdxuelregistermap */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregistermapstart(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregistermapstart(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxuelregistermapstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregisterraw(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *uel)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregisterraw(ValueCast(
      GXFILE_tgxfileobj,pgdx),uel);
    return result;
}  /* gdxuelregisterraw */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregisterrawstart(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregisterrawstart(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxuelregisterrawstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregisterstr(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *uel,
    SYSTEM_integer *uelnr)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregisterstr(ValueCast(
      GXFILE_tgxfileobj,pgdx),uel,uelnr);
    return result;
}  /* gdxuelregisterstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxuelregisterstrstart(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxuelregisterstrstart(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxuelregisterstrstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxumfinduel(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *uel,
    SYSTEM_integer *uelnr,
    SYSTEM_integer *uelmap)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxumfinduel(ValueCast(
      GXFILE_tgxfileobj,pgdx),uel,uelnr,uelmap);
    return result;
}  /* gdxumfinduel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxumuelget(
    SYSTEM_pointer pgdx,
    SYSTEM_integer uelnr,
    SYSTEM_ansichar *uel,
    SYSTEM_integer *uelmap)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxumuelget(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelnr,uel,uelmap);
    return result;
}  /* gdxumuelget */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxumuelinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer *uelcnt,
    SYSTEM_integer *highmap)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxumuelinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelcnt,highmap);
    return result;
}  /* gdxumuelinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxgetdomainelements(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer dimpos,
    SYSTEM_integer filternr,
    GXDEFS_tdomainindexproc dp,
    SYSTEM_integer *nrelem,
    SYSTEM_pointer uptr)
{
    SYSTEM_integer result;

    (ValueCast(GXFILE_tgxfileobj,pgdx))->
      GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_callbyref = 
      SYSTEM_false;
    result = GXFILE_tgxfileobj_DOT_gdxgetdomainelements(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,dimpos,filternr,dp,nrelem,uptr);
    return result;
}  /* gdxgetdomainelements */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    gdxcurrentdim(
    SYSTEM_pointer pgdx)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxcurrentdim(ValueCast(
      GXFILE_tgxfileobj,pgdx));
    return result;
}  /* gdxcurrentdim */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL gdxrenameuel(
    SYSTEM_pointer pgdx,
    const SYSTEM_ansichar *oldname,
    const SYSTEM_ansichar *newname)
{
    SYSTEM_integer result;

    result = GXFILE_tgxfileobj_DOT_gdxrenameuel(ValueCast(
      GXFILE_tgxfileobj,pgdx),oldname,newname);
    return result;
}  /* gdxrenameuel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    fgdxdatareadrawfastfilt(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    const SYSTEM_shortstring *uelfilterstr,
    GXDEFS_tdatastorefiltproc_f dp)
{
    SYSTEM_integer result;
    GXDEFS_tdatastorefiltproc local_dp;

    (ValueCast(GXFILE_tgxfileobj,pgdx))->
      GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_callbyref = 
      SYSTEM_true;
    PointerCast(SYSTEM_pointer,&local_dp) = ValueCast(SYSTEM_pointer,
      dp);
    result = GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,uelfilterstr,local_dp);
    return result;
}  /* fgdxdatareadrawfastfilt */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    fgdxgetdomainelements(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer dimpos,
    SYSTEM_integer filternr,
    GXDEFS_tdomainindexproc_f dp,
    SYSTEM_integer *nrelem,
    SYSTEM_pointer uptr)
{
    SYSTEM_integer result;
    GXDEFS_tdomainindexproc local_dp;

    (ValueCast(GXFILE_tgxfileobj,pgdx))->
      GXFILE_tgxfileobj_DOT_gdxgetdomainelements_dp_callbyref = 
      SYSTEM_true;
    PointerCast(SYSTEM_pointer,&local_dp) = ValueCast(SYSTEM_pointer,
      dp);
    result = GXFILE_tgxfileobj_DOT_gdxgetdomainelements(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,dimpos,filternr,local_dp,nrelem,
      uptr);
    return result;
}  /* fgdxgetdomainelements */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxacronymadd(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar aname,
    SYSTEM_P3_pansichar txt,
    SYSTEM_integer aindx)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxacronymadd(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,aname),
        PCHUTIL_pchartostr(_t2,255,txt),aindx);
    }
    return result;
}  /* cgdxacronymadd */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxacronymgetinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    SYSTEM_P3_pansichar aname,
    SYSTEM_P3_pansichar txt,
    SYSTEM_integer *aindx)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_aname;
    SYSTEM_shortstring local_txt;

    result = GXFILE_tgxfileobj_DOT_gdxacronymgetinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),n,local_aname,local_txt,aindx);
    SYSUTILS_P3_strpcopy(aname,local_aname);
    SYSUTILS_P3_strpcopy(txt,local_txt);
    return result;
}  /* cgdxacronymgetinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxacronymname(
    SYSTEM_pointer pgdx,
    SYSTEM_double v,
    SYSTEM_P3_pansichar aname)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_aname;

    result = GXFILE_tgxfileobj_DOT_gdxacronymname(ValueCast(
      GXFILE_tgxfileobj,pgdx),v,local_aname);
    SYSUTILS_P3_strpcopy(aname,local_aname);
    return result;
}  /* cgdxacronymname */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxacronymsetinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    SYSTEM_P3_pansichar aname,
    SYSTEM_P3_pansichar txt,
    SYSTEM_integer aindx)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxacronymsetinfo(ValueCast(
        GXFILE_tgxfileobj,pgdx),n,PCHUTIL_pchartostr(_t1,255,aname),
        PCHUTIL_pchartostr(_t2,255,txt),aindx);
    }
    return result;
}  /* cgdxacronymsetinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cgdxaddalias(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar id1,
    SYSTEM_P3_pansichar id2)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxaddalias(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,id1),
        PCHUTIL_pchartostr(_t2,255,id2));
    }
    return result;
}  /* cgdxaddalias */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxaddsettext(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar txt,
    SYSTEM_integer *txtnr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxaddsettext(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,txt),
        txtnr);
    }
    return result;
}  /* cgdxaddsettext */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatareadrawfastfilt(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_P3_ppointerarray uelfilterstr,
    GXDEFS_tdatastorefiltproc dp)
{
    SYSTEM_integer result;
    SYSTEM_integer uelfilterstr_i;
    GXDEFS_tgdxstrindex uelfilterstr_s;

    (ValueCast(GXFILE_tgxfileobj,pgdx))->
      GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt_dp_callbyref = 
      SYSTEM_false;
    { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,pgdx))->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((uelfilterstr_i = 1) <=  _stop) do {
        PCHUTIL_pchartostr(uelfilterstr_s[uelfilterstr_i - 1],255,ValueCast(
          SYSTEM_P3_pansichar,(*uelfilterstr)[uelfilterstr_i - 1]));
      } while (uelfilterstr_i++ !=  _stop);

    }
    result = GXFILE_tgxfileobj_DOT_gdxdatareadrawfastfilt(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,uelfilterstr_s,dp);
    return result;
}  /* cgdxdatareadrawfastfilt */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatareadslice(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_ppointerarray uelfilterstr,
    SYSTEM_integer *dimen,
    GXDEFS_tdatastoreproc dp)
{
    SYSTEM_integer result;
    SYSTEM_integer uelfilterstr_i;
    GXDEFS_tgdxstrindex uelfilterstr_s;

    { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,pgdx))->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((uelfilterstr_i = 1) <=  _stop) do {
        PCHUTIL_pchartostr(uelfilterstr_s[uelfilterstr_i - 1],255,ValueCast(
          SYSTEM_P3_pansichar,(*uelfilterstr)[uelfilterstr_i - 1]));
      } while (uelfilterstr_i++ !=  _stop);

    }
    result = GXFILE_tgxfileobj_DOT_gdxdatareadslice(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelfilterstr_s,dimen,dp);
    return result;
}  /* cgdxdatareadslice */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatareadstr(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_ppointerarray keystr,
    SYSTEM_double *values,
    SYSTEM_integer *dimfrst)
{
    SYSTEM_integer result;
    SYSTEM_integer keystr_i;
    GXDEFS_tgdxstrindex keystr_s;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadstr(ValueCast(
      GXFILE_tgxfileobj,pgdx),keystr_s,values,dimfrst);
    if (result != 0) 
      { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,
          pgdx))->GXFILE_tgxfileobj_DOT_fcurrentdim;
        if ((keystr_i = 1) <=  _stop) do {
          SYSUTILS_P3_strpcopy(ValueCast(SYSTEM_P3_pansichar,(*keystr)[
            keystr_i - 1]),keystr_s[keystr_i - 1]);
        } while (keystr_i++ !=  _stop);

      }
    return result;
}  /* cgdxdatareadstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatasliceuels(
    SYSTEM_pointer pgdx,
    const SYSTEM_integer *slicekeyint,
    SYSTEM_P3_ppointerarray keystr)
{
    SYSTEM_integer result;
    SYSTEM_integer keystr_i;
    GXDEFS_tgdxstrindex keystr_s;

    result = GXFILE_tgxfileobj_DOT_gdxdatasliceuels(ValueCast(
      GXFILE_tgxfileobj,pgdx),slicekeyint,keystr_s);
    if (result != 0) 
      { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,
          pgdx))->GXFILE_tgxfileobj_DOT_fcurrentdim;
        if ((keystr_i = 1) <=  _stop) do {
          SYSUTILS_P3_strpcopy(ValueCast(SYSTEM_P3_pansichar,(*keystr)[
            keystr_i - 1]),keystr_s[keystr_i - 1]);
        } while (keystr_i++ !=  _stop);

      }
    return result;
}  /* cgdxdatasliceuels */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatawritemapstart(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar syid,
    SYSTEM_P3_pansichar expltxt,
    SYSTEM_integer dimen,
    SYSTEM_integer typ,
    SYSTEM_integer userinfo)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxdatawritemapstart(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,syid),
        PCHUTIL_pchartostr(_t2,255,expltxt),dimen,typ,userinfo);
    }
    return result;
}  /* cgdxdatawritemapstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatawriterawstart(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar syid,
    SYSTEM_P3_pansichar expltxt,
    SYSTEM_integer dimen,
    SYSTEM_integer typ,
    SYSTEM_integer userinfo)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxdatawriterawstart(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,syid),
        PCHUTIL_pchartostr(_t2,255,expltxt),dimen,typ,userinfo);
    }
    return result;
}  /* cgdxdatawriterawstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatawritestr(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_ppointerarray keystr,
    const SYSTEM_double *values)
{
    SYSTEM_integer result;
    SYSTEM_integer keystr_i;
    GXDEFS_tgdxstrindex keystr_s;

    { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,pgdx))->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((keystr_i = 1) <=  _stop) do {
        PCHUTIL_pchartostr(keystr_s[keystr_i - 1],255,ValueCast(
          SYSTEM_P3_pansichar,(*keystr)[keystr_i - 1]));
      } while (keystr_i++ !=  _stop);

    }
    result = GXFILE_tgxfileobj_DOT_gdxdatawritestr(ValueCast(
      GXFILE_tgxfileobj,pgdx),keystr_s,values);
    return result;
}  /* cgdxdatawritestr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxdatawritestrstart(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar syid,
    SYSTEM_P3_pansichar expltxt,
    SYSTEM_integer dimen,
    SYSTEM_integer typ,
    SYSTEM_integer userinfo)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxdatawritestrstart(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,syid),
        PCHUTIL_pchartostr(_t2,255,expltxt),dimen,typ,userinfo);
    }
    return result;
}  /* cgdxdatawritestrstart */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxgetdllversion(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar v)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_v;

    result = GXFILE_tgxfileobj_DOT_gdxgetdllversion(ValueCast(
      GXFILE_tgxfileobj,pgdx),local_v);
    SYSUTILS_P3_strpcopy(v,local_v);
    return result;
}  /* cgdxgetdllversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cgdxerrorstr(
    SYSTEM_pointer pgdx,
    SYSTEM_integer errnr,
    SYSTEM_P3_pansichar errmsg)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_errmsg;

    result = GXFILE_tgxfileobj_DOT_gdxerrorstr(ValueCast(
      GXFILE_tgxfileobj,pgdx),errnr,local_errmsg);
    SYSUTILS_P3_strpcopy(errmsg,local_errmsg);
    return result;
}  /* cgdxerrorstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxfileversion(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar filestr,
    SYSTEM_P3_pansichar producestr)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_filestr;
    SYSTEM_shortstring local_producestr;

    result = GXFILE_tgxfileobj_DOT_gdxfileversion(ValueCast(
      GXFILE_tgxfileobj,pgdx),local_filestr,local_producestr);
    SYSUTILS_P3_strpcopy(filestr,local_filestr);
    SYSUTILS_P3_strpcopy(producestr,local_producestr);
    return result;
}  /* cgdxfileversion */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxfindsymbol(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar syid,
    SYSTEM_integer *synr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxfindsymbol(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,syid),
        synr);
    }
    return result;
}  /* cgdxfindsymbol */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxgetelemtext(
    SYSTEM_pointer pgdx,
    SYSTEM_integer txtnr,
    SYSTEM_P3_pansichar txt,
    SYSTEM_integer *node)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_txt;

    result = GXFILE_tgxfileobj_DOT_gdxgetelemtext(ValueCast(
      GXFILE_tgxfileobj,pgdx),txtnr,local_txt,node);
    SYSUTILS_P3_strpcopy(txt,local_txt);
    return result;
}  /* cgdxgetelemtext */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cgdxgetuel(
    SYSTEM_pointer pgdx,
    SYSTEM_integer uelnr,
    SYSTEM_P3_pansichar uel)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_uel;

    result = GXFILE_tgxfileobj_DOT_gdxgetuel(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelnr,local_uel);
    SYSUTILS_P3_strpcopy(uel,local_uel);
    return result;
}  /* cgdxgetuel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxopenappend(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar filename,
    SYSTEM_P3_pansichar producer,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxopenappend(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,filename),
        PCHUTIL_pchartostr(_t2,255,producer),errnr);
    }
    return result;
}  /* cgdxopenappend */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cgdxopenread(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar filename,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxopenread(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,filename),
        errnr);
    }
    return result;
}  /* cgdxopenread */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxopenwrite(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar filename,
    SYSTEM_P3_pansichar producer,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxopenwrite(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,filename),
        PCHUTIL_pchartostr(_t2,255,producer),errnr);
    }
    return result;
}  /* cgdxopenwrite */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxopenwriteex(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar filename,
    SYSTEM_P3_pansichar producer,
    SYSTEM_integer compr,
    SYSTEM_integer *errnr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxopenwriteex(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,filename),
        PCHUTIL_pchartostr(_t2,255,producer),compr,errnr);
    }
    return result;
}  /* cgdxopenwriteex */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsettracelevel(
    SYSTEM_pointer pgdx,
    SYSTEM_integer n,
    SYSTEM_P3_pansichar s)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxsettracelevel(ValueCast(
        GXFILE_tgxfileobj,pgdx),n,PCHUTIL_pchartostr(_t1,255,s));
    }
    return result;
}  /* cgdxsettracelevel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymboladdcomment(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_P3_pansichar txt)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxsymboladdcomment(ValueCast(
        GXFILE_tgxfileobj,pgdx),synr,PCHUTIL_pchartostr(_t1,255,txt));
    }
    return result;
}  /* cgdxsymboladdcomment */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymbolgetcomment(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer n,
    SYSTEM_P3_pansichar txt)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_txt;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolgetcomment(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,n,local_txt);
    SYSUTILS_P3_strpcopy(txt,local_txt);
    return result;
}  /* cgdxsymbolgetcomment */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymbolgetdomainx(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_P3_ppointerarray domainids)
{
    SYSTEM_integer result;
    SYSTEM_integer domainids_i;
    GXDEFS_tgdxstrindex domainids_s;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolgetdomainx(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,domainids_s);
    if (result != 0) 
      { register SYSTEM_int32 _stop = 
          GXFILE_tgxfileobj_DOT_gdxsymboldim(ValueCast(
          GXFILE_tgxfileobj,pgdx),synr);
        if ((domainids_i = 1) <=  _stop) do {
          SYSUTILS_P3_strpcopy(ValueCast(SYSTEM_P3_pansichar,(*
            domainids)[domainids_i - 1]),domainids_s[domainids_i - 1]);
        } while (domainids_i++ !=  _stop);

      }
    return result;
}  /* cgdxsymbolgetdomainx */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymbolinfo(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_P3_pansichar syid,
    SYSTEM_integer *dimen,
    SYSTEM_integer *typ)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_syid;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolinfo(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,local_syid,dimen,typ);
    SYSUTILS_P3_strpcopy(syid,local_syid);
    return result;
}  /* cgdxsymbolinfo */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymbolinfox(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_integer *reccnt,
    SYSTEM_integer *userinfo,
    SYSTEM_P3_pansichar expltxt)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_expltxt;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolinfox(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,reccnt,userinfo,local_expltxt);
    SYSUTILS_P3_strpcopy(expltxt,local_expltxt);
    return result;
}  /* cgdxsymbolinfox */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymbolsetdomain(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_ppointerarray domainids)
{
    SYSTEM_integer result;
    SYSTEM_integer domainids_i;
    GXDEFS_tgdxstrindex domainids_s;

    { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,pgdx))->
        GXFILE_tgxfileobj_DOT_fcurrentdim;
      if ((domainids_i = 1) <=  _stop) do {
        PCHUTIL_pchartostr(domainids_s[domainids_i - 1],255,ValueCast(
          SYSTEM_P3_pansichar,(*domainids)[domainids_i - 1]));
      } while (domainids_i++ !=  _stop);

    }
    result = GXFILE_tgxfileobj_DOT_gdxsymbolsetdomain(ValueCast(
      GXFILE_tgxfileobj,pgdx),domainids_s);
    return result;
}  /* cgdxsymbolsetdomain */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxsymbolsetdomainx(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_P3_ppointerarray domainids)
{
    SYSTEM_integer result;
    SYSTEM_integer domainids_i;
    GXDEFS_tgdxstrindex domainids_s;

    { register SYSTEM_int32 _stop = GXFILE_tgxfileobj_DOT_gdxsymboldim(ValueCast(
        GXFILE_tgxfileobj,pgdx),synr);
      if ((domainids_i = 1) <=  _stop) do {
        PCHUTIL_pchartostr(domainids_s[domainids_i - 1],255,ValueCast(
          SYSTEM_P3_pansichar,(*domainids)[domainids_i - 1]));
      } while (domainids_i++ !=  _stop);

    }
    result = GXFILE_tgxfileobj_DOT_gdxsymbolsetdomainx(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,domainids_s);
    return result;
}  /* cgdxsymbolsetdomainx */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxuelregistermap(
    SYSTEM_pointer pgdx,
    SYSTEM_integer umap,
    SYSTEM_P3_pansichar uel)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxuelregistermap(ValueCast(
        GXFILE_tgxfileobj,pgdx),umap,PCHUTIL_pchartostr(_t1,255,uel));
    }
    return result;
}  /* cgdxuelregistermap */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxuelregisterraw(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar uel)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxuelregisterraw(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,uel));
    }
    return result;
}  /* cgdxuelregisterraw */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxuelregisterstr(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar uel,
    SYSTEM_integer *uelnr)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxuelregisterstr(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,uel),
        uelnr);
    }
    return result;
}  /* cgdxuelregisterstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxumfinduel(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar uel,
    SYSTEM_integer *uelnr,
    SYSTEM_integer *uelmap)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;

      result = GXFILE_tgxfileobj_DOT_gdxumfinduel(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,uel),
        uelnr,uelmap);
    }
    return result;
}  /* cgdxumfinduel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL cgdxumuelget(
    SYSTEM_pointer pgdx,
    SYSTEM_integer uelnr,
    SYSTEM_P3_pansichar uel,
    SYSTEM_integer *uelmap)
{
    SYSTEM_integer result;
    SYSTEM_shortstring local_uel;

    result = GXFILE_tgxfileobj_DOT_gdxumuelget(ValueCast(
      GXFILE_tgxfileobj,pgdx),uelnr,local_uel,uelmap);
    SYSUTILS_P3_strpcopy(uel,local_uel);
    return result;
}  /* cgdxumuelget */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    cgdxrenameuel(
    SYSTEM_pointer pgdx,
    SYSTEM_P3_pansichar oldname,
    SYSTEM_P3_pansichar newname)
{
    SYSTEM_integer result;

    {
      SYSTEM_shortstring _t1;
      SYSTEM_shortstring _t2;

      result = GXFILE_tgxfileobj_DOT_gdxrenameuel(ValueCast(
        GXFILE_tgxfileobj,pgdx),PCHUTIL_pchartostr(_t1,255,oldname),
        PCHUTIL_pchartostr(_t2,255,newname));
    }
    return result;
}  /* cgdxrenameuel */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    bgdxdatareadstr(
    SYSTEM_pointer pgdx,
    SYSTEM_shortstring *keystr,
    SYSTEM_double *values,
    SYSTEM_integer *dimfrst)
{
    SYSTEM_integer result;
    SYSTEM_integer keystr_i, keystr_j, keystr_len;
    GXDEFS_tgdxstrindex keystr_s;

    result = GXFILE_tgxfileobj_DOT_gdxdatareadstr(ValueCast(
      GXFILE_tgxfileobj,pgdx),keystr_s,values,dimfrst);
    if (result != 0) 
      { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,
          pgdx))->GXFILE_tgxfileobj_DOT_fcurrentdim;
        if ((keystr_i = 1) <=  _stop) do {
          keystr_len = ValueCast(SYSTEM_int32,keystr_s[keystr_i - 1][0]);
          { register SYSTEM_int32 _stop = keystr_len;
            if ((keystr_j = 1) <=  _stop) do {
              keystr[keystr_i - 1][keystr_j - 1] = keystr_s[
                keystr_i - 1][keystr_j];
            } while (keystr_j++ !=  _stop);

          }
        
        } while (keystr_i++ !=  _stop);

      }
    return result;
}  /* bgdxdatareadstr */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    bgdxdatasliceuels(
    SYSTEM_pointer pgdx,
    const SYSTEM_integer *slicekeyint,
    SYSTEM_shortstring *keystr)
{
    SYSTEM_integer result;
    SYSTEM_integer keystr_i, keystr_j, keystr_len;
    GXDEFS_tgdxstrindex keystr_s;

    result = GXFILE_tgxfileobj_DOT_gdxdatasliceuels(ValueCast(
      GXFILE_tgxfileobj,pgdx),slicekeyint,keystr_s);
    if (result != 0) 
      { register SYSTEM_int32 _stop = (ValueCast(GXFILE_tgxfileobj,
          pgdx))->GXFILE_tgxfileobj_DOT_fcurrentdim;
        if ((keystr_i = 1) <=  _stop) do {
          keystr_len = ValueCast(SYSTEM_int32,keystr_s[keystr_i - 1][0]);
          { register SYSTEM_int32 _stop = keystr_len;
            if ((keystr_j = 1) <=  _stop) do {
              keystr[keystr_i - 1][keystr_j - 1] = keystr_s[
                keystr_i - 1][keystr_j];
            } while (keystr_j++ !=  _stop);

          }
        
        } while (keystr_i++ !=  _stop);

      }
    return result;
}  /* bgdxdatasliceuels */

Extern_C _P3_DllExport Function(SYSTEM_integer )  STDCALL 
    bgdxsymbolgetdomainx(
    SYSTEM_pointer pgdx,
    SYSTEM_integer synr,
    SYSTEM_shortstring *domainids)
{
    SYSTEM_integer result;
    SYSTEM_integer domainids_i, domainids_j, domainids_len;
    GXDEFS_tgdxstrindex domainids_s;

    result = GXFILE_tgxfileobj_DOT_gdxsymbolgetdomainx(ValueCast(
      GXFILE_tgxfileobj,pgdx),synr,domainids_s);
    if (result != 0) 
      { register SYSTEM_int32 _stop = 
          GXFILE_tgxfileobj_DOT_gdxsymboldim(ValueCast(
          GXFILE_tgxfileobj,pgdx),synr);
        if ((domainids_i = 1) <=  _stop) do {
          domainids_len = ValueCast(SYSTEM_int32,domainids_s[
            domainids_i - 1][0]);
          { register SYSTEM_int32 _stop = domainids_len;
            if ((domainids_j = 1) <=  _stop) do {
              domainids[domainids_i - 1][domainids_j - 1] = 
                domainids_s[domainids_i - 1][domainids_j];
            } while (domainids_j++ !=  _stop);

          }
        
        } while (domainids_i++ !=  _stop);

      }
    return result;
}  /* bgdxsymbolgetdomainx */

/* Library gdxdclib */

/* THIS IS THE DLL EXIT ROUTINE */
Extern_C void _P3_DllFini()
{ if (--SYSTEM_dll_refcount != 0) return;
  _P3_Finalizing = 1;
  _P3_DLL_UNWIND();
  _Final_Module_gxfile();
  _Final_Module_runner();
  _Final_Module_gdlaudit();
  _Final_Module_paldoorg();
  _Final_Module_gmsglobx();
  _Final_Module_datastorage();
  _Final_Module_gmsheapnew();
  _Final_Module_strhash();
  _Final_Module_gmsstrm();
  _Final_Module_xcompress();
  _Final_Module_gmslibname();
  _Final_Module_clibtypes();
  _Final_Module_gmsdata();
  _Final_Module_gmsobj();
  _Final_Module_strutilx();
  _Final_Module_gmsglob();
  _Final_Module_gxdefs();
  _Final_Module_pchutil();
  _Final_Module_gmsgen();
  _Final_Module_gmsspecs();
  _Final_Module_idglobal_p3();
  _Final_Module_p3threads();
  _Final_Module_p3ieeefp();
  _Final_Module_p3process();
  _Final_Module_p3utils();
  _Final_Module_p3library();
  _Final_Module_math_p3();
  _Final_Module_sysutils_p3();
  _Final_Module_exceptions();
  _Final_Module_p3private();
  _Final_Module_system_p3();
  _Final_Module_p3platform();
} /* _P3_DllFini */

/* THIS IS THE MAIN DLL ENTRY ROUTINE */ 
Extern_C SYSTEM_integer _P3_DllInit()
{
  if (SYSTEM_dll_refcount++ != 0) return 0;
  _P3_DLL_INIT();

  _Init_Module_p3platform();
  _Init_Module_system_p3();
  _Init_Module_p3private();
  _Init_Module_exceptions();
  _Init_Module_sysutils_p3();
  _Init_Module_math_p3();
  _Init_Module_p3library();
  _Init_Module_p3utils();
  _Init_Module_p3process();
  _Init_Module_p3ieeefp();
  _Init_Module_p3threads();
  _Init_Module_idglobal_p3();
  _Init_Module_gmsspecs();
  _Init_Module_gmsgen();
  _Init_Module_pchutil();
  _Init_Module_gxdefs();
  _Init_Module_gmsglob();
  _Init_Module_strutilx();
  _Init_Module_gmsobj();
  _Init_Module_gmsdata();
  _Init_Module_clibtypes();
  _Init_Module_gmslibname();
  _Init_Module_xcompress();
  _Init_Module_gmsstrm();
  _Init_Module_strhash();
  _Init_Module_gmsheapnew();
  _Init_Module_datastorage();
  _Init_Module_gmsglobx();
  _Init_Module_paldoorg();
  _Init_Module_gdlaudit();
  _Init_Module_runner();
  _Init_Module_gxfile();


  return SYSTEM_exitcode;
}  /* gdxdclib */

