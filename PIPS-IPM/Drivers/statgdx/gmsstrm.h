#ifndef _P3___gmsstrm___H
#define _P3___gmsstrm___H

cnstdef {GMSSTRM_sofrombeginning = 0};
cnstdef {GMSSTRM_sofromcurrent = 1};
cnstdef {GMSSTRM_sofromend = 2};
cnstdef {GMSSTRM_fmcreate = 65535};
cnstdef {GMSSTRM_fmopenread = 0};
cnstdef {GMSSTRM_fmopenwrite = 1};
cnstdef {GMSSTRM_fmopenreadwrite = 2};
extern SYSTEM_integer GMSSTRM_buffersize;
cnstdef {GMSSTRM_strmerrornoerror = 0};
cnstdef {GMSSTRM_strmerrorioresult = 1};
cnstdef {GMSSTRM_strmerrorgamsheader = 2};
cnstdef {GMSSTRM_strmerrornopassword = 3};
cnstdef {GMSSTRM_strmerrorintegrity = 4};
cnstdef {GMSSTRM_strmerrorzlib = 5};
typedef SYSTEM_byte GMSSTRM_rwtype; /* Anonymous */ enum{GMSSTRM_rw_byte,GMSSTRM_rw_bool,GMSSTRM_rw_char,GMSSTRM_rw_word,
  GMSSTRM_rw_integer,GMSSTRM_rw_int64,GMSSTRM_rw_double,
  GMSSTRM_rw_string,GMSSTRM_rw_pchar,GMSSTRM_rw_pstring};
typedef _P3STR_7 _arr_0GMSSTRM[10];
extern _arr_0GMSSTRM GMSSTRM_rwtypetext;
typedef struct GMSSTRM_txstream_OD_S* GMSSTRM_txstream; /* sy_class */
typedef struct GMSSTRM_txstream_OD_S {  /* Objects of 'txstream' */
  SYSTEM_classreference_t CD;  /* = &GMSSTRM_txstream_CD */
} GMSSTRM_txstream_OD;


Procedure GMSSTRM_txstream_DOT_parwrite(
  GMSSTRM_txstream self,
  GMSSTRM_rwtype t);

Procedure GMSSTRM_txstream_DOT_parcheck(
  GMSSTRM_txstream self,
  GMSSTRM_rwtype t);

Prototype Function(SYSTEM_int64 ) (*GMSSTRM_txstream_DOT_getposition_T)(
  GMSSTRM_txstream self);

Prototype Procedure (*GMSSTRM_txstream_DOT_setposition_T)(
  GMSSTRM_txstream self,
  SYSTEM_int64 p);

Prototype Function(SYSTEM_int64 ) (*GMSSTRM_txstream_DOT_getsize_T)(
  GMSSTRM_txstream self);

Prototype Function(SYSTEM_longword ) (*GMSSTRM_txstream_DOT_read_T)(
  GMSSTRM_txstream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count);

Prototype Function(SYSTEM_longword ) (*GMSSTRM_txstream_DOT_write_T)(
  GMSSTRM_txstream self,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword count);

Procedure GMSSTRM_txstream_DOT_writestring(
  GMSSTRM_txstream self,
  const SYSTEM_ansichar *s);

Procedure GMSSTRM_txstream_DOT_writepstring(
  GMSSTRM_txstream self,
  SYSTEM_P3_pshortstring ps);

Procedure GMSSTRM_txstream_DOT_writedouble(
  GMSSTRM_txstream self,
  SYSTEM_double x);

Procedure GMSSTRM_txstream_DOT_writeinteger(
  GMSSTRM_txstream self,
  SYSTEM_integer n);

Procedure GMSSTRM_txstream_DOT_writeint64(
  GMSSTRM_txstream self,
  SYSTEM_int64 n);

Procedure GMSSTRM_txstream_DOT_writebyte(
  GMSSTRM_txstream self,
  SYSTEM_byte b);

Procedure GMSSTRM_txstream_DOT_writeword(
  GMSSTRM_txstream self,
  SYSTEM_word w);

Procedure GMSSTRM_txstream_DOT_writebool(
  GMSSTRM_txstream self,
  SYSTEM_boolean b);

Procedure GMSSTRM_txstream_DOT_writechar(
  GMSSTRM_txstream self,
  SYSTEM_ansichar c);

Procedure GMSSTRM_txstream_DOT_writepchar(
  GMSSTRM_txstream self,
  SYSTEM_P3_pansichar p,
  SYSTEM_integer l);

Function(SYSTEM_ansichar *) GMSSTRM_txstream_DOT_readstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSSTRM_txstream self);

Prototype Function(SYSTEM_double ) (*GMSSTRM_txstream_DOT_readdouble_T)(
  GMSSTRM_txstream self);

Function(SYSTEM_double ) GMSSTRM_txstream_DOT_readdouble(
  GMSSTRM_txstream self);

Prototype Function(SYSTEM_integer ) (*
  GMSSTRM_txstream_DOT_readinteger_T)(
  GMSSTRM_txstream self);

Function(SYSTEM_integer ) GMSSTRM_txstream_DOT_readinteger(
  GMSSTRM_txstream self);

Function(SYSTEM_byte ) GMSSTRM_txstream_DOT_readbyte(
  GMSSTRM_txstream self);

Prototype Function(SYSTEM_word ) (*GMSSTRM_txstream_DOT_readword_T)(
  GMSSTRM_txstream self);

Function(SYSTEM_word ) GMSSTRM_txstream_DOT_readword(
  GMSSTRM_txstream self);

Prototype Function(SYSTEM_int64 ) (*GMSSTRM_txstream_DOT_readint64_T)(
  GMSSTRM_txstream self);

Function(SYSTEM_int64 ) GMSSTRM_txstream_DOT_readint64(
  GMSSTRM_txstream self);

Function(SYSTEM_boolean ) GMSSTRM_txstream_DOT_readbool(
  GMSSTRM_txstream self);

Function(SYSTEM_ansichar ) GMSSTRM_txstream_DOT_readchar(
  GMSSTRM_txstream self);

Procedure GMSSTRM_txstream_DOT_readpchar(
  GMSSTRM_txstream self,
  SYSTEM_P3_pansichar *p,
  SYSTEM_integer *l);

Procedure GMSSTRM_txstream_DOT_readpstring(
  GMSSTRM_txstream self,
  SYSTEM_P3_pshortstring *ps);
extern void * const GMSSTRM_txstream_VT[];
extern const SYSTEM_classdescriptor_t GMSSTRM_txstream_CD;


typedef struct GMSSTRM_txfilestream_OD_S* GMSSTRM_txfilestream; /* sy_class */
typedef struct GMSSTRM_txfilestream_OD_S {  /* Objects of 'txfilestream' */
  SYSTEM_classreference_t CD;  /* = &GMSSTRM_txfilestream_CD */
  P3UTILS_tp3filehandle GMSSTRM_txfilestream_DOT_fs;
  SYSTEM_boolean GMSSTRM_txfilestream_DOT_fileisopen;
  SYSTEM_int64 GMSSTRM_txfilestream_DOT_physposition;
  SYSTEM_shortstring GMSSTRM_txfilestream_DOT_fpassword;
  SYSTEM_shortstring GMSSTRM_txfilestream_DOT_ffilename;
  SYSTEM_integer GMSSTRM_txfilestream_DOT_flastioresult;
} GMSSTRM_txfilestream_OD;


Procedure GMSSTRM_txfilestream_DOT_setlastioresult(
  GMSSTRM_txfilestream self,
  SYSTEM_integer v);

Function(SYSTEM_integer ) GMSSTRM_txfilestream_DOT_getlastioresult(
  GMSSTRM_txfilestream self);

Procedure GMSSTRM_txfilestream_DOT_setpassword(
  GMSSTRM_txfilestream self,
  const SYSTEM_ansichar *s);

Function(SYSTEM_boolean ) GMSSTRM_txfilestream_DOT_getusespassword(
  GMSSTRM_txfilestream self);

Function(SYSTEM_ansichar *) GMSSTRM_txfilestream_DOT_randstring(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSSTRM_txfilestream self,
  SYSTEM_integer l);

Function(SYSTEM_int64 ) GMSSTRM_txfilestream_DOT_getsize(
  GMSSTRM_txfilestream self);

Function(SYSTEM_int64 ) GMSSTRM_txfilestream_DOT_getposition(
  GMSSTRM_txfilestream self);

Procedure GMSSTRM_txfilestream_DOT_setposition(
  GMSSTRM_txfilestream self,
  SYSTEM_int64 p);

Constructor(GMSSTRM_txfilestream ) GMSSTRM_txfilestream_DOT_create(
  GMSSTRM_txfilestream self,
  const SYSTEM_ansichar *afilename,
  SYSTEM_word amode);

Destructor(GMSSTRM_txfilestream ) GMSSTRM_txfilestream_DOT_destroy(
  GMSSTRM_txfilestream self);

Procedure GMSSTRM_txfilestream_DOT_applypassword(
  GMSSTRM_txfilestream self,
  GMSGEN_pbytedataarray pr,
  GMSGEN_pbytedataarray pw,
  SYSTEM_integer len,
  SYSTEM_int64 offs);

Function(SYSTEM_longword ) GMSSTRM_txfilestream_DOT_read(
  GMSSTRM_txfilestream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count);

Function(SYSTEM_longword ) GMSSTRM_txfilestream_DOT_write(
  GMSSTRM_txfilestream self,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword count);
extern void * const GMSSTRM_txfilestream_VT[];
extern const SYSTEM_classdescriptor_t GMSSTRM_txfilestream_CD;


typedef struct GMSSTRM_tcompressheader_S {
  SYSTEM_byte cxtyp,cxb1,cxb2;
} GMSSTRM_tcompressheader;

typedef struct GMSSTRM_tcompressbuffer_S *GMSSTRM_pcompressbuffer;
typedef struct GMSSTRM_tcompressbuffer_S {
  GMSSTRM_tcompressheader cxheader;
  SYSTEM_byte cxdata;
} GMSSTRM_tcompressbuffer;

typedef struct GMSSTRM_tbufferedfilestream_OD_S* 
  GMSSTRM_tbufferedfilestream; /* sy_class */
typedef struct GMSSTRM_tbufferedfilestream_OD_S {  /* Objects of 'tbufferedfilestream' */
  SYSTEM_classreference_t CD;  /* = &GMSSTRM_tbufferedfilestream_CD */
  P3UTILS_tp3filehandle GMSSTRM_txfilestream_DOT_fs;
  SYSTEM_boolean GMSSTRM_txfilestream_DOT_fileisopen;
  SYSTEM_int64 GMSSTRM_txfilestream_DOT_physposition;
  SYSTEM_shortstring GMSSTRM_txfilestream_DOT_fpassword;
  SYSTEM_shortstring GMSSTRM_txfilestream_DOT_ffilename;
  SYSTEM_integer GMSSTRM_txfilestream_DOT_flastioresult;
  GMSGEN_pbytedataarray GMSSTRM_tbufferedfilestream_DOT_bufptr;
  GMSSTRM_pcompressbuffer GMSSTRM_tbufferedfilestream_DOT_cbufptr;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_bufsize;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_cbufsize;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_nrloaded;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_nrread;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_nrwritten;
  SYSTEM_boolean GMSSTRM_tbufferedfilestream_DOT_fcompress;
  SYSTEM_boolean GMSSTRM_tbufferedfilestream_DOT_fcancompress;
  SYSTEM_shortstring GMSSTRM_tbufferedfilestream_DOT_floadpath;
} GMSSTRM_tbufferedfilestream_OD;


Function(SYSTEM_boolean ) GMSSTRM_tbufferedfilestream_DOT_fillbuffer(
  GMSSTRM_tbufferedfilestream self);

Procedure GMSSTRM_tbufferedfilestream_DOT_setcompression(
  GMSSTRM_tbufferedfilestream self,
  SYSTEM_boolean v);

Function(SYSTEM_int64 ) GMSSTRM_tbufferedfilestream_DOT_getposition(
  GMSSTRM_tbufferedfilestream self);

Procedure GMSSTRM_tbufferedfilestream_DOT_setposition(
  GMSSTRM_tbufferedfilestream self,
  SYSTEM_int64 p);

Function(SYSTEM_int64 ) GMSSTRM_tbufferedfilestream_DOT_getsize(
  GMSSTRM_tbufferedfilestream self);

Constructor(GMSSTRM_tbufferedfilestream ) 
  GMSSTRM_tbufferedfilestream_DOT_create(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode);

Constructor(GMSSTRM_tbufferedfilestream ) 
  GMSSTRM_tbufferedfilestream_DOT_createwithpath(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode,
  const SYSTEM_ansichar *loadpath);

Destructor(GMSSTRM_tbufferedfilestream ) 
  GMSSTRM_tbufferedfilestream_DOT_destroy(
  GMSSTRM_tbufferedfilestream self);

Function(SYSTEM_boolean ) GMSSTRM_tbufferedfilestream_DOT_flushbuffer(
  GMSSTRM_tbufferedfilestream self);

Function(SYSTEM_longword ) GMSSTRM_tbufferedfilestream_DOT_read(
  GMSSTRM_tbufferedfilestream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count);

Function(SYSTEM_ansichar ) 
  GMSSTRM_tbufferedfilestream_DOT_readcharacter(
  GMSSTRM_tbufferedfilestream self);

Function(SYSTEM_longword ) GMSSTRM_tbufferedfilestream_DOT_write(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_untyped *buffer,
  SYSTEM_longword count);

Function(SYSTEM_boolean ) GMSSTRM_tbufferedfilestream_DOT_iseof(
  GMSSTRM_tbufferedfilestream self);

Function(SYSTEM_ansichar *) 
  GMSSTRM_tbufferedfilestream_DOT_getloadpath(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  GMSSTRM_tbufferedfilestream self);

Procedure GMSSTRM_tbufferedfilestream_DOT_setloadpath(
  GMSSTRM_tbufferedfilestream self,
  const SYSTEM_ansichar *s);
extern void * const GMSSTRM_tbufferedfilestream_VT[];
extern const SYSTEM_classdescriptor_t GMSSTRM_tbufferedfilestream_CD;


typedef struct GMSSTRM_tmibufferedstream_OD_S* 
  GMSSTRM_tmibufferedstream; /* sy_class */
typedef struct GMSSTRM_tmibufferedstream_OD_S {  /* Objects of 'tmibufferedstream' */
  SYSTEM_classreference_t CD;  /* = &GMSSTRM_tmibufferedstream_CD */
  P3UTILS_tp3filehandle GMSSTRM_txfilestream_DOT_fs;
  SYSTEM_boolean GMSSTRM_txfilestream_DOT_fileisopen;
  SYSTEM_int64 GMSSTRM_txfilestream_DOT_physposition;
  SYSTEM_shortstring GMSSTRM_txfilestream_DOT_fpassword;
  SYSTEM_shortstring GMSSTRM_txfilestream_DOT_ffilename;
  SYSTEM_integer GMSSTRM_txfilestream_DOT_flastioresult;
  GMSGEN_pbytedataarray GMSSTRM_tbufferedfilestream_DOT_bufptr;
  GMSSTRM_pcompressbuffer GMSSTRM_tbufferedfilestream_DOT_cbufptr;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_bufsize;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_cbufsize;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_nrloaded;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_nrread;
  SYSTEM_longword GMSSTRM_tbufferedfilestream_DOT_nrwritten;
  SYSTEM_boolean GMSSTRM_tbufferedfilestream_DOT_fcompress;
  SYSTEM_boolean GMSSTRM_tbufferedfilestream_DOT_fcancompress;
  SYSTEM_shortstring GMSSTRM_tbufferedfilestream_DOT_floadpath;
  SYSTEM_byte GMSSTRM_tmibufferedstream_DOT_order_word;
  SYSTEM_byte GMSSTRM_tmibufferedstream_DOT_order_integer;
  SYSTEM_byte GMSSTRM_tmibufferedstream_DOT_order_double;
  SYSTEM_byte GMSSTRM_tmibufferedstream_DOT_size_word;
  SYSTEM_byte GMSSTRM_tmibufferedstream_DOT_size_integer;
  SYSTEM_byte GMSSTRM_tmibufferedstream_DOT_size_double;
  SYSTEM_boolean GMSSTRM_tmibufferedstream_DOT_normalorder;
} GMSSTRM_tmibufferedstream_OD;


Procedure GMSSTRM_tmibufferedstream_DOT_determinebyteorder(
  GMSSTRM_tmibufferedstream self);

Constructor(GMSSTRM_tmibufferedstream ) 
  GMSSTRM_tmibufferedstream_DOT_create(
  GMSSTRM_tmibufferedstream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode);

Constructor(GMSSTRM_tmibufferedstream ) 
  GMSSTRM_tmibufferedstream_DOT_createwithpath(
  GMSSTRM_tmibufferedstream self,
  const SYSTEM_ansichar *filename,
  SYSTEM_word mode,
  const SYSTEM_ansichar *loadpath);

Procedure GMSSTRM_tmibufferedstream_DOT_reversebytes(
  GMSSTRM_tmibufferedstream self,
  SYSTEM_pointer psrc,
  SYSTEM_pointer pdest,
  SYSTEM_integer sz);

Function(SYSTEM_integer ) GMSSTRM_tmibufferedstream_DOT_goodbyteorder(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_double ) GMSSTRM_tmibufferedstream_DOT_readdouble(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_integer ) GMSSTRM_tmibufferedstream_DOT_readinteger(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_word ) GMSSTRM_tmibufferedstream_DOT_readword(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_int64 ) GMSSTRM_tmibufferedstream_DOT_readint64(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_boolean ) GMSSTRM_tmibufferedstream_DOT_wordsneedflip(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_boolean ) GMSSTRM_tmibufferedstream_DOT_intsneedflip(
  GMSSTRM_tmibufferedstream self);

Procedure GMSSTRM_tmibufferedstream_DOT_writegmsinteger(
  GMSSTRM_tmibufferedstream self,
  SYSTEM_integer n);

Procedure GMSSTRM_tmibufferedstream_DOT_writegmsdouble(
  GMSSTRM_tmibufferedstream self,
  SYSTEM_double d);

Function(SYSTEM_integer ) GMSSTRM_tmibufferedstream_DOT_readgmsinteger(
  GMSSTRM_tmibufferedstream self);

Function(SYSTEM_double ) GMSSTRM_tmibufferedstream_DOT_readgmsdouble(
  GMSSTRM_tmibufferedstream self);
extern void * const GMSSTRM_tmibufferedstream_VT[];
extern const SYSTEM_classdescriptor_t GMSSTRM_tmibufferedstream_CD;


typedef struct GMSSTRM_tgzipinputstream_OD_S* GMSSTRM_tgzipinputstream; /* sy_class */
typedef struct GMSSTRM_tgzipinputstream_OD_S {  /* Objects of 'tgzipinputstream' */
  SYSTEM_classreference_t CD;  /* = &GMSSTRM_tgzipinputstream_CD */
  XCOMPRESS_pgzfile GMSSTRM_tgzipinputstream_DOT_pgz;
  GMSGEN_pbytedataarray GMSSTRM_tgzipinputstream_DOT_bufptr;
  SYSTEM_longword GMSSTRM_tgzipinputstream_DOT_bufsize;
  SYSTEM_longword GMSSTRM_tgzipinputstream_DOT_cbufsize;
  SYSTEM_longword GMSSTRM_tgzipinputstream_DOT_nrloaded;
  SYSTEM_longword GMSSTRM_tgzipinputstream_DOT_nrread;
} GMSSTRM_tgzipinputstream_OD;


Constructor(GMSSTRM_tgzipinputstream ) 
  GMSSTRM_tgzipinputstream_DOT_create(
  GMSSTRM_tgzipinputstream self,
  const SYSTEM_ansichar *fn,
  SYSTEM_ansichar *errmsg);

Destructor(GMSSTRM_tgzipinputstream ) 
  GMSSTRM_tgzipinputstream_DOT_destroy(
  GMSSTRM_tgzipinputstream self);

Function(SYSTEM_longword ) GMSSTRM_tgzipinputstream_DOT_read(
  GMSSTRM_tgzipinputstream self,
  SYSTEM_untyped *buffer,
  SYSTEM_longword count);

Procedure GMSSTRM_tgzipinputstream_DOT_readline(
  GMSSTRM_tgzipinputstream self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer *len,
  SYSTEM_integer maxinp,
  SYSTEM_ansichar *lastchar);
extern void * const GMSSTRM_tgzipinputstream_VT[];
extern const SYSTEM_classdescriptor_t GMSSTRM_tgzipinputstream_CD;


typedef SYSTEM_byte GMSSTRM_tfilesignature; /* Anonymous */ enum{GMSSTRM_fsign_text,GMSSTRM_fsign_blocktext,GMSSTRM_fsign_gzip};
typedef SYSTEM_byte _enm_1GMSSTRM; /* Anonymous */ enum{GMSSTRM_fm_read,GMSSTRM_fm_write};
typedef struct GMSSTRM_tbinarytextfileio_OD_S* 
  GMSSTRM_tbinarytextfileio; /* sy_class */
typedef struct GMSSTRM_tbinarytextfileio_OD_S {  /* Objects of 'tbinarytextfileio' */
  SYSTEM_classreference_t CD;  /* = &GMSSTRM_tbinarytextfileio_CD */
  GMSSTRM_tbufferedfilestream GMSSTRM_tbinarytextfileio_DOT_fs;
  GMSSTRM_tgzipinputstream GMSSTRM_tbinarytextfileio_DOT_gzfs;
  _enm_1GMSSTRM GMSSTRM_tbinarytextfileio_DOT_frw;
  GMSSTRM_tfilesignature GMSSTRM_tbinarytextfileio_DOT_ffilesignature;
  SYSTEM_byte GMSSTRM_tbinarytextfileio_DOT_fmajorversionread,GMSSTRM_tbinarytextfileio_DOT_fminorversionread;
  SYSTEM_int64 GMSSTRM_tbinarytextfileio_DOT_frewindpoint;
} GMSSTRM_tbinarytextfileio_OD;


Function(SYSTEM_integer ) 
  GMSSTRM_tbinarytextfileio_DOT_getlastioresult(
  GMSSTRM_tbinarytextfileio self);

Constructor(GMSSTRM_tbinarytextfileio ) 
  GMSSTRM_tbinarytextfileio_DOT_openforread(
  GMSSTRM_tbinarytextfileio self,
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *password,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg);

Constructor(GMSSTRM_tbinarytextfileio ) 
  GMSSTRM_tbinarytextfileio_DOT_openforwrite(
  GMSSTRM_tbinarytextfileio self,
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *producer,
  const SYSTEM_ansichar *password,
  GMSSTRM_tfilesignature signature,
  SYSTEM_boolean comp,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg);

Destructor(GMSSTRM_tbinarytextfileio ) 
  GMSSTRM_tbinarytextfileio_DOT_destroy(
  GMSSTRM_tbinarytextfileio self);

Function(SYSTEM_integer ) GMSSTRM_tbinarytextfileio_DOT_read(
  GMSSTRM_tbinarytextfileio self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer count);

Function(SYSTEM_ansichar ) GMSSTRM_tbinarytextfileio_DOT_readcharacter(
  GMSSTRM_tbinarytextfileio self);

Procedure GMSSTRM_tbinarytextfileio_DOT_readline(
  GMSSTRM_tbinarytextfileio self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer *len,
  SYSTEM_integer maxinp,
  SYSTEM_ansichar *lastchar);

Function(SYSTEM_integer ) GMSSTRM_tbinarytextfileio_DOT_write(
  GMSSTRM_tbinarytextfileio self,
  SYSTEM_untyped *buffer,
  SYSTEM_integer count);

Function(SYSTEM_boolean ) GMSSTRM_tbinarytextfileio_DOT_usespassword(
  GMSSTRM_tbinarytextfileio self);

Procedure GMSSTRM_tbinarytextfileio_DOT_rewind(
  GMSSTRM_tbinarytextfileio self);
extern void * const GMSSTRM_tbinarytextfileio_VT[];
extern const SYSTEM_classdescriptor_t GMSSTRM_tbinarytextfileio_CD;



Procedure GMSSTRM_uncompresstextfile(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *fo,
  const SYSTEM_ansichar *password,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg);

Procedure GMSSTRM_compresstextfile(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *fo,
  const SYSTEM_ansichar *password,
  SYSTEM_boolean comp,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg);

Procedure GMSSTRM_compressfromstdin(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *password,
  SYSTEM_boolean comp,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg);

Procedure GMSSTRM_uncompresstostdout(
  const SYSTEM_ansichar *fn,
  const SYSTEM_ansichar *password,
  SYSTEM_integer *errnr,
  SYSTEM_ansichar *errmsg);

extern void _Init_Module_gmsstrm(void);
extern void _Final_Module_gmsstrm(void);

#endif /* ! defined _P3___gmsstrm___H */
