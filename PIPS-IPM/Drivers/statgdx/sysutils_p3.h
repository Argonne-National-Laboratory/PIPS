#ifndef _P3___sysutils_p3___H
#define _P3___sysutils_p3___H

/**** C code included from sysutils_p3.pas(45:1): 3 lines ****/
#if defined(_WIN32)
#  define snprintf _snprintf
#endif
cnstdef {SYSUTILS_P3_fareadonly = 1};
cnstdef {SYSUTILS_P3_fahidden = 2};
cnstdef {SYSUTILS_P3_fasysfile = 4};
cnstdef {SYSUTILS_P3_favolumeid = 8};
cnstdef {SYSUTILS_P3_fadirectory = 16};
cnstdef {SYSUTILS_P3_faarchive = 32};
cnstdef {SYSUTILS_P3_fasymlink = 64};
cnstdef {SYSUTILS_P3_faanyfile = 63};
cnstdef {SYSUTILS_P3_hoursperday = 24};
cnstdef {SYSUTILS_P3_minsperhour = 60};
cnstdef {SYSUTILS_P3_secspermin = 60};
cnstdef {SYSUTILS_P3_msecspersec = 1000};
cnstdef {SYSUTILS_P3_minsperday = 1440};
cnstdef {SYSUTILS_P3_secsperday = 86400};
cnstdef {SYSUTILS_P3_msecsperday = 86400000};
cnstdef {SYSUTILS_P3_datedelta = 693594};
cnstdef {SYSUTILS_P3_unixdatedelta = 25569};
cnstdef {SYSUTILS_P3_max_path = 260};
typedef SYSTEM_uint8 _sub_1SYSUTILS_P3;
typedef SYSTEM_byte _arr_0SYSUTILS_P3[2];
typedef struct SYSUTILS_P3_wordrec_S {
  union{
    struct{
      SYSTEM_byte lo,hi;
    } _c1;
    struct{
      _arr_0SYSUTILS_P3 bytes;
    } _c2;
  } _u;
} SYSUTILS_P3_wordrec;

typedef SYSTEM_uint8 _sub_3SYSUTILS_P3;
typedef SYSTEM_byte _arr_2SYSUTILS_P3[4];
typedef struct SYSUTILS_P3_longrec_S {
  union{
    struct{
      SYSTEM_word lo,hi;
    } _c1;
    struct{
      _arr_2SYSUTILS_P3 bytes;
    } _c2;
  } _u;
} SYSUTILS_P3_longrec;

/* (skipped fwd ref decl:) typedef SYSUTILS_P3_tbytearray *pbytearray; */
typedef SYSTEM_uint16 _sub_4SYSUTILS_P3;
typedef SYSTEM_byte SYSUTILS_P3_tbytearray[32768];
typedef SYSUTILS_P3_tbytearray *SYSUTILS_P3_pbytearray; /* previous fwd ptr */
/* (skipped fwd ref decl:) typedef SYSUTILS_P3_twordarray *pwordarray; */
typedef SYSTEM_uint16 _sub_5SYSUTILS_P3;
typedef SYSTEM_word SYSUTILS_P3_twordarray[16384];
typedef SYSUTILS_P3_twordarray *SYSUTILS_P3_pwordarray; /* previous fwd ptr */

Prototype Procedure (*SYSUTILS_P3_tprocedure)(void);

typedef SYSTEM_shortstring SYSUTILS_P3_tfilename;
typedef SYSTEM_longword SYSUTILS_P3_dword;
typedef SYSTEM_pointer SYSUTILS_P3_searchhandle_t;
typedef SYSTEM_longword SYSUTILS_P3_t__u_long;
typedef SYSUTILS_P3_t__u_long SYSUTILS_P3_t__ino_t;
typedef SYSTEM_longint SYSUTILS_P3_t__off_t;
typedef SYSTEM_uint8 _sub_7SYSUTILS_P3;
typedef SYSTEM_ansichar _arr_6SYSUTILS_P3[256];
typedef struct SYSUTILS_P3_dirent_S {
  SYSUTILS_P3_t__ino_t d_ino;
  SYSUTILS_P3_t__off_t d_off;
  SYSTEM_word d_reclen;
  SYSTEM_byte d_type;
  _arr_6SYSUTILS_P3 d_name;
} SYSUTILS_P3_dirent;

typedef SYSUTILS_P3_dirent SYSUTILS_P3_tdirent;
typedef SYSUTILS_P3_tdirent *SYSUTILS_P3_pdirent;
typedef SYSUTILS_P3_pdirent *SYSUTILS_P3_ppdirent;
typedef SYSUTILS_P3_t__u_long SYSUTILS_P3_mode_t;
typedef struct SYSUTILS_P3_tfiletime_S {
  SYSUTILS_P3_dword dwlowdatetime;
  SYSUTILS_P3_dword dwhighdatetime;
} SYSUTILS_P3_tfiletime;

typedef SYSTEM_uint16 _sub_9SYSUTILS_P3;
typedef SYSTEM_ansichar _arr_8SYSUTILS_P3[260];
typedef SYSTEM_uint8 _sub_11SYSUTILS_P3;
typedef SYSTEM_ansichar _arr_10SYSUTILS_P3[14];
typedef struct SYSUTILS_P3_twin32finddata_S {
  SYSUTILS_P3_dword dwfileattributes;
  SYSUTILS_P3_tfiletime ftcreationtime;
  SYSUTILS_P3_tfiletime ftlastaccesstime;
  SYSUTILS_P3_tfiletime ftlastwritetime;
  SYSUTILS_P3_dword nfilesizehigh;
  SYSUTILS_P3_dword nfilesizelow;
  SYSUTILS_P3_dword dwreserved0;
  SYSUTILS_P3_dword dwreserved1;
  _arr_8SYSUTILS_P3 cfilename;
  _arr_10SYSUTILS_P3 calternatefilename;
} SYSUTILS_P3_twin32finddata;

typedef struct SYSUTILS_P3_tsearchrec_S {
  SYSTEM_integer time;
  SYSTEM_integer size;
  SYSTEM_integer attr;
  SYSUTILS_P3_tfilename name;
  SYSTEM_integer excludeattr;
  SYSUTILS_P3_searchhandle_t findhandle;
  SYSUTILS_P3_twin32finddata finddata;
  SYSTEM_shortstring pathonly;
  SYSTEM_shortstring pattern;
  SYSUTILS_P3_mode_t mode;
} SYSUTILS_P3_tsearchrec;

typedef struct SYSUTILS_P3_ttimestamp_S {
  SYSTEM_integer time;
  SYSTEM_integer date;
} SYSUTILS_P3_ttimestamp;

extern SYSTEM_ansichar SYSUTILS_P3_pathdelim, SYSUTILS_P3_drivedelim, 
  SYSUTILS_P3_pathsep;

Function(SYSTEM_pointer ) SYSUTILS_P3_allocmem(
  SYSTEM_cardinal sz);

Function(SYSTEM_ansichar *) SYSUTILS_P3_uppercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) SYSUTILS_P3_lowercase(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) SYSUTILS_P3_comparestr(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_integer ) SYSUTILS_P3_comparetext(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_boolean ) SYSUTILS_P3_sametext(
  const SYSTEM_ansichar *s1,
  const SYSTEM_ansichar *s2);

Function(SYSTEM_ansichar *) SYSUTILS_P3_trim(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) SYSUTILS_P3_trimleft(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) SYSUTILS_P3_trimright(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) SYSUTILS_P3_inttostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 n);

Function(SYSTEM_ansichar *) SYSUTILS_P3_inttohex(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_int64 v,
  SYSTEM_integer d);

Function(SYSTEM_integer ) SYSUTILS_P3_strtoint(
  const SYSTEM_ansichar *s);

Function(SYSTEM_int64 ) SYSUTILS_P3_strtoint64(
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) SYSUTILS_P3_fileage(
  const SYSTEM_ansichar *filename);

Function(SYSTEM_boolean ) SYSUTILS_P3_fileexists(
  const SYSTEM_ansichar *filename);

Function(SYSTEM_boolean ) SYSUTILS_P3_directoryexists(
  const SYSTEM_ansichar *directory);

Function(SYSTEM_integer ) SYSUTILS_P3_findfirst(
  const SYSTEM_ansichar *path,
  SYSTEM_integer attr,
  SYSUTILS_P3_tsearchrec *f);

Function(SYSTEM_integer ) SYSUTILS_P3_findnext(
  SYSUTILS_P3_tsearchrec *f);

Procedure SYSUTILS_P3_findclose(
  SYSUTILS_P3_tsearchrec *f);

Function(SYSTEM_boolean ) SYSUTILS_P3_deletefile(
  const SYSTEM_ansichar *filename);

Function(SYSTEM_boolean ) SYSUTILS_P3_renamefile(
  const SYSTEM_ansichar *oldname,
  const SYSTEM_ansichar *newname);

Function(SYSTEM_ansichar *) SYSUTILS_P3_changefileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename,
  const SYSTEM_ansichar *extension);

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfilepath(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename);

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfiledir(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename);

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfilename(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename);

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractfileext(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename);

Function(SYSTEM_ansichar *) SYSUTILS_P3_extractshortpathname(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *filename);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_filedatetodatetime(
  SYSTEM_integer filedate);

Function(SYSTEM_integer ) SYSUTILS_P3_datetimetofiledate(
  SYSTEM_P3_tdatetime datetime);

Function(SYSTEM_ansichar *) SYSUTILS_P3_getcurrentdir(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret);

Function(SYSTEM_boolean ) SYSUTILS_P3_setcurrentdir(
  const SYSTEM_ansichar *dir);

Function(SYSTEM_boolean ) SYSUTILS_P3_createdir(
  const SYSTEM_ansichar *dir);

Function(SYSTEM_boolean ) SYSUTILS_P3_removedir(
  const SYSTEM_ansichar *dir);

Function(SYSTEM_cardinal ) SYSUTILS_P3_strlen(
  SYSTEM_P3_pansichar str);

Function(SYSTEM_P3_pansichar ) SYSUTILS_P3_strpcopy(
  SYSTEM_P3_pansichar dest,
  const SYSTEM_ansichar *src);

Function(SYSTEM_ansichar *) SYSUTILS_P3_floattostr(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_double v);

Function(SYSUTILS_P3_ttimestamp *) SYSUTILS_P3_datetimetotimestamp(
  SYSUTILS_P3_ttimestamp *result,
  SYSTEM_P3_tdatetime datetime);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_timestamptodatetime(
  const SYSUTILS_P3_ttimestamp *timestamp);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_encodedate(
  SYSTEM_word year,
  SYSTEM_word month,
  SYSTEM_word day);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_encodetime(
  SYSTEM_word hour,
  SYSTEM_word _min,
  SYSTEM_word sec,
  SYSTEM_word msec);

Function(SYSTEM_boolean ) SYSUTILS_P3_tryencodedate(
  SYSTEM_word year,
  SYSTEM_word month,
  SYSTEM_word day,
  SYSTEM_P3_tdatetime *date);

Function(SYSTEM_boolean ) SYSUTILS_P3_tryencodetime(
  SYSTEM_word hour,
  SYSTEM_word _min,
  SYSTEM_word sec,
  SYSTEM_word msec,
  SYSTEM_P3_tdatetime *time);

Procedure SYSUTILS_P3_decodedate(
  SYSTEM_P3_tdatetime datetime,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day);

Function(SYSTEM_boolean ) SYSUTILS_P3_decodedatefully(
  SYSTEM_P3_tdatetime datetime,
  SYSTEM_word *year,
  SYSTEM_word *month,
  SYSTEM_word *day,
  SYSTEM_word *dow);

Procedure SYSUTILS_P3_decodetime(
  SYSTEM_P3_tdatetime datetime,
  SYSTEM_word *hour,
  SYSTEM_word *_min,
  SYSTEM_word *sec,
  SYSTEM_word *msec);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_date(void);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_time(void);

Function(SYSTEM_P3_tdatetime ) SYSUTILS_P3_now(void);

Function(SYSTEM_boolean ) SYSUTILS_P3_isleapyear(
  SYSTEM_word year);
typedef SYSTEM_uint8 _sub_12SYSUTILS_P3;
typedef SYSTEM_word SYSUTILS_P3_tdaytable[12];
typedef SYSUTILS_P3_tdaytable *SYSUTILS_P3_pdaytable;
typedef SYSUTILS_P3_tdaytable _arr_13SYSUTILS_P3[2];
extern _arr_13SYSUTILS_P3 SYSUTILS_P3_monthdays;

Function(SYSTEM_ansichar *) SYSUTILS_P3_syserrormessage(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  SYSTEM_integer errorcode);

Procedure SYSUTILS_P3_sleep(
  SYSTEM_cardinal milliseconds);

Function(SYSTEM_ansichar *) SYSUTILS_P3_includetrailingpathdelimiter(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_ansichar *) SYSUTILS_P3_excludetrailingpathdelimiter(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *s);

Function(SYSTEM_integer ) SYSUTILS_P3_lastdelimiter(
  const SYSTEM_ansichar *delimiters,
  const SYSTEM_ansichar *s);

Procedure SYSUTILS_P3_freeandnil(
  SYSTEM_untyped *obj);

Function(SYSTEM_ansichar *) SYSUTILS_P3_getenvironmentvariable(
  SYSTEM_ansichar *result,
  SYSTEM_uint8 _len_ret,
  const SYSTEM_ansichar *name);

extern void _Init_Module_sysutils_p3(void);
extern void _Final_Module_sysutils_p3(void);

#endif /* ! defined _P3___sysutils_p3___H */
