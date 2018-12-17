#ifndef _P3___gmsgen___H
#define _P3___gmsgen___H

typedef _P3SET_255 GMSGEN_tcharset;
typedef SYSTEM_uint32 GMSGEN_tbigindex;
typedef SYSTEM_ansichar GMSGEN_tansichararray[10000001];
typedef GMSGEN_tansichararray *GMSGEN_pansichararray;
typedef SYSTEM_double GMSGEN_doublearray[10000001];
typedef GMSGEN_doublearray *GMSGEN_pdoublearray;
typedef SYSTEM_uint32 _sub_0GMSGEN;
typedef SYSTEM_double GMSGEN_doublearrayone[10000000];
typedef GMSGEN_doublearrayone *GMSGEN_pdoublearrayone;
typedef SYSTEM_textfile *GMSGEN_ptextfile;
typedef SYSTEM_longint GMSGEN_longintarray[10000001];
typedef GMSGEN_longintarray *GMSGEN_plongintarray;
typedef SYSTEM_byte GMSGEN_tbytedataarray[10000001];
typedef GMSGEN_tbytedataarray *GMSGEN_pbytedataarray;
typedef SYSTEM_uint32 _sub_1GMSGEN;
typedef SYSTEM_longint GMSGEN_longintarrayone[10000000];
typedef GMSGEN_longintarrayone *GMSGEN_plongintarrayone;
typedef SYSTEM_uint32 _sub_2GMSGEN;
typedef SYSTEM_boolean GMSGEN_tbooleanarrayone[10000000];
typedef GMSGEN_tbooleanarrayone *GMSGEN_pbooleanarrayone;
typedef SYSTEM_uint32 _sub_3GMSGEN;
typedef SYSTEM_byte GMSGEN_tbytearrayone[10000000];
typedef GMSGEN_tbytearrayone *GMSGEN_pbytearrayone;
typedef SYSTEM_byte GMSGEN_tfileaction; /* Anonymous */ enum{GMSGEN_forread,GMSGEN_forwrite,GMSGEN_forappend};
typedef SYSTEM_uint32 _sub_4GMSGEN;
typedef SYSTEM_integer GMSGEN_tintegerarrayone[10000000];
typedef GMSGEN_tintegerarrayone *GMSGEN_pintegerarrayone;

extern void _Init_Module_gmsgen(void);
extern void _Final_Module_gmsgen(void);

#endif /* ! defined _P3___gmsgen___H */
