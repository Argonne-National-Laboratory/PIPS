#ifndef _P3___exceptions___H
#define _P3___exceptions___H

typedef struct EXCEPTIONS_eabort_OD_S* EXCEPTIONS_eabort; /* sy_class */
typedef struct EXCEPTIONS_eabort_OD_S {  /* Objects of 'eabort' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eabort_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eabort_OD;

extern void * const EXCEPTIONS_eabort_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eabort_CD;


typedef struct EXCEPTIONS_eabstracterror_OD_S* 
  EXCEPTIONS_eabstracterror; /* sy_class */
typedef struct EXCEPTIONS_eabstracterror_OD_S {  /* Objects of 'eabstracterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eabstracterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eabstracterror_OD;

extern void * const EXCEPTIONS_eabstracterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eabstracterror_CD;


typedef struct EXCEPTIONS_eexternal_OD_S* EXCEPTIONS_eexternal; /* sy_class */
typedef struct EXCEPTIONS_eexternal_OD_S {  /* Objects of 'eexternal' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eexternal_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eexternal_OD;

extern void * const EXCEPTIONS_eexternal_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eexternal_CD;


typedef struct EXCEPTIONS_eaccessviolation_OD_S* 
  EXCEPTIONS_eaccessviolation; /* sy_class */
typedef struct EXCEPTIONS_eaccessviolation_OD_S {  /* Objects of 'eaccessviolation' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eaccessviolation_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eaccessviolation_OD;

extern void * const EXCEPTIONS_eaccessviolation_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eaccessviolation_CD;


typedef struct EXCEPTIONS_ematherror_OD_S* EXCEPTIONS_ematherror; /* sy_class */
typedef struct EXCEPTIONS_ematherror_OD_S {  /* Objects of 'ematherror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_ematherror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_ematherror_OD;

extern void * const EXCEPTIONS_ematherror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_ematherror_CD;


typedef struct EXCEPTIONS_einvalidargument_OD_S* 
  EXCEPTIONS_einvalidargument; /* sy_class */
typedef struct EXCEPTIONS_einvalidargument_OD_S {  /* Objects of 'einvalidargument' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_einvalidargument_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_einvalidargument_OD;

extern void * const EXCEPTIONS_einvalidargument_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidargument_CD;


typedef struct EXCEPTIONS_einvalidop_OD_S* EXCEPTIONS_einvalidop; /* sy_class */
typedef struct EXCEPTIONS_einvalidop_OD_S {  /* Objects of 'einvalidop' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_einvalidop_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_einvalidop_OD;

extern void * const EXCEPTIONS_einvalidop_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidop_CD;


typedef struct EXCEPTIONS_eoverflow_OD_S* EXCEPTIONS_eoverflow; /* sy_class */
typedef struct EXCEPTIONS_eoverflow_OD_S {  /* Objects of 'eoverflow' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eoverflow_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eoverflow_OD;

extern void * const EXCEPTIONS_eoverflow_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eoverflow_CD;


typedef struct EXCEPTIONS_eunderflow_OD_S* EXCEPTIONS_eunderflow; /* sy_class */
typedef struct EXCEPTIONS_eunderflow_OD_S {  /* Objects of 'eunderflow' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eunderflow_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eunderflow_OD;

extern void * const EXCEPTIONS_eunderflow_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eunderflow_CD;


typedef struct EXCEPTIONS_ezerodivide_OD_S* EXCEPTIONS_ezerodivide; /* sy_class */
typedef struct EXCEPTIONS_ezerodivide_OD_S {  /* Objects of 'ezerodivide' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_ezerodivide_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_ezerodivide_OD;

extern void * const EXCEPTIONS_ezerodivide_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_ezerodivide_CD;


typedef struct EXCEPTIONS_econtrolc_OD_S* EXCEPTIONS_econtrolc; /* sy_class */
typedef struct EXCEPTIONS_econtrolc_OD_S {  /* Objects of 'econtrolc' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_econtrolc_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_econtrolc_OD;

extern void * const EXCEPTIONS_econtrolc_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_econtrolc_CD;


typedef struct EXCEPTIONS_einterror_OD_S* EXCEPTIONS_einterror; /* sy_class */
typedef struct EXCEPTIONS_einterror_OD_S {  /* Objects of 'einterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_einterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_einterror_OD;

extern void * const EXCEPTIONS_einterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_einterror_CD;


typedef struct EXCEPTIONS_edivbyzero_OD_S* EXCEPTIONS_edivbyzero; /* sy_class */
typedef struct EXCEPTIONS_edivbyzero_OD_S {  /* Objects of 'edivbyzero' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_edivbyzero_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_edivbyzero_OD;

extern void * const EXCEPTIONS_edivbyzero_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_edivbyzero_CD;


typedef struct EXCEPTIONS_erangeerror_OD_S* EXCEPTIONS_erangeerror; /* sy_class */
typedef struct EXCEPTIONS_erangeerror_OD_S {  /* Objects of 'erangeerror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_erangeerror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_erangeerror_OD;

extern void * const EXCEPTIONS_erangeerror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_erangeerror_CD;


typedef struct EXCEPTIONS_eintoverflow_OD_S* EXCEPTIONS_eintoverflow; /* sy_class */
typedef struct EXCEPTIONS_eintoverflow_OD_S {  /* Objects of 'eintoverflow' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eintoverflow_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eintoverflow_OD;

extern void * const EXCEPTIONS_eintoverflow_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eintoverflow_CD;


typedef struct EXCEPTIONS_estackoverflow_OD_S* 
  EXCEPTIONS_estackoverflow; /* sy_class */
typedef struct EXCEPTIONS_estackoverflow_OD_S {  /* Objects of 'estackoverflow' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_estackoverflow_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_estackoverflow_OD;

extern void * const EXCEPTIONS_estackoverflow_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_estackoverflow_CD;


typedef struct EXCEPTIONS_eassertionfailed_OD_S* 
  EXCEPTIONS_eassertionfailed; /* sy_class */
typedef struct EXCEPTIONS_eassertionfailed_OD_S {  /* Objects of 'eassertionfailed' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eassertionfailed_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eassertionfailed_OD;

extern void * const EXCEPTIONS_eassertionfailed_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eassertionfailed_CD;


typedef struct EXCEPTIONS_econverterror_OD_S* EXCEPTIONS_econverterror; /* sy_class */
typedef struct EXCEPTIONS_econverterror_OD_S {  /* Objects of 'econverterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_econverterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_econverterror_OD;

extern void * const EXCEPTIONS_econverterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_econverterror_CD;


typedef struct EXCEPTIONS_eexternalexception_OD_S* 
  EXCEPTIONS_eexternalexception; /* sy_class */
typedef struct EXCEPTIONS_eexternalexception_OD_S {  /* Objects of 'eexternalexception' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eexternalexception_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eexternalexception_OD;

extern void * const EXCEPTIONS_eexternalexception_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eexternalexception_CD;


typedef struct EXCEPTIONS_eheapexception_OD_S* 
  EXCEPTIONS_eheapexception; /* sy_class */
typedef struct EXCEPTIONS_eheapexception_OD_S {  /* Objects of 'eheapexception' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eheapexception_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eheapexception_OD;

extern void * const EXCEPTIONS_eheapexception_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eheapexception_CD;


typedef struct EXCEPTIONS_eoutofmemory_OD_S* EXCEPTIONS_eoutofmemory; /* sy_class */
typedef struct EXCEPTIONS_eoutofmemory_OD_S {  /* Objects of 'eoutofmemory' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eoutofmemory_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eoutofmemory_OD;

extern void * const EXCEPTIONS_eoutofmemory_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eoutofmemory_CD;


typedef struct EXCEPTIONS_einvalidpointer_OD_S* 
  EXCEPTIONS_einvalidpointer; /* sy_class */
typedef struct EXCEPTIONS_einvalidpointer_OD_S {  /* Objects of 'einvalidpointer' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_einvalidpointer_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_einvalidpointer_OD;

extern void * const EXCEPTIONS_einvalidpointer_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidpointer_CD;


typedef struct EXCEPTIONS_einouterror_OD_S* EXCEPTIONS_einouterror; /* sy_class */
typedef struct EXCEPTIONS_einouterror_OD_S {  /* Objects of 'einouterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_einouterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_einouterror_OD;

extern void * const EXCEPTIONS_einouterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_einouterror_CD;


typedef struct EXCEPTIONS_eintfcasterror_OD_S* 
  EXCEPTIONS_eintfcasterror; /* sy_class */
typedef struct EXCEPTIONS_eintfcasterror_OD_S {  /* Objects of 'eintfcasterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eintfcasterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eintfcasterror_OD;

extern void * const EXCEPTIONS_eintfcasterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eintfcasterror_CD;


typedef struct EXCEPTIONS_einvalidcast_OD_S* EXCEPTIONS_einvalidcast; /* sy_class */
typedef struct EXCEPTIONS_einvalidcast_OD_S {  /* Objects of 'einvalidcast' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_einvalidcast_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_einvalidcast_OD;

extern void * const EXCEPTIONS_einvalidcast_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidcast_CD;


typedef struct EXCEPTIONS_eoserror_OD_S* EXCEPTIONS_eoserror; /* sy_class */
typedef struct EXCEPTIONS_eoserror_OD_S {  /* Objects of 'eoserror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_eoserror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_eoserror_OD;

extern void * const EXCEPTIONS_eoserror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_eoserror_CD;


typedef struct EXCEPTIONS_epackageerror_OD_S* EXCEPTIONS_epackageerror; /* sy_class */
typedef struct EXCEPTIONS_epackageerror_OD_S {  /* Objects of 'epackageerror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_epackageerror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_epackageerror_OD;

extern void * const EXCEPTIONS_epackageerror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_epackageerror_CD;


typedef struct EXCEPTIONS_esafecallexception_OD_S* 
  EXCEPTIONS_esafecallexception; /* sy_class */
typedef struct EXCEPTIONS_esafecallexception_OD_S {  /* Objects of 'esafecallexception' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_esafecallexception_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_esafecallexception_OD;

extern void * const EXCEPTIONS_esafecallexception_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_esafecallexception_CD;


typedef struct EXCEPTIONS_evarianterror_OD_S* EXCEPTIONS_evarianterror; /* sy_class */
typedef struct EXCEPTIONS_evarianterror_OD_S {  /* Objects of 'evarianterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_evarianterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_evarianterror_OD;

extern void * const EXCEPTIONS_evarianterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_evarianterror_CD;


typedef struct EXCEPTIONS_estreamerror_OD_S* EXCEPTIONS_estreamerror; /* sy_class */
typedef struct EXCEPTIONS_estreamerror_OD_S {  /* Objects of 'estreamerror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_estreamerror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_estreamerror_OD;

extern void * const EXCEPTIONS_estreamerror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_estreamerror_CD;


typedef struct EXCEPTIONS_efcreateerror_OD_S* EXCEPTIONS_efcreateerror; /* sy_class */
typedef struct EXCEPTIONS_efcreateerror_OD_S {  /* Objects of 'efcreateerror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_efcreateerror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_efcreateerror_OD;

extern void * const EXCEPTIONS_efcreateerror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_efcreateerror_CD;


typedef struct EXCEPTIONS_efopenerror_OD_S* EXCEPTIONS_efopenerror; /* sy_class */
typedef struct EXCEPTIONS_efopenerror_OD_S {  /* Objects of 'efopenerror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_efopenerror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_efopenerror_OD;

extern void * const EXCEPTIONS_efopenerror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_efopenerror_CD;


typedef struct EXCEPTIONS_elisterror_OD_S* EXCEPTIONS_elisterror; /* sy_class */
typedef struct EXCEPTIONS_elisterror_OD_S {  /* Objects of 'elisterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_elisterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_elisterror_OD;

extern void * const EXCEPTIONS_elisterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_elisterror_CD;


typedef struct EXCEPTIONS_epropertyconverterror_OD_S* 
  EXCEPTIONS_epropertyconverterror; /* sy_class */
typedef struct EXCEPTIONS_epropertyconverterror_OD_S {  /* Objects of 'epropertyconverterror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_epropertyconverterror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_epropertyconverterror_OD;

extern void * const EXCEPTIONS_epropertyconverterror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_epropertyconverterror_CD;


typedef struct EXCEPTIONS_epropertyerror_OD_S* 
  EXCEPTIONS_epropertyerror; /* sy_class */
typedef struct EXCEPTIONS_epropertyerror_OD_S {  /* Objects of 'epropertyerror' */
  SYSTEM_classreference_t CD;  /* = &EXCEPTIONS_epropertyerror_CD */
  SYSTEM_shortstring SYSTEM_exception_DOT_message;
} EXCEPTIONS_epropertyerror_OD;

extern void * const EXCEPTIONS_epropertyerror_VT[];
extern const SYSTEM_classdescriptor_t EXCEPTIONS_epropertyerror_CD;


typedef SYSTEM_byte EXCEPTIONS_texceptions;
/* Anonymous */
enum{EXCEPTIONS_cnoexception,EXCEPTIONS_cexception,
     EXCEPTIONS_cexternal,EXCEPTIONS_caccessviolation,EXCEPTIONS_cabort,
     EXCEPTIONS_cmatherror,EXCEPTIONS_ccontrolc,EXCEPTIONS_cinouterror,
     EXCEPTIONS_cconverterror,EXCEPTIONS_cinvalidcast,EXCEPTIONS_cinterror,
     EXCEPTIONS_crangeerror,EXCEPTIONS_cassertionfailed,
     EXCEPTIONS_cheapexception,EXCEPTIONS_coutofmemory,
     EXCEPTIONS_cabstracterror};

Function(SYSTEM_tobject ) EXCEPTIONS_create_exception_by_code(
  EXCEPTIONS_texceptions code,
  const SYSTEM_ansichar *msg);
/**** C code included from exceptions.pas(89:1): 6 lines ****/
/* All the exception stuff from p3io.c which is called from apps using
   exceptions (try statements) must be prototyped in exceptions.h: */

extern void _P3_Create_Exception(EXCEPTIONS_texceptions code, char* CMsg);
extern void _P3_Free_Exception();
extern void _P3_Std_Exception_Handler(SYSTEM_tobject obj);
typedef SYSTEM_byte EXCEPTIONS_p3_signal; /* Anonymous */ enum{EXCEPTIONS_p3_sigdummy};

extern void _Init_Module_exceptions(void);
extern void _Final_Module_exceptions(void);

#endif /* ! defined _P3___exceptions___H */
