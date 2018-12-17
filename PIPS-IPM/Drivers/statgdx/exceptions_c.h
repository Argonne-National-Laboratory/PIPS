#include "p3io.h"
#include "exceptions.h"


void * const EXCEPTIONS_eabort_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eabort' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eabort_CD = {
  _P3str1("\006eabort"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eabort_OD), EXCEPTIONS_eabort_VT, NULL};


void * const EXCEPTIONS_eabstracterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eabstracterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eabstracterror_CD = {
  _P3str1("\016eabstracterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eabstracterror_OD), EXCEPTIONS_eabstracterror_VT, NULL};


void * const EXCEPTIONS_eexternal_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eexternal' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eexternal_CD = {
  _P3str1("\011eexternal"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eexternal_OD), EXCEPTIONS_eexternal_VT, NULL};


void * const EXCEPTIONS_eaccessviolation_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eaccessviolation' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eaccessviolation_CD = {
  _P3str1("\020eaccessviolation"),
  &EXCEPTIONS_eexternal_CD, NULL, 0,
  sizeof(EXCEPTIONS_eaccessviolation_OD),
  EXCEPTIONS_eaccessviolation_VT, NULL};


void * const EXCEPTIONS_ematherror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'ematherror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_ematherror_CD = {
  _P3str1("\012ematherror"),
  &EXCEPTIONS_eexternal_CD, NULL, 0,
  sizeof(EXCEPTIONS_ematherror_OD), EXCEPTIONS_ematherror_VT, NULL};


void * const EXCEPTIONS_einvalidargument_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'einvalidargument' */
const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidargument_CD = {
  _P3str1("\020einvalidargument"),
  &EXCEPTIONS_ematherror_CD, NULL, 0,
  sizeof(EXCEPTIONS_einvalidargument_OD),
  EXCEPTIONS_einvalidargument_VT, NULL};


void * const EXCEPTIONS_einvalidop_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'einvalidop' */
const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidop_CD = {
  _P3str1("\012einvalidop"),
  &EXCEPTIONS_ematherror_CD, NULL, 0,
  sizeof(EXCEPTIONS_einvalidop_OD), EXCEPTIONS_einvalidop_VT, NULL};


void * const EXCEPTIONS_eoverflow_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eoverflow' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eoverflow_CD = {
  _P3str1("\011eoverflow"),
  &EXCEPTIONS_ematherror_CD, NULL, 0,
  sizeof(EXCEPTIONS_eoverflow_OD), EXCEPTIONS_eoverflow_VT, NULL};


void * const EXCEPTIONS_eunderflow_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eunderflow' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eunderflow_CD = {
  _P3str1("\012eunderflow"),
  &EXCEPTIONS_ematherror_CD, NULL, 0,
  sizeof(EXCEPTIONS_eunderflow_OD), EXCEPTIONS_eunderflow_VT, NULL};


void * const EXCEPTIONS_ezerodivide_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'ezerodivide' */
const SYSTEM_classdescriptor_t EXCEPTIONS_ezerodivide_CD = {
  _P3str1("\013ezerodivide"),
  &EXCEPTIONS_ematherror_CD, NULL, 0,
  sizeof(EXCEPTIONS_ezerodivide_OD), EXCEPTIONS_ezerodivide_VT, NULL};


void * const EXCEPTIONS_econtrolc_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'econtrolc' */
const SYSTEM_classdescriptor_t EXCEPTIONS_econtrolc_CD = {
  _P3str1("\011econtrolc"),
  &EXCEPTIONS_eexternal_CD, NULL, 0,
  sizeof(EXCEPTIONS_econtrolc_OD), EXCEPTIONS_econtrolc_VT, NULL};


void * const EXCEPTIONS_einterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'einterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_einterror_CD = {
  _P3str1("\011einterror"),
  &EXCEPTIONS_eexternal_CD, NULL, 0,
  sizeof(EXCEPTIONS_einterror_OD), EXCEPTIONS_einterror_VT, NULL};


void * const EXCEPTIONS_edivbyzero_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'edivbyzero' */
const SYSTEM_classdescriptor_t EXCEPTIONS_edivbyzero_CD = {
  _P3str1("\012edivbyzero"),
  &EXCEPTIONS_einterror_CD, NULL, 0,
  sizeof(EXCEPTIONS_edivbyzero_OD), EXCEPTIONS_edivbyzero_VT, NULL};


void * const EXCEPTIONS_erangeerror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'erangeerror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_erangeerror_CD = {
  _P3str1("\013erangeerror"),
  &EXCEPTIONS_einterror_CD, NULL, 0,
  sizeof(EXCEPTIONS_erangeerror_OD), EXCEPTIONS_erangeerror_VT, NULL};


void * const EXCEPTIONS_eintoverflow_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eintoverflow' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eintoverflow_CD = {
  _P3str1("\014eintoverflow"),
  &EXCEPTIONS_einterror_CD, NULL, 0,
  sizeof(EXCEPTIONS_eintoverflow_OD), EXCEPTIONS_eintoverflow_VT, NULL};


void * const EXCEPTIONS_estackoverflow_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'estackoverflow' */
const SYSTEM_classdescriptor_t EXCEPTIONS_estackoverflow_CD = {
  _P3str1("\016estackoverflow"),
  &EXCEPTIONS_eexternal_CD, NULL, 0,
  sizeof(EXCEPTIONS_estackoverflow_OD), EXCEPTIONS_estackoverflow_VT, NULL};


void * const EXCEPTIONS_eassertionfailed_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eassertionfailed' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eassertionfailed_CD = {
  _P3str1("\020eassertionfailed"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eassertionfailed_OD),
  EXCEPTIONS_eassertionfailed_VT, NULL};


void * const EXCEPTIONS_econverterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'econverterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_econverterror_CD = {
  _P3str1("\015econverterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_econverterror_OD), EXCEPTIONS_econverterror_VT, NULL};


void * const EXCEPTIONS_eexternalexception_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eexternalexception' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eexternalexception_CD = {
  _P3str1("\022eexternalexception"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eexternalexception_OD),
  EXCEPTIONS_eexternalexception_VT, NULL};


void * const EXCEPTIONS_eheapexception_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eheapexception' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eheapexception_CD = {
  _P3str1("\016eheapexception"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eheapexception_OD), EXCEPTIONS_eheapexception_VT, NULL};


void * const EXCEPTIONS_eoutofmemory_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eoutofmemory' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eoutofmemory_CD = {
  _P3str1("\014eoutofmemory"),
  &EXCEPTIONS_eheapexception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eoutofmemory_OD), EXCEPTIONS_eoutofmemory_VT, NULL};


void * const EXCEPTIONS_einvalidpointer_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'einvalidpointer' */
const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidpointer_CD = {
  _P3str1("\017einvalidpointer"),
  &EXCEPTIONS_eheapexception_CD, NULL, 0,
  sizeof(EXCEPTIONS_einvalidpointer_OD), EXCEPTIONS_einvalidpointer_VT, NULL};


void * const EXCEPTIONS_einouterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'einouterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_einouterror_CD = {
  _P3str1("\013einouterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_einouterror_OD), EXCEPTIONS_einouterror_VT, NULL};


void * const EXCEPTIONS_eintfcasterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eintfcasterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eintfcasterror_CD = {
  _P3str1("\016eintfcasterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eintfcasterror_OD), EXCEPTIONS_eintfcasterror_VT, NULL};


void * const EXCEPTIONS_einvalidcast_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'einvalidcast' */
const SYSTEM_classdescriptor_t EXCEPTIONS_einvalidcast_CD = {
  _P3str1("\014einvalidcast"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_einvalidcast_OD), EXCEPTIONS_einvalidcast_VT, NULL};


void * const EXCEPTIONS_eoserror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'eoserror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_eoserror_CD = {
  _P3str1("\010eoserror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_eoserror_OD), EXCEPTIONS_eoserror_VT, NULL};


void * const EXCEPTIONS_epackageerror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'epackageerror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_epackageerror_CD = {
  _P3str1("\015epackageerror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_epackageerror_OD), EXCEPTIONS_epackageerror_VT, NULL};


void * const EXCEPTIONS_esafecallexception_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'esafecallexception' */
const SYSTEM_classdescriptor_t EXCEPTIONS_esafecallexception_CD = {
  _P3str1("\022esafecallexception"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_esafecallexception_OD),
  EXCEPTIONS_esafecallexception_VT, NULL};


void * const EXCEPTIONS_evarianterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'evarianterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_evarianterror_CD = {
  _P3str1("\015evarianterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_evarianterror_OD), EXCEPTIONS_evarianterror_VT, NULL};


void * const EXCEPTIONS_estreamerror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'estreamerror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_estreamerror_CD = {
  _P3str1("\014estreamerror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_estreamerror_OD), EXCEPTIONS_estreamerror_VT, NULL};


void * const EXCEPTIONS_efcreateerror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'efcreateerror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_efcreateerror_CD = {
  _P3str1("\015efcreateerror"),
  &EXCEPTIONS_estreamerror_CD, NULL, 0,
  sizeof(EXCEPTIONS_efcreateerror_OD), EXCEPTIONS_efcreateerror_VT, NULL};


void * const EXCEPTIONS_efopenerror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'efopenerror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_efopenerror_CD = {
  _P3str1("\013efopenerror"),
  &EXCEPTIONS_estreamerror_CD, NULL, 0,
  sizeof(EXCEPTIONS_efopenerror_OD), EXCEPTIONS_efopenerror_VT, NULL};


void * const EXCEPTIONS_elisterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'elisterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_elisterror_CD = {
  _P3str1("\012elisterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_elisterror_OD), EXCEPTIONS_elisterror_VT, NULL};


void * const EXCEPTIONS_epropertyconverterror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'epropertyconverterror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_epropertyconverterror_CD = {
  _P3str1("\025epropertyconverterror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_epropertyconverterror_OD),
  EXCEPTIONS_epropertyconverterror_VT, NULL};


void * const EXCEPTIONS_epropertyerror_VT[] = {(void*)&
  SYSTEM_tobject_DOT_destroy};

/* Class descriptor for 'epropertyerror' */
const SYSTEM_classdescriptor_t EXCEPTIONS_epropertyerror_CD = {
  _P3str1("\016epropertyerror"),
  &SYSTEM_exception_CD, NULL, 0,
  sizeof(EXCEPTIONS_epropertyerror_OD), EXCEPTIONS_epropertyerror_VT, NULL};


Function(SYSTEM_tobject ) EXCEPTIONS_create_exception_by_code(
  EXCEPTIONS_texceptions code,
  const SYSTEM_ansichar *msg)
{
  SYSTEM_tobject result;

  switch (ValueCast(EXCEPTIONS_texceptions,code)) {
  case EXCEPTIONS_cexception:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&SYSTEM_exception_CD)),msg));
    break;
  case EXCEPTIONS_cexternal:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_eexternal_CD)),
        msg));
    break;
  case EXCEPTIONS_caccessviolation:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&
        EXCEPTIONS_eaccessviolation_CD)),msg));
    break;
  case EXCEPTIONS_cabort:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_eabort_CD)),msg));
    break;
  case EXCEPTIONS_cmatherror:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_ematherror_CD)),
        msg));
    break;
  case EXCEPTIONS_ccontrolc:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_econtrolc_CD)),
        msg));
    break;
  case EXCEPTIONS_cinouterror:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_einouterror_CD)),
        msg));
    break;
  case EXCEPTIONS_cassertionfailed:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&
        EXCEPTIONS_eassertionfailed_CD)),msg));
    break;
  case EXCEPTIONS_cconverterror:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_econverterror_CD)),
        msg));
    break;
  case EXCEPTIONS_cinvalidcast:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_einvalidcast_CD)),
        msg));
    break;
  case EXCEPTIONS_cinterror:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_einterror_CD)),
        msg));
    break;
  case EXCEPTIONS_crangeerror:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_erangeerror_CD)),
        msg));
    break;
  case EXCEPTIONS_cheapexception:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_eheapexception_CD)),
        msg));
    break;
  case EXCEPTIONS_coutofmemory:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_eoutofmemory_CD)),
        msg));
    break;
  case EXCEPTIONS_cabstracterror:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&EXCEPTIONS_eabstracterror_CD)),
        msg));
    break;
  default:
    result = ValueCast(SYSTEM_tobject,SYSTEM_exception_DOT_create(ValueCast(
        SYSTEM_exception,_P3alloc_object(&SYSTEM_exception_CD)),msg));
  } /* switch */
  return result;
}  /* create_exception_by_code */

/* unit exceptions */
void _Init_Module_exceptions(void)
{
} /* _Init_Module_exceptions */

void _Final_Module_exceptions(void)
{
} /* _Final_Module_exceptions */
