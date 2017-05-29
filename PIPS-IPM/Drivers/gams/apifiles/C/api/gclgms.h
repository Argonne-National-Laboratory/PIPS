/* global constants (symbol dimensions, symbol layout, etc.)
 * that might ultimately come from global
 * created: Steve Dirkse, July 2007
 */

#if ! defined(_GCLGMS_H_)
#     define  _GCLGMS_H_

#define GLOBAL_MAX_INDEX_DIM       20
#define GLOBAL_UEL_IDENT_SIZE      64  /* implies len of 63 */
#define ITERLIM_INFINITY           2000000000
#define GMS_MAX_SOLVERS            150 /* maximum number of solver that can be stored */

#define GMS_MAX_INDEX_DIM      GLOBAL_MAX_INDEX_DIM
#define GMS_UEL_IDENT_SIZE     GLOBAL_UEL_IDENT_SIZE  /* implies len of 63 */

#define GMS_SSSIZE            256 /* short string size */

#define GMS_VARTYPE_UNKNOWN  0
#define GMS_VARTYPE_BINARY   1
#define GMS_VARTYPE_INTEGER  2
#define GMS_VARTYPE_POSITIVE 3
#define GMS_VARTYPE_NEGATIVE 4
#define GMS_VARTYPE_FREE     5
#define GMS_VARTYPE_SOS1     6
#define GMS_VARTYPE_SOS2     7
#define GMS_VARTYPE_SEMICONT 8
#define GMS_VARTYPE_SEMIINT  9
#define GMS_VARTYPE_MAX      10

/* base value used by GDX to store equation types:
 * returned as userInfo value by gdxSymbolInfoX */
#define GMS_EQU_USERINFO_BASE 53

#define GMS_EQUTYPE_E        0
#define GMS_EQUTYPE_G        1
#define GMS_EQUTYPE_L        2
#define GMS_EQUTYPE_N        3
#define GMS_EQUTYPE_X        4
#define GMS_EQUTYPE_C        5
#define GMS_EQUTYPE_MAX      6

#define GMS_VAL_LEVEL    0
#define GMS_VAL_MARGINAL 1
#define GMS_VAL_LOWER    2
#define GMS_VAL_UPPER    3
#define GMS_VAL_SCALE    4
#define GMS_VAL_MAX      5

enum gdxSpecValue {
  sv_valund  = 0,
  sv_valna   = 1,
  sv_valpin  = 2,
  sv_valmin  = 3,
  sv_valeps  = 4,
  sv_normal  = 5,
  sv_acronym = 6  };

#define GMS_SVIDX_UNDEF  0
#define GMS_SVIDX_NA     1
#define GMS_SVIDX_PINF   2
#define GMS_SVIDX_MINF   3
#define GMS_SVIDX_EPS    4
#define GMS_SVIDX_NORMAL 5
#define GMS_SVIDX_ACR    6
#define GMS_SVIDX_MAX    7

enum gdxSyType {
  dt_set   = 0,
  dt_par   = 1,
  dt_var   = 2,
  dt_equ   = 3,
  dt_alias = 4  };

#define GMS_DT_SET       0
#define GMS_DT_PAR       1
#define GMS_DT_VAR       2
#define GMS_DT_EQU       3
#define GMS_DT_ALIAS     4
#define GMS_DT_MAX       5

#define GMS_SV_UNDEF  1.0E300   /* undefined                */
#define GMS_SV_NA     2.0E300   /* not available/applicable */
#define GMS_SV_PINF   3.0E300   /* plus infinity            */
#define GMS_SV_MINF   4.0E300   /* minus infinity           */
#define GMS_SV_EPS    5.0E300   /* epsilon                  */
#define GMS_SV_ACR   10.0E300   /* potential/real acronym   */
#define GMS_SV_NAINT 2100000000 /* not available/applicable for integers */

#define STAT_OK    0
#define STAT_NOPT  1
#define STAT_INFES 2
#define STAT_UNBND 3
#define STAT_EVAL  4
#define STAT_UNKNW 5
#define STAT_REDEF 6
#define STAT_DEPND 7
#define STAT_REDIR 8
#define STAT_MAX   9

#define SS_MAX 14
#define MS_MAX 20

#if defined(__cplusplus)
extern "C" {
#endif

typedef double gdxSVals_t[GMS_SVIDX_MAX];
typedef double gdxValues_t[GMS_VAL_MAX];
typedef char gdxStrIndex_t[GMS_MAX_INDEX_DIM][GMS_SSSIZE];
typedef char *gdxStrIndexPtrs_t[GMS_MAX_INDEX_DIM];
typedef int gdxUelIndex_t[GMS_MAX_INDEX_DIM];
#define GDXSTRINDEXPTRS_INIT(idx,idxPtrs) do {int GDXSTRINDEXPTRS_i; for (GDXSTRINDEXPTRS_i=0; GDXSTRINDEXPTRS_i < GMS_MAX_INDEX_DIM;  GDXSTRINDEXPTRS_i++) (idxPtrs)[GDXSTRINDEXPTRS_i] = (idx)[GDXSTRINDEXPTRS_i];} while (0)

extern const char *gmsGdxTypeText[GMS_DT_MAX];
extern const char *gmsVarTypeText[GMS_VARTYPE_MAX];
extern const char *gmsValTypeText[GMS_VAL_MAX];
extern const char *gmsSVText[GMS_SVIDX_MAX];

extern const double gmsDefRecVar[GMS_VARTYPE_MAX][GMS_VAL_MAX];
extern const double gmsDefRecEqu[GMS_EQUTYPE_MAX][GMS_VAL_MAX];

extern const char * const rcStat[STAT_MAX];
extern const char * const solveStatusTxt[SS_MAX];
extern const char * const modelStatusTxt[MS_MAX];

int
gmsFixEquType (int userInfo);
int
gmsFixVarType (int userInfo);
#if defined(__cplusplus)
}
#endif

#endif /* if ! defined(_GCLGMS_H_) */
