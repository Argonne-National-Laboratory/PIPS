#if ! defined(_GDXSTATIC_H_)
#     define  _GDXSTATIC_H_

#include "gclgms.h"

struct gdxRec;
typedef struct gdxRec *gdxHandle_t;

#define gdxAcronymCount gdxacronymcount
#define gdxAcronymGetInfo cgdxacronymgetinfo
#define gdxAcronymName cgdxacronymname
#define gdxAddSetText cgdxaddsettext
#define gdxClose gdxclose
#define gdxDataReadDone gdxdatareaddone
#define gdxDataReadRaw gdxdatareadraw
#define gdxDataReadRawStart gdxdatareadrawstart
#define gdxDataReadStr cgdxdatareadstr
#define gdxDataReadStrStart gdxdatareadstrstart
#define gdxDataWriteDone gdxdatawritedone
#define gdxDataWriteRaw gdxdatawriteraw
#define gdxDataWriteRawStart cgdxdatawriterawstart
#define gdxDataWriteStr cgdxdatawritestr
#define gdxDataWriteStrStart cgdxdatawritestrstart
#define gdxErrorStr cgdxerrorstr
#define gdxFileVersion cgdxfileversion
#define gdxFileVersion cgdxfileversion
#define gdxFindSymbol cgdxfindsymbol
#define gdxGetElemText cgdxgetelemtext
#define gdxGetLastError gdxgetlasterror
#define gdxGetDLLVersion cgdxgetdllversion
#define gdxMapValue gdxmapvalue
#define gdxOpenRead cgdxopenread
#define gdxOpenWrite cgdxopenwrite
#define gdxSymbIndxMaxLength gdxsymbindxmaxlength
#define gdxSymbMaxLength gdxsymbmaxlength
#define gdxSymbolGetComment cgdxsymbolgetcomment
#define gdxSymbolGetDomain gdxsymbolgetdomain
#define gdxSymbolGetDomainX cgdxsymbolgetdomainx
#define gdxSymbolDim gdxsymboldim
#define gdxSymbolInfo cgdxsymbolinfo
#define gdxSymbolInfoX cgdxsymbolinfox
#define gdxSystemInfo gdxsysteminfo
#define gdxUELRegisterDone gdxuelregisterdone
#define gdxUELRegisterRawStart gdxuelregisterrawstart
#define gdxUELRegisterRaw cgdxuelregisterraw
#define gdxUMUelGet cgdxumuelget

#if defined(__cplusplus)
extern "C" {
#endif

int gdxCreate    (gdxHandle_t *pgdx, char *msgBuf, int msgBufLen);
int gdxFree      (gdxHandle_t *pgdx);
int _P3_DllInit();

int gdxAcronymCount (gdxHandle_t pgdx);
int gdxAcronymGetInfo (gdxHandle_t pgdx, int N, char *AName, char *Txt, int *AIndx);
int gdxAcronymName (gdxHandle_t pgdx, double V, char *AName);
int gdxAddSetText (gdxHandle_t pgdx, const char *Txt, int *TxtNr);
int gdxClose (gdxHandle_t pgdx);
int gdxDataReadDone (gdxHandle_t pgdx);
int gdxDataReadRaw (gdxHandle_t pgdx, int KeyInt[], double Values[], int *DimFrst);
int gdxDataReadRawStart (gdxHandle_t pgdx, int SyNr, int *NrRecs);
int gdxDataReadStr (gdxHandle_t pgdx, char *KeyStr[], double Values[], int *DimFrst);
int gdxDataReadStrStart (gdxHandle_t pgdx, int SyNr, int *NrRecs);
int gdxDataWriteDone (gdxHandle_t pgdx);
int gdxDataWriteRaw (gdxHandle_t pgdx, const int KeyInt[], const double Values[]);
int gdxDataWriteRawStart (gdxHandle_t pgdx, const char *SyId, const char *ExplTxt, int Dimen, int Typ, int UserInfo);
int gdxDataWriteStr (gdxHandle_t pgdx, const char *KeyStr[], const double Values[]);
int gdxDataWriteStrStart (gdxHandle_t pgdx, const char *SyId, const char *ExplTxt, int Dimen, int Typ, int UserInfo);
int gdxGetDLLVersion (gdxHandle_t pgdx, char *V);
int gdxErrorStr (gdxHandle_t pgdx, int ErrNr, char *ErrMsg);
int gdxFileVersion (gdxHandle_t pgdx, char *FileStr, char *ProduceStr);
int gdxFindSymbol (gdxHandle_t pgdx, const char *SyId, int *SyNr);
int gdxGetElemText (gdxHandle_t pgdx, int TxtNr, char *Txt, int *Node);
int gdxGetLastError (gdxHandle_t pgdx);
int gdxMapValue (gdxHandle_t pgdx, double D, int *sv);
int gdxOpenRead (gdxHandle_t pgdx, const char *FileName, int *ErrNr);
int gdxOpenWrite (gdxHandle_t pgdx, const char *FileName, const char *Producer, int *ErrNr);
int gdxSymbIndxMaxLength (gdxHandle_t pgdx, int SyNr, int LengthInfo[]);
int gdxSymbMaxLength (gdxHandle_t pgdx);
int gdxSymbolGetComment (gdxHandle_t pgdx, int SyNr, int N, char *Txt);
int gdxSymbolGetDomain (gdxHandle_t pgdx, int SyNr, int DomainSyNrs[]);
int gdxSymbolGetDomainX (gdxHandle_t pgdx, int SyNr, char *DomainIDs[]);
int gdxSymbolDim (gdxHandle_t pgdx, int SyNr);
int gdxSymbolInfo (gdxHandle_t pgdx, int SyNr, char *SyId, int *Dimen, int *Typ);
int gdxSymbolInfoX (gdxHandle_t pgdx, int SyNr, int *RecCnt, int *UserInfo, char *ExplTxt);
int gdxSystemInfo (gdxHandle_t pgdx, int *SyCnt, int *UelCnt);
int gdxUELRegisterDone (gdxHandle_t pgdx);
int gdxUELRegisterRaw (gdxHandle_t pgdx, const char *Uel);
int gdxUELRegisterRawStart (gdxHandle_t pgdx);
int gdxUMUelGet (gdxHandle_t pgdx, int UelNr, char *Uel, int *UelMap);

#if defined(__cplusplus)
}
#endif

#endif /* #if ! defined(_GDXSTATIC_H_) */
