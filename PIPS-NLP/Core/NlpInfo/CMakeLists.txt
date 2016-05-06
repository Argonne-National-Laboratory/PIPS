include_directories(../..)

add_library(updateFromAMPL
  getAmplFunction.C NlpInfoAMPL.C
)

add_library(sInfoQPLP
  sNlpInfoFIX.C NlpInfoFIX.C NlpInfo.C
)

add_library(StochInfoFromAMPL
   getAmplFunctionNew.C sNlpInfoFromNL.C
)

add_library(structStochInfo
	StructJuMPsInfo.C
)

add_library(updateFromCallBack
   NlpInfoCallBack.C
)

