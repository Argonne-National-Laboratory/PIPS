

add_library(updateFromAMPL
  getAmplFunction.C NlpInfoAMPL.C
)

add_library(sInfoQPLP
  sNlpInfoFIX.C
)

add_library(StochInfoFromAMPL
   getAmplFunctionNew.C sNlpInfoFromNL.C
)

add_library(updateFromCallBack
   NlpInfoCallBack.C
)

