include_directories(../)

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel1 parmodel1.cpp)
target_link_libraries(parmodel1
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel parmodel.cpp)
target_link_libraries(parmodel
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel2 parmodel2.cpp)
target_link_libraries(parmodel2
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel3 parmodel3.cpp)
target_link_libraries(parmodel3
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel9 parmodel9.cpp)
target_link_libraries(parmodel9
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel9_yk parmodel9_yk.cpp)
target_link_libraries(parmodel9_yk
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parleastsquares parleastsquares.cpp)
target_link_libraries(parleastsquares
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel_lp parmodel_eqieq_lp.cpp)
target_link_libraries(parmodel_lp
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )


# a testing example for PIPSNLP parallel structure interface
add_executable(parmodel_1STonly parmodel_1STonly.cpp)
target_link_libraries(parmodel_1STonly
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack  
  LU_solvers
  nlpsparse 
  nlpdense
  pipsnlpFLSStoch 
  nlpstoch 
  nlpstochla 
  nlpgensparse 
  nlpbase 
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

# Test for problem with linking constraints 
add_executable(parmodelLink parmodelLink.cpp)
target_link_libraries(parmodelLink
  parpipsnlp
  reducedSolverSpecial
  ${WHOLE_ARCHIVE}
  updateFromCallBack
  LU_solvers
  nlpsparse
  nlpdense
  pipsnlpFLSStoch
  nlpstoch
  nlpstochla
  nlpgensparse
  nlpbase
  genStructparallel
  structStochInfo
  global_var
  ${COIN_LIBS}
  ${MA57_LIBRARY}
  ${MA27_LIBRARY}
  ${UMF_ALLINONE_LIBRARY}
  ${METIS_LIBRARY}
  ${PARDISO_LIBRARY}
  ${NO_WHOLE_ARCHIVE}
  ${MATH_LIBS} )

