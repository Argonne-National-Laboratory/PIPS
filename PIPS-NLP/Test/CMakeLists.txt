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
