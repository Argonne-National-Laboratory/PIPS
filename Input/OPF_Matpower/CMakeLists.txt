


add_library(dcopffactory dcopflow.cpp)

add_library(opfcinput ps.cpp)

add_library(graphpart graphPart.cpp)

add_library(dcopfGP ps.cpp dcopflow.cpp graphPart.cpp)


#message(STATUS "------ ${METIS_LIBRARY}")
add_executable(generateGraph tempMain.cpp)
  target_link_libraries(generateGraph 
     opfcinput graphpart
     ${METIS_LIBRARY} ${MATH_LIBS})
     
