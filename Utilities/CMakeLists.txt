
set(EXECUTABLE_OUTPUT_PATH "${PROJECT_BINARY_DIR}")


if(HAVE_SML)

add_executable(dumpSmlModel dumpSmlModel.cpp)
target_link_libraries(dumpSmlModel ${SML_LIBRARY} ${AMPL_LIBRARY} ${DL_LIBRARY} ${FL_LIBRARY})

endif(HAVE_SML)



