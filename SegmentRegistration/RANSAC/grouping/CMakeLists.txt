#
# (C) 2010 Groupsac

# Library 
IF (WIN32)
  SET(LAPACK_BLAS_LIB ${LAPACK_BLAS_LIB} 
    ${CMAKE_SOURCE_DIR}/thirdParty/windowsLib/libf2c.lib
    ${CMAKE_SOURCE_DIR}/thirdParty/windowsLib/lapack.lib
    ${CMAKE_SOURCE_DIR}/thirdParty/windowsLib/blas.lib)
ELSE (WIN32)
  SET(LAPACK_BLAS_LIB ${LAPACK_BLAS_LIB} lapack blas)
ENDIF (WIN32)

ADD_EXECUTABLE(mainTestGrouping_GroupByFlows testGroupByFlows.cpp)
add_custom_target(mainTestGrouping_GroupByFlows.run mainTestGrouping_GroupByFlows)
TARGET_LINK_LIBRARIES(mainTestGrouping_GroupByFlows "CppUnitLite" ${LAPACK_BLAS_LIB})
add_test(mainTestGrouping_GroupByFlows mainTestGrouping_GroupByFlows)

