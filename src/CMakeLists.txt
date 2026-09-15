
add_library (trustregion
  prima_common/consts.F90  
  prima_common/debug.F90
  prima_common/huge.F90
  prima_common/inf.F90
  prima_common/infnan.F90
  prima_common/infos.F90
  prima_common/linalg.F90
  prima_common/memory.F90
  prima_common/powalg.F90
  prima_common/string.F90
  prima_common/univar.F90
  getact.F90
  qdec.F90
  trsapp.F90
  trsbox.F90
  trslin.F90
  trsunc.F90
  trustregion.F90
)

set (PRIMA_REAL_PRECISION "64" CACHE STRING "real precision")
set (PRIMA_INTEGER_KIND "0" CACHE STRING "integer kind")
target_compile_definitions (trustregion PUBLIC "PRIMA_REAL_PRECISION=${PRIMA_REAL_PRECISION}")
target_compile_definitions (trustregion PUBLIC "PRIMA_INTEGER_KIND=${PRIMA_INTEGER_KIND}")
set_target_properties(trustregion PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_PREPROCESS ON
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories (trustregion PRIVATE prima_common)
#target_include_directories (trustregion PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/modules)
target_include_directories (trustregion PUBLIC
  $<INSTALL_INTERFACE:include/trustregion/modules>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>)

message(WARNING "module dir ${CMAKE_CURRENT_BINARY_DIR}/modules")
message(WARNING "install include dir ${CMAKE_INSTALL_INCLUDEDIR}/trustregion")

install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/trustregion)
