cmake_minimum_required(VERSION 3.17.2)
#project(${SKBUILD_PROJECT_NAME} LANGUAGES C Fortran)
project(trustregion LANGUAGES C Fortran)

find_package(
  Python
  COMPONENTS Interpreter Development.Module NumPy
  REQUIRED)

message("python executable is ${Python_EXECUTABLE}")

# F2PY headers
execute_process(
  COMMAND "${Python_EXECUTABLE}" -c
          "import numpy.f2py; print(numpy.f2py.get_include())"
  OUTPUT_VARIABLE F2PY_INCLUDE_DIR
  OUTPUT_STRIP_TRAILING_WHITESPACE)

message("f2py dir is ${F2PY_INCLUDE_DIR}")

add_library(fortranobject OBJECT "${F2PY_INCLUDE_DIR}/fortranobject.c")
target_link_libraries(fortranobject PUBLIC Python::NumPy)
target_include_directories(fortranobject PUBLIC "${F2PY_INCLUDE_DIR}")
set_property(TARGET fortranobject PROPERTY POSITION_INDEPENDENT_CODE ON)

file(GLOB F90_SOURCE_FILES 
  "${CMAKE_CURRENT_SOURCE_DIR}/../../src/*.F90" 
  "${CMAKE_CURRENT_SOURCE_DIR}/../../src/prima_common/*.F90"
)

#add_custom_command(
#  OUTPUT _trsmodule.c _trs-f2pywrappers.f
#  DEPENDS ${F90_SOURCE_FILES}
#  VERBATIM
#  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
#          "${CMAKE_CURRENT_SOURCE_DIR}/../../src/trustregion.F90" -m _trs --lower)

add_custom_command(
  OUTPUT _trsmodule.c _trs-f2pywrappers.f
  DEPENDS ${F90_SOURCE_FILES}
  VERBATIM
  COMMAND "${Python_EXECUTABLE}" -m numpy.f2py
          -c "${F90_SOURCE_FILES}" -m _trs --lower)

python_add_library(_trs MODULE "${CMAKE_CURRENT_BINARY_DIR}/_trsmodule.c"
                   "${F90_SOURCE_FILES}" WITH_SOABI)
target_link_libraries(_trs PRIVATE fortranobject)

install(TARGETS _trs DESTINATION trustregion/)

# Taken from PRIMA:
# The following are not added to the component because scikit-build-core automatically
# detects files in python/<package_name> and adds them to the wheel. These commands
# are here in case one is building via cmake directly and not via scikit-build-core
file(GLOB SUPPORTING_PY_FILES "${CMAKE_CURRENT_SOURCE_DIR}/trustregion/*.py")
install (FILES ${SUPPORTING_PY_FILES} DESTINATION trustregion)
