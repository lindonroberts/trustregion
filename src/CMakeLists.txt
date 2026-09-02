
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
  trustregion.F90
)

set (PRIMA_REAL_PRECISION "64" CACHE STRING "real precision")
set (PRIMA_INTEGER_KIND "0" CACHE STRING "integer kind")
target_compile_definitions (trustregion PUBLIC "PRIMA_REAL_PRECISION=${PRIMA_REAL_PRECISION}")
target_compile_definitions (trustregion PUBLIC "PRIMA_INTEGER_KIND=${PRIMA_INTEGER_KIND}")
set_target_properties(trustregion PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  Fortran_PREPROCESS ON
  Fortran_MODULE_DIRECTORY mod)
target_include_directories (trustregion PRIVATE prima_common)
target_include_directories (trustregion PUBLIC
  $<INSTALL_INTERFACE:include/prima/mod>
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/mod>)

# because of calcfc_internal nested subroutine in cobylb.f90
if (HAVE_WARN_EXECSTACK)
  target_link_options (trustregionf PUBLIC "-Wl,--no-warn-execstack")
endif ()

# On Windows, when building shared libs, visible functions must be explicitly listed.
# Unlike C there are no standard keywords to do this in Fortran but only compiler-specific pragmas:
# !GCC$ attributes dllexport :: bobyqa_c
# !DEC$ attributes dllexport :: bobyqa_c
# the downside is that gfortran will raise -Wattributes warnings, and Intel compilers will complain
# that the directives are not standard conforming. See https://github.com/libprima/prima/issues/70.
# Another option is to use .def files referencing exported symbols:
# https://learn.microsoft.com/en-us/cpp/build/exporting-from-a-dll-using-def-files?view=msvc-170
# Mangling is different between GNU and Intel/LLVM compiler families so we use a dedicated file for each.
# Symbol names can be seen in object files with the dumpbin tool (or objdump) to write a new .def:
# dumpbin /symbols cobyla.f90.obj
# 017 00000000 UNDEF  no_type       External     | COBYLA_MOD_mp_COBYLA
# x86_64-w64-mingw32-objdump -t cobyla.f90.obj
# [  4](sec  1)(fl 0x00)(ty   20)(scl   2) (nx 0) 0x000000000000218b __cobyla_mod_MOD_cobyla
if (WIN32 AND BUILD_SHARED_LIBS)
  if (EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/trustregion-${CMAKE_Fortran_COMPILER_ID}.def)
    target_sources (trustregion PRIVATE trustregion-${CMAKE_Fortran_COMPILER_ID}.def)
  else ()
    message (STATUS "Assuming Intel symbol mangling")
    target_sources (trustregion PRIVATE trustregion-Intel.def)
  endif ()
endif ()

if (WIN32)
  set_target_properties(trustregion PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
endif()

install (DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/trustregion)
