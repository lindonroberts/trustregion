cmake_minimum_required (VERSION 3.15)

SET(TRUSTREGION_VERSION_MAJOR 1)
SET(TRUSTREGION_VERSION_MINOR 0)
SET(TRUSTREGION_VERSION_PATCH 0)  
SET(TRUSTREGION_VERSION ${TRUSTREGION_VERSION_MAJOR}.${TRUSTREGION_VERSION_MINOR}.${TRUSTREGION_VERSION_PATCH})

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (trustregion 
         VERSION ${TRUSTREGION_VERSION_MAJOR}.${TRUSTREGION_VERSION_MINOR}
         LANGUAGES C Fortran)

option (BUILD_SHARED_LIBS "shared/static" ON)

option(BUILD_TESTING "Build tests" ON)

include(GNUInstallDirs)

add_subdirectory(src)
add_subdirectory(interfaces/c)

option (TRUSTREGION_ENABLE_PYTHON "Python binding" OFF)
if(TRUSTREGION_ENABLE_PYTHON)
    enable_language(CXX)
    add_subdirectory(interfaces/python)
endif()

if(BUILD_TESTING)
    add_subdirectory(examples/fortran)
    add_subdirectory(examples/c)
endif()

install(
    TARGETS trustregion trustregion_c_api
    EXPORT trustregion-targets
    INCLUDES DESTINATION include
)