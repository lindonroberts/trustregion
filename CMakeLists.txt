cmake_minimum_required (VERSION 3.13)

SET(TRUSTREGION_VERSION_MAJOR 0)
SET(TRUSTREGION_VERSION_MINOR 1)
SET(TRUSTREGION_VERSION_PATCH 0)  
SET(TRUSTREGION_VERSION ${TRUSTREGION_VERSION_MAJOR}.${TRUSTREGION_VERSION_MINOR}.${TRUSTREGION_VERSION_PATCH})

set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type")

project (trustregion 
         VERSION ${TRUSTREGION_VERSION_MAJOR}.${TRUSTREGION_VERSION_MINOR}
         LANGUAGES C Fortran)

option (BUILD_SHARED_LIBS "shared/static" ON)

option(BUILD_TESTING "Build tests" ON)

add_subdirectory(src)
add_subdirectory(interfaces/c)

if(BUILD_TESTING)
    add_subdirectory(examples/fortran)
    add_subdirectory(examples/c)
endif()
