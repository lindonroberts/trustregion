
# Find the module development requirements (requires FindPython from 3.17 or
# scikit-build-core's built-in backport)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module NumPy)
find_package(pybind11 CONFIG REQUIRED)

# Variables TRUSTREGION_INCLUDE_DIR and TRUSTREGION_LINK_DIR can be used to set the locations
# of trustregion.h and the main trustregion library using pip like 
#     $ pip install . --config-settings=cmake.define.TRUSTREGION_INCLUDE_DIR=/path/to/trustregion
# If trustregion is installed in a standard system location this should not be needed
# If the variables are not set by pip, they point to the local build locations
if(NOT TRUSTREGION_INCLUDE_DIR)
	set(TRUSTREGION_INCLUDE_DIR "${CMAKE_SOURCE_DIR}/interfaces/c")
endif()
#if(NOT TRUSTREGION_C_API_LINK_DIR)
#	set(TRUSTREGION_C_API_LINK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../build/interfaces/c")
#endif()
#if(NOT TRUSTREGION_LINK_DIR)
#	set(TRUSTREGION_LINK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../build/src")
#endif()

# Add a library using FindPython's tooling (pybind11 also provides a helper like this)
#python_add_library(_trustregion MODULE src/pytrustregion.cpp WITH_SOABI)
pybind11_add_module(_trustregion "${CMAKE_CURRENT_SOURCE_DIR}/src/pytrustregion.cpp")
target_include_directories(_trustregion PRIVATE ${TRUSTREGION_INCLUDE_DIR})
#target_link_libraries(_trustregion PRIVATE pybind11::headers trustregion_c_api trustregion)
target_link_libraries(_trustregion PRIVATE trustregion_c_api trustregion)
#target_link_directories(_trustregion PRIVATE ${TRUSTREGION_C_API_LINK_DIR} ${TRUSTREGION_LINK_DIR})

# This is passing in the version as a define just as an example
target_compile_definitions(_trustregion PRIVATE VERSION_INFO=${PROJECT_VERSION})

# The install directory is the output (wheel) directory
install(TARGETS _trustregion DESTINATION trustregion/ COMPONENT trustregion_python_c_ext)

# Install supplementary .py files such as __init__.py
file(GLOB SUPPORTING_PY_FILES "${CMAKE_CURRENT_SOURCE_DIR}/src/*.py")

message(WARNING "extra files ${SUPPORTING_PY_FILES}")

# On Windows, the compiled _trustregion module (.pyd file) must have the 
# DLL file in the same directory (not just %PATH%), so manually copy
# it here
#if(WIN32)
#    list(APPEND SUPPORTING_PY_FILES "${TRUSTREGION_C_API_LINK_DIR}/trustregion_c_api.dll")
#    list(APPEND SUPPORTING_PY_FILES "${TRUSTREGION_LINK_DIR}/trustregion.dll")
#else()
#    list(APPEND SUPPORTING_PY_FILES "${TRUSTREGION_C_API_LINK_DIR}/libtrustregion_c_api.so")
#    list(APPEND SUPPORTING_PY_FILES "${TRUSTREGION_LINK_DIR}/libtrustregion.so")
#endif()

install (FILES ${SUPPORTING_PY_FILES} DESTINATION trustregion/)
