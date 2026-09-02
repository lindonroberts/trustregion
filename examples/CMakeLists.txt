
FILE( GLOB EXAMPLES *.F90 )

FOREACH( example_file ${EXAMPLES} )
    # Cut off the file extension and directory path
    GET_FILENAME_COMPONENT( example_exec ${example_file} NAME_WE )
	
    ADD_EXECUTABLE( ${example_exec} ${example_file} )
	
	TARGET_LINK_LIBRARIES( ${example_exec} PRIVATE trustregion )

ENDFOREACH( example_file ${EXAMPLES} )

#add_executable (ex_trsapp ex_trsapp.F90)
#target_link_libraries (ex_trsapp PRIVATE trustregion)
#target_include_directories (ex_trsapp PRIVATE ${CMAKE_BINARY_DIR}/src)
#set_target_properties(ex_trsapp PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/examples/${name}_${number}/mod)
#if (WIN32)
#  set_target_properties(ex_trsapp PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)
#endif()
#
#install (TARGETS ex_trsapp)
