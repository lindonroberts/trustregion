
FILE( GLOB EXAMPLES *.c )

FOREACH( example_file ${EXAMPLES} )
    # Cut off the file extension and directory path
    GET_FILENAME_COMPONENT( example_exec ${example_file} NAME_WE )
	
    ADD_EXECUTABLE( ${example_exec} ${example_file} )
	
	TARGET_LINK_LIBRARIES( ${example_exec} PRIVATE trustregion_c_api )

ENDFOREACH( example_file ${EXAMPLES} )
