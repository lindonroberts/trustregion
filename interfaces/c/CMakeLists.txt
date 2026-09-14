add_library(trustregion_c_api
    c_trustregion.F90
)

set(HEADERFILE trustregion.h)

set_target_properties(trustregion_c_api PROPERTIES
    PUBLIC_HEADER ${HEADERFILE}
)

target_link_libraries(trustregion_c_api
    PRIVATE
        trustregion
)

target_include_directories(trustregion_c_api
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
)

install(TARGETS trustregion_c_api)
install(FILES ${HEADERFILE} DESTINATION include)