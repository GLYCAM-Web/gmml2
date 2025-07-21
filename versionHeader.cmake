file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/version.txt" VERSION_VAR)
configure_file(${CMAKE_CURRENT_LIST_DIR}/include/version.h.in ${CMAKE_CURRENT_LIST_DIR}/include/version.h @ONLY)
