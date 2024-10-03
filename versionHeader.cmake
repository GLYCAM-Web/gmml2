file(STRINGS "${CMAKE_CURRENT_LIST_DIR}/version.txt" VERSION_VAR)
configure_file(${CMAKE_CURRENT_LIST_DIR}/includes/version.h.in ${CMAKE_CURRENT_LIST_DIR}/includes/version.h @ONLY)
