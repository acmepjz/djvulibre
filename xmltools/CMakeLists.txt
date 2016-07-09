include_directories( ../libdjvu 
                     ..)

# ******************************************************************************
# build xml tools
add_executable(djvutoxml djvutoxml.cpp)
target_link_libraries(djvutoxml djvulibre)

add_executable(djvuxmlparser djvuxmlparser.cpp)
target_link_libraries(djvuxmlparser djvulibre)

install(TARGETS djvutoxml ${INSTALL_TARGETS_DEFAULT_ARGS})
install(TARGETS djvuxmlparser ${INSTALL_TARGETS_DEFAULT_ARGS})

# ******************************************************************************
# configure manpage
configure_file(djvuxml.1.cmake ${CMAKE_CURRENT_BINARY_DIR}/djvuxml.1)
# install man pages
file(GLOB manpages *.1)
set(manpages ${manpages} ${CMAKE_CURRENT_BINARY_DIR}/djvuxml.1)

install(FILES ${manpages} DESTINATION ${MAN_INSTALL_DIR}/man1)

# ******************************************************************************
# install dtd's
file(GLOB docdefs ../share/pubtest/*.dtd)
install(FILES ${docdefs} DESTINATION ${DATA_INSTALL_DIR}/djvu/pubtest)
