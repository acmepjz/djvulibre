file(READ "configure.ac" CONFIGURE_AC)
string(REGEX MATCH "AC_INIT\\(\\[minidjvu\\],\\[[0-9.]+\\]" MDJVU_VERSION ${CONFIGURE_AC})
string(REGEX REPLACE "^AC_INIT\\(\\[minidjvu\\],\\[([0-9.]+)\\]$" "\\1" MDJVU_VERSION ${MDJVU_VERSION})
message(STATUS "The version of minidjvu read from configure.ac is: ${MDJVU_VERSION}")

add_definitions(-DMDJVU_VERSION="${MDJVU_VERSION}")

include_directories( . 
                     ./include
					 ${ADDITIONAL_INCLUDES} )

File (GLOB SOURCES ./src/*/*.c ./src/*/*.cpp)

add_library(libminidjvu STATIC ${SOURCES})

# ******************************************************************************
# build minidjvu
add_executable(minidjvu tools/minidjvu.c)
target_link_libraries(minidjvu libminidjvu ${ADDITIONAL_LIBRARIES})

install(TARGETS minidjvu ${INSTALL_TARGETS_DEFAULT_ARGS})

# ******************************************************************************
# install man pages
set(manpages doc/minidjvu.1)

install(FILES ${manpages} DESTINATION ${MAN_INSTALL_DIR}/man1)
