project(djvulibre)

cmake_minimum_required(VERSION 2.6)

# ******************************************************************************
# add some predefined variables
set(LIB_SUFFIX "" CACHE STRING "Define suffix of directory name (32/64)" )  

set(BIN_INSTALL_DIR     bin)
set(INCLUDE_INSTALL_DIR include)
set(LIB_INSTALL_DIR     lib${LIB_SUFFIX})
set(MAN_INSTALL_DIR     share/man)
set(DATA_INSTALL_DIR    share)
set(PLUGIN_INSTALL_DIR  plugins)

option(BUILD_TOOLS "build the different tools" ON)
option(BUILD_XMLTOOLS "build the xml tools" OFF)

set(INSTALL_TARGETS_DEFAULT_ARGS RUNTIME DESTINATION ${BIN_INSTALL_DIR}
                                 LIBRARY DESTINATION ${LIB_INSTALL_DIR}
                                 ARCHIVE DESTINATION ${LIB_INSTALL_DIR})

if(MSVC)
    add_definitions(-D_CRT_SECURE_NO_WARNINGS)
endif(MSVC)

if(MINGW)
    add_definitions(-Wno-deprecated)
endif(MINGW)

if(MINGW)
    add_definitions(-DHAVE_STDINT_H=1)
endif(MINGW)
# ******************************************************************************
# cmake handles the dependencies itself
add_subdirectory(libdjvu)
if(BUILD_TOOLS)
    add_subdirectory(tools)
endif()
if(BUILD_XMLTOOLS)
    add_subdirectory(xmltools)
endif()
