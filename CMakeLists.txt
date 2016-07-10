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

include(CheckIncludeFiles)
check_include_files(sys/mman.h HAVE_SYS_MMAN_H)
include(CheckFunctionExists)
check_function_exists(mmap HAVE_MMAP)
if(HAVE_MMAP AND HAVE_SYS_MMAN_H)
	add_definitions(-DHAS_MEMMAP=1)
endif()

option(BUILD_TOOLS "build the different tools" ON)
option(BUILD_XMLTOOLS "build the xml tools" OFF)

find_package(JPEG)
find_package(TIFF)

if(JPEG_FOUND)
	option(HAVE_JPEG "Define if you have the IJG JPEG library." ON)
else()
	option(HAVE_JPEG "Define if you have the IJG JPEG library." OFF)
endif()

if(TIFF_FOUND)
	option(HAVE_TIFF "Define if you have libtiff." ON)
else()
	option(HAVE_TIFF "Define if you have libtiff." OFF)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT AND NOT WIN32)
	option(HAVE_PTHREAD "Pthread support" ON)
else()
	option(HAVE_PTHREAD "Pthread support" OFF)
endif()

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

add_definitions(-DHAVE_NAMESPACES=1)

set(ADDITIONAL_INCLUDES )
set(ADDITIONAL_LIBRARIES )

if(CMAKE_USE_PTHREADS_INIT AND HAVE_PTHREAD)
    add_definitions(-DHAVE_PTHREAD=1)
    set(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()

if(JPEG_FOUND AND HAVE_JPEG)
    add_definitions(-DHAVE_JPEG=1 -DNEED_JPEG_DECODER=1)
	set(ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES} ${JPEG_INCLUDE_DIR})
    set(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${JPEG_LIBRARIES})
endif()

if(TIFF_FOUND AND HAVE_TIFF)
    add_definitions(-DHAVE_TIFF=1)
	set(ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES} ${TIFF_INCLUDE_DIR})
    set(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${TIFF_LIBRARIES})
endif()

# ******************************************************************************
# cmake handles the dependencies itself
add_subdirectory(libdjvu)
if(BUILD_TOOLS)
    add_subdirectory(tools)
endif()
if(BUILD_XMLTOOLS)
    add_subdirectory(xmltools)
endif()
