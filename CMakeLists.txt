project(djvulibre)

cmake_minimum_required(VERSION 2.6)

Set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

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
option(BUILD_MINIDJVU "build the minidjvu tools" OFF)
option(BUILD_TSV2DSED "build the tsv2dsed tools" OFF)

find_package(JPEG)
find_package(TIFF)
find_package(ZLIB)
find_package(FreeImage)
find_package(Leptonica)

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

if(ZLIB_FOUND)
	option(HAVE_ZLIB "Define if you have zlib." ON)
else()
	option(HAVE_ZLIB "Define if you have zlib." OFF)
endif()

if(FreeImage_FOUND)
	option(HAVE_FREEIMAGE "Define if you have FreeImage." ON)
else()
	option(HAVE_FREEIMAGE "Define if you have FreeImage." OFF)
endif()

if(Leptonica_FOUND)
	option(HAVE_LEPT "Define if you want to use Leptonica." ON)
else()
	option(HAVE_LEPT "Define if you want to use Leptonica." OFF)
endif()

set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads)
if(CMAKE_USE_PTHREADS_INIT AND NOT WIN32)
	option(HAVE_PTHREAD "Pthread support" ON)
else()
	option(HAVE_PTHREAD "Pthread support" OFF)
endif()

file(READ "configure.ac" CONFIGURE_AC)
string(REGEX MATCH "m4_define\\(\\[v_maj\\],[ ]*\\[[0-9]+\\]\\)" V_MAJ ${CONFIGURE_AC})
string(REGEX REPLACE "^.*\\[([0-9]+)\\].*$" "\\1" V_MAJ ${V_MAJ})
string(REGEX MATCH "m4_define\\(\\[v_min\\],[ ]*\\[[0-9]+\\]\\)" V_MIN ${CONFIGURE_AC})
string(REGEX REPLACE "^.*\\[([0-9]+)\\].*$" "\\1" V_MIN ${V_MIN})
string(REGEX MATCH "m4_define\\(\\[v_mic\\],[ ]*\\[[0-9]+\\]\\)" V_MIC ${CONFIGURE_AC})
string(REGEX REPLACE "^.*\\[([0-9]+)\\].*$" "\\1" V_MIC ${V_MIC})
set(DJVULIBRE_VERSION "${V_MAJ}.${V_MIN}.${V_MIC}")
message(STATUS "The version of djvulibre read from configure.ac is: ${DJVULIBRE_VERSION}")

find_package(git)
if(GIT_FOUND)
	exec_program(${GIT_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR} ARGS "describe"
		OUTPUT_VARIABLE VERSION RETURN_VALUE GIT_RETURN_VALUE)
	if(GIT_RETURN_VALUE STREQUAL "0")
		set(DJVULIBRE_VERSION ${VERSION})
		message(STATUS "The version of djvulibre read from git is: ${DJVULIBRE_VERSION}")
	endif()
endif()
add_definitions(-DDJVULIBRE_VERSION="${DJVULIBRE_VERSION}")

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

if(ZLIB_FOUND AND HAVE_ZLIB)
    add_definitions(-DHAVE_ZLIB=1)
	set(ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES} ${ZLIB_INCLUDE_DIRS})
    set(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${ZLIB_LIBRARIES})
endif()

if(FreeImage_FOUND AND HAVE_FREEIMAGE)
    add_definitions(-DHAVE_FREEIMAGE=1)
	set(ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES} ${FreeImage_INCLUDE_DIRS})
    set(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${FreeImage_LIBRARIES})
endif()

if(Leptonica_FOUND AND HAVE_LEPT)
    add_definitions(-DHAVE_LEPT=1)
	set(ADDITIONAL_INCLUDES ${ADDITIONAL_INCLUDES} ${Leptonica_INCLUDE_DIRS})
    set(ADDITIONAL_LIBRARIES ${ADDITIONAL_LIBRARIES} ${Leptonica_LIBRARIES})
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
if(BUILD_MINIDJVU)
	add_subdirectory(minidjvu)
endif()
if(BUILD_TSV2DSED)
	add_subdirectory(tsv2dsed)
endif()
