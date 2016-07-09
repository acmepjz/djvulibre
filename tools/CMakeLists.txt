# libdjvu tools

include_directories( ../libdjvu 
                     ..)

# ******************************************************************************
# add a convenience library even though cmake doesn't really support this
set(JB2_SRC jb2tune.cpp 
            jb2cmp/classify.cpp 
            jb2cmp/cuts.cpp 
            jb2cmp/frames.cpp 
            jb2cmp/patterns.cpp)

add_library(jb2 STATIC ${JB2_SRC})

# ******************************************************************************
# add a macro for all tools
macro(add_djvutool name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} djvulibre)
    install(TARGETS ${name} ${INSTALL_TARGETS_DEFAULT_ARGS})
endmacro(add_djvutool name)

add_djvutool(bzz bzz.cpp)
add_djvutool(c44 c44.cpp)
add_djvutool(cjb2 cjb2.cpp)
add_djvutool(cpaldjvu cpaldjvu.cpp)
add_djvutool(csepdjvu csepdjvu.cpp)
add_djvutool(ddjvu ddjvu.cpp tiff2pdf.c)
add_djvutool(djvm djvm.cpp)
add_djvutool(djvmcvt djvmcvt.cpp)
add_djvutool(djvudump djvudump.cpp)
add_djvutool(djvups djvups.cpp)
add_djvutool(djvuextract djvuextract.cpp)
add_djvutool(djvumake djvumake.cpp)
add_djvutool(djvused djvused.cpp)
add_djvutool(djvutxt djvutxt.cpp)
add_djvutool(djvuserve djvuserve.cpp)

target_link_libraries(cjb2 jb2 djvulibre)
target_link_libraries(cpaldjvu jb2 djvulibre)
target_link_libraries(csepdjvu jb2 djvulibre)

# ******************************************************************************
# install man pages
file(GLOB manpages *.1)
install(FILES ${manpages} DESTINATION ${MAN_INSTALL_DIR}/man1)

# ******************************************************************************
# install shell scripts
if(NOT WIN32)
install(FILES djvudigital any2djvu DESTINATION ${BIN_INSTALL_DIR})
endif(NOT WIN32)
