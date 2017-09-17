# ******************************************************************************
# building djvu library
set(DJVU_SRC
        DjVuGlobal.cpp 
        DjVuGlobalMemory.cpp 
        GOS.cpp 
        GException.cpp 
        GContainer.cpp
        GString.cpp 
        GThreads.cpp 
        GRect.cpp 
        atomic.cpp 
        GSmartPointer.cpp 
        MMX.cpp
        Arrays.cpp 
        GBitmap.cpp 
        GPixmap.cpp 
        GScaler.cpp 
        DjVuPalette.cpp
        ByteStream.cpp 
        IFFByteStream.cpp 
        BSByteStream.cpp 
        BSEncodeByteStream.cpp
        ZPCodec.cpp 
        JB2Image.cpp 
        JB2EncodeCodec.cpp 
        IW44Image.cpp
        IW44EncodeCodec.cpp 
        MMRDecoder.cpp 
        JPEGDecoder.cpp 
        DataPool.cpp 
        GURL.cpp
        DjVuFileCache.cpp 
        DjVuPort.cpp 
        GMapAreas.cpp 
        DjVuAnno.cpp 
        GIFFManager.cpp
        DjVuText.cpp 
        DjVuInfo.cpp 
        DjVuFile.cpp 
        DjVuImage.cpp 
        DjVuDocument.cpp
        DjVmDir.cpp 
        DjVmDoc.cpp 
        DjVmDir0.cpp 
        DjVmNav.cpp 
        DjVuNavDir.cpp
        DjVuErrorList.cpp 
        DjVuDocEditor.cpp 
        DjVuToPS.cpp 
        DjVuDumpHelper.cpp
        DjVuMessageLite.cpp 
        DjVuMessage.cpp 
        GUnicode.cpp 
        UnicodeByteStream.cpp
        XMLParser.cpp 
        XMLTags.cpp 
        miniexp.cpp 
        ddjvuapi.cpp 
        debug.cpp
)

include_directories(${ADDITIONAL_INCLUDES})

add_definitions(-DDJVUAPI_EXPORT -DMINILISPAPI_EXPORT)

add_library(djvulibre SHARED ${DJVU_SRC})

target_link_libraries(djvulibre ${ADDITIONAL_LIBRARIES})

install(TARGETS djvulibre ${INSTALL_TARGETS_DEFAULT_ARGS})

# ******************************************************************************
# installing djvu headers
install(FILES ddjvuapi.h miniexp.h DESTINATION ${INCLUDE_INSTALL_DIR}/libdjvu)

# ******************************************************************************
# install djvu data files
file(GLOB_RECURSE datafiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/../share/*.xml)

foreach(filename ${datafiles})
    get_filename_component(installname ${filename} PATH)
    string(REPLACE "../share/djvu/" "" installname ${installname})
    install(FILES ${filename} DESTINATION ${DATA_INSTALL_DIR}/${installname})
endforeach(filename ${datafiles})
