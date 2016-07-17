# ******************************************************************************
# build tsv2dsed
add_executable(tsv2dsed tsv2dsed.cpp)
target_link_libraries(tsv2dsed )

install(TARGETS tsv2dsed ${INSTALL_TARGETS_DEFAULT_ARGS})
