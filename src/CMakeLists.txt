project(ngstrefftz)

cmake_minimum_required(VERSION 3.1)

include_directories(../dep)

find_package(NGSolve CONFIG REQUIRED
    HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
    )

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

install (FILES
        ../py/__init__.py 
        ../py/DGeq.py 
        ../py/prodmesh.py 
        DESTINATION ngstrefftz
        COMPONENT ngstrefftz
        )

add_ngsolve_python_module(ngstrefftz 
    python_trefftz.cpp
    diffopmapped.hpp
    scalarmappedfe.cpp
    #trefftzheatfe.cpp
    trefftzfespace.cpp
    monomialfespace.cpp
    specialcoefficientfunction.cpp
    ../dep/tents.cpp
    twavetents.cpp
    meshtentslab.cpp
    )

#add_compile_options(PRIVATE -Wall -Wextra -pedantic)
#add_compile_options(-g -O0)

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS ngstrefftz DESTINATION ngstrefftz)
