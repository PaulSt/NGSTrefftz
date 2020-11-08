project(trefftzngs)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
    HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
    )

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

install (FILES
        python/__init__.py 
        python/DGeq.py 
        python/prodmesh.py 
        DESTINATION trefftzngs
        COMPONENT trefftzngs
        )

add_ngsolve_python_module(trefftzngs trefftzngspy.cpp
    trefftz/helpers.cpp
    trefftz/trefftzwavefe.cpp
    trefftz/qtrefftzwavefe.cpp
    trefftz/monomialfe.cpp
    trefftz/diffopmapped.hpp
    trefftz/trefftzfespace.cpp
    trefftz/monomialfespace.cpp
    trefftz/scalarmappedfe.cpp
    trefftz/specialcoefficientfunction.cpp
    tents/tents.cpp
    trefftz/evolvetent.cpp
    trefftz/meshtentslab.cpp
    )

#add_compile_options(PRIVATE -Wall -Wextra -pedantic)
#add_compile_options(-g -O0)

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS trefftzngs DESTINATION trefftzngs)
