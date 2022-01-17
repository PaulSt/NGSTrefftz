project(ngstrefftz)

cmake_minimum_required(VERSION 3.1)

include_directories(../external_dependencies/ngstents/src)

find_package(NGSolve CONFIG REQUIRED
    HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
    )

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


if (NOT TARGET _pytents)
add_ngsolve_python_module(_pytents
  ../external_dependencies/ngstents/src/python_tents.cpp
  ../external_dependencies/ngstents/src/tents.cpp
  )
install(TARGETS _pytents DESTINATION ngstents)
endif()

add_ngsolve_python_module(_trefftz 
    python_trefftz.cpp
    diffopmapped.hpp
    scalarmappedfe.cpp
    planewavefe.cpp
    #trefftzheatfe.cpp
    trefftzfespace.cpp
    monomialfespace.cpp
    specialcoefficientfunction.cpp
    twavetents.cpp
    svdtrefftz.cpp
    #airy.cpp
    meshtentslab.cpp
    )
target_link_libraries(_trefftz PRIVATE _pytents)

#add_compile_options(PRIVATE -Wall -Wextra -pedantic)
#add_compile_options(-g -O0)

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS _trefftz DESTINATION ngstrefftz)
install (FILES
        __init__.py 
        DGeq.py 
        prodmesh.py 
        DESTINATION ngstrefftz
        #COMPONENT ngstrefftz
        )

include(CTest)
add_test(NAME svdtest COMMAND python3 ../.github/tests/svdt.py -v)
add_test(NAME trefftz COMMAND python3 ../.github/tests/trefftz.py -v)
add_test(NAME tents COMMAND python3 ../.github/tests/tents.py -v)
