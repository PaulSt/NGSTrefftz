project(ngstrefftz)
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(EXT_TTEST "Add helpful functions for testing" OFF)
option(USE_LAPACK "build ngstrefftz with LAPACK" ON)


find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. $ENV{NETGENDIR}/../Resources/CMake /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
    )
# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)


if (NOT TARGET _pytents)
include_directories(../external_dependencies/ngstents/src)
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
    trefftzfespace.cpp
    specialcoefficientfunction.cpp
    twavetents.cpp
    embtrefftz.cpp
    )
if(EXT_TTEST)
    add_compile_definitions(EXT_TTEST)
    target_sources(_trefftz PRIVATE
    monomialfespace.cpp 
    mesh1dtents.cpp
    #airy.cpp #for testing, requires boost
    )
endif()
target_link_libraries(_trefftz PRIVATE _pytents)
#add_compile_options(PRIVATE -Wall -Wextra -pedantic)
#add_compile_options(-g -O0)

if (USE_LAPACK) # currently needed for SVD in embTrefftz, TODO: switch to ngsolve interface..
    find_package(BLAS)
    find_package(LAPACK)
    if(LAPACK_FOUND AND BLAS_FOUND)
       set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
    target_link_libraries(_trefftz PRIVATE ${lapackblas_libraries})
endif (USE_LAPACK)


message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS _trefftz DESTINATION ngstrefftz)
install (FILES
        __init__.py 
        DESTINATION ngstrefftz
        #COMPONENT ngstrefftz
        )

include(CTest)
add_test(NAME embtrefftz COMMAND python3 ../.github/tests/embt.py -v)
add_test(NAME trefftz COMMAND python3 ../.github/tests/trefftz.py -v)
#add_test(NAME tents COMMAND python3 ../.github/tests/tents.py -v)
