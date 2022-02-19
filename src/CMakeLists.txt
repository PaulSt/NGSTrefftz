########################################################
project(ngstrefftz)
########################################################

option(USE_LAPACK "build ngstrefftz with LAPACK" ON)

### Requirements
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()
cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. $ENV{NETGENDIR}/../Resources/CMake /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
    )

### check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(NGSTREFFTZ_INSTALL_DIR_PYTHON ${CMAKE_INSTALL_PREFIX}/${NGSOLVE_INSTALL_DIR_PYTHON}/ngstrefftz)

### add ngstents dependency if needed
if (NOT TARGET _pytents)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../external_dependencies/ngstents/src)
    set(NGSTENTS_INSTALL_DIR_PYTHON ${CMAKE_INSTALL_PREFIX}/${NGSOLVE_INSTALL_DIR_PYTHON}/ngstents)
    add_ngsolve_python_module(_pytents
      ${CMAKE_CURRENT_SOURCE_DIR}/../external_dependencies/ngstents/src/python_tents.cpp
      ${CMAKE_CURRENT_SOURCE_DIR}/../external_dependencies/ngstents/src/tents.cpp
      )
    set_target_properties(_pytents PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ngstents")
    install(TARGETS _pytents DESTINATION ${NGSTENTS_INSTALL_DIR_PYTHON})
endif()

### ngstrefftz
add_ngsolve_python_module(_trefftz 
    python_trefftz.cpp
    diffopmapped.hpp
    scalarmappedfe.cpp
    planewavefe.cpp
    trefftzfespace.cpp
    specialcoefficientfunction.cpp
    twavetents.cpp
    embtrefftz.cpp
    monomialfespace.cpp 
    mesh1dtents.cpp
    #airy.cpp #for testing, requires boost
    )
target_link_libraries(_trefftz PRIVATE _pytents)
set_target_properties(_trefftz PROPERTIES
    #ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/ngstrefftz"
    #RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
)



#add_compile_options(PRIVATE -Wall -Wextra -pedantic)
#add_compile_options(-g -O0)

### link lapack if found
if (USE_LAPACK) # currently needed for SVD in embTrefftz, TODO: switch to ngsolve interface..
    find_package(BLAS)
    find_package(LAPACK)
    if(LAPACK_FOUND AND BLAS_FOUND)
       set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
    endif()
    target_link_libraries(_trefftz PRIVATE ${lapackblas_libraries})
endif (USE_LAPACK)

### install
message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS _trefftz DESTINATION ${NGSTREFFTZ_INSTALL_DIR_PYTHON})
install (FILES
        __init__.py 
        DESTINATION ${NGSTREFFTZ_INSTALL_DIR_PYTHON}
        #COMPONENT ngstrefftz
        )

### tests
include(CTest)
file(COPY 
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/embt.py
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/trefftz.py
    ${CMAKE_CURRENT_SOURCE_DIR}/../test/tents.py
     DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/Testing)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/ngstrefftz)
add_test(NAME embtrefftz COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/Testing/embt.py -v)
add_test(NAME trefftz COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/Testing/trefftz.py -v)
add_test(NAME tents COMMAND python3 ${CMAKE_CURRENT_BINARY_DIR}/Testing/tents.py -v)
