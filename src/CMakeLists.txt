cmake_minimum_required(VERSION 3.1)
########################################################
project(ngstrefftz)
########################################################
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(BoldCyan    "${Esc}[1;36m")
endif()

### Requirements
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(NGSolve CONFIG REQUIRED
    HINTS $ENV{NETGENDIR}/.. $ENV{NETGENDIR}/../Resources/CMake /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen $ENV{NGSolve_DIR}
    )

### check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
if (NOT DEFINED PY_INSTALL_DIR)
    set(PY_INSTALL_DIR ${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "subdir of install dir for python modules" FORCE)
endif()
set(NGSTREFFTZ_INSTALL_DIR_PYTHON ${PY_INSTALL_DIR}/ngstrefftz)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} CACHE PATH "build dir")


function(add_ngsolve_python_module2 target)
    if(NOT CMAKE_BUILD_TYPE)
      message(STATUS "Setting build type to NGSolve build type: ${NGSOLVE_CMAKE_BUILD_TYPE}")
      set(CMAKE_BUILD_TYPE ${NGSOLVE_CMAKE_BUILD_TYPE} PARENT_SCOPE)
    endif()
    add_library(${target} SHARED ${ARGN})

    #find_package (Python3 COMPONENTS Interpreter Development)
    #target_include_directories(${target} PRIVATE ${Python3_INCLUDE_DIRS})
    find_package(PythonInterp 3 REQUIRED)
    find_package(PythonLibs 3 REQUIRED)
    target_include_directories(${target} PRIVATE ${PYTHON_INCLUDE_DIR})

    if(NETGEN_BUILD_FOR_CONDA)
        message(STATUS "netgen build for conda")
        if(APPLE)
            target_link_options(${target} PUBLIC -undefined dynamic_lookup)
        else(APPLE)
        #target_link_libraries(${target} PUBLIC ${Python3_LIBRARIES})
        target_link_libraries(${target} PUBLIC ${PYTHON_LIBRARY})
        endif(APPLE)
    endif(NETGEN_BUILD_FOR_CONDA)

    set_target_properties(${target} PROPERTIES PREFIX "" CXX_STANDARD 17)

    #if(NETGEN_USE_GUI AND WIN32)
    target_link_libraries(${target} PUBLIC ngsolve)
    #else(NETGEN_USE_GUI AND WIN32)
    #  target_link_libraries(${target} PUBLIC solve)
    #endif(NETGEN_USE_GUI AND WIN32)

    if(WIN32)
      set_target_properties( ${target} PROPERTIES SUFFIX ".pyd" )
    else(WIN32)
      set_target_properties(${target} PROPERTIES SUFFIX ".so")
    endif(WIN32)

    set_target_properties(${target} PROPERTIES INSTALL_RPATH "${NGSOLVE_LIBRARY_DIR}")
endfunction()

### add ngstents dependency if needed
if (NOT TARGET _pytents)
    include_directories(${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src)
    set(NGSTENTS_INSTALL_DIR_PYTHON ${PY_INSTALL_DIR}/ngstents)
    add_ngsolve_python_module2(_pytents
      ${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src/python_tents.cpp
      ${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src/tents.cpp
      )
    target_compile_definitions(_pytents PRIVATE NGSTENT_EXPORTS)
    set_target_properties(_pytents PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents")
    install(TARGETS _pytents DESTINATION ${CMAKE_INSTALL_PREFIX}/${NGSTENTS_INSTALL_DIR_PYTHON})
endif()

### ngstrefftz
add_ngsolve_python_module2(_trefftz 
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
#target_compile_options(_trefftz PRIVATE "$<$<CONFIG:DEBUG>:-W;-Wall;-Wextra;-Wpedantic;-Wno-unused-parameter;-g;>")
target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_EXPORTS)
target_link_libraries(_trefftz PRIVATE _pytents)
set_target_properties(_trefftz PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz"
)
file(COPY ${CMAKE_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz)

### link lapack if found
if (USE_LAPACK) # currently needed for SVD in embTrefftz, TODO: switch to ngsolve interface..
    #find_package(BLAS)
    #find_package(LAPACK)
    #if(LAPACK_FOUND AND BLAS_FOUND)
       #set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
       #target_link_libraries(_trefftz PRIVATE ${lapackblas_libraries})
       #target_compile_definitions(_trefftz PRIVATE ET_LAPACK)
    #endif()
    target_link_libraries(_trefftz INTERFACE "${LAPACK_LIBRARIES}")
    #target_link_libraries(_trefftz INTERFACE ngs_lapack)
endif (USE_LAPACK)
#if(NOT LAPACK_FOUND OR NOT BLAS_FOUND)
    #target_compile_definitions(_trefftz PRIVATE -DET_LAPACK)
#endif()

### install
message("${BoldCyan}With 'make install' the package will be installed to: ${CMAKE_INSTALL_PREFIX}
Make sure to add it to the python path: 
-->   export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}:$PYTHONPATH")
install(TARGETS _trefftz DESTINATION ${CMAKE_INSTALL_PREFIX}/${NGSTREFFTZ_INSTALL_DIR_PYTHON}/)
install(FILES
        __init__.py 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${NGSTREFFTZ_INSTALL_DIR_PYTHON}
        #COMPONENT ngstrefftz
       )

### tests 
include(CTest)
file(COPY 
    ${CMAKE_SOURCE_DIR}/../test/embt.py
    ${CMAKE_SOURCE_DIR}/../test/trefftz.py
    ${CMAKE_SOURCE_DIR}/../test/tents.py
     DESTINATION ${CMAKE_BINARY_DIR}/Testing)
add_test(NAME embtrefftz COMMAND python3 -m doctest ${CMAKE_BINARY_DIR}/Testing/embt.py)
add_test(NAME trefftz COMMAND python3 -m doctest ${CMAKE_BINARY_DIR}/Testing/trefftz.py)
add_test(NAME tents COMMAND python3 -m doctest ${CMAKE_BINARY_DIR}/Testing/tents.py)
#WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
set_tests_properties(embtrefftz trefftz tents
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:$ENV{PYTHONPATH}")
