cmake_minimum_required(VERSION 3.1)
########################################################
project(ngstrefftz)
########################################################
if(NOT WIN32)
  string(ASCII 27 Esc)
  set(ColourReset "${Esc}[m")
  set(Red         "${Esc}[31m")
  set(BoldCyan    "${Esc}[1;36m")
  set(BoldWhite   "${Esc}[1;37m")
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

function(add_ngsolve_python_module2 target)
    if(NOT CMAKE_BUILD_TYPE)
      message(STATUS "Setting build type to NGSolve build type: ${NGSOLVE_CMAKE_BUILD_TYPE}")
      set(CMAKE_BUILD_TYPE ${NGSOLVE_CMAKE_BUILD_TYPE} PARENT_SCOPE)
    endif()
    add_library(${target} SHARED ${ARGN})

    if(${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
        find_package(Python3 REQUIRED COMPONENTS Development.Module)
        if(NOT EMSCRIPTEN)
            find_package(Python3 COMPONENTS Interpreter Development.Embed)
        endif()
    else()
      find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
    endif()
    target_include_directories(${target} PRIVATE ${Python3_INCLUDE_DIRS})

    if(NETGEN_BUILD_FOR_CONDA AND NOT WIN32)
        if(APPLE)
            target_link_options(${target} PUBLIC -undefined dynamic_lookup)
        endif(APPLE)
    else(NETGEN_BUILD_FOR_CONDA AND NOT WIN32)
        target_link_libraries(${target} PUBLIC ${PYTHON_LIBRARY})
    endif(NETGEN_BUILD_FOR_CONDA AND NOT WIN32)

    set_target_properties(${target} PROPERTIES PREFIX "" CXX_STANDARD 17)
    target_link_libraries(${target} PUBLIC ngsolve)

    if(WIN32)
      set_target_properties( ${target} PROPERTIES SUFFIX ".pyd" )
    else(WIN32)
      set_target_properties(${target} PROPERTIES SUFFIX ".so")
    endif(WIN32)

    set_target_properties(${target} PROPERTIES INSTALL_RPATH "${NGSOLVE_LIBRARY_DIR}")
endfunction()

### check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(PY_INSTALL_DIR ${NGSOLVE_INSTALL_DIR_PYTHON} CACHE STRING "subdir of install dir for python modules")
set(NGSTREFFTZ_INSTALL_DIR_PYTHON ${PY_INSTALL_DIR}/ngstrefftz)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} CACHE PATH "build dir")

### add ngstents dependency if needed
if (NOT TARGET _pytents)
    include_directories(${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src)
    set(NGSTENTS_INSTALL_DIR_PYTHON ${PY_INSTALL_DIR}/ngstents)
    add_ngsolve_python_module2(_pytents
      ${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src/python_tents.cpp
      ${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src/tents.cpp
      )
    if(WIN32)
        target_link_libraries(_pytents PUBLIC ${PYTHON_LIBRARY})
    endif(WIN32)
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
    specialintegrator.cpp
    twavetents.cpp
    embtrefftz.cpp
    monomialfespace.cpp 
    mesh1dtents.cpp
    #airy.cpp #for testing, requires boost
    )
if(WIN32)
    target_link_libraries(_trefftz PUBLIC ${PYTHON_LIBRARY})
endif(WIN32)
#target_compile_options(_trefftz PRIVATE "$<$<CONFIG:DEBUG>:-W;-Wall;-Wextra;-Wpedantic;-Wno-unused-parameter;-g;>")
#set_target_properties(_trefftz PROPERTIES COMPILE_DEFINITIONS "NGSTREFFTZ_EXPORTS")
set(NGSTREFFTZ_USE_LAPACK ${NGSOLVE_USE_LAPACK} CACHE BOOL "Try to use lapack from ngsolve")
if(NGSTREFFTZ_USE_LAPACK AND NOT WIN32)
    target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_USE_LAPACK)
endif()
if(USE_MKL AND NGSOLVE_USE_MKL)
    find_package(MKL REQUIRED)
    if(APPLE OR WIN32)
      set(MKL_RT_NAME mkl_rt.2)
    else()
      set(MKL_RT_NAME libmkl_rt.so.2)
    endif()
    find_library(MKLRT ${MKL_RT_NAME} HINTS ${MKL_ROOT}/lib)
    message("Linking with MKL_RT=${MKLRT}")
    target_link_libraries(_trefftz PUBLIC ${MKLRT})
endif()
target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_EXPORTS)
target_link_libraries(_trefftz PRIVATE _pytents)
set_target_properties(_trefftz PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz"
)
file(COPY ${CMAKE_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz)

### install
message("${BoldCyan}
------------------------------------------------------------------------
  Build type: ${CMAKE_BUILD_TYPE}
  Compiler: ${CMAKE_CXX_COMPILER}
  Flags: ${flags}

  Enabled functionality from ngsolve:
    LAPACK: ............ ${NGSTREFFTZ_USE_LAPACK} 
    MKL: ............... ${NGSOLVE_USE_MKL}

With 'make install' the package will be installed to: ${CMAKE_INSTALL_PREFIX}
Make sure to add it to the python path: 
-->   export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}:$PYTHONPATH
------------------------------------------------------------------------
\n${ColourReset}")

install(TARGETS _trefftz DESTINATION ${CMAKE_INSTALL_PREFIX}/${NGSTREFFTZ_INSTALL_DIR_PYTHON}/)
install(FILES
        __init__.py 
        DESTINATION ${CMAKE_INSTALL_PREFIX}/${NGSTREFFTZ_INSTALL_DIR_PYTHON}
        #COMPONENT ngstrefftz
       )

### tests 
include(CTest)
file(COPY 
    ${CMAKE_SOURCE_DIR}/../test/dg.py
    ${CMAKE_SOURCE_DIR}/../test/embt.py
    ${CMAKE_SOURCE_DIR}/../test/trefftz.py
    ${CMAKE_SOURCE_DIR}/../test/tents.py
    DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing)
add_test(NAME embtrefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/embt.py)
add_test(NAME trefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/trefftz.py)
add_test(NAME tents COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/tents.py)
#WORKING_DIRECTORY ${ CMAKE_CURRENT_SOURCE_DIR }
set_tests_properties(embtrefftz trefftz tents
    PROPERTIES ENVIRONMENT "PYTHONPATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON}:$ENV{PYTHONPATH}")

# make uninstall
add_custom_target("uninstall" COMMENT "Uninstall installed files")
add_custom_command(
    TARGET "uninstall"
    POST_BUILD
    COMMENT "Uninstall files with install_manifest.txt"
    COMMAND xargs rm -vf < install_manifest.txt || echo Nothing in
            install_manifest.txt to be uninstalled!
)
