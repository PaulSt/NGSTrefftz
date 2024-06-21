cmake_minimum_required(VERSION 3.16)
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

file(DOWNLOAD https://raw.githubusercontent.com/NGSolve/ngsolve-addon-template/main/ngsolve_addon.cmake ${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake)
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_BINARY_DIR}/cmake_modules")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} CACHE PATH "build dir")

### ngstrefftz
set(addon_name _trefftz)
project(_trefftz)
include_directories(${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src)
include(${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake)
add_ngsolve_addon( _trefftz 
    ${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src/tents.cpp
    src/python_trefftz.cpp
    src/diffopmapped.hpp
    src/scalarmappedfe.cpp
    src/planewavefe.cpp
    src/trefftzfespace.cpp
    src/specialcoefficientfunction.cpp
    src/specialintegrator.cpp
    src/twavetents.cpp
    src/embtrefftz.cpp
    src/monomialfespace.cpp 
    src/mesh1dtents.cpp
    src/condensedg.cpp
    src/pufe.cpp
    src/pufespace.cpp
    #src/airy.cpp #for testing, requires boost
  )
target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_EXPORTS)
target_compile_definitions(_trefftz PRIVATE NGSTENT_EXPORTS)

### link LAPACK or MKL if possible
set(NGSTREFFTZ_USE_MKL ${NGSOLVE_USE_MKL} CACHE BOOL "Try to use mkl")
set(NGSTREFFTZ_USE_LAPACK ${NGSOLVE_USE_LAPACK} CACHE BOOL "Try to use lapack from ngsolve")
if (NGSTREFFTZ_USE_MKL)
    file(DOWNLOAD 
    https://raw.githubusercontent.com/NGSolve/ngsolve/master/cmake/cmake_modules/FindMKL.cmake 
    "${CMAKE_BINARY_DIR}/cmake_modules/FindMKL.cmake")

    #file(COPY ${CMAKE_SOURCE_DIR}/FindMKL.cmake DESTINATION 
    #"${CMAKE_BINARY_DIR}/cmake_modules/FindMKL.cmake")

    #message(STATUS "MKL_ROOT=${MKL_ROOT}")
    #message(STATUS "NGSOLVED=${NGSolve_DIR}")


    #set(PIP_MKL_ROOT FALSE CACHE BOOL "Set to true if MKL_ROOT is set by pip")
    #message(STATUS "PIP_MKL_ROOT=${PIP_MKL_ROOT}")
    #if(PIP_MKL_ROOT AND LINUX)
        #set(MKL_ROOT "${NGSolve_DIR}/../../.." CACHE PATH "MKL root directory" )
        #set(MKL_LIBRARY "${MKL_ROOT}/lib/libmkl_rt.so.2" CACHE PATH "MKL library" )
        #set(MKL_INCLUDE_DIR "${MKL_ROOT}/include" CACHE PATH "MKL include directory" )
        #message(STATUS "set MKL_LIBRARY=${MKL_LIBRARY}")
    #elseif(PIP_MKL_ROOT AND WIN32)
        #set(MKL_ROOT "${NGSolve_DIR}/../../../.." CACHE PATH "MKL root directory" )
        #set(MKL_LIBRARY "${MKL_ROOT}/Library/lib/mkl_rt.lib" CACHE PATH "MKL library" )
        #set(MKL_INCLUDE_DIR "${MKL_ROOT}/Library/include" CACHE PATH "MKL include directory" )
    #endif()

    set(MKL_ARCH "intel64")
    if(NOT MKL_INTERFACE)
        set(MKL_INTERFACE "lp64")
    endif(NOT MKL_INTERFACE)
    set(MKL_MULTI_THREADED ON CACHE BOOL "Use threaded MKL libs")
    if(NOT MKL_MULTI_THREADED)
      set(MKL_THREADING sequential)
    endif()
    set(MKL_STATIC OFF CACHE BOOL "Link static MKL")
    set(MKL_SDL ON CACHE BOOL "Link single dynamic MKL lib")
    if(NOT MKL_LINK)
      if(MKL_STATIC)
        set(MKL_LINK static)
      elseif(MKL_SDL)
        set(MKL_LINK sdl)
      else()
        set(MKL_LINK dynamic)
      endif()
    endif()

    file(WRITE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/findmklconf.py 
        "from contextlib import suppress;import importlib.metadata; \
        [print( f.locate().parent.resolve().absolute().as_posix()) \
        if f.match('MKLConfig.cmake') else 0 for f in importlib.metadata.files('mkl-devel')]")
    execute_process(COMMAND ${Python3_EXECUTABLE} ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/findmklconf.py
        ERROR_QUIET OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE MKL_CONF_FROM_PY)
    if(MKL_CONF_FROM_PY)
        message(STATUS "MKL_CONF_FROM_PY=${MKL_CONF_FROM_PY}")
        message(STATUS "pylibdir=${python3_library_dir}")
        set(MKL_ROOT "${MKL_CONF_FROM_PY}/../../.." CACHE PATH "MKL root directory" )
        set(MKL_LIBRARY "${MKL_ROOT}/lib/libmkl_rt.so.2" CACHE PATH "MKL library" )
        set(MKL_INCLUDE_DIR "${MKL_ROOT}/include" CACHE PATH "MKL include directory" )
    endif()
    #set(MKL_VERSION_H "${MKL_ROOT}/include/mkl_version.h")

    #set(d lib/cmake/mkl)
    #find_package(MKL CONFIG REQUIRED HINTS ${MKL_CONF} ${MKL_ROOT}/${d} PATHS /opt/intel/oneapi/mkl/latest/${d} ~/.local/${d})

    unset(USE_MPI CACHE)
    find_package(MKL)

    #if(USE_MUMPS)
        ## include scalapack
        #set( LAPACK_LIBRARIES "${MKL_LIBRARIES}")
        #set( ENABLE_SCALAPACK ON)
    #else(USE_MUMPS)
        #set( LAPACK_LIBRARIES "${MKL_MINIMAL_LIBRARIES}")
    #endif(USE_MUMPS)

    if(MKL_FOUND)
        message(STATUS "Found MKL_MINIMAL_LIBRARIES=${MKL_MINIMAL_LIBRARIES}")
        set( LAPACK_LIBRARIES "${MKL_MINIMAL_LIBRARIES}")
    else()
        message("${Red}MKL not found!${ColourReset}")
    endif()
endif (NGSTREFFTZ_USE_MKL)
if(NGSTREFFTZ_USE_LAPACK AND NOT LAPACK_LIBRARIES)
    message(STATUS "Looking for LAPACK")
    find_package(LAPACK)
endif(NGSTREFFTZ_USE_LAPACK AND NOT LAPACK_LIBRARIES)
if(LAPACK_LIBRARIES)
    #target_link_libraries(_trefftz PRIVATE $<BUILD_INTERFACE:ngs_lapack>)
    message(STATUS "Linking with LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
    target_link_libraries(_trefftz PRIVATE ${LAPACK_LIBRARIES})
    target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_USE_LAPACK)
else(LAPACK_LIBRARIES)
    message("${Red}LAPACK not found, restricted usage${ColourReset}")
endif(LAPACK_LIBRARIES)

### install
set_target_properties(_trefftz PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz")
file(COPY ${CMAKE_SOURCE_DIR}/src/__init__.py DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz)
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/ngstrefftz)
install(TARGETS _trefftz DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})
#ngsolve_generate_stub_files(${addon_name})

message("${BoldCyan}
------------------------------------------------------------------------
  Build type: ${CMAKE_BUILD_TYPE}
  Compiler: ${CMAKE_CXX_COMPILER}
  Flags: ${flags}

  Enabled functionality:
    LAPACK: ............ ${NGSTREFFTZ_USE_LAPACK} 
    MKL: ............... ${NGSTREFFTZ_USE_MKL}

With 'make install' the package will be installed to: ${CMAKE_INSTALL_PREFIX}
Make sure to add it to the python path: 
-->   export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${python3_library_dir}:$PYTHONPATH
------------------------------------------------------------------------
\n${ColourReset}")

### tests 
include(CTest)
file(COPY 
    ${CMAKE_SOURCE_DIR}/test/dg.py
    ${CMAKE_SOURCE_DIR}/test/embt.py
    ${CMAKE_SOURCE_DIR}/test/trefftz.py
    ${CMAKE_SOURCE_DIR}/test/tents.py
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
