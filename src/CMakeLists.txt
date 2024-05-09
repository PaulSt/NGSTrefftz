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
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find NGSolve and Netgen using python
if(CMAKE_VERSION VERSION_LESS "3.18")
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
else()
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
endif()
set(Netgen_DIR "" CACHE PATH "Path to directory containing NetgenConfig.cmake")
set(NGSolve_DIR "" CACHE PATH "Path to directory containing NGSolveConfig.cmake")
execute_process(COMMAND ${Python3_EXECUTABLE} -m netgen.config OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE Netgen_DIR)
execute_process(COMMAND ${Python3_EXECUTABLE} -m ngsolve.config OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE NGSolve_DIR)
find_package(NGSolve CONFIG REQUIRED)

### check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir

#set(INSTDIR_FROM_NGSOLVE ${CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT} CACHE BOOL "Try to get install dir from NGSolve")
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    execute_process(COMMAND ${Python3_EXECUTABLE} -m site --user-site OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE install_dir RESULT_VARIABLE ret)
    if (NOT ret EQUAL 0)
        # user site directory is disabled (are we in a virtual environment?)
        set(install_dir ${Python3_SITEARCH})
    endif()
    message("Install dir from NGSolve: ${install_dir}")
    set(CMAKE_INSTALL_PREFIX ${install_dir} CACHE PATH "Install dir" FORCE)
else(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    #set(PY_INSTALL_DIR ${NGSOLVE_INSTALL_DIR_PYTHON} CACHE STRING "subdir of install dir for python modules")
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import os.path, sysconfig;print(os.path.relpath(sysconfig.get_path('platlib'), sysconfig.get_path('data')))" OUTPUT_VARIABLE PY_INSTALL_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
    message("python packages install dir: ${PY_INSTALL_DIR}")
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} CACHE PATH "build dir")

### install message
set(NGSTREFFTZ_USE_MKL ${NGSOLVE_USE_MKL} CACHE BOOL "Try to use mkl")
set(NGSTREFFTZ_USE_LAPACK ${NGSOLVE_USE_LAPACK} CACHE BOOL "Try to use lapack from ngsolve")
set(INSTALL_TENTS ON CACHE BOOL "Install ngstents")
message("${BoldCyan}
------------------------------------------------------------------------
  Build type: ${CMAKE_BUILD_TYPE}
  Compiler: ${CMAKE_CXX_COMPILER}
  Flags: ${flags}

  Enabled functionality:
    LAPACK: ............ ${NGSTREFFTZ_USE_LAPACK} 
    MKL: ............... ${NGSTREFFTZ_USE_MKL}
    ngstents: .......... ${INSTALL_TENTS}

With 'make install' the package will be installed to: ${CMAKE_INSTALL_PREFIX}
Make sure to add it to the python path: 
-->   export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}:$PYTHONPATH
------------------------------------------------------------------------
\n${ColourReset}")

### add ngstents dependency if needed
if (INSTALL_TENTS)# AND NOT TARGET _pytents)
    include_directories(${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src)
    add_library( _pytents SHARED
      ${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src/python_tents.cpp
      ${CMAKE_SOURCE_DIR}/../external_dependencies/ngstents/src/tents.cpp
      )
    target_link_libraries(_pytents PUBLIC ngsolve Python3::Module)
    set_target_properties(_pytents PROPERTIES PREFIX "" CXX_STANDARD 17)
    if(WIN32)
      set_target_properties(_pytents PROPERTIES SUFFIX ".pyd" )
    else(WIN32)
      set_target_properties(_pytents PROPERTIES SUFFIX ".so")
    endif(WIN32)
    target_compile_definitions(_pytents PRIVATE NGSTENT_EXPORTS)
    set_target_properties(_pytents PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents")
    file(WRITE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents/__init__.py "from ._pytents import TentSlab, Tent")
    install(TARGETS _pytents DESTINATION ${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}/ngstents)
    install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents/__init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}/ngstents)
endif(INSTALL_TENTS)

### ngstrefftz
add_library( _trefftz SHARED
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
    condensedg.cpp
    #airy.cpp #for testing, requires boost
  )
target_link_libraries(_trefftz PUBLIC ngsolve Python3::Module)
target_link_libraries(_trefftz PRIVATE _pytents)
set_target_properties(_trefftz PROPERTIES PREFIX "" CXX_STANDARD 17)
if(WIN32)
  set_target_properties(_trefftz PROPERTIES SUFFIX ".pyd" )
else(WIN32)
  set_target_properties(_trefftz PROPERTIES SUFFIX ".so")
endif(WIN32)
target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_EXPORTS)
#target_compile_options(_trefftz PRIVATE "$<$<CONFIG:DEBUG>:-W;-Wall;-Wextra;-Wpedantic;-Wno-unused-parameter;-g;>")
#set_target_properties(_trefftz PROPERTIES COMPILE_DEFINITIONS "NGSTREFFTZ_EXPORTS")

### link LAPACK or MKL if possible
if (NGSTREFFTZ_USE_MKL)
    #list(GET NGSOLVE_INCLUDE_DIRS 0 NGSOLVE_MKL_INCLUDE_DIRS)
    file(DOWNLOAD https://raw.githubusercontent.com/NGSolve/ngsolve/master/cmake/cmake_modules/FindMKL.cmake ${CMAKE_BINARY_DIR}/cmake_modules/FindMKL.cmake)
    set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH}" "${CMAKE_BINARY_DIR}/cmake_modules")

    set(MKL_MULTI_THREADED ON CACHE BOOL "Use threaded MKL libs")
    if(NOT MKL_MULTI_THREADED)
      set(MKL_THREADING sequential)
    endif()
    set(MKL_STATIC OFF CACHE BOOL "Link static MKL")
    set(MKL_SDL ON CACHE BOOL "Link single dynamic MKL lib")
    set(MKL_ARCH "intel64")
    message("MKL_ARCH=${MKL_ARCH}")
    if(NOT MKL_LINK)
      if(MKL_STATIC)
        set(MKL_LINK static)
      elseif(MKL_SDL AND NOT USE_MPI)
        set(MKL_LINK sdl)
      else()
        set(MKL_LINK dynamic)
      endif()
    endif()
    find_package(MKL) # HINTS ${MKL_ROOT} ${CMAKE_INSTALL_PREFIX}/../../cmake/mkl ${Python3_ROOT_DIR}/lib/cmake/mkl)
    #message("MKL minimal: ${MKL_MINIMAL_LIBRARIES}")
    #message("MKL libraries: ${MKL_LIBRARIES}")
    #message("MKL link flags: ${MKL_LINK_FLAGS}")

    #if(MKL_FOUND)
    if(MKL_MINIMAL_LIBRARIES)
        set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
        message("Setting LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
        target_include_directories(_trefftz PRIVATE "${MKL_INCLUDE_DIRS}")
    endif()
    if(MKL_LIBRARY)
        message("MKL_LIBRARY=${MKL_LIBRARY}")
        #if(WIN32)
            #set(LAPACK_LIBRARIES "${MKL_MINIMAL_LIBRARIES}")
        #else()
            #set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
        #endif()
    #if(NOT MKL_LIBRARIES)
        #if(APPLE OR WIN32)
          #set(MKL_RT_NAME mkl_rt.2)
        #else()
          #set(MKL_RT_NAME libmkl_rt.so.2)
        #endif()
        get_filename_component(MKL_RT_NAME ${MKL_LIBRARY} NAME)
        message("MKL_RT_NAME=${MKL_RT_NAME}")
        find_library(MKLRT ${MKL_RT_NAME} HINTS ${MKL_ROOT}/lib ${MKL_LIBRARY})
        message("MKLRT=${MKLRT}")
        if(MKLRT)
            set(LAPACK_LIBRARIES ${MKLRT})
            message("MKLRT setting LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
        endif()
    #else()
        #message("${Red}MKL not found, switching off NGSTREFFTZ_USE_LAPACK!${ColourReset}")
        #set(NGSTREFFTZ_USE_LAPACK OFF)
    endif()
endif (NGSTREFFTZ_USE_MKL)
if(NGSTREFFTZ_USE_LAPACK AND NOT LAPACK_LIBRARIES)
    find_package(LAPACK)
endif(NGSTREFFTZ_USE_LAPACK AND NOT LAPACK_LIBRARIES)
if(LAPACK_LIBRARIES)
    #target_link_libraries(_trefftz PRIVATE $<BUILD_INTERFACE:ngs_lapack>)
    message("Linking with LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
    target_link_libraries(_trefftz PRIVATE ${LAPACK_LIBRARIES})
    target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_USE_LAPACK)
else(LAPACK_LIBRARIES)
    message("${Red}LAPACK not found, restricted usage${ColourReset}")
endif(LAPACK_LIBRARIES)

### install
set_target_properties(_trefftz PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz")
file(COPY ${CMAKE_SOURCE_DIR}/__init__.py DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz)
install(TARGETS _trefftz DESTINATION ${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}/ngstrefftz)
install(FILES __init__.py DESTINATION ${CMAKE_INSTALL_PREFIX}/${PY_INSTALL_DIR}/ngstrefftz) #COMPONENT ngstrefftz

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
