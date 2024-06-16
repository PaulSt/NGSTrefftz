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

### add ngstents dependency if needed
set(INSTALL_TENTS ON CACHE BOOL "Install ngstents")
if (INSTALL_TENTS)# AND NOT TARGET _pytents)
    include_directories(${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src)
    set(addon_name _pytents)
    project(_pytents)
    include(${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake)
    add_ngsolve_addon( _pytents 
      ${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src/python_tents.cpp
      ${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src/tents.cpp
      )
    target_compile_definitions(_pytents PRIVATE NGSTENT_EXPORTS)
    set_target_properties(_pytents PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents")
    file(WRITE ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents/__init__.py "from ._pytents import TentSlab, Tent")
    set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/ngstents)
    install(TARGETS _pytents DESTINATION ${install_dir})
    install(FILES ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstents/__init__.py DESTINATION ${install_dir})
endif(INSTALL_TENTS)

### ngstrefftz
set(addon_name _trefftz)
project(_trefftz)
include(${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake)
add_ngsolve_addon( _trefftz 
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
    #src/airy.cpp #for testing, requires boost
  )
target_link_libraries(_trefftz PRIVATE _pytents)
target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_EXPORTS)

### link LAPACK or MKL if possible
set(NGSTREFFTZ_USE_MKL ${NGSOLVE_USE_MKL} CACHE BOOL "Try to use mkl")
set(NGSTREFFTZ_USE_LAPACK ${NGSOLVE_USE_LAPACK} CACHE BOOL "Try to use lapack from ngsolve")
if (NGSTREFFTZ_USE_MKL)
    #list(GET NGSOLVE_INCLUDE_DIRS 0 NGSOLVE_MKL_INCLUDE_DIRS)
    file(DOWNLOAD https://raw.githubusercontent.com/NGSolve/ngsolve/master/cmake/cmake_modules/FindMKL.cmake ${CMAKE_BINARY_DIR}/cmake_modules/FindMKL.cmake)

    set(MKL_MULTI_THREADED ON CACHE BOOL "Use threaded MKL libs")
    if(NOT MKL_MULTI_THREADED)
      set(MKL_THREADING sequential)
    endif()
    set(MKL_STATIC OFF CACHE BOOL "Link static MKL")
    set(MKL_SDL ON CACHE BOOL "Link single dynamic MKL lib")
    set(MKL_ARCH "intel64")
    if(NOT MKL_LINK)
      if(MKL_STATIC)
        set(MKL_LINK static)
      elseif(MKL_SDL AND NOT USE_MPI)
        set(MKL_LINK sdl)
      else()
        set(MKL_LINK dynamic)
      endif()
    endif()
    execute_process(COMMAND ${Python3_EXECUTABLE} -c "import sys; print(sys.prefix)"
        OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE rootdir)
    if(WIN32)
        set(MKL_ROOT ${rootdir})
        set(MKL_LIBRARY ${root_dir}/Library/lib/mkl_rt.lib)
        set(MKL_INCLUDE_DIR ${root_dir}/Library/include)
    elseif(UNIX AND NOT APPLE)
        set(MKL_ROOT ${rootdir})
        set(MKL_LIBRARY ${root_dir}/lib/libmkl_rt.so.2)
        set(MKL_INCLUDE_DIR ${root_dir}/include)
    endif()
    find_package(MKL) # HINTS ${MKL_ROOT} ${CMAKE_INSTALL_PREFIX}/../../cmake/mkl ${Python3_ROOT_DIR}/lib/cmake/mkl)
    #message("MKL minimal: ${MKL_MINIMAL_LIBRARIES}")
    #message("MKL libraries: ${MKL_LIBRARIES}")
    #message("MKL link flags: ${MKL_LINK_FLAGS}")

    #if(MKL_FOUND)
    if(MKL_MINIMAL_LIBRARIES)
        set(LAPACK_LIBRARIES "${MKL_LIBRARIES}")
        message(STATUS "Setting LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
        target_include_directories(_trefftz PRIVATE "${MKL_INCLUDE_DIRS}")
    endif()
    if(MKL_LIBRARY)
        message(STATUS "MKL_LIBRARY=${MKL_LIBRARY}")
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
        message(STATUS "MKL_RT_NAME=${MKL_RT_NAME}")
        find_library(MKLRT ${MKL_RT_NAME} HINTS ${MKL_ROOT}/lib ${MKL_LIBRARY})
        message(STATUS "MKLRT=${MKLRT}")
        if(MKLRT)
            set(LAPACK_LIBRARIES ${MKLRT})
            message(STATUS "MKLRT setting LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
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
    message(STATUS "Linking with LAPACK_LIBRARIES=${LAPACK_LIBRARIES}")
    target_link_libraries(_trefftz PRIVATE ${LAPACK_LIBRARIES})
    target_compile_definitions(_trefftz PRIVATE NGSTREFFTZ_USE_LAPACK)
else(LAPACK_LIBRARIES)
    message("${Red}LAPACK not found, restricted usage${ColourReset}")
endif(LAPACK_LIBRARIES)

### doc 
# first we can indicate the documentation build as an option and set it to ON by default
option(BUILD_DOC "Build documentation" ON)

# check if Doxygen is installed
find_package(Doxygen)
if (DOXYGEN_FOUND)
    # set input and output files
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/../docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)

    # request to configure the file
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    message("Doxygen build started")

    # note the option ALL which allows to build the docs together with the application
    add_custom_target( doc_doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
else (DOXYGEN_FOUND)
  message("Doxygen need to be installed to generate the doxygen documentation")
endif (DOXYGEN_FOUND)

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
    ngstents: .......... ${INSTALL_TENTS}

With 'make install' the package will be installed to: ${CMAKE_INSTALL_PREFIX}
Make sure to add it to the python path: 
-->   export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${ADDON_INSTALL_DIR_PYTHON}:$PYTHONPATH
------------------------------------------------------------------------
\n${ColourReset}")

### tests 
include(CTest)
file(COPY 
    ${CMAKE_SOURCE_DIR}/../test/dg.py
    ${CMAKE_SOURCE_DIR}/../test/embt.py
    ${CMAKE_SOURCE_DIR}/../test/trefftz.py
    ${CMAKE_SOURCE_DIR}/../test/tents.py
    ${CMAKE_SOURCE_DIR}/../test/constraint_trefftz.py
    DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing)
add_test(NAME embtrefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/embt.py)
add_test(NAME trefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/trefftz.py)
add_test(NAME constraint_trefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/constraint_trefftz.py)
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
