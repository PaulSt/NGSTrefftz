cmake_minimum_required(VERSION 3.18)
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

function(download_if_possible url file)
  get_filename_component(dir "${file}" DIRECTORY)
  file(MAKE_DIRECTORY "${dir}")
  file(DOWNLOAD "${url}" "${file}.tmp" STATUS status SHOW_PROGRESS TIMEOUT 10)
  list(GET status 0 code)
  if(NOT code EQUAL 0)
    if(EXISTS "${file}")
      message(WARNING "offline: using existing ${file}")
    else()
      message(FATAL_ERROR "failed to download ${file}")
    endif()
  else()
    file(RENAME "${file}.tmp" "${file}")
  endif()
endfunction()

set(ngsolve_addon_commit_hash refs/heads/main)
download_if_possible (
    "https://raw.githubusercontent.com/NGSolve/ngsolve-addon-template/${ngsolve_addon_commit_hash}/ngsolve_addon.cmake"
    "${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_BINARY_DIR}/cmake_modules)

# -----------------------------------------------------------------------------
# ngsolve-targets.cmake may reference MKL::MKL before NGSolveConfig.cmake has imported it
if(NOT TARGET MKL::MKL)
    set(MKL_INTERFACE "lp64" CACHE STRING "MKL integer interface" FORCE)
    set(MKL_LINK "sdl" CACHE STRING "MKL linking model" FORCE)
    set(_NGSTREFFTZ_MKL_HINTS)
    if(MKL_ROOT)
        list(APPEND _NGSTREFFTZ_MKL_HINTS "${MKL_ROOT}/lib/cmake/mkl")
    endif()
    if(DEFINED ENV{MKLROOT} AND NOT "$ENV{MKLROOT}" STREQUAL "")
        list(APPEND _NGSTREFFTZ_MKL_HINTS "$ENV{MKLROOT}/lib/cmake/mkl")
    endif()
    list(APPEND _NGSTREFFTZ_MKL_HINTS
        "/opt/intel/oneapi/mkl/latest/lib/cmake/mkl"
        "$ENV{HOME}/.local/lib/cmake/mkl"
    )
    find_package(MKL CONFIG QUIET HINTS ${_NGSTREFFTZ_MKL_HINTS})
endif()

# Older NGSolve wheels used the PyPI mkl-devel package.
if(NOT TARGET MKL::MKL)
    find_package(Python3 QUIET COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        execute_process(
            COMMAND "${Python3_EXECUTABLE}" -c
                "from importlib.metadata import files; fs=files('mkl-devel') or []; print(next((str(p.locate().parent.resolve()) for p in fs if p.name == 'MKLConfig.cmake'), ''))"
            RESULT_VARIABLE _NGSTREFFTZ_MKL_QUERY_RESULT
            OUTPUT_VARIABLE _NGSTREFFTZ_MKL_CONFIG_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET
        )
        if(_NGSTREFFTZ_MKL_QUERY_RESULT EQUAL 0
           AND _NGSTREFFTZ_MKL_CONFIG_DIR
           AND EXISTS "${_NGSTREFFTZ_MKL_CONFIG_DIR}/MKLConfig.cmake")
            find_package(MKL CONFIG QUIET
                HINTS "${_NGSTREFFTZ_MKL_CONFIG_DIR}"
                NO_DEFAULT_PATH
            )
        endif()
    endif()
endif()

if(TARGET MKL::MKL)
    message(STATUS "Preloaded MKL::MKL for the NGSolve CMake package")
endif()

unset(_NGSTREFFTZ_MKL_HINTS)
unset(_NGSTREFFTZ_MKL_QUERY_RESULT)
unset(_NGSTREFFTZ_MKL_CONFIG_DIR)

function(ngstrefftz_relocate_ngsolve_openblas)
    if(NOT TARGET ngs_lapack OR APPLE)
        return()
    endif()

    get_target_property(_ngs_lapack_libs ngs_lapack INTERFACE_LINK_LIBRARIES)
    if(NOT _ngs_lapack_libs)
        return()
    endif()

    set(_ngs_openblas_paths_are_stale FALSE)
    foreach(_lib IN LISTS _ngs_lapack_libs)
        if(IS_ABSOLUTE "${_lib}"
           AND "${_lib}" MATCHES "ngsolve_openblas"
           AND NOT EXISTS "${_lib}")
            set(_ngs_openblas_paths_are_stale TRUE)
        endif()
    endforeach()

    if(NOT _ngs_openblas_paths_are_stale)
        return()
    endif()

    find_package(Python3 REQUIRED COMPONENTS Interpreter)

    execute_process(
        COMMAND "${Python3_EXECUTABLE}" -c
            "import importlib.util, pathlib, sys; s=importlib.util.find_spec('ngsolve_openblas'); p=pathlib.Path(next(iter(s.submodule_search_locations))) if s and s.submodule_search_locations else None; libs=([*p.glob('libopenblas*.so*'), *p.glob('libgfortran*.so*')] if p and sys.platform.startswith('linux') else ([p/'openblas.lib'] if p and sys.platform.startswith('win') else [])); print(';'.join(str(x.resolve()) for x in libs if x.exists()))"
        RESULT_VARIABLE _ngs_openblas_query_result
        OUTPUT_VARIABLE _ngs_openblas_libraries
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_VARIABLE _ngs_openblas_query_error
    )

    if(NOT _ngs_openblas_query_result EQUAL 0 OR NOT _ngs_openblas_libraries)
        message(FATAL_ERROR
            "NGSolve's ngs_lapack target contains stale ngsolve_openblas paths, "
            "but the active Python environment does not provide usable OpenBLAS "
            "libraries. Python error: ${_ngs_openblas_query_error}")
    endif()

    message(STATUS
        "Relocating NGSolve OpenBLAS link interface to: ${_ngs_openblas_libraries}")
    set_property(TARGET ngs_lapack PROPERTY
        INTERFACE_LINK_LIBRARIES "${_ngs_openblas_libraries}")
endfunction()
# -----------------------------------------------------------------------------

### ngstrefftz

set(addon_name ngstrefftz)
project(ngstrefftz)
include_directories(SYSTEM ${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src) # SYSTEM to suppress warnings from headers
include(${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake)
ngstrefftz_relocate_ngsolve_openblas()
add_ngsolve_addon(ngstrefftz
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
    src/boxintegral.cpp
    src/tp0fespace.cpp
    #src/airy.cpp #for testing, requires boost
  )
target_compile_definitions(ngstrefftz PRIVATE NGSTREFFTZ_EXPORTS)
target_compile_definitions(ngstrefftz PRIVATE NGSTENT_EXPORTS)

# enable compiler warinigs, when building with g++ or clang++
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(ngstrefftz PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wno-vla # ngsolve uses VLAs anyway
    )
    set_source_files_properties( 
        ${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src/tents.cpp 
        PROPERTIES COMPILE_FLAGS "-w" #suppress all warnings
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(ngstrefftz PRIVATE
        -Wall
        -Wextra
        -pedantic
        -Wtype-limits
        -Wuninitialized
    )
    set_source_files_properties( 
        ${CMAKE_SOURCE_DIR}/external_dependencies/ngstents/src/tents.cpp 
        PROPERTIES COMPILE_FLAGS "-w" #suppress all warnings
    )
endif()

### user options
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} CACHE PATH "build dir")
set(NGSTREFFTZ_USE_GTEST OFF CACHE BOOL "Try to use gtest")
set(NGSTREFFTZ_USE_LAPACK ${NGSOLVE_USE_LAPACK} CACHE BOOL
    "Use the BLAS/LAPACK backend exported by NGSolve")
set(NGSTREFFTZ_USE_MKL ${NGSOLVE_USE_MKL})

### use the LAPACK backend exported by NGSolve hopefully...
if(NGSTREFFTZ_USE_LAPACK)
    if(NOT NGSOLVE_USE_LAPACK)
        message(FATAL_ERROR "NGSTREFFTZ_USE_LAPACK=ON, but NGSolve was built without LAPACK")
    elseif(NOT TARGET ngs_lapack)
        message(FATAL_ERROR "NGSolve reports LAPACK support, but did not export the ngs_lapack target")
    endif()
    message(STATUS "Linking NGSTrefftz against the NGSolve ngs_lapack target")
    target_link_libraries(ngstrefftz PRIVATE ngs_lapack)
    target_compile_definitions(ngstrefftz PRIVATE NGSTREFFTZ_USE_LAPACK)
else()
    message(WARNING "LAPACK disabled for NGSTrefftz, restricted usage!")
endif()

### doc
option(NGSTREFFTZ_BUILD_DOC "Build documentation" OFF)
if (NGSTREFFTZ_BUILD_DOC)
    find_package(Doxygen REQUIRED)
    set(DOXYGEN_IN ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in)
    set(DOXYGEN_OUT ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile)
    configure_file(${DOXYGEN_IN} ${DOXYGEN_OUT} @ONLY)
    add_custom_target( doxygen
        COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYGEN_OUT}
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen"
        VERBATIM )
endif (NGSTREFFTZ_BUILD_DOC)

### install
set_target_properties(ngstrefftz PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz)
file(COPY ${CMAKE_SOURCE_DIR}/src/__init__.py DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/ngstrefftz)
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/ngstrefftz)
install(TARGETS ngstrefftz DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})
ngsolve_generate_stub_files(${addon_name})

message("${BoldCyan}
------------------------------------------------------------------------
  Build type: ${CMAKE_BUILD_TYPE}
  Compiler: ${CMAKE_CXX_COMPILER}
  Flags: ${flags}

  Enabled functionality:
    LAPACK: ............ ${NGSTREFFTZ_USE_LAPACK}
    MKL: ............... ${NGSTREFFTZ_USE_MKL}
    Doxygen: ........... ${NGSTREFFTZ_BUILD_DOC}
    GTest: ............. ${NGSTREFFTZ_USE_GTEST}

With 'make install' the package will be installed to: ${CMAKE_INSTALL_PREFIX}
Make sure to add it to the python path:
-->   export PYTHONPATH=${CMAKE_INSTALL_PREFIX}/${python3_library_dir}:$PYTHONPATH
------------------------------------------------------------------------
\n${ColourReset}")

### tests
enable_testing()
# set up python integration tests
file(COPY
    ${CMAKE_SOURCE_DIR}/test/dg.py
    ${CMAKE_SOURCE_DIR}/test/embt.py
    ${CMAKE_SOURCE_DIR}/test/trefftz.py
    ${CMAKE_SOURCE_DIR}/test/tents.py
    ${CMAKE_SOURCE_DIR}/test/conforming_trefftz.py
    ${CMAKE_SOURCE_DIR}/test/boxint.py
    DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing)
add_test(NAME embtrefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/embt.py)
add_test(NAME conforming_trefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/conforming_trefftz.py)
add_test(NAME trefftz COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/trefftz.py)
add_test(NAME tents COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/tents.py)
add_test(NAME box COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/boxint.py)
#WORKING_DIRECTORY ${ CMAKE_CURRENT_SOURCE_DIR }
set_tests_properties(embtrefftz trefftz tents
    PROPERTIES ENVIRONMENT PYTHONPATH=${CMAKE_LIBRARY_OUTPUT_DIRECTORY}:${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON}:$ENV{PYTHONPATH})

if(NGSTREFFTZ_USE_GTEST)
    find_package(GTest REQUIRED)
    include(GoogleTest)
    if(CMAKE_VERSION VERSION_GREATER_EQUAL "3.18")
        find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Embed)
    endif()
    add_executable(trefftz_gtest ${CMAKE_SOURCE_DIR}/test/embtrefftz_gtest.cpp)

    target_include_directories(trefftz_gtest BEFORE PRIVATE
        $<TARGET_PROPERTY:ngsolve,INTERFACE_INCLUDE_DIRECTORIES>
        ${CMAKE_SOURCE_DIR}/src
    )
    target_link_libraries(trefftz_gtest PRIVATE
        GTest::gtest_main
        Python3::Python
        ngsolve
        ngstrefftz
    )
    if(TARGET ngs_lapack)
        target_link_libraries(trefftz_gtest PRIVATE ngs_lapack)
    endif()
    set_target_properties(trefftz_gtest PROPERTIES
        PREFIX ""
        CXX_STANDARD 20
        CXX_STANDARD_REQUIRED ON
    )
    gtest_add_tests(TARGET trefftz_gtest)
endif()

# make uninstall
add_custom_target("uninstall" COMMENT "Uninstall installed files")
add_custom_command(
    TARGET "uninstall"
    POST_BUILD
    COMMENT "Uninstall files with install_manifest.txt"
    COMMAND xargs rm -vf < install_manifest.txt || echo Nothing in
            install_manifest.txt to be uninstalled!
)
