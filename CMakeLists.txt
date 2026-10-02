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

### ngstrefftz

set(addon_name ngstrefftz)
project(ngstrefftz)
include(${CMAKE_BINARY_DIR}/cmake_modules/ngsolve_addon.cmake)
set(NGSTREFFTZ_USE_NGSTENTS AUTO CACHE STRING "Build with ngstents support")
set_property(CACHE NGSTREFFTZ_USE_NGSTENTS PROPERTY STRINGS AUTO ON OFF)
if(NOT NGSTREFFTZ_USE_NGSTENTS STREQUAL "OFF")
  if(NOT ngstents_DIR)
    execute_process(COMMAND ${Python3_EXECUTABLE} -m ngstents.config
      OUTPUT_VARIABLE ngstents_DIR OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  endif()
  find_package(ngstents CONFIG QUIET)
  if(ngstents_FOUND)
    set(NGSTREFFTZ_USE_NGSTENTS ON)
  elseif(NGSTREFFTZ_USE_NGSTENTS STREQUAL "ON")
    message(FATAL_ERROR "NGSTREFFTZ_USE_NGSTENTS=ON but ngstents was not found")
  else()
    set(NGSTREFFTZ_USE_NGSTENTS OFF)
    message(WARNING "ngstents not found, building without tent-pitching support")
  endif()
endif()
add_ngsolve_addon(ngstrefftz
    src/python_trefftz.cpp
    src/diffopmapped.hpp
    src/scalarmappedfe.cpp
    src/planewavefe.cpp
    src/trefftzfespace.cpp
    src/specialcoefficientfunction.cpp
    src/specialintegrator.cpp
    src/embtrefftz.cpp
    src/lapack_loader.cpp
    src/monomialfespace.cpp
    src/condensedg.cpp
    src/pufe.cpp
    src/pufespace.cpp
    src/boxintegral.cpp
    src/tp0fespace.cpp
    #src/airy.cpp #for testing, requires boost
  )
target_compile_definitions(ngstrefftz PRIVATE NGSTREFFTZ_EXPORTS)
if(NGSTREFFTZ_USE_NGSTENTS)
  target_sources(ngstrefftz PRIVATE src/twavetents.cpp src/mesh1dtents.cpp)
  target_compile_definitions(ngstrefftz PRIVATE NGSTREFFTZ_USE_NGSTENTS)
  target_link_libraries(ngstrefftz PRIVATE ngstents::core)
  ngsolve_addon_set_relative_rpath(ngstrefftz "../ngstents")
endif()
# add dynamic loader library to link libraries (on linux)
if(CMAKE_DL_LIBS)
    target_link_libraries(ngstrefftz PRIVATE ${CMAKE_DL_LIBS})
endif()

# enable compiler warinigs, when building with g++ or clang++
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    target_compile_options(ngstrefftz PRIVATE
        -Wall
        -Wextra
        -Wpedantic
        -Wno-vla # ngsolve uses VLAs anyway
    )
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(ngstrefftz PRIVATE
        -Wall
        -Wextra
        -pedantic
        -Wtype-limits
        -Wuninitialized
    )
endif()

### user options
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} CACHE PATH "build dir")
set(NGSTREFFTZ_USE_GTEST OFF CACHE BOOL "Try to use gtest")

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
    Doxygen: ........... ${NGSTREFFTZ_BUILD_DOC}
    GTest: ............. ${NGSTREFFTZ_USE_GTEST}
    ngstents: .......... ${NGSTREFFTZ_USE_NGSTENTS}

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
set(tests_with_pythonpath embtrefftz trefftz)
if(NGSTREFFTZ_USE_NGSTENTS)
  add_test(NAME tents COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/tents.py)
  list(APPEND tests_with_pythonpath tents)
endif()
add_test(NAME box COMMAND python3 -m doctest ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/Testing/boxint.py)
#WORKING_DIRECTORY ${ CMAKE_CURRENT_SOURCE_DIR }
set_tests_properties(${tests_with_pythonpath}
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
