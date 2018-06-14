project(trefftzngs)

cmake_minimum_required(VERSION 3.1)

find_package(NGSolve CONFIG REQUIRED
  HINTS $ENV{NETGENDIR}/.. /opt/netgen/ /Applications/Netgen.app/Contents/Resources/CMake C:/netgen
)

# check if CMAKE_INSTALL_PREFIX is set by user, if not install in NGSolve python dir
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
  set(CMAKE_INSTALL_PREFIX ${NGSOLVE_INSTALL_DIR}/${NGSOLVE_INSTALL_DIR_PYTHON} CACHE PATH "Install dir" FORCE)
endif(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)

add_ngsolve_python_module(trefftzngs trefftzngspy.cpp
	trefftz/helpers.cpp
	trefftz/trefftzelement.cpp
	trefftz/diffopmapped.hpp
	trefftz/trefftzfespace.cpp
	trefftz/scalarmappedelement.cpp
  trefftz/clipcoefficientfunction.cpp
  trefftz/trefftzcoefficient.hpp
  trefftz/evolvetent.hpp
  )

message("With 'make install' the python package will be installed to: ${CMAKE_INSTALL_PREFIX}")
install(TARGETS trefftzngs DESTINATION .)
