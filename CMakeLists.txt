#----------------------------------------------------------
#
cmake_minimum_required(VERSION 3.15)
PROJECT(fada)
#
#----------------------------------------------------------
SET(CMAKE_BUILD_TYPE "Release")
# SET(CMAKE_BUILD_TYPE "Debug")
#
SET(CMAKE_CXX_STANDARD 14)
if (CMAKE_BUILD_TYPE STREQUAL "Debug")
  # add_definitions(-DARMA_EXTRA_DEBUG)
	add_definitions(-DDEBUG)
elseif (CMAKE_BUILD_TYPE STREQUAL "Release")
  SET(CMAKE_C_FLAGS "-O3")
  SET(CMAKE_CXX_FLAGS "-O3")
  add_definitions(-DARMA_NO_DEBUG)
endif()
ADD_DEFINITIONS(-DARMA_USE_HDF5)

#----------------------------------------------------------
SET(PACKAGES
  HDF5
  Armadillo
  superlu
  amgcl
)
SET(LINK_LIBRARIES "")
SET(INCLUDE_DIRECTORIES "")
FOREACH(PACKAGE ${PACKAGES})
  FIND_PACKAGE(${PACKAGE} REQUIRED)
  string(TOUPPER ${PACKAGE} PACKAGE_UPPER)
  list(APPEND LINK_LIBRARIES ${${PACKAGE_UPPER}_LIBRARIES})
  list(APPEND INCLUDE_DIRECTORIES ${${PACKAGE_UPPER}_INCLUDE_DIRS})
ENDFOREACH(PACKAGE)
#----------------------------------------------------------
# packages without "FindXXX.cmake"
SET(PACKAGES
  umfpack
)
FOREACH(PACKAGE ${PACKAGES})
  string(TOUPPER ${PACKAGE} PACKAGE_UPPER)
  FIND_LIBRARY(UMFPACK_LIBRARIES NAMES ${PACKAGE} lib${PACKAGE} ${PACKAGE_UPPER}
                 PATH_SUFFIXES ${CMAKE_LIBRARY_ARCHITECTURE}
                 HINTS ${CMAKE_INSTALL_PREFIX}/External/lib)
  if(NOT ${PACKAGE_UPPER}_LIBRARIES)
      message(FATAL_ERROR "*** ${PACKAGE} library not find (consider installing suite-sparse)")
  endif()
  list(APPEND LINK_LIBRARIES ${${PACKAGE_UPPER}_LIBRARIES})
ENDFOREACH(PACKAGE)

MESSAGE("{LINK_LIBRARIES} ${LINK_LIBRARIES}")
MESSAGE("{INCLUDE_DIRECTORIES} ${INCLUDE_DIRECTORIES}")
#
#----------------------------------------------------------
AUX_SOURCE_DIRECTORY(Fada SRC)
ADD_LIBRARY(fadalib ${SRC})
TARGET_INCLUDE_DIRECTORIES(fadalib PRIVATE ${INCLUDE_DIRECTORIES})
#
ADD_EXECUTABLE(fada fadamain.cpp)
TARGET_INCLUDE_DIRECTORIES(fada PUBLIC ${INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(fada fadalib  ${LINK_LIBRARIES})
#
ADD_EXECUTABLE(testarma testarma.cpp)
TARGET_INCLUDE_DIRECTORIES(testarma PUBLIC ${INCLUDE_DIRECTORIES})
TARGET_LINK_LIBRARIES(testarma fadalib  ${LINK_LIBRARIES})
#
#----------------------------------------------------------
AUX_SOURCE_DIRECTORY(pybind SRCPY)
FIND_PACKAGE(carma CONFIG REQUIRED)
pybind11_add_module(pyfada  MODULE ${SRCPY} ${SRC})
TARGET_LINK_LIBRARIES(
    pyfada PUBLIC
    ${LINK_LIBRARIES}
    carma::carma
)
#----------------------------------------------------------
#-- copy .py files
#
FILE(GLOB PYBINFILES ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
FOREACH(PYBINFILE ${PYBINFILES})
  get_filename_component(FILE ${PYBINFILE} NAME)
  configure_file(${PYBINFILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILE} COPYONLY)
ENDFOREACH(PYBINFILE)
