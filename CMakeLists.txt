#
#-----------------------------#
#
cmake_minimum_required(VERSION 3.15)
PROJECT(fada)
#
#-----------------------------#
#
SET(CMAKE_BUILD_TYPE "Release")
#
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_C_FLAGS "-O3")
SET(CMAKE_CXX_FLAGS "-O3")
##
#
find_package(Armadillo REQUIRED)
message("ARMADILLO_INCLUDE_DIRS ${ARMADILLO_INCLUDE_DIRS}")

find_library(SUPERLU_LIBRARY NAMES "superlu")
message("SUPERLU_LIBRARY ${SUPERLU_LIBRARY}")
FIND_LIBRARY(UMFPACK_LIBRARIES NAMES umfpack libumfpack UMFPACK libUMFPACK
               PATH_SUFFIXES ${CMAKE_LIBRARY_ARCHITECTURE}
               HINTS ${CMAKE_INSTALL_PREFIX}/External/lib)
if(NOT UMFPACK_LIBRARIES)
    message(FATAL_ERROR "*** umfpack library not find (consider installing suite-sparse)")
endif()

INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
AUX_SOURCE_DIRECTORY(Fada SRC)
#
ADD_LIBRARY(fadalib ${SRC})
#
ADD_EXECUTABLE(fada fadamain.cpp)
TARGET_LINK_LIBRARIES(fada fadalib  ${UMFPACK_LIBRARIES} ${ARMADILLO_LIBRARIES} )
#
ADD_EXECUTABLE(testarma testarma.cpp)
TARGET_LINK_LIBRARIES(testarma fadalib  ${UMFPACK_LIBRARIES} ${ARMADILLO_LIBRARIES})
#
find_package(pybind11 REQUIRED)
message("pybind11 ${ARMADILLO_INCLUDE_DIRS}")
AUX_SOURCE_DIRECTORY(pybind SRCPY)
pybind11_add_module(pyfada  ${SRCPY} ${SRC})
target_link_libraries(pyfada PRIVATE pybind11::module ${UMFPACK_LIBRARIES} ${ARMADILLO_LIBRARIES})
#
#-- copy .py files
#
FILE(GLOB PYBINFILES ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
FOREACH(PYBINFILE ${PYBINFILES})
  get_filename_component(FILE ${PYBINFILE} NAME)
  configure_file(${PYBINFILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILE} COPYONLY)
ENDFOREACH(PYBINFILE)
