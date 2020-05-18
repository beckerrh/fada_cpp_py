#
#-----------------------------#
#
cmake_minimum_required(VERSION 3.17)
PROJECT(fada)
#
#-----------------------------#
#
#
MESSAGE("....... bonjour ...... fada ....")
#
SET(CMAKE_C_FLAGS "-O3")
SET(CMAKE_CXX_FLAGS "-O3")
##
#
find_package(Armadillo REQUIRED)
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
AUX_SOURCE_DIRECTORY(Fada SRC)
#
ADD_LIBRARY(fadalib ${SRC})
#
ADD_EXECUTABLE(fada main.cpp)
TARGET_LINK_LIBRARIES(fada fadalib ${ARMADILLO_LIBRARIES})
#
ADD_EXECUTABLE(testarma testarma.cpp)
TARGET_LINK_LIBRARIES(testarma fadalib ${ARMADILLO_LIBRARIES})
#
find_package(pybind11 REQUIRED)
pybind11_add_module(pyfada  pybind/operator.cpp ${SRC})
target_link_libraries(pyfada PRIVATE pybind11::module fadalib ${ARMADILLO_LIBRARIES})
#
FILE(GLOB PYBINFILES ${CMAKE_CURRENT_SOURCE_DIR}/*.py)
FOREACH(PYBINFILE ${PYBINFILES})
  get_filename_component(FILE ${PYBINFILE} NAME)
  configure_file(${PYBINFILE} ${CMAKE_CURRENT_BINARY_DIR}/${FILE} COPYONLY)
ENDFOREACH(PYBINFILE)
