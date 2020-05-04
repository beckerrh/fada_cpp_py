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
ADD_EXECUTABLE(fada ${SRC} main.cpp)
TARGET_LINK_LIBRARIES(fada  ${ARMADILLO_LIBRARIES})
#
ADD_EXECUTABLE(testarma ${SRC} testarma.cpp)
TARGET_LINK_LIBRARIES(testarma  ${ARMADILLO_LIBRARIES})
#
find_package(pybind11 REQUIRED)
pybind11_add_module(pyfada  pybind/operator.cpp ${SRC})
target_link_libraries(pyfada PRIVATE pybind11::module ${ARMADILLO_LIBRARIES})
#
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/test.py ${CMAKE_CURRENT_BINARY_DIR}/test.py COPYONLY)
