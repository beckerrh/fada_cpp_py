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
