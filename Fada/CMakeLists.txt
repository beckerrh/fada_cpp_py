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

#
AUX_SOURCE_DIRECTORY(. SRC)
#
INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
ADD_EXECUTABLE(fada ${SRC})
#
TARGET_LINK_LIBRARIES(fada  ${ARMADILLO_LIBRARIES})
