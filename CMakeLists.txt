cmake_minimum_required(VERSION 2.8)
project(pasi)

# compiler options

add_definitions(-std=c++11)
add_definitions(-Wall)


# CoDiPack
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/codi/include" )

# boost/functional
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/boost/functional/include" )

# Eigen
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/eigen" )

# VTK
find_package(VTK 7 REQUIRED)
include(${VTK_USE_FILE})
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

add_subdirectory( src )

file(COPY data DESTINATION ${CMAKE_BINARY_DIR} )
