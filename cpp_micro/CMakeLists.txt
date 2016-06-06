project(pasi)
cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)
add_definitions(-Wall)

# VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules)


# CoDiPack
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/codi/include" )

# boost/functional
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/functional/include" )

# Eigen
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/eigen3" )

# open babel
# include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/openbabel/include" )
find_package(OpenBabel2 REQUIRED)
include_directories(${OPENBABEL2_INCLUDE_DIR})


set( SRC
    main.cpp
    ParticleSystem.cpp
    ParticleRenderer.cpp
    ConfigManager.cpp
    Grid.cpp
    #FunctionPlotter.cpp
)

set( HEADERS
    ParticleSystem.hpp
    ParticleRenderer.hpp
    ConfigManager.hpp
    Grid.hpp
    #FunctionPlotter.hpp
)

add_executable(${PROJECT_NAME} ${SRC} ${HEADERS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} ${OPENBABEL2_LIBRARIES})

