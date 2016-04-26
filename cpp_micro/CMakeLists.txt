project(pasi)
cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)
add_definitions(-Wall)

# VTK
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

#CoDiPack
include_directories( "${CMAKE_CURRENT_SOURCE_DIR}/codi/include" )

set( SRC
    main.cpp
    CParticleSystem.cpp
    CParticleRenderer.cpp
)

set( HEADERS
    CParticleSystem.hpp
    CParticleRenderer.hpp
)

add_executable(${PROJECT_NAME} ${SRC} ${HEADERS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

