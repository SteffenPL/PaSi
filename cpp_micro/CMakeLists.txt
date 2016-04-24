project(pasi)
cmake_minimum_required(VERSION 2.8)

add_definitions(-std=c++11)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

set( SRC
    main.cpp
)

add_executable(${PROJECT_NAME} ${SRC})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})

