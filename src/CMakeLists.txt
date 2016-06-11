cmake_minimum_required(VERSION 2.8)


set( SRC
    main.cpp
    ConfigManager.cpp
    ForceField.cpp
    State.cpp
    MoleculeViewer.cpp
)

set( HEADERS
    ConfigManager.hpp
    ForceField.hpp
    State.hpp
    MoleculeViewer.hpp
)

add_executable(${PROJECT_NAME} ${SRC} ${HEADERS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} )


