cmake_minimum_required(VERSION 2.8)


set( SRC
    ConfigManager.cpp
    ForceField.cpp
    MoleculeState.cpp
    MoleculeViewer.cpp
    Domain.cpp
)

set( HEADERS
    ConfigManager.hpp
    ForceField.hpp
    MoleculeState.hpp
    MoleculeViewer.hpp
    Domain.hpp
)

set( SIM_SRC
    main.cpp )

add_executable(${PROJECT_NAME} ${SIM_SRC} ${SRC} ${HEADERS})
target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES} )


# molecule editor

# Find includes in corresponding build directories
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

# Find the QtWidgets library
find_package(Qt5Widgets)

set( GUI_SRC
    main_editor.cpp
    MoleculeEditor.cpp )

set( GUI_HEADERS
    MoleculeEditor.hpp )

qt5_wrap_ui(
    MOC_UI
    MoleculeEditor.ui)


add_executable(${PROJECT_NAME}_editor ${GUI_SRC} ${SRC} ${GUI_HEADERS} ${HEADERS} ${MOC_UI} )
target_link_libraries(${PROJECT_NAME}_editor ${VTK_LIBRARIES} Qt5::Widgets Qt5::Core Qt5::Gui )


execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory  ${CMAKE_SOURCE_DIR}/data
                                                            ${CMAKE_BINARY_DIR}/data)
