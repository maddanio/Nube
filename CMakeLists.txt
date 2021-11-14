cmake_minimum_required(VERSION 3.10)
project(Nube)
add_executable(
    Nube
    src/Axis.cpp
    src/Camera.cpp
    src/Canvas.cpp
    src/Cloud.cpp
    src/Cumulus.cpp
    src/Exception.cpp
    src/FluidCPU.cpp
    src/Lsystem.cpp
    src/Main.cpp
    src/Metacloud.cpp
    src/Model.cpp
    src/Morph.cpp
    src/Mountain.cpp
    src/OBJLoader.cpp
    src/PrecomputeCPU.cpp
    src/Shader.cpp
    src/Tool.cpp    
)
target_include_directories(
    Nube
    PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/dependencies
    /opt/homebrew/include
)
set_property(TARGET Nube PROPERTY CXX_STANDARD 17)
find_package(PkgConfig REQUIRED)
pkg_check_modules(
    GLFW
    REQUIRED
    IMPORTED_TARGET
    glfw3
)
pkg_check_modules(
    GLM
    REQUIRED
    IMPORTED_TARGET
    glm
)
find_package(OpenGL REQUIRED)
target_link_libraries(Nube PRIVATE OpenGL::GL PkgConfig::GLFW PkgConfig::GLM)
target_compile_options(Nube PRIVATE -Wno-deprecated-declarationsdeprecated-declarations)
