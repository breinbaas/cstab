SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
cmake_minimum_required(VERSION 3.1)
project(pstab 
    VERSION 0.1
    DESCRIPTION "Cpp implementation of slope stability calculations"
    LANGUAGES CXX
)
add_library(clipper2lib STATIC clipper2/clipper.engine.cpp clipper2/clipper.engine.h)
target_include_directories(clipper2lib PUBLIC clipper2)
target_compile_features(clipper2lib PUBLIC cxx_std_11)
add_executable(pstab main.cpp)
target_link_libraries(pstab PUBLIC clipper2lib)