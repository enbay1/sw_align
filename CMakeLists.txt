cmake_minimum_required(VERSION 2.8)

project(sw_align)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 11)
SET(CMAKE_CXX_FLAGS "-O3")

file(GLOB_RECURSE SOURCES src/*.cpp src/*.c ext/src/*.c ext/src/*.cpp)
file(GLOB_RECURSE HEADERS src/*.hpp src/*.h ext/include/*.h ext/include/*.hpp)
add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
