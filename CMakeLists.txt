cmake_minimum_required(VERSION 3.8)
project(bcpt)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES bcpt.cpp)
add_executable(bcpt ${SOURCE_FILES})