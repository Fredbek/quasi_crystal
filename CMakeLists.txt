cmake_minimum_required(VERSION 3.22.1)


project(MonteCarlo)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -DNDEBUG -fomit-frame-pointer -march=native")
add_executable(${PROJECT_NAME} main.cpp)
