cmake_minimum_required(VERSION 3.6)
project(laplace_equation)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  ${OpenMP_C_FLAGS}")
endif()

set(SOURCE_FILES main.cpp)
add_executable(laplace_equation ${SOURCE_FILES})