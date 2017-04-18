cmake_minimum_required(VERSION 3.6)
set(target laplace_equation_openmp)


set(SOURCE_FILES main.cpp)
add_executable(${target} ${SOURCE_FILES})