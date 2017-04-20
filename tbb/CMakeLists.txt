set(target laplace_equation_tbb)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

find_package(TBB REQUIRED)
if(TBB_FOUND)
    include_directories(${TBB_INCLUDE_DIRS})
endif()


set(SOURCE_FILES main.cpp)
add_executable(${target} ${SOURCE_FILES})
target_link_libraries(${target} ${TBB_LIBRARIES})
target_compile_definitions(${target} PRIVATE ${TBB_DEFINITIONS})