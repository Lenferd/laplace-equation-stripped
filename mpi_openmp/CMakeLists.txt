set(target laplace_equation_mpi)

find_package( MPI REQUIRED )
include_directories( ${MPI_INCLUDE_PATH} )


set(SOURCE_FILES main.cpp)
add_executable(${target} ${SOURCE_FILES})
target_link_libraries(${target} ${MPI_LIBRARIES})