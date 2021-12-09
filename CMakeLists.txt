cmake_minimum_required(VERSION 3.19)
project(openMP)

set(CMAKE_C_STANDARD 99)
FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(Matlab)
matlab_add_mex(
        NAME openMp
        SRC smw.c arrayProduct.c
)

add_executable(piProgramFalseSharingWithPadding piProgramFalseSharingWithPadding.c)
add_executable(piProgramFalseSharingWithCritical piProgramFalseSharingWithCritical.c)
add_executable(piProgramWithPragma piProgramWithPragma.c)
add_executable(piProgramWithForAndReduction piProgramWithForAndReduction.c)
add_executable(matrix matrix.c)
add_executable(SerialPiProgram SerialPiProgram.c)
add_executable(matrix_mul_openmp matrix_mul_openmp.c)
add_executable(modifiedSmw modifiedSmw.c)
add_executable(modifiedSmwOpt6Switch smw6switch/modifiedSmwOpt6Switch.c)
add_executable(modifiedSmwOpt18Switch smw18switch/modifiedSmwOpt18Switch.c)
add_executable(test test.c)