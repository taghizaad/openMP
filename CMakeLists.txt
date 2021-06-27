cmake_minimum_required(VERSION 3.19)
project(openMP C)

set(CMAKE_C_STANDARD 99)

add_executable(piProgram piProgram.c)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
    message("OPENMP FOUND")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
add_executable(main main.c)
add_executable(piProgramFalseSharingWithPadding piProgramFalseSharingWithPadding.c)
add_executable(piProgramFalseSharingWithCritical piProgramFalseSharingWithCritical.c)