cmake_minimum_required (VERSION 3.10.2)
project (ECDLP-Diem)

set (CMAKE_CXX_STANDARD 17)
set (CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Debug)
endif()

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
set (CMAKE_CXX_FLAGS_DEBUG "-fsanitize=undefined -g")
set (CMAKE_CXX_FLAGS_RELEASE "-march=native -O3")

include_directories ("include")

enable_testing ()

add_executable (ff-basic src/finitefield.cpp test/ff-basic.cpp)
add_executable (poly-basic test/poly-basic.cpp)
add_executable (fibonacci test/fibonacci.cpp)
add_executable (semaev test/semaev.cpp)
add_executable (determinant test/determinant.cpp)
add_executable (weil src/finitefield.cpp test/weil.cpp)
add_executable (grobner src/finitefield.cpp test/grobner.cpp)

add_test (ff-basic-test ff-basic)
add_test (poly-basic-test poly-basic)
add_test (semaev-test semaev)
add_test (determinant-test determinant)
add_test (weil-test weil)
add_test (grobner-test grobner)
add_test (fibonacci-test fibonacci)

