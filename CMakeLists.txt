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

add_executable (ff-basic test/ff-basic.cpp)

add_test (maintest ff-basic)

